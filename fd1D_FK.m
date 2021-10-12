function flag_exit=fd1D_FK(num_meth, outfreq, L, totalt, dx, sgm, fname_out, dt)
%
% 1D Finite element and Finite difference program for 
% the monodomain equation. It uses Operator Splitting 
% with Implicit/Explicit Euler for the diffusion equation and 
% Explicit Euler for the Ionic Currents
%
% The funtion receives as an input parameter the numerical scheme to be 
% used (num_meth)
%  0 finite differences (explicit scheme)
%  1 finite differences (implicit scheme)
%  2 finite differences (Crack-Nicholson scheme)
% 
% The second parameter referes to the frequency at which 
% the potential will be saved in the output database (outfreq)
%
% The code outputs file output.mat with the following variables:
%
%  x  : A vector of length nnd, the total number of nodes in the model, 
%       with position of the nodes (discretized cable)
% Vsol: A N by (nnd+1) matrix with the AP along the cable. The first column
%       has the time, and N is the total number of frames that have been
%       saved during the simulation.
% Jion: A N by (nnd+1) matrix with the transmembrane ionic current along
%       the cable. The first column has the time, and N is the total number
%       of frames that have been saved during the simulation.
%
% The file is read by program postprocess.m that perform the postprocessing
%
% The discretized system is written, at increment k, as follows
%
%  B1*V^(k+1) = V^* + B2*V^(k)
%
% where B1 and B2 are matrices
%
% B1 = Id     for Explicit scheme, num_meth=0
% B1 = Id-A   for Implicit scheme, num_meth=1
% B1 = Id-A/2 for Cranck-Nicholson scheme, num_meth=2
%
% B2 = A     for Explicit scheme, num_meth=0
% B2 = 0     for Implicit scheme, num_meth=1
% B2 = A/2   for Cranck-Nicholson scheme, num_meth=2
%
% Id is the nnd*nnd identity matrix, and A is the finite difference matrix
% associated with the second order difference operator for the diffusion
% term of the cable equation
%

%
% ----------------------------------------------
%   MODIFIABLE PARAMETERS
% ----------------------------------------------
% General Data
%
% L=4.0;           % total length, cm
% dx=0.01;         % mesh size, cm
% sgm=0.001; %1.0e-3;      % normalized tissue conductivity, cm^2/ms
Cm=1.0;          % normalized cell membrane capacitance, dimensionless
% fname_out='output_dx01.mat'; % output filename
%
% Time integration
%
% dt=0.02;         % stable time step to integrate ionic model
% totalt = 400;    % Total simulation time, ms 
%
% Stimulation
%
dtstim=1.0;      % Stimulation duration, ms
tbegin=0.0;      % Begining of stimulus, ms
node_st=1;       % Node for stimulation
Stim_curr=-60;  % Stimulation current, -
BCL=2000;        % basic cycle length, ms
% ----------------------------------------------
%   END OF MODIFIABLE PARAMETERS
% ----------------------------------------------
% if(num_meth==0)                   % Time step is automatically set with
%  dt=min([0.02 0.45*Cm*dx^2/sgm]); % the explicit method to avoid numerical
%  outfreq = fix(0.02*outfreq/dt);  % instabilities. Sampling rate is changed
% end                               % accordingly
ntmp = fix((totalt)/dt);
flag_stm=1.0;
Istim=0.0;
invCm = 1.0/Cm;
print_sim_char(num_meth,dt,dtstim,tbegin,totalt,BCL,Stim_curr)
%
% ----------------------------------------------
%
% Space discretization
%
if(num_meth<=2)
  x=0.0:dx:L;                        % Nodal coordinates
  nnd = length(x);
  fprintf('The problem has %6d nodes\n',nnd);
%
% Assembling A matrix
%
  fprintf('Assembling finite difference matrices ...\n');
  [B1,B2]=assembling(num_meth,sgm,Cm,nnd,dx,dt);
else
  x=0.0;
  nnd = 1;
  fprintf('The problem has %6d nodes\n',nnd);
end
%
% ----------------------------------------------
%
% Dimensioning Database for postprocessing
%
% Forst column is time
%
N = length(1:outfreq:ntmp);
Vsol=zeros(N,nnd+1);         % Potential at each node
Jion=zeros(N,nnd+1);         % Ionic current at each node
contout=1;
%
% ----------------------------------------------
%
% initializing ionic models
%
p=mod_param();       % loading model parameters
V=initial_cond(nnd); % initializing cell models
%
fprintf('Initiating temporal integration\n');
t=0.0;
tend=tbegin+dtstim;
%
xsol=zeros(nnd,1);
Istim=zeros(nnd,1);

%
% Writing to Output file
%
Vsol(contout,:) = [t V.V'];
Jion(contout,:) = [t V.Itot'];
contout=contout+1;
%
tic
for i=1:ntmp
%
% Integrating cellular model in the nodes to obtain V* (Step 1)
%
  if (t>=tbegin && t<=tend)
    Istim(node_st)=Stim_curr;
    flag_stm=0.0;
  else
    Istim(node_st)=0.0;
    if(flag_stm==0.0)
      tbegin=tbegin+BCL;
      tend=tbegin+dtstim;
      fprintf('Next Stimulation at: %10.2f ms\n',tbegin);
      flag_stm=1.0;
    end
  end
  V=currents(V,p,dt,Istim);
  xsol = V.V - invCm*V.Itot*dt;      
%
% xsol continains the value of V*
% Solving for the propagation of electric signal (Step 2)
%
  switch(num_meth)
    case 0 % Explicit Finite differences
      V.V = xsol + B2*V.V;
    case 1 % Implicit Finite differences
      [V.V,flag]=pcg(B1,xsol,1.0e-8,100,[],[],xsol);
    case 2
      [V.V,flag]=pcg(B1,xsol+B2*V.V,1.0e-8,100,[],[],xsol);
  end
  xsol=V.V;
  t=t+dt;
  if (mod(i,500)==0)
    fprintf('Step: %10d Time: %15.4f Vmax: %15.4f\n',i,t,max(xsol));
    if((max(xsol)<0.4)&&(i>500)&&(i<=1000))
       fprintf('No electric activity is registered...\n');
       fprintf('Increase the stimulation current and rerun.\n');
       break;
    end
  end
  if (mod(i,outfreq)==0)
    Vsol(contout,:) = [t V.V'];
    Jion(contout,:) = [t V.Itot'];
    contout = contout + 1;
  end
end
toc
save(fname_out,'x','Vsol','Jion');
flag_exit=1;
end
%
%-------------------------------------------
%
function p=mod_param()
%
% This function returns the model parameters
% for the Bueno et al. Model
%
%
% Constants
%

      p.Vm=0.3; p.Vp=0.13; p.Vq=0.006; p.Vr=0.006;
%
%     GATE H (v)
%
      p.TAUHm1=60.0; p.TAUHm2=1150.0; p.TAUHp=1.4506;
%
%     GATE W (w)
%
      p.TAUWm1=60.0; p.TAUWm2=15.0; p.TAUWp=200.0; p.TAUWinf=0.07;
      p.Vw = 0.03; p.Wfstinf = 0.94; p.KW = 65.0;
%
%     GATE S (s)
%
      p.TAUS1=2.7342; p.TAUS2=16.0; p.KS=2.0994; p.Vs=0.9087;
%
%     CURRENT Ifi
%
      p.gfi=1.0/0.11; p.Vv=1.55;
%
%     CURRENT Iso
%
      p.TAUSo1=30.0181; p.TAUSo2=0.9957; p.KSo=2.0458; 
      p.Vso=0.65; p.TAUo1=400.00; p.TAUo2=6.0; p.Vo=0.0;
%
%     CURRENT Isi
%
      p.gsi=1.0/1.8875;
end

%-------------------------------------------

function V=initial_cond(nnd)
%
% This function sets the initial conditions
% for the Bueno et al. Model
%
% V is the main structure which contains the
% membrane potential, gate variable, and Total 
% current
%
V.V=zeros(nnd,1); 
V.Hh=ones(nnd,1); 
V.Ww=ones(nnd,1); 
V.Ss=zeros(nnd,1);
V.Itot=zeros(nnd,1);
end

%-------------------------------------------
%
function V=currents(V,par,dt,Istim)
%
% This function computes current and gate variable
% for the 4 variable Fenton and Karma Model
%
    m=(V.V>par.Vm);
    p=(V.V>par.Vp);
    q=(V.V>par.Vq);
    r=(V.V>par.Vr);
%
%     Gate Hh
%  
    Hinf=(V.V<par.Vq);
    TAUHm = (1.0-q)*par.TAUHm1 + q*par.TAUHm2;
    DHh=((1.0-m).*(Hinf-V.Hh))./TAUHm - m.*V.Hh/par.TAUHp;
%
%     Gate Ww 
%
    Winf = (1.0-r).*(1.0-V.V/par.TAUWinf) + r*par.Wfstinf;
    TAUWm = par.TAUWm1 + ...
          0.5*(par.TAUWm2-par.TAUWm1)*(1.0+tanh(par.KW*(V.V-par.Vw)));
    DWw = ((1.0-p).*(Winf-V.Ww))./TAUWm - p.*V.Ww/par.TAUWp;
%
%     Gate Ss 
%
    TAUS = (1.0-p)*par.TAUS1 + p*par.TAUS2;
    DSs = (0.5*(1.0+tanh(par.KS*(V.V-par.Vs)))-V.Ss)./TAUS;

    V.Hh= V.Hh+DHh*dt;
    V.Ww= V.Ww+DWw*dt;
    V.Ss= V.Ss+DSs*dt;
%
% Currents
%
    Ifi = zeros(length(V.V),1);
    Id1 = V.V>par.Vm;
    Ifi(Id1) = -(par.gfi*V.Hh(Id1).*(V.V(Id1)-par.Vm)).*(par.Vv-V.V(Id1));

    Isi = zeros(length(V.V),1);
    Id1 = (V.V>par.Vp);
    Isi(Id1) = -par.gsi*(V.Ww(Id1).*V.Ss(Id1));
    
    TauSo = ones(length(V.V),1)*par.TAUo2;
    Id1 = V.V>par.Vp;
    Id2 = V.V<par.Vr;
    TauSo(Id1)=par.TAUSo1 + ...
              (par.TAUSo2-par.TAUSo1)*0.5*(1.0+tanh(par.KSo*(V.V(Id1)-par.Vso)));
    TauSo(Id2) = par.TAUo1;
    Iso = (V.V-par.Vo)./TauSo;
    Iso(Id1) = 1./TauSo(Id1);
%  
    V.Itot = Ifi+Iso+Isi+Istim;
end
%
%--------------------------------------------------------------------------
%
% FINITE DIFFERENCE FUNCTIONS
%
%-------------------------------------------

function [B1,B2]=assembling(num_meth,sgm,Cm,nnd,dx,dt)
%
% Function for assembling the A matrix for the FD scheme
%
%  (flag*A+I)V^(k+1) = V^* + (1+flag)*A*V^k
%
%   flag: 0 (Explicit Euler)
%         1 (Implicit Euler)
%        0.5(Crank-Nicolson)
%
% num_meth: Numerical scheme
% sgm     : Conductivity
%  Cm     : Capacitance
% nnd     : Number of nodes
% dx      : element size
% dt      : time step
%
A=sparse(zeros(nnd));
scale = sgm*dt/(Cm*dx^2);
fprintf('Assembling matrix A\n');
for ii=1:nnd
  switch ii
    case 1
        A(ii,ii) = -2*scale;
        A(ii,ii+1) = 2*scale;
    case nnd
        A(ii,ii) = -2*scale;
        A(ii,ii-1) = 2*scale;
    otherwise
        A(ii,ii) = -2*scale;
        A(ii,ii+1) = scale;
        A(ii,ii-1) = scale;
  end
end
switch num_meth
    case 0
        B1 = [];
        B2 = A;
    case 1
        B1 = sparse(eye(nnd))-A;
        B2 = [];
    case 2
        B1 = sparse(eye(nnd))-0.5*A;
        B2 = 0.5*A;
    otherwise
        B1=[];
        B2=[];
end
fprintf('Assembling finished\n');
end


function print_sim_char(num_meth,dt,dtstim,tbegin,totalt,BCL,Istim)
%
%
 fprintf('Simulation Data:\n');
 fprintf('Solving Bueno et al. minimum model\n');
 switch num_meth
     case 0
       fprintf('FD with Explicit Euler scheme\n');
       fprintf('Time step calculated to avoid numerical instabilities\n');
     case 1
       fprintf('FD with Implicit Euler scheme\n');
     case 2
       fprintf('FD with Crank-Nicolson scheme\n');
     otherwise
       fprintf('Integrating AP model only\n');         
 end
 fprintf('Time increment       : %12.4e msec\n',dt);
 fprintf('Stim. duration       : %12.2f msec\n',dtstim);
 fprintf('First stimulation    : %12.2f msec\n',tbegin);
 fprintf('Total simulation time: %12.2f msec\n',totalt);
 fprintf('Basic Cycle Length   : %12.2f msec\n',BCL);
 fprintf('Stimulation current  : %12.2f mA\n',Istim);
end
