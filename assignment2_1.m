%%% EXERCISE 1
%%% Relationship between the conduction velocity and the conductance.
% Problem setting:
% 	Cable length
L= 4; %cm
% 	Total Simulation Time
totalt = 400; %ms
% 	Space discretization
dx = [0.01 0.02 0.025]; %cm
% 	Time step
dt = 0.02; % ms
% 	Output Frecuenc
outfreq = 10;
% 	Numerical Method: Implicit Euler, Option 1
num_meth = 1;
% Run simulations for different values of the conductance, 
sigma=[0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003];
sigmatxt = {'0005', '001', '0015', '002','0025', '003'};
dxtxt = {'02','01','025'};
% and plot the value of the conductance versus the Conduction velocity.
% 
% 1) What can you say about the relationship between these two variables
%    (Hint: Use a loglog plot)
% 2) Repeat the same procedure for dx=0.01cm and dx=0.025cm.
%  Comment on the results?. 
% 
% Call the program a follows:
% >> flag_exit=fd1D_FK(1, 10);
% 
% Note: Remmember to update the name of the output file if you want to 
% keep the result of all simulations

%% SIMULATION
for i = 1:length(sigma)
    for j = 1:length(dx)
        fname_out = strcat("output_dx",string(dxtxt{j}),"_sgm",string(sigmatxt{i}));
        flag_exit = fd1D_FK(num_meth, outfreq, L, totalt, dx(j), sigma(i), fname_out, dt);
        postprocess(fname_out);
    end
end
%% PLOTTING
figure;
subplot(131);
CVs=zeros(1,length(sigma));
for i = 1:length(sigma)
   load(strcat("output_dx",string(dxtxt{1}),"_sgm",string(sigmatxt{i})));
   CVs(i) = CV;
end
interpreterlatex
loglog(sigma,CVs,'.-'); grid on;
title("\textbf{dx = 0.01cm}");
xlabel("Conductance $[cm^2/ms]$");
ylabel("Conduction Velocity $[cm/s]$");

subplot(132);
CVs=zeros(1,length(sigma));
for i = 1:length(sigma)
   load(strcat("output_dx",string(dxtxt{2}),"_sgm",string(sigmatxt{i})));
   CVs(i) = CV;
end
loglog(sigma,CVs,'.-'); grid on;
title("\textbf{dx = 0.02cm}");
xlabel("Conductance $[cm^2/ms]$");
ylabel("Conduction Velocity $[cm/s]$");

subplot(133);
CVs=zeros(1,length(sigma));
for i = 1:length(sigma)
   load(strcat("output_dx",string(dxtxt{3}),"_sgm",string(sigmatxt{i})));
   CVs(i) = CV;
end
loglog(sigma,CVs,'.-'); grid on;
title("\textbf{dx = 0.025cm}");
xlabel("Conductance $[cm^2/ms]$");
ylabel("Conduction Velocity $[cm/s]$");




