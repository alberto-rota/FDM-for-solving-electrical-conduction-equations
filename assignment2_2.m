%%% EXERCISE 2
%%%Exercise 2:  Calibrating conductance to achieve a given conduction velocity
% Using the result from the previous exercise for the relationship between 
% the Conduction Velocity and the tissue conductance, propose a method to 
% identify the required value of the conductance in order to achieve the 
% target conduction velocity given in the following.
% 
% Problem setting:
% 	Target velocity: 48 cm/s
targetcv = 48;
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

figure;
CVs=zeros(1,length(sigma));
for i = 1:length(sigma)
   load(strcat("output_dx",string(dxtxt{2}),"_sgm",string(sigmatxt{i})));
   CVs(i) = CV;
end
plot(sigma,CVs,'.-'); grid on;
title("\textbf{dx = 0.02cm}");
xlabel("Conductance $[cm^2/ms]$");
ylabel("Conduction Velocity $[cm/s]$");

% 1.	Describe the methodology used to solve the problem and report the
% value of the conductance. In case you propose an iterative method, report
% a table with the intermediate values showing the convergence to the solution.
% Note: Achieving the exact target velocity might not be possible in some 
% cases. Consider valid a result within less than 0.2% of the target value.
% Linear interpolation = Graph Lookup
targetsigma = interp1(CVs,sigma, targetcv);
hold on
scatter(targetsigma, targetcv, 'r');

% 2.	Calculate the maximum time step associated with the PDE for the 
% case of using an Explicit Euler integration scheme.

