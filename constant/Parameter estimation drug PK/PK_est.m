% =========================================================================
% PK_est.m
%
% PURPOSE (paper-guided):
%   Estimate pharmacokinetic (PK) parameters for EE and DNG used in the
%   integrated endocrine–PK model described in the paper.
%
%   The estimated parameters are later fixed and used in:
%     - main_model.m (forward simulation)
%     - circadian timing comparisons
%
% MODEL ASSUMPTIONS (as in the paper):
%   - Oral dosing
%   - Linear PK
%   - Multi-exponential disposition (closed-form solution)
%   - No circadian modulation is estimated here (handled separately)
%
% OUTPUT:
%   pkE : estimated EE PK parameters
%   pkP : estimated DNG PK parameters
% =========================================================================

function PK_est

% -------------------------------------------------------------------------
% [1] Top-level PK parameter estimation routine
%     Uses nonlinear least squares / optimization to fit concentration–time
%     data with a closed-form PK expression.
% -------------------------------------------------------------------------
parameterset = [80; 2; 40000; 1; 1];         % input initial guess


opts = optimset('TolX', 1e-12, 'Maxiter', 10000, 'MaxFunEvals', 10000);

% Run unconstrained optimization to estimate PK parameters.
% Parameters are optimized to minimize the mismatch between:
%   - observed concentration–time data
%   - model-predicted concentrations
[optimal_par, errorval, exitflag, ~] = fminsearch(@fminmerged, parameterset,opts);

parameterset = optimal_par.^2      % PK parameters

end


function [ls,sol] = fminmerged(parameterset)
ka = parameterset(1)^2;
k21 = parameterset(2)^2;
Vc = parameterset(3)^2;
alpha1 = parameterset(4)^2;
beta1 = parameterset(5)^2;
F =  0.90;         % 0.65 for EE

parameters = zeros(84,1);
parameters(57:84) = round(7/24,2);
parameters(1:21) = 3e4;
parameters(1+28:21+28) = 2e6;

x    = [0:0.01:28];

tend = 28;
tstep = 0.01;
xt = 0:tstep:tend; 

totalind   = ((tend)/tstep) + 1;

CE2 = zeros((tend), totalind);
CP4 = zeros((tend), totalind);
 

for la = 1:tend
    ind = round((la -1 + parameters(2*(tend)+la))*100 + 1);                 % index corresponding to dosing time
        N = (k21 - ka)/( (alpha1 - ka)*(beta1 - ka) );
        L = (k21 - alpha1)/( (ka - alpha1)*(beta1 - alpha1) );
        M = (k21 - beta1)/( (ka - beta1)*(alpha1 - beta1) );

        % Closed-form PK solution (sum of exponentials) for oral dosing.
        % This corresponds to the analytical PK expression used in the paper.
        CE2(la, ind:totalind)  =  ((ka*F*parameters(la))/Vc) * ( N*exp(-ka*(xt(1:(totalind +1-ind)))) + L*exp(-alpha1*(xt(1:(totalind +1-ind)))) + M*exp(-beta1*(xt(1:(totalind +1-ind)))) );        
        % Closed-form PK solution (sum of exponentials) for oral dosing.
        % This corresponds to the analytical PK expression used in the paper.
        CP4(la, ind:totalind)  =  ((ka*F*parameters(la + tend))/Vc) * ( N*exp(-ka*(xt(1:(totalind +1-ind)))) + L*exp(-alpha1*(xt(1:(totalind +1-ind)))) + M*exp(-beta1*(xt(1:(totalind +1-ind)))) );          
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalCE2 = zeros(1, totalind);
totalCP4 = zeros(1, totalind);

for su = 1:(tend)
    totalCE2 = totalCE2 + CE2(su, 1:totalind);
    totalCP4 = totalCP4 + CP4(su, 1:totalind);
end

totalCE2B = totalCE2;
totalCP4B = totalCP4;

% Load observed PK concentration–time data from the paper / supplement.
% Data typically consist of:
%   - sampling times
%   - measured plasma concentrations
load P4_PK_data.mat P4_PK_data             % load PK data

time_PK_data = P4_PK_data(:,1);            % use E2_PK_data for EE
PK_data      = P4_PK_data(:,2);

% Time points of blood sampling converted to time index
rnd1 = round((6.29+time_PK_data(1))*100+1);
rnd2 = round((6.29+time_PK_data(2))*100+1);
rnd3 = round((6.29+time_PK_data(3))*100+1);
rnd4 = round((6.29+time_PK_data(4))*100+1);
rnd5 = round((6.29+time_PK_data(5))*100+1);
rnd6 = round((6.29+time_PK_data(6))*100+1);
rnd7 = round((6.29+time_PK_data(7))*100+1);
rnd8 = round((6.29+time_PK_data(8))*100+1);
rnd9 = round((6.29+time_PK_data(9))*100+1);
rnd10 = round((6.29+time_PK_data(10))*100+1);
rnd11 = round((6.29+time_PK_data(11))*100+1);
rnd12 = round((6.29+time_PK_data(12))*100+1);
rnd13 = round((6.29+time_PK_data(13))*100+1);
rnd14 = round((6.29+time_PK_data(14))*100+1);
rnd15 = round((6.29+time_PK_data(15))*100+1);
rnd16 = round((6.29+time_PK_data(16))*100+1);
rnd17 = round((6.29+time_PK_data(17))*100+1);
rnd18 = round((6.29+time_PK_data(18))*100+1);

% Simulated DNG concentration
Con = [totalCP4B(rnd1); totalCP4B(rnd2); totalCP4B(rnd3); totalCP4B(rnd4); totalCP4B(rnd5); 
    totalCP4B(rnd6); totalCP4B(rnd7); totalCP4B(rnd8); totalCP4B(rnd9); totalCP4B(rnd10); 
    totalCP4B(rnd11); totalCP4B(rnd12); totalCP4B(rnd13); totalCP4B(rnd14); totalCP4B(rnd15); 
    totalCP4B(rnd16); totalCP4B(rnd17); totalCP4B(rnd18)];

% use for EE
% Con = [totalCE2B(rnd1); totalCE2B(rnd2); totalCE2B(rnd3); totalCE2B(rnd4); totalCE2B(rnd5); 
%     totalCE2B(rnd6); totalCE2B(rnd7); totalCE2B(rnd8); totalCE2B(rnd9); totalCE2B(rnd10); 
%     totalCE2B(rnd11); totalCE2B(rnd12); totalCE2B(rnd13); totalCE2B(rnd14); totalCE2B(rnd15); 
%     totalCE2B(rnd16); totalCE2B(rnd17); totalCE2B(rnd18)];

ls = (PK_data-Con)'*(PK_data-Con);


end

