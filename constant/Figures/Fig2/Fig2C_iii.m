% ======================================================================== 
% Fig2C_iii
%
% PURPOSE:
% Simulates pharmacokinetic (PK) profiles of ethinylestradiol (EE) and dienogest (DNG)
% over a 28-day cycle with possibly variable intake times, using a two-compartment
% PK model and overlays simulated profiles against observed clinical data.
% ========================================================================

function Fig2C_iii

% -------------------------------------------------------------------------
% Fixed parameters for daily dose values and intake times
% -------------------------------------------------------------------------
parameters = zeros(84,1);
parameters(1:21) =  3e4;		% Daily EE dose (ng), days 1–21
parameters(1+28:21+28) =  2e6;		% Daily DNG dose (ng), days 1–21
parameters(57:84) = round(7/24,3);	% Time of intake in fraction of day 

% -------------------------------------------------------------------------
% Time settings
% -------------------------------------------------------------------------
tend = 28;			% Total simulation duration (days)
tstep = 0.01;			% Time step (days)
xt = 0:tstep:tend; 		% Time vector (days)

% -------------------------------------------------------------------------
% Two-compartment PK parameters for EE
% -------------------------------------------------------------------------  
pkE = 1.0e+04 * [0.002107051969967		% EE absorption rate
   0.000461122242436				% EE peripheral ↔ central exchange
   7.456837524294440				% EE central volume
   0.002107082859459				% EE fast disposition
   0.000131895505995];				% EE slow disposition

kaEE = pkE(1);
k21EE = pkE(2);
VcEE = pkE(3);
alpha1EE = pkE(4);
beta1EE = pkE(5);
FEE =  0.65;					% EE bioavailability

% -------------------------------------------------------------------------
% PK parameters for Dienogest (DNG)
% -------------------------------------------------------------------------
pkP = 1.0e+04 * [0.003731412943223
   0.000930569740760
   2.207272982252732
   0.001484183803122
   0.000174898270145];                                    %INPUT

kaDNG = pkP(1);					% DNG absorption rate
k21DNG = pkP(2);				% DNG peripheral ↔ central exchange
VcDNG = pkP(3);					% DNG central volume
alpha1DNG = pkP(4);				% DNG central volume
beta1DNG = pkP(5);				% DNG central volume
FDNG =  0.90;					% DNG bioavailability


totalind   = ((tend)/tstep) + 1;		% Total number of simulation points

% Preallocate matrices to store daily concentration profiles
CE2 = zeros((tend), totalind);		% EE concentrations
CP4 = zeros((tend), totalind);		% DNG concentrations
 
% -------------------------------------------------------------------------
% Loop over each dosing day, apply dosing and compute concentration profile
% -------------------------------------------------------------------------
for la = 1:tend
     % Determine index in time vector where dose is taken
     ind = round((la -1 + parameters(2*(tend)+la))*100 + 1);                 % Offset by intake time

	% EE PK model coefficients
        NEE = (k21EE - kaEE)/( (alpha1EE - kaEE)*(beta1EE - kaEE) );
        LEE = (k21EE - alpha1EE)/( (kaEE - alpha1EE)*(beta1EE - alpha1EE) );
        MEE = (k21EE - beta1EE)/( (kaEE - beta1EE)*(alpha1EE - beta1EE) );

        CE2(la, ind:totalind)  =  ((kaEE*FEE*parameters(la))/VcEE) * ( NEE*exp(-kaEE*(xt(1:(totalind +1-ind)))) + LEE*exp(-alpha1EE*(xt(1:(totalind +1-ind)))) + MEE*exp(-beta1EE*(xt(1:(totalind +1-ind)))) );        


	% DNG PK model coefficients
        NDNG = (k21DNG - kaDNG)/( (alpha1DNG - kaDNG)*(beta1DNG - kaDNG) );
        LDNG = (k21DNG - alpha1DNG)/( (kaDNG - alpha1DNG)*(beta1DNG - alpha1DNG) );
        MDNG = (k21DNG - beta1DNG)/( (kaDNG - beta1DNG)*(alpha1DNG - beta1DNG) );

        CP4(la, ind:totalind)  =  ((kaDNG*FDNG*parameters(la + tend))/VcDNG) * ( NDNG*exp(-kaDNG*(xt(1:(totalind +1-ind)))) + LDNG*exp(-alpha1DNG*(xt(1:(totalind +1-ind)))) + MDNG*exp(-beta1DNG*(xt(1:(totalind +1-ind)))) );        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Sum concentrations from multiple days to get total profile
% -------------------------------------------------------------------------
totalCE2B = zeros(1, totalind);		
totalCP4B = zeros(1, totalind);

for su = 1:(tend)
    totalCE2B = totalCE2B + CE2(su, 1:totalind);
    totalCP4B = totalCP4B + CP4(su, 1:totalind);
end

% -------------------------------------------------------------------------
% Plot EE profile with observed data overlay
% -------------------------------------------------------------------------
figure(1)
%------------------------------------------------------------
% Time index conversion for observed EE data (x_E2)
%------------------------------------------------------------
% The vector `xt` is the time axis of the simulation, defined in days
% from t = 0 to t = 28, with steps of 0.01 day (i.e., 100 points per day).
% Thus, each entry in `xt` corresponds to a specific time index:
%     index = round(time_in_days * 100 + 1)
%
% Each time point is mapped to its corresponding index in `xt`
% using the formula above. These indices are used to extract the time
% values `xt(rnd_i)` for plotting observed EE data against the simulation.
rnd1 = round((6.29)*100+1);
rnd2 = round((6.29+0.34/24)*100+1);
rnd3 = round((6.29+0.68/24)*100+1);
rnd4 = round((6.29+1.04/24)*100+1);
rnd5 = round((6.29+1.51/24)*100+1);
rnd6 = round((6.29+1.98/24)*100+1);
rnd7 = round((6.29+2.56/24)*100+1);
rnd8 = round((6.29+2.94/24)*100+1);
rnd9 = round((6.29+3.51/24)*100+1);
rnd10 = round((6.29+3.89/24)*100+1);
rnd11 = round((6.29+4.93/24)*100+1);
rnd12 = round((6.29+5.96/24)*100+1);
rnd13 = round((6.29+7.00/24)*100+1);
rnd14 = round((6.29+8.03/24)*100+1);
rnd15 = round((6.29+8.97/24)*100+1);
rnd16 = round((6.29+9.91/24)*100+1);
rnd17 = round((6.29+11.98/24)*100+1);
rnd18 = round((6.29+23.99/24)*100+1);

% Sample timepoints for observed EE data
x_E2 = [xt(rnd1); xt(rnd2); xt(rnd3); xt(rnd4); xt(rnd5); xt(rnd6); xt(rnd7); xt(rnd8); xt(rnd9);
    xt(rnd10); xt(rnd11); xt(rnd12); xt(rnd13); xt(rnd14); xt(rnd15); xt(rnd16); xt(rnd17); xt(rnd18)];

% Observed EE concentration data
y_E2 = [18.44; 50.81; 111.03; 126.84; 122.34; 116.70; 100.90; 84.35; 82.47; 76.84; 62.55;
    52.41; 45.27; 41.15; 37.78; 36.66; 34.44; 21.45]./1000;

plot(xt, totalCE2B, 'b', 'LineWidth', 1)
hold on
plot(x_E2, y_E2, 'ob', 'MarkerSize', 5,'linewidth', 3)
legend 'Simulation' 'Data'
set(gca,'FontSize',20)
xlim([6.29 7.29]);
xticks([6.29 7.29]);
xticklabels([0 24]);
ylim([0 0.15]);
yticks([0 0.15]);
xlabel('$t$ after dose [hours]','Interpreter','latex')
ylabel('$EE$ [ng/mL]','Interpreter','latex')

figure(2)
%------------------------------------------------------------
% Each observed DNG sampling time is first
% converted to days and then mapped to its corresponding index in `xt`
% using the formula: index = round(time_in_days * 100 + 1). 
%
%These indices are used to extract the time
% values `xt(rnd_i)` for plotting the observed DNG data against the
% simulated concentration–time profile.
%------------------------------------------------------------
rnd1 = round((6.29+0.08/24)*100+1);
rnd2 = round((6.29+0.33/24)*100+1);
rnd3 = round((6.29+0.66/24)*100+1);
rnd4 = round((6.29+1.03/24)*100+1);
rnd5 = round((6.29+1.49/24)*100+1);
rnd6 = round((6.29+2.05/24)*100+1);
rnd7 = round((6.29+2.52/24)*100+1);
rnd8 = round((6.29+2.98/24)*100+1);
rnd9 = round((6.29+3.55/24)*100+1);
rnd10 = round((6.29+4.01/24)*100+1);
rnd11 = round((6.29+4.94/24)*100+1);
rnd12 = round((6.29+6.06/24)*100+1);
rnd13 = round((6.29+7.09/24)*100+1);
rnd14 = round((6.29+8.02/24)*100+1);
rnd15 = round((6.29+9.04/24)*100+1);
rnd16 = round((6.29+9.97/24)*100+1);
rnd17 = round((6.29+12.01/24)*100+1);
rnd18 = round((6.29+23.98/24)*100+1);

% Sample timepoints for observed DNG data
x_P4 = [xt(rnd1); xt(rnd2); xt(rnd3); xt(rnd4); xt(rnd5); xt(rnd6); xt(rnd7); xt(rnd8); xt(rnd9);
    xt(rnd10); xt(rnd11); xt(rnd12); xt(rnd13); xt(rnd14); xt(rnd15); xt(rnd16); xt(rnd17); xt(rnd18)];

% Observed DNG concentration data
y_P4 = [10.77; 32.14; 60.80; 63.75; 64.80; 63.07; 58.73; 55.26; 51.96; 48.84; 45.20;
    40.69; 36.53; 33.58; 29.59; 27.69; 25.97; 11.46];
plot(xt, totalCP4B, 'm', 'LineWidth', 1)
hold on
plot(x_P4, y_P4, 'om', 'MarkerSize', 5,'linewidth', 3)
set(gca,'FontSize',20)
legend 'Simulation' 'Data'
set(gca,'FontSize',20)
xlim([6.29 7.29]);
xticks([6.29 7.29]);
xticklabels([0 24]);
ylim([0 65]);
yticks([0 60]);
yticklabels([0 60]);
xlabel('$t$ after dose [hours]','Interpreter','latex')
ylabel('$DNG$ [ng/mL]','Interpreter','latex')


end

