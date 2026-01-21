% ========================================================================= 
% Fig4.m
%
% PURPOSE:
% Compare hormone profiles and ovarian dynamics between two different
% dosing times (morning vs evening) for EE + DNG contraceptive regimens 
% when only LH circadian rhythm exists.
% Specifically, this script simulates a 3-cycle (3*28-day) treatment schedule
% with constant EE and DNG doses, taken either at 11:00 or 22:00.
%
% The model integrates:
% (1) Dosing schedule for 3 cycles (with fixed intake time per case)
% (2) Delay differential menstrual cycle model (dde23)
% (3) Closed-form PK calculations for EE and DNG per dose
% (4) Circadian modulation of LH production
% (5) Endogenous hormone construction from ovarian states
% (6) Plotting of outputs for comparison: EE, DNG, RPLH, LH, P4, etc.
%
% CODE STRUCTURE:
% A) Set up dose vector and intake time for each case
% B) Solve DDE model (solve_mod)
% C) Compute hormone outputs and PK profiles
% D) Store results per intake time
% E) Plot comparisons across key variables
%
% NOTES:
% - This script is paired with the main_model.m endocrine system.
% - Uses hard-coded parameters and fixed values, all explained below.
% =========================================================================

function Fig4

% -------------------------------------------------------------------------
% [A] Initialize parameter vector and simulation arrays
% - parameters vector is 3×3×28 = 252 elements:
% * First 84: EE daily doses (ng)
% * Next 84: DNG daily doses (ng)
% * Final 84: Intake times per day (fraction of day)
% - Each row of *_res arrays will store outputs for one dosing time (11:00 or 22:00)
% -------------------------------------------------------------------------
parameters = zeros(3*3*28,1);			% dose + intake vector
E2_t0 = 1; E2_tf = 21; 				% EE given on days 1–21 per cycle
P4_t0 = E2_t0 + 3*28;				% DNG block starts at index 85
P4_tf = P4_t0 + E2_tf - E2_t0;			% Also 21-day window

% Preallocate result storage matrices
P4res    = zeros(2, 8401);			% endogenous P4
RPLHcirc = zeros(2, 8401);			% LH synthesis circadian component
EE_con   = zeros(2, 8401);			% EE concentration
DNG_con   = zeros(2, 8401);			% DNG concentration
RPLHres  = zeros(2, 8401);			% LH in reserve pool/pituitary
LHres    = zeros(2, 8401);			% LH blood concentration
RcFres   = zeros(2, 8401);			% Recruited follicular mass
Lut3res  = zeros(2, 8401);			% Luteal stage 3 mass
Lut4res  = zeros(2, 8401);			% Luteal stage 4 mass

time = [11,22];					% Test 11:00 vs 22:00 dosing

for i = 1:2
  oras = time(i);                               % Intake time (hours)
  parameters(1+6*28:9*28) = round(oras/24,2);	% Store intake as fraction of day

% -----------------------------------------------------------------------
% Fill daily EE and DNG doses (ng) for 3 cycles
% Classic 21/7 regimen: dose from day 1–21, rest for 7 days
% Each dose multiplied by 1000 to convert µg to ng
% -----------------------------------------------------------------------
dose_EE = 28.2*1000;                                      % Input EE dose in ng

   parameters(E2_t0:E2_tf) = dose_EE;   		   % EE doses for treatment cycle 1
   parameters(E2_t0+28:E2_tf+28) = dose_EE; 		   % EE doses for treatment cycle 2
   parameters(E2_t0+2*28:E2_tf+2*28) = dose_EE;   	   % EE doses for treatment cycle 3
   
dose_DNG = 320*1000;                                       % Input DNG dose in ng
   parameters(P4_t0:P4_tf) = dose_DNG;			   % DNG doses for treatment cycle 1
   parameters(P4_t0+28:P4_tf+28) = dose_DNG;		   % DNG doses for treatment cycle 2
   parameters(P4_t0+2*28:P4_tf+2*28) = dose_DNG;	   % DNG doses for treatment cycle 3
  
% -----------------------------------------------------------------------
% Circadian parameters: amplitude and phase
% par1–4, par7-8: E2, P4, FSH modulations (set to 0 here)
% par5–6: LH modulation
% -----------------------------------------------------------------------
parameterssquared = [0; 0.0915; 0; 0.3487; 0.2242; 0.5730; 0; 0.5480];

% -----------------------------------------------------------------------
% Coefficients for constructing ovarian hormone concentrations
% from follicular and luteal compartment masses:
%
%   E2(t) = e0 + e1*GrF(t) + e2*DomF(t) + e3*Lut4(t)
%   P4(t) = p1*Lut3(t) + p2*Lut4(t)
%   Inh(t) = h0 + h1*DomF(t) + h2*Lut2(t) + h3*Lut3(t);
%
% The constants e0–e3, p1–p2, and h0-h3 were obtained from prior parameter
% calibration (see S1 Table 1) and map tissue masses to circulating
% hormone concentrations.
% -----------------------------------------------------------------------
e0 = 57.60/1000;         	% basal E2
e1 = 0.0269/1000;        	% E2 contribution per unit growing follicle mass
e2 = 0.4196/1000;        	% E2 contribution per unit dominant follicle mass
e3 = 0.4923/1000;        	% E2 contribution per unit luteal stage 4 mass
p1 = 0.0032;             	% P4 contribution per unit luteal stage 3 mass
p2 = 0.1188;             	% P4 contribution per unit luteal stage 4 mass
h0 = 0.6606;  			% basal Inh
h1 = 0.0193; 			% Inh contribution per unit dominant follicle mass
h2 = 0.0159; 			% Inh contribution per unit luteal stage 2 mass
h3 = 0.0119; 			% Inh contribution per unit luteal stage 3 mass

% -----------------------------------------------------------------------
% [B] Solve endocrine model for given intake time
% solve_mod() integrates delay differential system for 3 cycles (84 days)
% -----------------------------------------------------------------------
solutions = solve_mod(parameters)';  
RPLH = solutions(:,1)';		% LH in the pituitary
LH = solutions(:,2)';		% Circulating LH concentration
FSH = solutions(:,4)';		% Circulating FSH concentration
RcF = solutions(:,5)';		% Recruited follicular mass
GrF = solutions(:,6)';		% Growing follicular mass
DomF = solutions(:,7)';		% Dominant follicular mass
Lut2 = solutions(:,11)';	% Luteal stage 2 tissue mass
Lut3 = solutions(:,12)';	% Luteal stage 3 tissue mass
Lut4 = solutions(:,13)';	% Luteal stage 4 tissue mass

% Compute baseline (non‑circadian, non‑drug) endogenous steroid profiles
E2bef  = e0  + e1*GrF  + e2*DomF  + e3*Lut4;	% Estrogen
P4bef  = p1*Lut3  + p2*Lut4;			% Progesterone

% -----------------------------------------------------------------------
% [C] Closed‑form two‑compartment PK model for EE and DNG
%
% Time grid:
%   xt spans 0–84 days (3 cycles) with resolution 0.01 day (≈14.4 min),
%   so 100 samples per day are used.
% -----------------------------------------------------------------------
tend = 3*28;				% total simulation length in days (3 cycles)
tstep = 0.01;				% time step in days
xt = 0:tstep:tend;			% time vector
totalind   = ((tend)/tstep) + 1;  	% number of discrete time points

% PK parameter vector for EE:
%   ka     : first‑order absorption rate constant
%   k21    : inter‑compartmental transfer rate
%   Vc     : central volume of distribution
%   alpha1  : fast disposition (distribution) rate
%   beta1   : slow disposition (elimination) rate
pkE = 1.0e+04 * [0.002107051969967		
   0.000461122242436
   7.456837524294440
   0.002107082859459
   0.000131895505995];

kaEE = pkE(1);		% absorption rate of EE
k21EE = pkE(2);		% peripheral ↔ central exchange rate
VcEE = pkE(3);		% central volume of distribution
alpha1EE = pkE(4);	% fast disposition rate
beta1EE = pkE(5);	% slow disposition rate
FEE =  0.65;		% EE bioavailability

% PK parameter vector for DNG (same structure as EE)
pkP = 1.0e+04 * [0.003731412943223
   0.000930569740760
   2.207272982252732
   0.001484183803122
   0.000174898270145];                                    

kaDNG = pkP(1);		% absorption rate of EE
k21DNG = pkP(2);	% peripheral ↔ central exchange rate
VcDNG = pkP(3);		% central volume of distribution
alpha1DNG = pkP(4);	% fast disposition rate
beta1DNG = pkP(5);	% slow disposition rate
FDNG =  0.90; 		% DNG bioavailability

% Preallocate matrices storing the concentration contribution of each
% daily dose at all subsequent time points
CE2 = zeros((tend), totalind);		% EE concentration contributions
CP4 = zeros((tend), totalind);		% DNG concentration contributions
 
% -----------------------------------------------------------------------
% [D] Loop over each administered daily dose and compute its contribution
% to the total plasma concentration using the closed‑form solution of the
% two‑compartment oral PK model.
%
% Intake time (stored in parameters) is converted to a discrete index in xt:
%   index = round( (day_index − 1 + intake_time_in_days) * 100 + 1 )
% because xt has 100 points per day (Δt = 0.01 day).
% -----------------------------------------------------------------------
for la = 1:tend
    % Convert dosing time of day "la" to index in time vector xt
    ind = round((la -1 + parameters(2*(tend)+la))*100 + 1);                 

	% EE PK coefficients
        NEE = (k21EE - kaEE)/( (alpha1EE - kaEE)*(beta1EE - kaEE) );
        LEE = (k21EE - alpha1EE)/( (kaEE - alpha1EE)*(beta1EE - alpha1EE) );
        MEE = (k21EE - beta1EE)/( (kaEE - beta1EE)*(alpha1EE - beta1EE) );

	% EE concentration profile generated by dose "la" from its intake time onward
        CE2(la, ind:totalind)  =  ((kaEE*FEE*parameters(la))/VcEE) * ( NEE*exp(-kaEE*(xt(1:(totalind +1-ind)))) + LEE*exp(-alpha1EE*(xt(1:(totalind +1-ind)))) + MEE*exp(-beta1EE*(xt(1:(totalind +1-ind)))) );        
	
	% EE PK coefficients
        NDNG = (k21DNG - kaDNG)/( (alpha1DNG - kaDNG)*(beta1DNG - kaDNG) );
        LDNG = (k21DNG - alpha1DNG)/( (kaDNG - alpha1DNG)*(beta1DNG - alpha1DNG) );
        MDNG = (k21DNG - beta1DNG)/( (kaDNG - beta1DNG)*(alpha1DNG - beta1DNG) );

	% DNG concentration profile generated by dose "la"
        CP4(la, ind:totalind)  =  ((kaDNG*FDNG*parameters(la + tend))/VcDNG) * ( NDNG*exp(-kaDNG*(xt(1:(totalind +1-ind)))) + LDNG*exp(-alpha1DNG*(xt(1:(totalind +1-ind)))) + MDNG*exp(-beta1DNG*(xt(1:(totalind +1-ind)))) );        
end


% -----------------------------------------------------------------------
% [E] Superposition of all prior daily dose contributions
% Total plasma concentration at each time point is obtained by summing
% across rows of CE2 and CP4.
% -----------------------------------------------------------------------
totalCE2 = zeros(1, totalind);
totalCP4 = zeros(1, totalind);

for su = 1:(tend)
    totalCE2 = totalCE2 + CE2(su, 1:totalind);		% total EE concentration time series
    totalCP4 = totalCP4 + CP4(su, 1:totalind);		% total DNG concentration time series
end

% Circadian modulation term applied to LH reserve pool production
RPLHcirc(i,:) = 1 + parameterssquared(5)*cos((2*pi)*(xt-parameterssquared(6)));

% Store results for current intake time (morning or evening)
EE_con(i,:)   = totalCE2;		% EE concentration 
DNG_con(i,:)   = totalCP4;		% DNG concentration 
RPLHres(i,:)  = RPLH;			% pituitary LH 
LHres(i,:)    = LH;			% circulating LH 
RcFres(i,:)   = RcF;			% recruited follicular mass
Lut3res(i,:)  = Lut3;			% luteal stage 3 tissue mass
Lut4res(i,:)  = Lut4;			% luteal stage 4 tissue mass

% Add circadian modulation to estrogen and progesterone
E2     = E2bef     + (parameterssquared(1))*E2bef.*cos(2*pi*(xt - (parameterssquared(2))));  % endogenous E2
P4     = P4bef    + (parameterssquared(3))*P4bef.*cos(2*pi*(xt - (parameterssquared(4))));   % endogenous P4

% Compute Inhibin based on follicle/luteal states
Inh           = h0 + h1*DomF + h2*Lut2 + h3*Lut3;

% Store final progesterone profile for current simulation iteration
P4res(i,:) = P4;

end


% Plots EE profile zoomed in on day 28-29 (morning)
figure(1)
plot(xt, EE_con(1,:), 'Color', [1, 215/255, 0],'LineWidth', 1);  
set(gca,'FontSize',20)
xlim([28 29]);
xticks([28 28+11/24 29]);
xticklabels([0 11 24]);
xtickangle(90)
%ylim([0 0.14]);
yticks([0 0.14]);
xlabel('$t$ [hours]','Interpreter','latex')
ylabel('$EE$ [ng/mL]','Interpreter','latex')
legend 'morning' 

% Plots LH production circadian component zoomed in on day 28-29 (morning)
figure(2)
plot(xt, RPLHcirc(1,:), 'Color', [1, 215/255, 0],'LineWidth', 1);  
set(gca,'FontSize',20)
xlim([28 29]);
xticks([28 28+11/24 29]);
xticklabels([0 11 24]);
xtickangle(90)
%ylim([0.6 1.3]);
yticks([0.6 1.3]);
xlabel('$t$ [hours]','Interpreter','latex')
ylabel('$RP_{LH}$ circadian rhythm [IU/(Lday)]','Interpreter','latex')
legend 'morning' 

% Plots EE profile zoomed in on day 28-29 (evening)
figure(3)
plot(xt, EE_con(2,:), 'k','LineWidth', 1);
set(gca,'FontSize',20)
xlim([28 29]);
xticks([28 28+22/24 29]);
xticklabels([0 22 24]);
xtickangle(90)
%ylim([0 0.14]);
yticks([0 0.14]);
xlabel('$t$ [hours]','Interpreter','latex')
ylabel('$EE$ [ng/mL]','Interpreter','latex')
legend 'evening'

% Plots LH production circadian component zoomed in on day 28-29 (evening)
figure(4)
plot(xt, RPLHcirc(2,:), 'k','LineWidth', 1);
set(gca,'FontSize',20)
xlim([28 29]);
xticks([28 28+22/24 29]);
xticklabels([0 22 24]);
xtickangle(90)
%ylim([0.6 1.3]);
yticks([0.6 1.3]);
xlabel('$t$ [hours]','Interpreter','latex')
ylabel('$RP_{LH}$ circadian rhythm [IU/(Lday)]','Interpreter','latex')
legend 'evening'

% Plots RPLH dynamics over days 28–56 (second cycle) for morning and evening dosing
figure(5)
plot(xt, RPLHres(1,:), 'Color', [1, 215/255, 0],'LineWidth', 1);  
hold on
plot(xt, RPLHres(2,:), 'k','LineWidth', 1);
set(gca,'FontSize',20)
xlim([28 56]);
xticks([28 56]);
%ylim([0 3000]);
yticks([0 3000]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$RP_{LH}$ [IU]','Interpreter','latex')
legend 'morning' 'evening'

% Plots LH dynamics over days 28–56 (second cycle) for morning and evening dosing
figure(6)
plot(xt, LHres(1,:), 'Color', [1, 215/255, 0],'LineWidth', 1);  
hold on
plot(xt, LHres(2,:), 'k','LineWidth', 1);
set(gca,'FontSize',20)
xlim([28 56]);
xticks([28 56]);
%ylim([0 80]);
yticks([0 80]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$LH$ [IU/L]','Interpreter','latex')
legend 'morning' 'evening'

% Plots recruited follicular dynamics over days 28–56 (second cycle) for morning and evening dosing
figure(7)
plot(xt, RcFres(1,:), 'Color', [1, 215/255, 0],'LineWidth', 1);  
hold on
plot(xt, RcFres(2,:), 'k','LineWidth', 1);
set(gca,'FontSize',20)
xlim([28 56]);
xticks([28 56]);
%ylim([0 120]);
yticks([0 120]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$RcF$ [$\mu$g]','Interpreter','latex')
legend 'morning' 'evening'

% Plots luteal stage 3 tissue dynamics over days 28–56 (second cycle) for morning and evening dosing
figure(8)
plot(xt, Lut3res(1,:), 'Color', [1, 215/255, 0],'LineWidth', 1);  
hold on
plot(xt, Lut3res(2,:), 'k','LineWidth', 1);
set(gca,'FontSize',20)
xlim([28 56]);
xticks([28 56]);
%ylim([0 60]);
yticks([0 60]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$Lut_3$ [$\mu$g]','Interpreter','latex')
legend 'morning' 'evening'

% Plots luteal stage 4 tissue dynamics over days 28–56 (second cycle) for morning and evening dosing
figure(9)
plot(xt, Lut4res(1,:), 'Color', [1, 215/255, 0],'LineWidth', 1);  
hold on
plot(xt, Lut4res(2,:), 'k','LineWidth', 1);
set(gca,'FontSize',20)
xlim([28 56]);
xticks([28 56]);
%ylim([0 30]);
yticks([0 30]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$Lut_4$ [$\mu$g]','Interpreter','latex')
legend 'morning' 'evening'

% Plots endogenous P4 dynamics over days 28–56 (second cycle) for morning and evening dosing
figure(10)
plot(xt, P4res(1,:), 'Color', [1, 215/255, 0],'LineWidth', 1);  
hold on
plot(xt, P4res(2,:), 'k','LineWidth', 1);
hold on
yline(3, '--k', 'LineWidth', 0.5)
set(gca,'FontSize',20)
xlim([28 56]);
xticks([28 56]);
%ylim([0 4]);
yticks([0 3 4]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$P_4$ [ng/mL]','Interpreter','latex')
legend 'morning' 'evening'

% Plots DNG profile zoomed in on day 28-29 (morning)
figure(11)
plot(xt, DNG_con(1,:), 'Color', [1, 215/255, 0],'LineWidth', 1);  
set(gca,'FontSize',20)
xlim([28 29]);
xticks([28 28+11/24 29]);
xticklabels([0 11 24]);
xtickangle(90)
% ylim([0 0.14]);
% yticks([0 0.14]);
xlabel('$t$ [hours]','Interpreter','latex')
ylabel('$DNG$ [ng/mL]','Interpreter','latex')
legend 'morning' 

% Plots DNG profile zoomed in on day 28-29 (evening)
figure(12)
plot(xt, DNG_con(2,:), 'k','LineWidth', 1);
set(gca,'FontSize',20)
xlim([28 29]);
xticks([28 28+22/24 29]);
xticklabels([0 22 24]);
xtickangle(90)
% ylim([0 0.14]);
% yticks([0 0.14]);
xlabel('$t$ [hours]','Interpreter','latex')
ylabel('$DNG$ [ng/mL]','Interpreter','latex')
legend 'evening'

end

% -------------------------------------------------------------------------
% solve_mod(parameters)
%
% Runs the endocrine model as a delay differential equation (DDE) system:
%   - MATLAB solver: dde23
%   - Single discrete delay: dInh (Inhibin feedback delay)
%   - Returns solution evaluated on x = 0:0.01:3*28 (days)
%
% parameters input:
%   The same 3-block vector built in main_model:
%     (1) EE daily doses
%     (2) DNG daily doses
%     (3) daily intake times
% -------------------------------------------------------------------------
function [solution] = solve_mod(parameters)

x    = 0:0.01:3*28;
tdata = [0 3*28];

%Delay (days): inhibitory feedback uses delayed ovarian states
dInh =  1.5; 

% Initial state vector for DDE system (13 states)
% The state ordering matches model() below.
Init = [167.57; 11.81; 14.48; 11.41; 2.10; 4.12; 0.46; 1.06; 1.67; 4.16; 13.03; 16.48; 10.29];

u = parameters;

solution = dde23(@model, dInh, Init, tdata, [], u);
solution = deval(solution, x);

end


% -------------------------------------------------------------------------
% model(t, state, delay, u)
%
% Right-hand side of the DDE system:
%   - Pituitary hormone production/release (LH, FSH) with steroid feedback
%   - Ovarian follicle development and luteal progression
%   - Endogenous E2/P4 output from ovarian states
%   - Exogenous hormones (EE/DNG concentrations) from daily dosing schedule (u)
%   - Circadian forcing implemented via cosine terms in:
%       * LH production (par5, par6)
%       * FSH production (par7, par8)
%       * Endogenous E2/P4 outputs (par1..par4)
%
% Inputs:
%   t     : current time (days)
%   state : 13x1 current state vector
%   delay : delayed state(s) used for Inhibin feedback (dde23 provides this)
%   u     : 3*(3*28) vector encoding doses and intake times
% -------------------------------------------------------------------------
function dstate = model(t, state, delay, u)

% -----------------------------
% Fixed parameters (model fitted constants)
% -----------------------------
kLH =  0.9661; 
V0LH =  550.03; 
V1LH = 3329.19; 
KmLH =   136.05/1000;                                        
KiLHP =  6.78;            
cLHE =  0.0060*1000;  
cLHP =  1.98;    
VFSH = 294.90; 
kFSH = 14.59; 
cFSHE = 0.0151*1000^2;   
KiFSHInh = 16.83;              
cFSHP = 52.31;       
b = 0.0453;
c1 = 0.1036; 
c2 = 0.0577; 
c3 = 0.0170; 
c4 = 1.14; 
d1 = 0.7537; 
d2 = 0.6866; 
k1 = 0.6699; 
k2 = 0.6388; 
k3 = 0.9191; 
k4 = 1.88; 
alpha = 0.9505; 
gamma = 0.1615; 
w = 9.21;            
q = 5.11;   

% Endogenous hormone output coefficients (same as main_model)
e0 = 57.60/1000;     
e1 = 0.0269/1000;    
e2 = 0.4196/1000;    
e3 = 0.4923/1000;    
p1 = 0.0032;         
p2 = 0.1188;         
h0 = 0.6606;  
h1 = 0.0193; 
h2 = 0.0159; 
h3 = 0.0119; 
         
% Circadian forcing parameters (8 parameters)
additional_par =  [0; 0.0915; 0; 0.3487; 0.2242; 0.5730; 0; 0.5480];

% Map them to named scalars for readability in equations
par1 =         additional_par(1);	% E2 amplitude
par2 =         additional_par(2);	% E2 phase
par3 =         additional_par(3);	% P4 amplitude
par4 =         additional_par(4);	% P4 phase
par5 =         additional_par(5);	% LH production amplitude
par6 =         additional_par(6);	% LH production phase
par7 =         additional_par(7);	% FSH production amplitude
par8 =         additional_par(8);	% FSH production phase

% -----------------------------
% State vector unpacking (13 states)
% -----------------------------
RPLH            =  state(1); 	% reserve pool / releasable LH (pituitary)
LH              =  state(2); 	% circulating LH
RPFSH           =  state(3);	% reserve pool / releasable FSH (pituitary)
FSH             =  state(4); 	% circulating FSH
RcF             =  state(5);   	% recruited follicles
GrF             =  state(6);	% growing follicles
DomF            =  state(7);	% dominant follicle
Sc1             =  state(8);	% ovulation stage
Sc2             =  state(9);	% luteinization stage
Lut1            =  state(10);	% luteal stage 1
Lut2            =  state(11);	% luteal stage 2
Lut3            =  state(12);	% luteal stage 3
Lut4            =  state(13);	% luteal stage 4

% -----------------------------
% Exogenous hormone (drug concentration) computation at current time t
% -----------------------------
tend = 3*28;

% ConE2 and ConP4 store per-dose concentration at current time t
ConE2 = zeros(1,tend);
ConP4 = zeros(1,tend);

% PK parameter blocks (same as main_model)
pkE = 1.0e+04 * [0.002107051969967
   0.000461122242436
   7.456837524294440
   0.002107082859459
   0.000131895505995];

kaEE = pkE(1);
k21EE = pkE(2);
VcEE = pkE(3);
alpha1EE = pkE(4);
beta1EE = pkE(5);
FEE =  0.65;

pkP = 1.0e+04 * [0.003731412943223
   0.000930569740760
   2.207272982252732
   0.001484183803122
   0.000174898270145];                                     

kaDNG = pkP(1);
k21DNG = pkP(2);
VcDNG = pkP(3);
alpha1DNG = pkP(4);
beta1DNG = pkP(5);
FDNG =  0.90; 

% Loop over each day-dose and compute its contribution at time t
% If current time t is before that day's intake time, contribution is zero.
for la = 1:tend 
    if (0 <= t) && (t < (u(2*(tend)+la)+ la-1))
        ConE2(la) = 0;
        ConP4(la) = 0;
    elseif (u(2*tend+la)+ la-1) <= t
        
        NEE = (k21EE - kaEE)/( (alpha1EE - kaEE)*(beta1EE - kaEE) );
        LEE = (k21EE - alpha1EE)/( (kaEE - alpha1EE)*(beta1EE - alpha1EE) );
        MEE = (k21EE - beta1EE)/( (kaEE - beta1EE)*(alpha1EE - beta1EE) );

        ConE2(la)  =  ((kaEE*FEE*u(la))/VcEE) * ( NEE*exp(-kaEE*(t-(u(2*tend+1)+ la-1))) + LEE*exp(-alpha1EE*((t-(u(2*tend+1)+ la-1)))) + MEE*exp(-beta1EE*((t-(u(2*tend+1)+ la-1)))) );        

        NDNG = (k21DNG - kaDNG)/( (alpha1DNG - kaDNG)*(beta1DNG - kaDNG) );
        LDNG = (k21DNG - alpha1DNG)/( (kaDNG - alpha1DNG)*(beta1DNG - alpha1DNG) );
        MDNG = (k21DNG - beta1DNG)/( (kaDNG - beta1DNG)*(alpha1DNG - beta1DNG) );

        ConP4(la)  =  ((kaDNG*FDNG*u(la + tend))/VcDNG) * ( NDNG*exp(-kaDNG*((t-(u(2*tend+1)+ la-1)))) + LDNG*exp(-alpha1DNG*((t-(u(2*tend+1)+ la-1)))) + MDNG*exp(-beta1DNG*((t-(u(2*tend+1)+ la-1)))) );        
    else
    end
end

% Sum contributions across all daily doses to obtain total exogenous level at time t
    totalConE2 = sum(ConE2);
    totalConP4 = sum(ConP4);

% -------------------------------------------------------------------------
% Endogenous + exogenous hormones at time t
% - Endogenous E2/P4 from ovarian state variables, with cosine circadian forcing
% - Add scaled rug concentration terms
% -------------------------------------------------------------------------
E2              =    (e0 + e1*GrF + e2*DomF + e3*Lut4) + par1*(e0 + e1*GrF + e2*DomF + e3*Lut4)*cos(2*pi*(t-par2)) + 1.7*totalConE2;
P4              =    (p1*Lut3 + p2*Lut4) + par3*(p1*Lut3 + p2*Lut4)*cos(2*pi*(t-par4)) + 0.01*totalConP4;

% -------------------------------------------------------------------------
% Delayed inhibin feedback:
% Use delayed ovarian states (DomF, Lut2, Lut3) to compute InhdelayInh
% This delay implements physiological lag in negative feedback on FSH.
% -------------------------------------------------------------------------
statelag1       = delay(:,1);

DomFdelayInh     = statelag1(7);
Lut2delayInh    = statelag1(11);
Lut3delayInh    = statelag1(12);


InhdelayInh     = h0 + h1*DomFdelayInh + h2*Lut2delayInh + h3*Lut3delayInh;

% -------------------------------------------------------------------------
% Pituitary dynamics:
% - LH production: stimulated by E2 via Hill function, inhibited by P4
% - Circadian forcing: (1 + par5*cos(2*pi*(t-par6))) multiplies production term
% - Release term depends on RPLH pool, modulated by steroids (E2, P4)
% -------------------------------------------------------------------------
dRPLHdt     =   ((V0LH + (V1LH*E2^8/(KmLH^8 + E2^8)))/(1 + (P4/KiLHP))).*(1 + par5*cos((2*pi)*(t-par6))) - ((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2));

dLHdt       =   (1/2.5)*((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2)) - 14*LH;

% -------------------------------------------------------------------------
% FSH dynamics:
% - Production inhibited by delayed inhibin and P4, plus circadian forcing
% - Release depends on RPFSH pool and steroid modulation
% -------------------------------------------------------------------------
dRPFSHdt    =   (VFSH/(1 + (InhdelayInh/KiFSHInh) + (P4/w))).*(1 + par7*cos((2*pi)*(t-par8))) - ((kFSH*(1 + cFSHP*P4).*RPFSH)./(1+cFSHE*E2^2));

dFSHdt      =   (1/2.5)*((kFSH*(1 + cFSHP*P4)*RPFSH)/(1 + cFSHE*E2^2)) -  8.21*FSH;  

% -------------------------------------------------------------------------
% Ovarian dynamics:
% Recruitment -> growth -> dominance -> ovulation/luteal progression
% Parameterization follows the paper's endocrine module.
% -------------------------------------------------------------------------
dRcFdt      =   (b + c1*RcF)*(FSH/(1 + (P4/q))) - c2*(LH^alpha)*RcF;
dGrFdt      =    c2*(LH^alpha)*RcF - c3*LH*GrF; 
dDomFdt     =    c3*LH*GrF - c4*(LH^gamma)*DomF;
dSc1dt          =    c4.*(LH.^gamma).*DomF - d1.*Sc1;
dSc2dt          =    d1.*Sc1 - d2.*Sc2;
dLut1dt         =    d2.*Sc2 - k1.*Lut1;
dLut2dt         =    k1.*Lut1 - k2.*Lut2;
dLut3dt         =    k2.*Lut2 - k3.*Lut3;
dLut4dt         =    k3.*Lut3 - k4.*Lut4;

% Pack derivatives into 13x1 vector in the same order as state variables
dstate          = zeros(13,1);

dstate(1)       = dRPLHdt;
dstate(2)       = dLHdt;
dstate(3)       = dRPFSHdt;
dstate(4)       = dFSHdt;
dstate(5)       = dRcFdt;
dstate(6)       = dGrFdt;
dstate(7)       = dDomFdt;
dstate(8)       = dSc1dt;
dstate(9)       = dSc2dt;
dstate(10)      = dLut1dt;
dstate(11)      = dLut2dt;
dstate(12)      = dLut3dt;
dstate(13)      = dLut4dt;

end
