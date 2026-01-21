%% =========================================================================
% main_model.m
%
% PURPOSE:
%   This script implements the core mathematical model used in:
%   "Optimizing oral contraceptive timing: Daytime intake reduces
%    doses and enhances efficacy" (Gavina et al.)
%
%   The model integrates:
%     (1) Menstrual-cycle–dependent ovarian hormone dynamics
%     (2) Circadian modulation (implemented here as daily cosine forcing)
%     (3) Progesterone-based ovulation suppression criteria
%
% CODE STRUCTURE (high-level):
%   A) Build daily dosing schedule for 10 cycles (EE and DNG) + intake time
%   B) Solve the endocrine dynamics via a DDE model (dde23) -> LH, FSH, follicular masses
%   C) Compute endogenous E2/P4 from model states (baseline outputs)
%   D) Compute drug concentration (exogenous hormone) profiles for EE and DNG (closed-form)
%   E) Add circadian forcing to endogenous hormones + add scaled exogenous effects
%   F) Plot key outputs (FSH, LH, E2, P4 with threshold, Inh)
%
% NOTE:
%   - This file contains three functions in one .m file:
%       1) main_model (top-level runner)
%       2) solve_mod  (runs dde23 and returns solution trajectories)
%       3) model      (DDE right-hand side; endocrine + drug effects)
%   - For GitHub readers, the comments below explain the meaning of each block.
%
% =========================================================================

function main_model

% -------------------------------------------------------------------------
% [A] Build "parameters" vector encoding:
%   - Daily EE dose for each day of simulation (first 10*28 entries)
%   - Daily DNG dose for each day of simulation (next 10*28 entries)
%   - Daily intake time as a fraction of day (last 10*28 entries)
%
% The parameter vector length is 3*(10 cycles)*(28 days/cycle).
% Units:
%   - dose_* are entered in ng
%   - intake time is stored in "days" (i.e., oras/24)
% -------------------------------------------------------------------------
parameters = zeros(3*10*28,1);

% Input dosing time (clock hour in 0–24)
% Example: oras = 22 means 10 PM intake
oras = 22;                                      			               % input dosing time

% Store intake time for each day (as fraction of day, e.g., 22/24)
% parameters(1+2*10*28:3*10*28) corresponds to the 3rd block of the vector
parameters(1+2*10*28:3*10*28) = round(oras/24,2);

% Daily dose inputs (ng)
dose_EE = 28.2*1000;                   	                         % input EE dose in ng
dose_DNG = 320*1000;                  		                          % input DNG dose in ng
 
% Fill in daily doses across 10 cycles
% - EE doses go into parameters(1 : 10*28)
% - DNG doses go into parameters(10*28+1 : 2*10*28)
% Dosing is applied for days 1–21 in each 28-day cycle (classic 21/7 regimen)
for araw = 1:10  							                                       % number of cycles
        E2_t0 = 1;
        E2_tf = 21;
        % EE dosing on days 1–21 of each cycle
        parameters(E2_t0 + (araw-1)*28:E2_tf+ (araw-1)*28) =  dose_EE;                 
        
        % DNG block begins at offset 28*10 in the parameter vector
        P4_t0 = E2_t0 + 28*10;           
        P4_tf = P4_t0 + E2_tf - E2_t0;
        
        % DNG dosing on days 1–21 of each cycle
        parameters(P4_t0+ (araw-1)*28:P4_tf+ (araw-1)*28)  =  dose_DNG;               
   end   

% -------------------------------------------------------------------------
% Circadian-related parameter vector (8 parameters)
% Used later as:
%   par1..par4 : amplitude/phase for endogenous E2 and P4 cosine modulation
%   par5..par8 : amplitude/phase for pituitary production terms (LH/FSH) modulation
% The exact mapping is shown later inside model()
% -------------------------------------------------------------------------
parameterssquared = [0.0473; 0.0915; 0.0807; 0.3487; 0.2242; 0.5730; 0.1130; 0.5480];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Time grid (days)
% The endocrine model solution and plots use days as the independent variable.
% 0:0.01:10*28 corresponds to 10 cycles with step = 0.01 day (~14.4 minutes).
% -------------------------------------------------------------------------
x    = 0:0.01:10*28;


% -------------------------------------------------------------------------
% Coefficients for constructing endogenous E2 and P4 from ovarian state variables:
% E2bef = e0 + e1*GrF + e2*DomF + e3*Lut4
% P4bef = p1*Lut3 + p2*Lut4
%
% These represent baseline (non-circadian-forced, non-drug) hormone outputs
% from follicle/luteal compartments.
% -------------------------------------------------------------------------
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



% -------------------------------------------------------------------------
% [B] Solve the endocrine dynamics (DDE system) given the dosing/time vector
% solve_mod() returns a matrix of state trajectories over x (days)
% The solution ordering corresponds to the state vector inside model()
% -------------------------------------------------------------------------
solutions = solve_mod(parameters)';                                        

% Extract key states (index meanings documented again inside model()):
LH = solutions(:,2)';     % LH (IU/L)
FSH = solutions(:,4)';    % FSH (IU/L)
GrF = solutions(:,6)';    % Growing follicular mass
DomF = solutions(:,7)';   % Dominant follicular mass
Lut2 = solutions(:,11)';  % Luteal stage 2 tissue mass
Lut3 = solutions(:,12)';  % Luteal stage 3 tissue mass
Lut4 = solutions(:,13)';  % Luteal stage 4 tissue mass

% -------------------------------------------------------------------------
% [C] Baseline endogenous hormones from ovarian compartments
% These are "before" adding circadian forcing and exogenous drug contributions
% -------------------------------------------------------------------------
E2bef  = e0  + e1*GrF  + e2*DomF  + e3*Lut4;
P4bef  = p1*Lut3  + p2*Lut4;

% -------------------------------------------------------------------------
% [D] Drug concentration (exogenous hormone) profiles (EE and DNG)
% This block computes concentrations using a closed-form biexponential-like
% expression per daily dose, then sums them across all prior doses.
%
% IMPORTANT:
% - Time units in this PK block are still in "days" because xt spans 0..tend
%   with step 0.01 (consistent with endocrine time grid)
% - The exp() terms use rate constants ka/alpha/beta etc. already scaled
%   appropriately for the chosen time unit in the paper calibration.
% -------------------------------------------------------------------------
tend = 10*28;
tstep = 0.01;
xt = 0:tstep:tend; 
  
% PK parameter vector for EE (scaled by 1e4)
% Interpretation (typical 2-compartment oral model parameters):
%   ka   : absorption rate
%   k21  : inter-compartmental rate
%   Vc   : central volume of distribution
%   alpha/beta : disposition rates
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
FEE =  0.65;    % EE bioavailability

% PK parameter vector for DNG (scaled by 1e4)
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
FDNG =  0.90;  % DNG bioavailability


% totalind = number of time grid points in xt
totalind   = ((tend)/tstep) + 1;

% CE2(la, :) and CP4(la, :) store the concentration contribution from the la-th dose
% at all future time points (from its intake time onward).
CE2 = zeros((tend), totalind);
CP4 = zeros((tend), totalind);
 

for la = 1:tend
    % Convert daily intake time to an index in xt
    % Example: if intake time is 22/24 day, then ind shifts dose start
    ind = round((la -1 + parameters(2*(tend)+la))*100 + 1);                 % index corresponding to intake time
        % Coefficients for closed-form solution (EE)
        NEE = (k21EE - kaEE)/( (alpha1EE - kaEE)*(beta1EE - kaEE) );
        LEE = (k21EE - alpha1EE)/( (kaEE - alpha1EE)*(beta1EE - alpha1EE) );
        MEE = (k21EE - beta1EE)/( (kaEE - beta1EE)*(alpha1EE - beta1EE) );

        % Concentration contribution of EE dose la from time "ind" onward
        CE2(la, ind:totalind)  =  ((kaEE*FEE*parameters(la))/VcEE) * ( NEE*exp(-kaEE*(xt(1:(totalind +1-ind)))) + LEE*exp(-alpha1EE*(xt(1:(totalind +1-ind)))) + MEE*exp(-beta1EE*(xt(1:(totalind +1-ind)))) );        

        % Coefficients for closed-form solution (DNG)
        NDNG = (k21DNG - kaDNG)/( (alpha1DNG - kaDNG)*(beta1DNG - kaDNG) );
        LDNG = (k21DNG - alpha1DNG)/( (kaDNG - alpha1DNG)*(beta1DNG - alpha1DNG) );
        MDNG = (k21DNG - beta1DNG)/( (kaDNG - beta1DNG)*(alpha1DNG - beta1DNG) );

        % Concentration contribution of DNG dose la from time "ind" onward
        CP4(la, ind:totalind)  =  ((kaDNG*FDNG*parameters(la + tend))/VcDNG) * ( NDNG*exp(-kaDNG*(xt(1:(totalind +1-ind)))) + LDNG*exp(-alpha1DNG*(xt(1:(totalind +1-ind)))) + MDNG*exp(-beta1DNG*(xt(1:(totalind +1-ind)))) );        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Sum over all dose contributions to get total exogenous concentration profiles
% totalCE2, totalCP4 are 1 x totalind arrays (time series over xt)
% -------------------------------------------------------------------------
totalCE2 = zeros(1, totalind);
totalCP4 = zeros(1, totalind);

for su = 1:(tend)
    totalCE2 = totalCE2 + CE2(su, 1:totalind);
    totalCP4 = totalCP4 + CP4(su, 1:totalind);
end

% -------------------------------------------------------------------------
% [E] Final hormone signals:
%   - Endogenous component with circadian cosine forcing
%   - Plus scaled exogenous drug concentrations
%
% E2, P4 here are *endogenous only* (with circadian forcing)
% E2tot, P4tot are endogenous + exogenous (drug) combined signals
%
% NOTE:
%   The cosine term uses period 1 day because x is in days and the term is cos(2*pi*(x - phase))
% -------------------------------------------------------------------------

E2     = E2bef     + (parameterssquared(1))*E2bef.*cos(2*pi*(x - (parameterssquared(2))));  % endogenous E2
P4     = P4bef    + (parameterssquared(3))*P4bef.*cos(2*pi*(x - (parameterssquared(4))));   % endogenous P4

% Drug contributions are scaled (1.7 for EE, 0.01 for DNG) per the calibrated mapping
% in the paper's model/data fitting.

E2tot  = E2 + 1.7*totalCE2;
P4tot  = P4 + 0.01*totalCP4;

% Inhibin-like signal constructed from ovarian compartments (endogenous)
% This is used internally in the FSH production term (with delay in model())
Inh           = h0 + h1*DomF + h2*Lut2 + h3*Lut3;

% -------------------------------------------------------------------------
% [F] Figures: visualization of key outputs (first 3 cycles shown: 0–84 days)
% - FSH, LH: pituitary gonadotropins
% - E2, P4 : ovarian steroids (endogenous-only curves plotted here)
% - P4 plot includes yline(3): ovulation threshold reference (ng/mL)
% - Inh: inhibin
% -------------------------------------------------------------------------

figure(1)
plot(x, FSH, 'k','LineWidth', 1);  
set(gca,'FontSize',20)
xlim([0 84]);
xticks([0 28 56 84]);
xlabel('$t$ [days]', 'Interpreter','latex')
ylabel('$FSH$ [IU/L]', 'Interpreter','latex')

figure(2)
plot(x, LH, 'k','LineWidth', 1);  
set(gca,'FontSize',20)
xlim([0 84]);
xticks([0 28 56 84]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$LH$ [IU/L]','Interpreter','latex')
 
figure(3)
plot(x, E2, 'k','LineWidth', 1);  
set(gca,'FontSize',20)
xlim([0 84]);
xticks([0 28 56 84]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$E_2$ [ng/mL]','Interpreter','latex')
 
figure(4)
plot(x, P4, 'k','LineWidth', 1);  set(gca,'FontSize',20)
hold on
yline(3)
set(gca,'FontSize',20)
xlim([0 84]);
xticks([0 28 56 84]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$P_4$ [ng/mL]','Interpreter','latex')

figure(5)
plot(x, Inh, 'k','LineWidth', 1);  set(gca,'FontSize',20)
xlim([0 84]);
xticks([0 28 56 84]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$Inh$ [IU/mL]','Interpreter','latex')
 

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [solution] = solve_mod(parameters)
% -------------------------------------------------------------------------
% solve_mod(parameters)
%
% Runs the endocrine model as a delay differential equation (DDE) system:
%   - MATLAB solver: dde23
%   - Single discrete delay: dInh (Inhibin feedback delay)
%   - Returns solution evaluated on x = 0:0.01:10*28 (days)
%
% parameters input:
%   The same 3-block vector built in main_model:
%     (1) EE daily doses
%     (2) DNG daily doses
%     (3) daily intake times
% -------------------------------------------------------------------------
x    = 0:0.01:10*28;
tdata = [0 10*28];

% Delay (days): inhibitory feedback uses delayed ovarian states
dInh =  1.5; 

% Initial state vector for DDE system (13 states)
% The state ordering matches model() below.
Init = [167.57; 11.81; 14.48; 11.41; 2.10; 4.12; 0.46; 1.06; 1.67; 4.16; 13.03; 16.48; 10.29];

u = parameters;

solution = dde23(@model, dInh, Init, tdata, [], u);
solution = deval(solution, x);

end


function dstate = model(t, state, delay, u)
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
%   u     : 3*(10*28) vector encoding doses and intake times
% -------------------------------------------------------------------------

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

% Additional constants for follicle dynamics / feedback scaling
w = 9.21;            
q = 5.11;     

% Circadian forcing parameters (8 parameters)
additional_par =  [0.0473; 0.0915; 0.0807; 0.3487; 0.2242; 0.5730; 0.1130; 0.5480];


% Map them to named scalars for readability in equations
par1 =         additional_par(1);  % E2 amplitude
par2 =         additional_par(2);  % E2 phase
par3 =         additional_par(3);  % P4 amplitude
par4 =         additional_par(4);  % P4 phase
par5 =         additional_par(5);  % LH production amplitude
par6 =         additional_par(6);  % LH production phase
par7 =         additional_par(7);  % FSH production amplitude
par8 =         additional_par(8);  % FSH production phase


% -----------------------------
% State vector unpacking (13 states)
% -----------------------------
RPLH            =  state(1);   % reserve pool / releasable LH (pituitary)
LH              =  state(2);   % circulating LH
RPFSH           =  state(3);   % reserve pool / releasable FSH (pituitary)
FSH             =  state(4);   % circulating FSH
RcF             =  state(5);   % recruited follicles
GrF             =  state(6);   % growing follicles
DomF            =  state(7);   % dominant follicle
Sc1             =  state(8);   % ovulation stage
Sc2             =  state(9);   % luteinization stage
Lut1            =  state(10);  % luteal stage 1
Lut2            =  state(11);  % luteal stage 2
Lut3            =  state(12);  % luteal stage 3
Lut4            =  state(13);  % luteal stage 4

% -----------------------------
% Exogenous hormone (drug concentration) computation at current time t
% -----------------------------
tstart = 0;
tend = 10*28;
tend1 = tend + 28;

% ConE2 and ConP4 store per-dose concentration at current time t
ConE2 = zeros(1,tend1);
ConP4 = zeros(1,tend1);

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
        % EE closed-form coefficients
        NEE = (k21EE - kaEE)/( (alpha1EE - kaEE)*(beta1EE - kaEE) );
        LEE = (k21EE - alpha1EE)/( (kaEE - alpha1EE)*(beta1EE - alpha1EE) );
        MEE = (k21EE - beta1EE)/( (kaEE - beta1EE)*(alpha1EE - beta1EE) );
        
        % EE concentration contribution from dose la at time t
        ConE2(la)  =  ((kaEE*FEE*u(la))/VcEE) * ( NEE*exp(-kaEE*(t-(u(2*tend+1)+ la-1))) + LEE*exp(-alpha1EE*((t-(u(2*tend+1)+ la-1)))) + MEE*exp(-beta1EE*((t-(u(2*tend+1)+ la-1)))) );        
        
        % DNG closed-form coefficients
        NDNG = (k21DNG - kaDNG)/( (alpha1DNG - kaDNG)*(beta1DNG - kaDNG) );
        LDNG = (k21DNG - alpha1DNG)/( (kaDNG - alpha1DNG)*(beta1DNG - alpha1DNG) );
        MDNG = (k21DNG - beta1DNG)/( (kaDNG - beta1DNG)*(alpha1DNG - beta1DNG) );
        
        % DNG concentration contribution from dose la at time t
        ConP4(la)  =  ((kaDNG*FDNG*u(la + tend))/VcDNG) * ( NDNG*exp(-kaDNG*((t-(u(2*tend+1)+ la-1)))) + LDNG*exp(-alpha1DNG*((t-(u(2*tend+1)+ la-1)))) + MDNG*exp(-beta1DNG*((t-(u(2*tend+1)+ la-1)))) );        
    else
    end
end

% Sum contributions across all daily doses to obtain total exogenous level at time t
    totalConE2 = sum(ConE2);
    totalConP4 = sum(ConP4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
