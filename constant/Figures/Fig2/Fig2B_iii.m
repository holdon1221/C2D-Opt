% ========================================================================= 
% Fig2B_iii.m
%
% PURPOSE:
%   This script generates Fig. 2B(iii).
%
%   It compares:
%     (1) Target hormone profiles (FSH, LH, E2, P4) derived from 
%         experimental circadian data, and
%     (2) Simulated hormone trajectories from the circadian-coupled
%         HPO-axis DDE model using time-varying parameters.
%
% CODE STRUCTURE:
%   A) Define circadian amplitude/phase parameters for E2, P4, LH, FSH
%   B) Solve the DDE model with no circadian forcing (baseline)
%   C) Construct "target" profiles from cosine fits
%   D) Solve the model with full circadian forcing
%   E) Overlay target vs simulation results in plots
%
% NOTE:
%   - This file contains:
%       (1) Fig2B_iii (main function)
%       (2) solve_mod (DDE solver wrapper)
%       (3) model     (DDE right-hand side with circadian inputs)
%
% =========================================================================


function Fig2B_iii

clear all;
format long

% -------------------------------------------------------------------------
% [A] Circadian parameter vector (from optimization)
%   Each pair is [amplitude, acrophase (in days)] for:
%     [E2, P4, LH, FSH]
% -------------------------------------------------------------------------
parameters = [...
0.0473; 		% E2 circadian amplitude
0.0915; 		% E2 acrophase (fraction of day)
0.0807; 		% P4 amplitude
0.3487;			% P4 acrophase
0.2242; 		% LH amplitude (pituitary production)
0.5730; 		% LH acrophase
0.1130; 		% FSH amplitude (pituitary production)
0.5480];		% FSH acrophase

% Zero parameters used to simulate the baseline model without circadian forcing
parametersc = [0; 0; 0; 0; 0; 0; 0; 0]; 

% Time grid for simulation: 0–28 days in 0.01-day steps
x    = [0:0.01:28];

% -------------------------------------------------------------------------
% [B] Solve baseline (non-circadian) model
% -------------------------------------------------------------------------
solutionorig = solve_mod(parametersc)';		% returns [time x states]
LHorig  = solutionorig(:,2)';			% circulating LH
FSHorig = solutionorig(:,4)';			% circulating FSH
GrForig = solutionorig(:,6)';			% growing follicular mass
DomForig = solutionorig(:,7)';			% dominant follicular mass
Lut2orig = solutionorig(:,11)';			% luteal stage 2 tissue mass
Lut3orig = solutionorig(:,12)';			% luteal stage 3 tissue mass
Lut4orig = solutionorig(:,13)';			% luteal stage 4 tissue mass


% -------------------------------------------------------------------------
% [C] Construct "target" circadian hormone profiles (cosinor-based)
%   using baseline trajectories + cosine modulation from fitted amplitudes
%   and phases (from Fig. 2B(i)).
% -------------------------------------------------------------------------

% Constants for endogenous hormone reconstruction
e0 = 57.60;  		% baseline E2
e1 = 0.0269;		% E2 per unit GrF
e2 = 0.4196; 		% E2 per unit DomF
e3 = 0.4923;		% E2 per unit Lut4
p1 = 0.0032;		% P4 per unit Lut3
p2 = 0.1188;		% P4 per unit Lut4

% Reconstruct E2 and P4 without circadian forcing
E2orig  = e0  + e1*GrForig  + e2*DomForig  + e3*Lut4orig;
P4orig  = p1*Lut3orig  + p2*Lut4orig;

met(1) = 0.08026;               % E2 cosine fit amplitude
met(2) = 24.58/24;              % E2 cosine fit acrophase
met(3) = 0.1017;                % P4 cosine fit amplitude
met(4) = 8.46/24;               % P4 cosine fit acrophase
met(5) = 0.0531;                % LH cosine fit amplitude
met(6) = 18.27/24;              % LH cosine fit acrophase
met(7) = 0.0659;                % FSH cosine fit amplitude
met(8) = 16.46/24;              % FSH cosine fit acrophase

% Apply cosine modulation to create "target" circadian waveforms
E2circ    = E2orig    + met(1)*E2orig.*cos(2*pi*(x - met(2)));
P4circ    = P4orig    + met(3)*P4orig.*cos(2*pi*(x - met(4)));
LHcirc    = LHorig    + met(5)*LHorig.*cos(2*pi*(x - met(6)));
FSHcirc   = FSHorig   + met(7)*FSHorig.*cos(2*pi*(x - met(8)));

% -------------------------------------------------------------------------
% [D] Solve the model with circadian forcing (parameters ≠ 0)
% -------------------------------------------------------------------------
solutions = solve_mod(parameters)';
LH = solutions(:,2)';
FSH = solutions(:,4)';
GrF = solutions(:,6)';
DomF = solutions(:,7)';
Lut2 = solutions(:,11)';
Lut3 = solutions(:,12)';
Lut4 = solutions(:,13)';

e0 = 57.60;  
e1 = 0.0269;
e2 = 0.4196; 
e3 = 0.4923;
p1 = 0.0032;
p2 = 0.1188; 


E2bef  = e0  + e1*GrF  + e2*DomF  + e3*Lut4;
P4bef  = p1*Lut3  + p2*Lut4;

% Final endogenous E2/P4 after circadian modulation applied
E2    = E2bef   + (parameters(1))*E2bef.*cos((2*pi)*(x - parameters(2)));
P4    = P4bef   + (parameters(3))*P4bef.*cos((2*pi)*(x -  parameters(4)));

% -------------------------------------------------------------------------
% [E] Visualization: target vs simulated hormone profiles
%   Red = simulation, Black = target
% -------------------------------------------------------------------------
    figure(1)
    p = plot(x, FSHcirc, 'k', x, FSH, 'r', 'LineWidth',1);
    legend 'Target' 'Simulation'
    set(gca,'FontSize',20)
    xlim([0 28]);
    xticks([0 28]);
    ylim([0 20]);
    yticks([0 20]);
    xlabel('$t$ [days]','Interpreter','latex')
    ylabel('$FSH$ [IU/L]','Interpreter','latex')
    

    figure(2)
    p = plot(x, LHcirc, 'k', x, LH, 'r', 'LineWidth',1);
    legend 'Target' 'Simulation' 
    set(gca,'FontSize',20)
    xlim([0 28]);
    xticks([0 28]);
    ylim([0 150]);
    yticks([0 150]);
    xlabel('$t$ [days]','Interpreter','latex')
    ylabel('$LH$ [IU/L]','Interpreter','latex')

 
    figure(3)
    p = plot(x, E2circ, 'k', x, E2, 'r', 'LineWidth',1);
    legend 'Target' 'Simulation' 
    set(gca,'FontSize',20)
    xlim([0 28]);
    xticks([0 28]);
    ylim([0 300]);
    yticks([0 300]);
    yticklabels([0 0.3])
    xlabel('$t$ [days]','Interpreter','latex')
    ylabel('$E_2$ [ng/mL]','Interpreter','latex')
 
    figure(4)
    p = plot(x, P4circ, 'k', x, P4, 'r', 'LineWidth',1);
    legend 'Target' 'Simulation' 
    set(gca,'FontSize',20)
    xlim([0 28]);
    xticks([0 28]);
    ylim([0 20]);
    yticks([0 20]);
    xlabel('$t$ [days]','Interpreter','latex')
    ylabel('$P_4$ [ng/mL]','Interpreter','latex')
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:
%   Solves the endogenous HPO-axis model with optional circadian forcing
%   using MATLAB's delay differential equation solver `dde23`.
%
% INPUT:
%   parameterset (8x1 vector)
%     → Circadian modulation parameters:
%         param(1–2) = amplitude, phase of E2 forcing
%         param(3–4) = amplitude, phase of P4 forcing
%         param(5–6) = amplitude, phase of LH synthesis
%         param(7–8) = amplitude, phase of FSH synthesis
%
% OUTPUT:
%   solution (matrix)
%     → State trajectories evaluated on the time grid x = 0:0.01:28 days
%
% =========================================================================

function [solution] = solve_mod(parameterset)

% -------------------------------------------------------------------------
% Time vector: simulation over a single 28-day cycle
% Step size = 0.01 day (~14.4 minutes resolution)
% -------------------------------------------------------------------------
x    = 0:0.01:28;
tdata = [0 28];

% -------------------------------------------------------------------------
% Delay used in FSH production term (inhibin feedback delay)
% Value = 1.5 days based on prior model calibration
% -------------------------------------------------------------------------
dInh =  1.5;

% These values represent steady-state initial conditions for the 13 state variables
RPLH0 = 167.57;     		% releasable LH           
LH0 =   11.81;                  % circulating LH
RPFSH0 = 14.48;                 % releasable FSH
FSH0 =  11.41;                  % circulating FSH
ReF0 =  2.10;                   % recruited follicles
SeF0 =  4.12;                   % growing follicles
PrF0 =  0.46;          		% dominant follicles
Ov10 =  1.06;                   % ovulatory stage tissue
Ov20 =  1.67;                   % luteinization stage tissue
Lut10 = 4.16;                   % luteal stage 1 tissue
Lut20 = 13.03;                  % luteal stage 2 tissue
Lut30 = 16.48;                  % luteal stage 3 tissue
Lut40 =  10.29; 		% luteal stage 4 tissue

Init = [RPLH0; LH0; RPFSH0; FSH0; ReF0; SeF0; PrF0; 
        Ov10; Ov20; Lut10; Lut20; Lut30; Lut40];

% -------------------------------------------------------------------------
% Call dde23 to solve the model
% - @model: function handle to the DDE RHS
% - dInh: single delay for inhibin feedback
% - Init: history function = constant initial value for all t ≤ 0
% -------------------------------------------------------------------------
solution = dde23(@model, dInh, Init, tdata, [], parameterset);

% Evaluate the numerical solution on the uniform time grid
solution = deval(solution, x);	% output = [states x time] matrix

end


% =========================================================================
% model.m
%
% PURPOSE:
%   Defines the right-hand side of the DDE system
%   modeling the hormonal dynamics of the HPO axis.
%   Includes time-varying (circadian) modulation of hormone synthesis and feedback.
%
% STRUCTURE:
%   - Pituitary hormone production (RPLH, RPFSH → LH, FSH)
%   - Ovarian follicular development and luteal progression
%   - Steroid hormone (E2, P4) reconstruction from ovarian states
%   - Circadian forcing added via cosine terms in:
%       • Endogenous steroid synthesis
%       • LH and FSH synthesis rates
%
% INPUTS:
%   t      = current time (in days)
%   state  = 13x1 vector of state variables
%   delay  = delayed state values (used for inhibin feedback)
%   param  = 8x1 vector of circadian amplitude/phase values
%
% OUTPUT:
%   dstate = 13x1 vector of derivatives
%
% =========================================================================

function dstate = model(t, state, delay, param)

% -------------------------------------------------------------------------
% Unpack circadian modulation parameters from param vector
% These control amplitude and acrophase (peak timing) of:
%   - E2 and P4: param(1–4)
%   - LH and FSH: param(5–8)
% -------------------------------------------------------------------------
par1 =         param(1);  	% E2 amplitude
par2 =         param(2);	% E2 phase
par3 =         param(3);	% P4 amplitude
par4 =         param(4);	% P4 phase
par5 =         param(5);	% LH production amplitude
par6 =         param(6);	% LH phase
par7 =         param(7);	% FSH production amplitude
par8 =         param(8);	% FSH phase

% -------------------------------------------------------------------------
% Fixed model parameters (from prior calibration)
% Full descriptions, units, and physiological roles are detailed in
% Supplementary Section S1 Table 1 of the manuscript.
% -------------------------------------------------------------------------

% LH dynamics
kLH =  0.9661; 
V0LH =  550.03; 
V1LH = 3329.19; 
KmLH =   136.05;                                     
KiLHP =  6.78; 
cLHE =  0.0060; 
cLHP =  1.98; 

% FSH dynamics
VFSH = 294.90; 
kFSH = 14.59; 
cFSHE = 0.0151; 
KiFSHInh = 16.83;              
cFSHP = 52.31;

% Follicle recruitment and ovulation
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

% Steroid synthesis coefficients
e0 = 57.60;  
e1 = 0.0269;
e2 = 0.4196; 
e3 = 0.4923; 
p1 = 0.0032;
p2 = 0.1188; 

% Inhibin synthesis (proxy) coefficients
h0 = 0.6606;  
h1 = 0.0193; 
h2 = 0.0159; 
h3 = 0.0119; 

% Feedback scaling
w = 9.21;	% P4 → FSH inhibition
q = 5.11;       % P4 → follicle inhibition        

% -------------------------------------------------------------------------
% Unpack state variables
% -------------------------------------------------------------------------
RPLH        =  state(1); 
LH          =  state(2); 
RPFSH       =  state(3); 
FSH         =  state(4); 
RcF         =  state(5);   
GrF         =  state(6);
DomF        =  state(7);
Sc1         =  state(8);
Sc2         =  state(9);
Lut1        =  state(10);
Lut2        =  state(11);
Lut3        =  state(12);
Lut4        =  state(13); 

% -------------------------------------------------------------------------
% Reconstruct E2 and P4 from ovarian states + circadian modulation
% Endogenous steroid + time-dependent circadian modulation
% -------------------------------------------------------------------------
E2   =    (e0 + e1*GrF + e2*DomF + e3*Lut4) + par1*(e0 + e1*GrF + e2*DomF + e3*Lut4)*cos((2*pi)*(t -  par2));
P4   =    (p1*Lut3 + p2*Lut4) + par3*(p1*Lut3 + p2*Lut4)*cos((2*pi)*(t - par4));


% -------------------------------------------------------------------------
% Delayed inhibin feedback signal
% Computed using delayed values of DomF, Lut2, Lut3 from 1.5 days earlier
% -------------------------------------------------------------------------

statelag1       = delay(:,1);

DomFdelayInh     = statelag1(7);
Lut2delayInh    = statelag1(11);
Lut3delayInh    = statelag1(12);


InhdelayInh     = h0 + h1*DomFdelayInh + h2*Lut2delayInh + h3*Lut3delayInh;

% -------------------------------------------------------------------------
% Pituitary hormone dynamics
% Includes circadian forcing on production terms (via cosine)
% -------------------------------------------------------------------------

% RPLH synthesis → modulated by E2 stimulation and P4 inhibition

dRPLHdt     =   ((V0LH + (V1LH*E2^8/(KmLH^8 + E2^8)))/(1 + (P4/KiLHP))).*(1 + par5*cos((2*pi)*(t-par6))) - ...
                ((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2));

% LH release and clearance
dLHdt       =   (1/2.5)*((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2)) - 14*LH;

% RPFSH synthesis → inhibited by inhibin and P4, with circadian forcing
dRPFSHdt    =   (VFSH/(1 + (InhdelayInh/KiFSHInh) + (P4/w))).*(1 + par7*cos((2*pi)*(t-par8))) - ...
                 ((kFSH*(1 + cFSHP*P4).*RPFSH)./(1+cFSHE*E2^2));

% FSH release and clearance
dFSHdt      =   ((1/2.5)*((kFSH*(1 + cFSHP*P4)*RPFSH)/(1 + cFSHE*E2^2)) -  8.21*FSH); 

% -------------------------------------------------------------------------
% Ovarian follicle-luteal cascade
% Recruitment → Growth → Dominance → Ovulation → Luteal transformation
% -------------------------------------------------------------------------
dRcFdt      =   (b + c1*RcF)*(FSH/(1 + (P4/q))) - c2*(LH^alpha)*RcF;
dGrFdt      =    c2*(LH^alpha)*RcF - c3*LH*GrF; 
dDomFdt     =    c3*LH*GrF - c4*(LH^gamma)*DomF;
dSc1dt      =    c4*(LH^gamma)*DomF - d1*Sc1;
dSc2dt      =    d1*Sc1 - d2*Sc2;
dLut1dt     =    d2*Sc2 - k1*Lut1;
dLut2dt     =    k1*Lut1 - k2*Lut2;
dLut3dt     =    k2*Lut2 - k3*Lut3;
dLut4dt     =    k3*Lut3 - k4*Lut4;


% -------------------------------------------------------------------------
% Return derivatives as a column vector
% -------------------------------------------------------------------------
dstate      = zeros(13,1);

dstate(1)   = dRPLHdt;
dstate(2)   = dLHdt;
dstate(3)   = dRPFSHdt;
dstate(4)   = dFSHdt;
dstate(5)   = dRcFdt;
dstate(6)   = dGrFdt;
dstate(7)   = dDomFdt;
dstate(8)   = dSc1dt;
dstate(9)   = dSc2dt;
dstate(10)  = dLut1dt;
dstate(11)  = dLut2dt;
dstate(12)  = dLut3dt;
dstate(13)  = dLut4dt;

end
    



