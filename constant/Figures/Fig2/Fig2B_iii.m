%This script produces Fig 2B(iii), illustrating both target and simulated hormones with circadian modulation.
function Fig2B_iii

clear all;
format long

% Fitted circadian parameters for simulation: [E2_amp, E2_phase, P4_amp, P4_phase, LH_amp, LH_phase, FSH_amp, FSH_phase]
parameters = [0.0473; 0.0915; 0.0807; 0.3487; 0.2242; 0.5730; 0.1130; 0.5480];

% Zero modulation parameters for baseline hormone simulation
parametersc = [0; 0; 0; 0; 0; 0; 0; 0]; 
x    = [0:0.01:28];


% Generate target hormone profiles
solutionorig = solve_mod(parametersc)';
LHorig  = solutionorig(:,2)';
FSHorig = solutionorig(:,4)';
GrForig = solutionorig(:,6)';
DomForig = solutionorig(:,7)';
Lut2orig = solutionorig(:,11)';
Lut3orig = solutionorig(:,12)';
Lut4orig = solutionorig(:,13)';

% Parameters for E2 and P4 calculation
e0 = 57.60;  
e1 = 0.0269;
e2 = 0.4196; 
e3 = 0.4923;
p1 = 0.0032;
p2 = 0.1188; 

% Baseline hormone profiles
E2orig  = e0  + e1*GrForig  + e2*DomForig  + e3*Lut4orig;
P4orig  = p1*Lut3orig  + p2*Lut4orig;

% Parameters from cosine fitting of experimental data
met(1) = 0.08026;               % E2 cosine fit amplitude
met(2) = 24.58/24;              % E2 cosine fit acrophase
met(3) = 0.1017;                % P4 cosine fit amplitude
met(4) = 8.46/24;               % P4 cosine fit acrophase
met(5) = 0.0531;                % LH cosine fit amplitude
met(6) = 18.27/24;              % LH cosine fit acrophase
met(7) = 0.0659;                % FSH cosine fit amplitude
met(8) = 16.46/24;              % FSH cosine fit acrophase

% Target hormone profiles
E2circ    = E2orig    + met(1)*E2orig.*cos(2*pi*(x - met(2)));
P4circ    = P4orig    + met(3)*P4orig.*cos(2*pi*(x - met(4)));
LHcirc    = LHorig    + met(5)*LHorig.*cos(2*pi*(x - met(6)));
FSHcirc   = FSHorig   + met(7)*FSHorig.*cos(2*pi*(x - met(8)));

% Simulated profiles With optimized circadian parameters 
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

% Baseline hormone profiles for E2 and P4
E2bef  = e0  + e1*GrF  + e2*DomF  + e3*Lut4;
P4bef  = p1*Lut3  + p2*Lut4;

% E2 and P4 hormone profiles with circadian rhythm
E2    = E2bef   + (parameters(1))*E2bef.*cos((2*pi)*(x - parameters(2)));
P4    = P4bef   + (parameters(3))*P4bef.*cos((2*pi)*(x -  parameters(4)));


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

function [solution] = solve_mod(parameterset)
x    = 0:0.01:28;
tdata = [0 28];

dInh =  1.5;         % Delay for Inh feedback

% Initial values for 13 state variables
RPLH0 = 167.57;                  
LH0 =   11.81;                  
RPFSH0 = 14.48;                 
FSH0 =  11.41;                  
ReF0 =  2.10;                   
SeF0 =  4.12;                  
PrF0 =  0.46;          
Ov10 =  1.06;                   
Ov20 =  1.67;                  
Lut10 = 4.16;                  
Lut20 = 13.03;                  
Lut30 = 16.48;                 
Lut40 =  10.29; 

Init = [RPLH0; LH0; RPFSH0; FSH0; ReF0; SeF0; PrF0; 
        Ov10; Ov20; Lut10; Lut20; Lut30; Lut40];

solution = dde23(@model, dInh, Init, tdata, [], parameterset);
solution = deval(solution, x);

end


function dstate = model(t, state, delay, param)

% Unpack circadian rhythm parameters
par1 =         param(1); par2 =         param(2);    % parameters for E2 amp/phase
par3 =         param(3); par4 =         param(4);    % parameters for P4 amp/phase
par5 =         param(5); par6 =         param(6);    % parameters for LH amp/phase
par7 =         param(7); par8 =         param(8);    % parameters for FSH amp/phase

% Model parameters
kLH =  0.9661; 
V0LH =  550.03; 
V1LH = 3329.19; 
KmLH =   136.05;                                     
KiLHP =  6.78; 
cLHE =  0.0060; 
cLHP =  1.98; 
VFSH = 294.90; 
kFSH = 14.59; 
cFSHE = 0.0151; 
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
e0 = 57.60;  
e1 = 0.0269;
e2 = 0.4196; 
e3 = 0.4923; 
p1 = 0.0032;
p2 = 0.1188; 
h0 = 0.6606;  
h1 = 0.0193; 
h2 = 0.0159; 
h3 = 0.0119; 
w = 9.21;
q = 5.11;                

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



E2   =    (e0 + e1*GrF + e2*DomF + e3*Lut4) + par1*(e0 + e1*GrF + e2*DomF + e3*Lut4)*cos((2*pi)*(t -  par2));
P4   =    (p1*Lut3 + p2*Lut4) + par3*(p1*Lut3 + p2*Lut4)*cos((2*pi)*(t - par4));


statelag1       = delay(:,1);

DomFdelayInh     = statelag1(7);
Lut2delayInh    = statelag1(11);
Lut3delayInh    = statelag1(12);


InhdelayInh     = h0 + h1*DomFdelayInh + h2*Lut2delayInh + h3*Lut3delayInh;


dRPLHdt     =   ((V0LH + (V1LH*E2^8/(KmLH^8 + E2^8)))/(1 + (P4/KiLHP))).*(1 + par5*cos((2*pi)*(t-par6))) - ...
                ((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2));

dLHdt       =   (1/2.5)*((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2)) - 14*LH;

dRPFSHdt    =   (VFSH/(1 + (InhdelayInh/KiFSHInh) + (P4/w))).*(1 + par7*cos((2*pi)*(t-par8))) - ...
                 ((kFSH*(1 + cFSHP*P4).*RPFSH)./(1+cFSHE*E2^2));

dFSHdt      =   ((1/2.5)*((kFSH*(1 + cFSHP*P4)*RPFSH)/(1 + cFSHE*E2^2)) -  8.21*FSH); 

dRcFdt      =   (b + c1*RcF)*(FSH/(1 + (P4/q))) - c2*(LH^alpha)*RcF;
dGrFdt      =    c2*(LH^alpha)*RcF - c3*LH*GrF; 
dDomFdt     =    c3*LH*GrF - c4*(LH^gamma)*DomF;
dSc1dt      =    c4*(LH^gamma)*DomF - d1*Sc1;
dSc2dt      =    d1*Sc1 - d2*Sc2;
dLut1dt     =    d2*Sc2 - k1*Lut1;
dLut2dt     =    k1*Lut1 - k2*Lut2;
dLut3dt     =    k2*Lut2 - k3*Lut3;
dLut4dt     =    k3*Lut3 - k4*Lut4;



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
    



