% This code produces a parameter set for amplitude and acrophases of E2, P4, LH, and FSH circadian components of the menstrual cycle model.

function circ_par_est

clear all;
format long


tic

parameters = 0.05*ones(8,1);    % input initial guess

parameters = sqrt(parameters);  % to avoid negative parameters

%parameter arrangement = [E2amp; E2acro; P4amp; P4acro; LHamp; LHacro; FSHamp; FSHacro]; 
                                              
opts = optimset('MaxIter', 100000,'MaxFunEvals', 100000);
[optimal_par, errorval, exitflag, ~] = fminsearch(@fminmerged,parameters,opts);

parameters = optimal_par;
          
parameterssquared = parameters.^2		        % optimal parameter set

toc

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ls,sol] = fminmerged(parameters)

para1 = parameters(1).^2;
para2 = parameters(2).^2;
para3 = parameters(3).^2;
para4 = parameters(4).^2;

x    = [0:0.01:56];

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
       

solution = solve_mod(parameters)';
LH = solution(:,2)';
FSH = solution(:,4)';
GrF = solution(:,6)';
DomF = solution(:,7)';
Lut2 = solution(:,11)';
Lut3 = solution(:,12)';
Lut4 = solution(:,13)';

E2es  = e0  + e1*GrF  + e2*DomF  + e3*Lut4;
P4pr  = p1*Lut3  + p2*Lut4;
Inh   = h0 + h1*DomF + h2*Lut2 + h3*Lut3;

E2    = E2es + para1*E2es.*cos(2*pi*(x  - para2));
P4    = P4pr + para3*P4pr.*cos(2*pi*(x  - para4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parametersc = [0; 0; 0; 0; 0; 0; 0; 0];   

solutionold = solve_mod(parametersc)';
LHold  = solutionold(:,2)';
FSHold = solutionold(:,4)';
GrFold = solutionold(:,6)';
DomFold = solutionold(:,7)';
Lut2old = solutionold(:,11)';
Lut3old = solutionold(:,12)';
Lut4old = solutionold(:,13)';

E2old  = e0  + e1*GrFold  + e2*DomFold  + e3*Lut4old;
P4old  = p1*Lut3old  + p2*Lut4old;

met(1) = 0.08026;               % input E2 cosine fit amplitude
met(2) = 24.58/24;              % input E2 cosine fit acrophase
met(3) = 0.1017;                % input P4 cosine fit amplitude
met(4) = 8.46/24;               % input P4 cosine fit acrophase
met(5) = 0.0531;                % input LH cosine fit amplitude
met(6) = 18.27/24;              % input LH cosine fit acrophase
met(7) = 0.0659;                % input FSH cosine fit amplitude
met(8) = 16.46/24;              % input FSH cosine fit acrophase

E2circ2    = E2old    + met(1)*E2old.*cos(2*pi*(x - met(2)));
P4circ2    = P4old    + met(3)*P4old.*cos(2*pi*(x - met(4)));
LHcirc2    = LHold    + met(5)*LHold.*cos(2*pi*(x - met(6)));
FSHcirc2   = FSHold   + met(7)*FSHold.*cos(2*pi*(x - met(8)));

error1 = (E2 - E2circ2)./E2circ2;
error2 = (P4 - P4circ2)./P4circ2;
error3 = (LH - LHcirc2)./LHcirc2;
error4 = (FSH - FSHcirc2)./FSHcirc2;


ls1 = error1*error1';
ls2 = error2*error2';
ls3 = error3*error3';
ls4 = error4*error4';

ls = ls1 + ls2 + ls3 + ls4;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [solution] = solve_mod(parameterset)

x    = 0:0.01:56 ;
tdata = [0 56];

dInh =  1.5; 
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
param = parameterset.^2;

solution = dde23(@model, dInh, Init, tdata, [], param);
solution = deval(solution, x);

end



function dstate = model(t, state, delay, param)

par1 =         param(1);
par2 =         param(2);
par3 =         param(3);
par4 =         param(4);
par5 =         param(5);
par6 =         param(6);
par7 =         param(7);
par8 =         param(8);


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


dRPLHdt     =   ((V0LH + (V1LH*E2^8/(KmLH^8 + E2^8)))/(1 + (P4/KiLHP))).*(1 + par5*LH*cos((2*pi)*(t-par6))) - ...
                ((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2));

dLHdt       =   (1/2.5)*((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2)) - 14*LH;

dRPFSHdt    =   (VFSH/(1 + (InhdelayInh/KiFSHInh) + (P4/w))).*(1 + par7*FSH*cos((2*pi)*(t-par8))) - ...
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
    



