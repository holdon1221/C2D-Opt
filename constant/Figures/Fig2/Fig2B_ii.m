% This script draws Fig 2B(ii), depicting simulated hormone profiles without circadian rhythm.
function Fig2B_ii


x    = 0:0.001:28;

% Coefficients for E2 calculation based on follicle stages
e0 = 57.60/1000;         
e1 = 0.0269/1000;        
e2 = 0.4196/1000;        
e3 = 0.4923/1000;   

% Coefficients for P4 calculation
p1 = 0.0032;             
p2 = 0.1188;

% Coefficients for Inh calculation
h0 = 0.6606;  
h1 = 0.0193; 
h2 = 0.0159; 
h3 = 0.0119; 

tdata = [0 28];

% Time delay for Inh feedback (in days)
dInh =  1.5;

% Initial conditions for 13 state variables:
% [RPLH, LH, RPFSH, FSH, RcF, GrF, DomF, Sc1, Sc2, Lut1, Lut2, Lut3, Lut4]
Init = [167.57; 11.81; 14.48; 11.41; 2.10; 4.12; 0.46; 1.06; 1.67; 4.16; 13.03; 16.48; 10.29];

solution = dde23(@model, dInh, Init, tdata);
solutions = deval(solution, x);
LH = solutions(2,:)';
FSH = solutions(4,:)';
GrF = solutions(6,:)';
DomF = solutions(7,:)';
Lut2 = solutions(11,:)';
Lut3 = solutions(12,:)';
Lut4 = solutions(13,:)';


E2  = e0  + e1*GrF  + e2*DomF  + e3*Lut4;
P4  = p1*Lut3  + p2*Lut4;
Inh           = h0 + h1*DomF + h2*Lut2 + h3*Lut3;



figure(1)
plot(x, FSH, 'k','LineWidth', 1);  
set(gca,'FontSize',20)
xlim([0 max(x)]);
xticks([0 28]);
ylim([0 20]);
yticks([0 20]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$FSH$ [IU/L]','Interpreter','latex')


figure(2)
plot(x, LH, 'k','LineWidth', 1);  
set(gca,'FontSize',20)
xlim([0 max(x)]);
xticks([0 28]);
ylim([0 150]);
yticks([0 150]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$LH$ [IU/L]','Interpreter','latex')

figure(3)
plot(x, E2, 'k','LineWidth', 1);  
set(gca,'FontSize',20)
xlim([0 max(x)]);
xticks([0 28]);
ylim([0 0.3]);
yticks([0 0.3]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$E_2$ [ng/mL]','Interpreter','latex')

figure(4)
plot(x, P4, 'k','LineWidth', 1);  set(gca,'FontSize',20)
hold on
set(gca,'FontSize',20)
xlim([0 max(x)]);
xticks([0 28]);
yticks([0 20]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$P_4$ [ng/mL]','Interpreter','latex')


end


function dstate = model(t, state, delay, u)

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
w = 9.21;            
q = 5.11;            


RPLH            =  state(1); 
LH              =  state(2); 
RPFSH           =  state(3); 
FSH             =  state(4); 
RcF             =  state(5);   
GrF             =  state(6);
DomF            =  state(7);
Sc1             =  state(8);
Sc2             =  state(9);
Lut1            =  state(10);
Lut2            =  state(11);
Lut3            =  state(12);
Lut4            =  state(13);



E2              =    e0 + e1*GrF + e2*DomF + e3*Lut4;
P4              =    p1*Lut3 + p2*Lut4;


statelag1       = delay(:,1);

DomFdelayInh     = statelag1(7);
Lut2delayInh    = statelag1(11);
Lut3delayInh    = statelag1(12);


InhdelayInh     = h0 + h1*DomFdelayInh + h2*Lut2delayInh + h3*Lut3delayInh;

dRPLHdt     =   ((V0LH + (V1LH*E2^8/(KmLH^8 + E2^8)))/(1 + (P4/KiLHP))) - ((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2));

dLHdt       =   (1/2.5)*((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2)) - 14*LH;

dRPFSHdt    =   (VFSH/(1 + (InhdelayInh/KiFSHInh) + (P4/w))) - ((kFSH*(1 + cFSHP*P4).*RPFSH)./(1+cFSHE*E2^2));

dFSHdt      =   (1/2.5)*((kFSH*(1 + cFSHP*P4)*RPFSH)/(1 + cFSHE*E2^2)) -  8.21*FSH;  

dRcFdt      =   (b + c1*RcF)*(FSH/(1 + (P4/q))) - c2*(LH^alpha)*RcF;
dGrFdt      =    c2*(LH^alpha)*RcF - c3*LH*GrF; 
dDomFdt     =    c3*LH*GrF - c4*(LH^gamma)*DomF;
dSc1dt          =    c4.*(LH.^gamma).*DomF - d1.*Sc1;
dSc2dt          =    d1.*Sc1 - d2.*Sc2;
dLut1dt         =    d2.*Sc2 - k1.*Lut1;
dLut2dt         =    k1.*Lut1 - k2.*Lut2;
dLut3dt         =    k2.*Lut2 - k3.*Lut3;
dLut4dt         =    k3.*Lut3 - k4.*Lut4;


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
