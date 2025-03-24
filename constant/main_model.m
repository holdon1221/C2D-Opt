function main_model

parameters = zeros(9*28,1);
E2_t0 = 1; 
E2_tf = 21; 
P4_t0 = E2_t0 + 3*28;
P4_tf = P4_t0 + E2_tf - E2_t0;

oras = 10;                                      % input time
parameters(1+6*28:9*28) = round(oras/24,3);

dose_EE2 = 21*1000;                             % input EE dose in ng

   parameters(E2_t0:E2_tf) = dose_EE2;   
   parameters(E2_t0+28:E2_tf+28) = dose_EE2; 
   parameters(E2_t0+2*28:E2_tf+2*28) = dose_EE2; 
   
dose_DNG = 0*1000;                              % input DNG dose in ng
   parameters(P4_t0:P4_tf) = dose_DNG;
   parameters(P4_t0+28:P4_tf+28) = dose_DNG;
   parameters(P4_t0+2*28:P4_tf+2*28) = dose_DNG;
  
parameterssquared = [0.0301; 1.0987; 0.0773; 0.3329; 0.0111; 0.5541; 0.0100; 1.5843];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x    = 0:0.001:3*28;

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


solutions = solve_mod(parameters)';                                        
LH = solutions(:,2)';
FSH = solutions(:,4)';
GrF = solutions(:,6)';
DomF = solutions(:,7)';
Lut2 = solutions(:,11)';
Lut3 = solutions(:,12)';
Lut4 = solutions(:,13)';


E2bef  = e0  + e1*GrF  + e2*DomF  + e3*Lut4;
P4bef  = p1*Lut3  + p2*Lut4;


tend = 3*28;
tstep = 0.001;
xt = 0:tstep:tend; 
  
pkE = [0.65; 0; 0];                                     
pkP = [0.90; 0; 0];                                     

pkE(2) = 84.15;   
pkE(3) = 2.77; 
VDE2 =  2800*60;

pkP(2) = 79;
pkP(3) = 1.39; 
VDP4  = 40000; 

parame1(1:3*28)        =  parameters(1:tend);
parame1(3*28+1:6*28)   =  parameters(tend+1:2*tend);
parame1(6*28+1:9*28)   = parameters(2*tend+1:3*tend);

totalind   = ((tend)/tstep) + 1;

CE2 = zeros((tend), totalind);
CP4 = zeros((tend), totalind);
 

for la = 1:tend
    ind = round((la -1 + parame1(2*(tend)+la))*1000 + 1);                 % index corresponding to intake time
        CE2(la, ind:totalind) = (((pkE(1))*(pkE(2))*parame1(la))/(VDE2*((pkE(2))-(pkE(3))))).*((exp(-(pkE(3))*(xt(1:(totalind +1-ind)))))-(exp(-(pkE(2))*(xt(1:(totalind +1-ind))))));    
        CP4(la, ind:totalind)  =   (((pkP(1))*(pkP(2))*parame1(la+ (tend)))/(VDP4*((pkP(2))-(pkP(3))))).*((exp(-(pkP(3))*(xt(1:(totalind +1-ind)))))-(exp(-(pkP(2))*(xt(1:(totalind+1-ind))))));   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalCE2 = zeros(1, totalind);
totalCP4 = zeros(1, totalind);

for su = 1:(tend)
    totalCE2 = totalCE2 + CE2(su, 1:totalind);
    totalCP4 = totalCP4 + CP4(su, 1:totalind);
end


E2     = E2bef     + (parameterssquared(1))*E2bef.*cos(2*pi*(x - (parameterssquared(2))));  % endogenous E2
P4     = P4bef    + (parameterssquared(3))*P4bef.*cos(2*pi*(x - (parameterssquared(4))));   % endogenous P4

Inh           = h0 + h1*DomF + h2*Lut2 + h3*Lut3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


parame = parameters;


figure(1)
plot(x, FSH, 'k','LineWidth', 1);  
set(gca,'FontSize',20)
xlim([0 max(x)]);
xticks([0 28 56 84]);
xlabel('$t$ [days]', 'Interpreter','latex')
ylabel('$FSH$ [IU/L]', 'Interpreter','latex')
 
figure(2)
plot(x, LH, 'k','LineWidth', 1);  
set(gca,'FontSize',20)
xlim([0 max(x)]);
xticks([0 28 56 84]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$LH$ [IU/L]','Interpreter','latex')

figure(3)
plot(x, E2, 'k','LineWidth', 1);  
set(gca,'FontSize',20)
xlim([0 max(x)]);
xticks([0 28 56 84]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$E_2$ [ng/mL]','Interpreter','latex')

figure(4)
plot(x, P4, 'k','LineWidth', 1);  set(gca,'FontSize',20)
hold on
yline(3)
set(gca,'FontSize',20)
xlim([0 max(x)]);
xticks([0 28 56 84]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$P_4$ [ng/mL]','Interpreter','latex')

figure(5)
plot(x, Inh, 'k','LineWidth', 1);  set(gca,'FontSize',20)
xlim([0 max(x)]);
xticks([0 28 56 84]);
xlabel('$t$ [days]','Interpreter','latex')
ylabel('$Inh$ [IU/mL]','Interpreter','latex')


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [solution] = solve_mod(parameters)

x    = 0:0.001:3*28;
tdata = [0 3*28];

dInh =  1.5; 
Init = [167.57; 11.81; 14.48; 11.41; 2.10; 4.12; 0.46; 1.06; 1.67; 4.16; 13.03; 16.48; 10.29];
u = parameters;

solution = dde23(@model, dInh, Init, tdata, [], u);
solution = deval(solution, x);

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

additional_par =  [0.0301; 1.0987; 0.0773; 0.3329; 0.0111; 0.5541; 0.0100; 1.5843];


par1 =         additional_par(1);
par2 =         additional_par(2);
par3 =         additional_par(3);
par4 =         additional_par(4);
par5 =         additional_par(5);
par6 =         additional_par(6);
par7 =         additional_par(7);
par8 =         additional_par(8);


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

tstart = 0;
tend = 3*28;
tend1 = tend + 28;


ConE2 = zeros(1,tend1);
ConP4 = zeros(1,tend1);


pkE = [0.65; 0; 0];                                     
pkP = [0.90; 0; 0];                                    

 pkE(2) = 84.15;   
 pkE(3) = 2.77;  
 VdE2 =  2800*60;

pkP(2) = 79;
pkP(3) = 1.39; 
VdP4  = 40000; 


param1(1:3*28) = u(1:tend);
param1(3*28+1:6*28) = u(tend+1:2*tend);
param1(6*28+1:9*28) = u(2*tend+1:3*tend);  


for la = 1:tend 
    if (0 <= t) && (t < (param1(2*(tend)+la)+ la-1))
        ConE2(la) = 0;
        ConP4(la) = 0;
    else
        ConE2(la)  = (((pkE(1))*(pkE(2))*param1(la))/(VdE2*((pkE(2))-(pkE(3))))).*((exp(-(pkE(3))*(t-(param1(2*(tend)+la)+ la-1))))-...
                  (exp(-(pkE(2))*(t-(param1(2*(tend)+la)+ la-1)))));
        ConP4(la)  = (((pkP(1))*(pkP(2))*param1(la+ (tend)))/(VdP4*((pkP(2))-(pkP(3))))).*((exp(-(pkP(3))*(t-(param1(2*(tend)+la)+ la-1))))-...
                  (exp(-(pkP(2))*(t-(param1(2*(tend)+la)+ la-1)))));
    end
end


    totalConE2 = sum(ConE2);
    totalConP4 = sum(ConP4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


E2              =    (e0 + e1*GrF + e2*DomF + e3*Lut4) + par1*(e0 + e1*GrF + e2*DomF + e3*Lut4)*cos(2*pi*(t-par2)) + 1.7*totalConE2;
P4              =    (p1*Lut3 + p2*Lut4) + par3*(p1*Lut3 + p2*Lut4)*cos(2*pi*(t-par4)) + 0.01*totalConP4;


statelag1       = delay(:,1);

DomFdelayInh     = statelag1(7);
Lut2delayInh    = statelag1(11);
Lut3delayInh    = statelag1(12);


InhdelayInh     = h0 + h1*DomFdelayInh + h2*Lut2delayInh + h3*Lut3delayInh;

dRPLHdt     =   ((V0LH + (V1LH*E2^8/(KmLH^8 + E2^8)))/(1 + (P4/KiLHP))).*(1 + par5*LH*cos((2*pi)*(t-par6))) - ((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2));

dLHdt       =   (1/2.5)*((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2)) - 14*LH;

dRPFSHdt    =   (VFSH/(1 + (InhdelayInh/KiFSHInh) + (P4/w))).*(1 + par7*FSH*cos((2*pi)*(t-par8))) - ((kFSH*(1 + cFSHP*P4).*RPFSH)./(1+cFSHE*E2^2));

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