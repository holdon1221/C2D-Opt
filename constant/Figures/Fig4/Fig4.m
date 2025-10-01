function Fig4


parameters = zeros(3*3*28,1);
E2_t0 = 1; 
E2_tf = 21; 
P4_t0 = E2_t0 + 3*28;
P4_tf = P4_t0 + E2_tf - E2_t0;

P4res    = zeros(2, 8401);
RPLHcirc = zeros(2, 8401);
EE_con   = zeros(2, 8401);
DNG_con   = zeros(2, 8401);
RPLHres  = zeros(2, 8401);
LHres    = zeros(2, 8401);
RcFres   = zeros(2, 8401);
Lut3res  = zeros(2, 8401);
Lut4res  = zeros(2, 8401);

time = [11,22];

for i = 1:2
  oras = time(i);                                          % input time
  parameters(1+6*28:9*28) = round(oras/24,2);

dose_EE2 = 28.2*1000;                                      % input EE dose in ng

   parameters(E2_t0:E2_tf) = dose_EE2;   
   parameters(E2_t0+28:E2_tf+28) = dose_EE2; 
   parameters(E2_t0+2*28:E2_tf+2*28) = dose_EE2;   
   
dose_DNG = 320*1000;                                        % input DNG dose in ng
   parameters(P4_t0:P4_tf) = dose_DNG;
   parameters(P4_t0+28:P4_tf+28) = dose_DNG;
   parameters(P4_t0+2*28:P4_tf+2*28) = dose_DNG;
  
parameterssquared = [0; 0.0915; 0; 0.3487; 0.2242; 0.5730; 0; 0.5480];



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
RPLH = solutions(:,1)';
LH = solutions(:,2)';
FSH = solutions(:,4)';
RcF = solutions(:,5)';
GrF = solutions(:,6)';
DomF = solutions(:,7)';
Lut2 = solutions(:,11)';
Lut3 = solutions(:,12)';
Lut4 = solutions(:,13)';


E2bef  = e0  + e1*GrF  + e2*DomF  + e3*Lut4;
P4bef  = p1*Lut3  + p2*Lut4;


tend = 3*28;
tstep = 0.01;
xt = 0:tstep:tend; 
  
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

totalind   = ((tend)/tstep) + 1;

CE2 = zeros((tend), totalind);
CP4 = zeros((tend), totalind);
 

for la = 1:tend
    ind = round((la -1 + parameters(2*(tend)+la))*100 + 1);                 % index corresponding to intake time
        NEE = (k21EE - kaEE)/( (alpha1EE - kaEE)*(beta1EE - kaEE) );
        LEE = (k21EE - alpha1EE)/( (kaEE - alpha1EE)*(beta1EE - alpha1EE) );
        MEE = (k21EE - beta1EE)/( (kaEE - beta1EE)*(alpha1EE - beta1EE) );

        CE2(la, ind:totalind)  =  ((kaEE*FEE*parameters(la))/VcEE) * ( NEE*exp(-kaEE*(xt(1:(totalind +1-ind)))) + LEE*exp(-alpha1EE*(xt(1:(totalind +1-ind)))) + MEE*exp(-beta1EE*(xt(1:(totalind +1-ind)))) );        

        NDNG = (k21DNG - kaDNG)/( (alpha1DNG - kaDNG)*(beta1DNG - kaDNG) );
        LDNG = (k21DNG - alpha1DNG)/( (kaDNG - alpha1DNG)*(beta1DNG - alpha1DNG) );
        MDNG = (k21DNG - beta1DNG)/( (kaDNG - beta1DNG)*(alpha1DNG - beta1DNG) );

        CP4(la, ind:totalind)  =  ((kaDNG*FDNG*parameters(la + tend))/VcDNG) * ( NDNG*exp(-kaDNG*(xt(1:(totalind +1-ind)))) + LDNG*exp(-alpha1DNG*(xt(1:(totalind +1-ind)))) + MDNG*exp(-beta1DNG*(xt(1:(totalind +1-ind)))) );        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalCE2 = zeros(1, totalind);
totalCP4 = zeros(1, totalind);

for su = 1:(tend)
    totalCE2 = totalCE2 + CE2(su, 1:totalind);
    totalCP4 = totalCP4 + CP4(su, 1:totalind);
end


RPLHcirc(i,:) = 1 + parameterssquared(5)*cos((2*pi)*(xt-parameterssquared(6)));
EE_con(i,:)   = totalCE2;
DNG_con(i,:)   = totalCP4;
RPLHres(i,:)  = RPLH;
LHres(i,:)    = LH;
RcFres(i,:)   = RcF;
Lut3res(i,:)  = Lut3;
Lut4res(i,:)  = Lut4;


E2     = E2bef     + (parameterssquared(1))*E2bef.*cos(2*pi*(xt - (parameterssquared(2))));  % endogenous E2
P4     = P4bef    + (parameterssquared(3))*P4bef.*cos(2*pi*(xt - (parameterssquared(4))));   % endogenous P4

Inh           = h0 + h1*DomF + h2*Lut2 + h3*Lut3;

P4res(i,:) = P4;

end


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


function [solution] = solve_mod(parameters)

x    = 0:0.01:3*28;
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

additional_par =  [0; 0.0915; 0; 0.3487; 0.2242; 0.5730; 0; 0.5480];


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

tend = 3*28;

ConE2 = zeros(1,tend);
ConP4 = zeros(1,tend);

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

    totalConE2 = sum(ConE2);
    totalConP4 = sum(ConP4);


E2              =    (e0 + e1*GrF + e2*DomF + e3*Lut4) + par1*(e0 + e1*GrF + e2*DomF + e3*Lut4)*cos(2*pi*(t-par2)) + 1.7*totalConE2;
P4              =    (p1*Lut3 + p2*Lut4) + par3*(p1*Lut3 + p2*Lut4)*cos(2*pi*(t-par4)) + 0.01*totalConP4;


statelag1       = delay(:,1);

DomFdelayInh     = statelag1(7);
Lut2delayInh    = statelag1(11);
Lut3delayInh    = statelag1(12);


InhdelayInh     = h0 + h1*DomFdelayInh + h2*Lut2delayInh + h3*Lut3delayInh;

dRPLHdt     =   ((V0LH + (V1LH*E2^8/(KmLH^8 + E2^8)))/(1 + (P4/KiLHP))).*(1 + par5*cos((2*pi)*(t-par6))) - ((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2));

dLHdt       =   (1/2.5)*((kLH*(1 + cLHP*P4)*RPLH)/(1 + cLHE*E2)) - 14*LH;

dRPFSHdt    =   (VFSH/(1 + (InhdelayInh/KiFSHInh) + (P4/w))).*(1 + par7*cos((2*pi)*(t-par8))) - ((kFSH*(1 + cFSHP*P4).*RPFSH)./(1+cFSHE*E2^2));

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