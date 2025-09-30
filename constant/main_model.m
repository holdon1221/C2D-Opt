function main_model

parameters = zeros(3*10*28,1);


oras = 22;                                      			               % input dosing time
parameters(1+2*10*28:3*10*28) = round(oras/24,2);

dose_EE2 = 28.2*1000;                   	                         % input EE dose in ng

dose_DNG = 320*1000;                  		                          % input DNG dose in ng
 
   for araw = 1:10  							                                       % number of cycles
        E2_t0 = 1;
        E2_tf = 21;
        parameters(E2_t0 + (araw-1)*28:E2_tf+ (araw-1)*28) =  dose_EE2;                 

        P4_t0 = E2_t0 + 28*10;           
        P4_tf = P4_t0 + E2_tf - E2_t0;
        parameters(P4_t0+ (araw-1)*28:P4_tf+ (araw-1)*28)  =  dose_DNG;               
   end   


parameterssquared = [0.0473; 0.0915; 0.0807; 0.3487; 0.2242; 0.5730; 0.1130; 0.5480];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x    = 0:0.01:10*28;

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


tend = 10*28;
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

parame1(1:10*28)        =  parameters(1:tend);
parame1(10*28+1:2*10*28)   =  parameters(tend+1:2*tend);
parame1(2*10*28+1:3*10*28)   = parameters(2*tend+1:3*tend);

totalind   = ((tend)/tstep) + 1;

CE2 = zeros((tend), totalind);
CP4 = zeros((tend), totalind);
 

for la = 1:tend
    ind = round((la -1 + parame1(2*(tend)+la))*100 + 1);                 % index corresponding to intake time
        NEE = (k21EE - kaEE)/( (alpha1EE - kaEE)*(beta1EE - kaEE) );
        LEE = (k21EE - alpha1EE)/( (kaEE - alpha1EE)*(beta1EE - alpha1EE) );
        MEE = (k21EE - beta1EE)/( (kaEE - beta1EE)*(alpha1EE - beta1EE) );

        CE2(la, ind:totalind)  =  ((kaEE*FEE*parame1(la))/VcEE) * ( NEE*exp(-kaEE*(xt(1:(totalind +1-ind)))) + LEE*exp(-alpha1EE*(xt(1:(totalind +1-ind)))) + MEE*exp(-beta1EE*(xt(1:(totalind +1-ind)))) );        

        NDNG = (k21DNG - kaDNG)/( (alpha1DNG - kaDNG)*(beta1DNG - kaDNG) );
        LDNG = (k21DNG - alpha1DNG)/( (kaDNG - alpha1DNG)*(beta1DNG - alpha1DNG) );
        MDNG = (k21DNG - beta1DNG)/( (kaDNG - beta1DNG)*(alpha1DNG - beta1DNG) );

        CP4(la, ind:totalind)  =  ((kaDNG*FDNG*parame1(la + tend))/VcDNG) * ( NDNG*exp(-kaDNG*(xt(1:(totalind +1-ind)))) + LDNG*exp(-alpha1DNG*(xt(1:(totalind +1-ind)))) + MDNG*exp(-beta1DNG*(xt(1:(totalind +1-ind)))) );        
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

E2tot  = E2 + 1.7*totalCE2;
P4tot  = P4 + 0.01*totalCP4;


Inh           = h0 + h1*DomF + h2*Lut2 + h3*Lut3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

x    = 0:0.01:10*28;
tdata = [0 10*28];

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

additional_par =  [0.0473; 0.0915; 0.0807; 0.3487; 0.2242; 0.5730; 0.1130; 0.5480];


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
tend = 10*28;
tend1 = tend + 28;


ConE2 = zeros(1,tend1);
ConP4 = zeros(1,tend1);


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


param1(1:10*28) = u(1:tend);
param1(10*28+1:2*10*28) = u(tend+1:2*tend);
param1(2*10*28+1:3*10*28) = u(2*tend+1:3*tend);  


for la = 1:tend 
    if (0 <= t) && (t < (param1(2*(tend)+la)+ la-1))
        ConE2(la) = 0;
        ConP4(la) = 0;
   
   elseif (param1(2*tend+la)+ la-1) <= t
        
        NEE = (k21EE - kaEE)/( (alpha1EE - kaEE)*(beta1EE - kaEE) );
        LEE = (k21EE - alpha1EE)/( (kaEE - alpha1EE)*(beta1EE - alpha1EE) );
        MEE = (k21EE - beta1EE)/( (kaEE - beta1EE)*(alpha1EE - beta1EE) );

        ConE2(la)  =  ((kaEE*FEE*param1(la))/VcEE) * ( NEE*exp(-kaEE*(t-(param1(2*tend+1)+ la-1))) + LEE*exp(-alpha1EE*((t-(param1(2*tend+1)+ la-1)))) + MEE*exp(-beta1EE*((t-(param1(2*tend+1)+ la-1)))) );        

        NDNG = (k21DNG - kaDNG)/( (alpha1DNG - kaDNG)*(beta1DNG - kaDNG) );
        LDNG = (k21DNG - alpha1DNG)/( (kaDNG - alpha1DNG)*(beta1DNG - alpha1DNG) );
        MDNG = (k21DNG - beta1DNG)/( (kaDNG - beta1DNG)*(alpha1DNG - beta1DNG) );

        ConP4(la)  =  ((kaDNG*FDNG*param1(la + tend))/VcDNG) * ( NDNG*exp(-kaDNG*((t-(param1(2*tend+1)+ la-1)))) + LDNG*exp(-alpha1DNG*((t-(param1(2*tend+1)+ la-1)))) + MDNG*exp(-beta1DNG*((t-(param1(2*tend+1)+ la-1)))) );        
    else
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
