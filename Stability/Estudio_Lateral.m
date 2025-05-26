%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estudio_Lateral.m   Lateral stability analysis                     % 
%                                                                    %
%              Written by Sergio Esteban Roncero                     %
%              Updated by Sergio Esteban Roncero    Dec 2013         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Flight condition data for ALA VOLADORA                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=32.2;                                             % gravity constant
h=5000 ;                                               % altitude (ft)
ro=0.0020484;                                   % density (slugs/ft^3)
temp=500.8435;                                         % Temperature (R)
R=1716;                        % R-universal gas constant (ft^2/s^2*R)
gammaconstant=1.4;                                      % air constant
airspeed=sqrt(R*gammaconstant*temp);              % air speed (ft/sec)
u1=221;                                     % initial velocity (ft/sec)
mach=u1/airspeed;                                        % mach number
q1barw=.5*ro*(u1^2);                        % initial dynamic pressure
AoA=0;                                          % trim angleof atttack
iw=AoA;                                  % incidence angle of the wing
theta1=AoA*pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%%%%%   data for Cessna 182 - Airplane A Roskam volume 1          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reference geometry
Sw=174;
cbar=4.9;
span=36.0;

% flight condition
u1=220.1;
q1barw=49.6;

% mass data
w=2650;     % lb
m=w/32.2;   % lb
Ixx=948;    %slug*ft^2
Iyy= 1346;  %slug*ft^2
Izz=1967;   %slug*ft^2
Ixz=0;      %slug*ft^2

cl1 = 0.307;
cd1= 0.032;
ctx1= cd1;
cm1=0;
cmt1=-cm1;   % no tail

% lateral coefficients

% side-slip derivatives
cybeta=-.393;
clbeta=-.0923;
cnbeta=.0587;
cntbeta=0;

% `roll derivatives
clp=-.484;
cyp=-.075;
cnp=-.0278;

% yaw derivatives
cyr=.214;
clr=.0798;
cnr=-.0937;

% Lateral control moment derivatives %
cyda=0;
cydr=.187;
clda=.229;
cldr=0.147;
cnda=-.0216;
cndr=-.0645;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lateral Dimensional Stability Derivatives      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ybeta = (q1barw*Sw*cybeta)/m;
Yp = (q1barw*Sw*span*cyp)/(2*m*u1);
Yr = (q1barw*Sw*span*cyr)/(2*m*u1);
Yda = (q1barw*Sw*cyda)/(m);
Ydr = (q1barw*Sw*cydr)/(m);

Lbeta = (q1barw*Sw*span*clbeta)/Ixx;
Lp = (q1barw*Sw*span*span*clp)/(2*Ixx*u1);
Lr = (q1barw*Sw*span*span*clr)/(2*Ixx*u1);
Lda = (q1barw*Sw*span*clda)/(Ixx);
Ldr = (q1barw*Sw*span*cldr)/(Ixx);

Nbeta = (q1barw*Sw*span*cnbeta)/Izz;
Ntbeta = (q1barw*Sw*span*cntbeta)/Izz;
Np = (q1barw*Sw*span*span*cnp)/(2*Izz*u1);
Nr = (q1barw*Sw*span*span*cnr)/(2*Izz*u1);
Nda = (q1barw*Sw*span*cnda)/(Izz);
Ndr = (q1barw*Sw*span*cndr)/(Izz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms of the Matrix that finds the poles of the system response %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1=Ixz/Ixx;
B1=Ixz/Izz;

aal = (Ybeta)/u1;
bbl = Yp;
ccl = Yr-u1;
ddl = g*cos(theta1);
ddl2= 0;
eel = ((Lbeta + A1*(Nbeta+Ntbeta))/(1-A1*B1))/u1;
ffl = (Lp + A1*Np)/(1-A1*B1);
ggl = (Lr + A1*Nr)/(1-A1*B1);
hhl = 0;
hhl2 = 0;
iil = ((Nbeta + Ntbeta + B1*Lbeta)/(1-A1*B1))/u1;
jjl = (Np + B1*Lp)/(1-A1*B1);
kkl = (Nr + B1*Lr)/(1-A1*B1);
lll = 0;
lll2 = 0;
mml = 0;
nnl = 1;
ool = tan(theta1);
ppl = 0;
ppl2= 0;
fifth1 = 0;
fifth2 = 0;
fifth3 = 1/cos(theta1);
fifth4 = 0;
fifth5 = 0;
qql = Yda;
rrl = Ydr;
ssl = (Lda+A1*Nda)/(1-A1*B1);
ttl = (Ldr+A1*Ndr)/(1-A1*B1);
uul = (B1*Lda+Nda)/(1-A1*B1);
vvl = (B1*Ldr+Ndr)/(1-A1*B1);
wwl = 0;
xxl = 0;
yyl = 0;
zzl = 0;
Matrix_lat=[aal bbl ccl ddl ddl2; eel ffl ggl hhl hhl2; iil jjl kkl lll lll2; ...
    mml nnl ool ppl ppl2; fifth1 fifth2 fifth3 fifth4 fifth5];

Deflection_lat=[qql rrl; ssl ttl; uul vvl; wwl xxl; yyl zzl];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solve for the eigenvalues of the system response for the lateral      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vector,poleslat] = eig(Matrix_lat);
poleslat = eig(Matrix_lat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve for the aproximation of the eigenvalues for the lateral/direcional   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ductch Roll
w_n_d = sqrt(Nbeta + (1/u1)*(Ybeta*Nr - Nbeta*Yr))
xi_d = -(Nr + (Ybeta/u1))/(2*w_n_d)
% Poles Short Period
s1_d = - xi_d*w_n_d + i*(w_n_d*sqrt(1-xi_d^2))
s2_d = - xi_d*w_n_d - i*(w_n_d*sqrt(1-xi_d^2))
% Time to double or half
t_half_d = 0.693/(abs(xi_d)*w_n_d)
% Cycles to double or half
N_half_d = 0.110*w_n_d/(abs(-xi_d)*w_n_d)
% Logarithmic decrement
delta_d = -0.693/N_half_d

% Spiral Root
s3 = (Lbeta*Nr - Nbeta*Lr)/(Lbeta+Nbeta*A1)
Ts = - s3
% Spiral Root Stability Criterion
crit_spiral = (Lbeta*Nr - Nbeta*Lr)

% Rolling aproximation
s4 = Lp
Tr = -1/Lp


save Lateral.mat Matrix_lat Deflection_lat
