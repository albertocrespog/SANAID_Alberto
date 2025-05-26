%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estudio_Longitudinal.m   Longitudinal stability analysis for       % 
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
g=9.81                                           % gravity constant
h=5000;                                               % altitude (ft)
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
Sw=136;
cbar=5.4;
span=26.3;

% flight condition
u1=621;
q1barw=198;

% mass data
w=4000;     % lb
m=w/32.2;   % lb
Ixx=800;    %slug*ft^2
Iyy= 4800;  %slug*ft^2
Izz=5200;   %slug*ft^2
Ixz=200;      %slug*ft^2

cl1 = 0.307;
cd1= 0.032;
ctx1= cd1;
cm1=0;
cmt1=-cm1;   % no tail

% longitudinal coefficients
clu=0;
cdu=0;
cdalpha=.121;
ctxu=-3*ctx1;

% alpha derivatives
clalpha= 4.14;
clalphadot=1.7;
cmalpha = -.613;
cmalphadot=-7.27;

% q derivatives
clq = 3.9;
cmq = -12.4;

cmu=0;
cmtu=0;
cmtalpha=0;

% Longitudinal control and hinge moment derivatives %
cddeltae=0;
cldeltae=.43;
cmdeltae=-1.122;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Longitudinal Dimensional Stability Derivatives %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xu=(-q1barw*Sw*(cdu+(2*cd1)))/(m*u1);
Xtu=(q1barw*Sw*(ctxu+(2*ctx1)))/(m*u1);
Xalpha=(-q1barw*Sw*(cdalpha-cl1))/m;
Xdeltae=(-q1barw*Sw*cddeltae)/m;

Zu=(-q1barw*Sw*(clu+(2*cl1)))/(m*u1);
Zalpha=(-q1barw*Sw*(clalpha+cd1))/m;
Zalphadot=(-q1barw*Sw*clalphadot*cbar)/(2*m*u1);
Zq=(-q1barw*Sw*clq*cbar)/(2*m*u1);
Zdeltae=(-q1barw*Sw*cldeltae)/m;

Mu=(q1barw*Sw*cbar*(cmu+(2*cm1)))/(Iyy*u1);
Mtu=(q1barw*Sw*cbar*(cmtu+(2*cmt1)))/(Iyy*u1);
Malpha=(q1barw*Sw*cbar*cmalpha)/Iyy;
Malphadot=(q1barw*Sw*(cbar^2)*cmalphadot)/(2*Iyy*u1);
Mtalpha=(q1barw*Sw*cbar*cmtalpha)/Iyy;
Mq=(q1barw*Sw*(cbar^2)*cmq)/(2*Iyy*u1);
Mdeltae=(q1barw*Sw*cbar*cmdeltae)/Iyy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms of the Matrix that finds the poles of the system response %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aa = Xu + Xtu;
bb = Xalpha;
cc = 0;
dd = -g*cos(theta1);
ee = Zu/(u1 - Zalphadot);
ff = Zalpha/(u1 - Zalphadot);
gg = (Zq + u1)/(u1 - Zalphadot);
hh = (-g*sin(theta1))/(u1 - Zalphadot);
ii = Mu + Mtu + ((Zu*Malphadot)/(u1 - Zalphadot));
jj = Malpha + Mtalpha +((Zalpha*Malphadot)/(u1 - Zalphadot));
kk = Mq + (((Zq+u1)*Malphadot)/(u1 - Zalphadot));
ll = -(g*sin(theta1)*Malphadot)/(u1 - Zalphadot);
mm = 0;
nn = 0;
oo = 1;
pp = 0;
qq = Xdeltae;
rr = (Zdeltae/(u1 - Zalphadot));
ss = Mdeltae + ((Zdeltae*Malphadot)/(u1 - Zalphadot));
tt = 0;

Matrix_lon=[aa bb cc dd; ee ff gg hh; ii jj kk ll; mm nn oo pp];

Deflection_lon=[qq; rr; ss; tt];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solve for the eigenvalues of the system response for the longitudinal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vector,poleslon] = eig(Matrix_lon);
poleslon = eig(Matrix_lon)

omega_SP=sqrt(imag(poleslon(1,1))^2+real(poleslon(1,1))^2)
xi_SP=-real(poleslon(1,1))/omega_SP
T_SP=2*pi/(omega_SP*sqrt(1-xi_SP^2))
t_half_SP = 0.693/(abs(xi_SP)*omega_SP)
N_double_SP = 0.110*sqrt(1-xi_SP^2)/abs(xi_SP)
omega_SP_appox = sqrt((Zalpha*Mq/u1 - Malpha))
xi_SP_appox = -(Mq + Zalpha/u1 + Malphadot)/(2*omega_SP_appox)

omega_PH=sqrt(imag(poleslon(3,1))^2+real(poleslon(3,1))^2)
xi_PH=-real(poleslon(3,1))/omega_PH
T_PH=2*pi/(omega_PH*sqrt(1-xi_PH^2))
t_half_PH = 0.693/(abs(xi_PH)*omega_PH)
N_double_PH = 0.110*sqrt(1-xi_PH^2)/abs(xi_PH)
omega_PH_appox = sqrt(-g*Zu/u1)
xi_PH_appox = -Xu/(2*omega_PH_appox)

save Longitudinal.mat  Matrix_lon Deflection_lon AoA u1