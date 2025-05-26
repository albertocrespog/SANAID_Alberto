%% Defines mission properties
function [Data_ATM Performance] = Flight_Conditions_2020_v1(h,V)

%% Flight Conditons obtained using Standard Atmospheric Properties
Toffset = 0;
[rho,a,Temp,p,kvisc,ZorH]=stdatmo(h,Toffset,'SI');
% Atmospheric properties
%           rho:   Density            kg/m^3          slug/ft^3
%           a:     Speed of sound     m/s             ft/s
%           T:     Temperature        K              R
%           P:     Pressure           Pa              lbf/ft^2
%           nu:    Kinem. viscosity   m^2/s           ft^2/s
%           ZorH:  Height or altitude m               ft
Data_ATM.rho = rho;
Data_ATM.a = a;
Data_ATM.Temp = Temp;
Data_ATM.p = p;
Data_ATM.kvisc = kvisc;
Data_ATM.ZorH = ZorH;

Mach = V/a; % Mach number
% Saves data regarding flight conditions
Performance.h = h;
Performance.Mach = Mach;
Performance.Temp = Temp;
Performance.rho = rho;
Performance.p = p;
Performance.a = a;
Performance.V = V;
Performance.q_inf = 0.5*rho*V^2;