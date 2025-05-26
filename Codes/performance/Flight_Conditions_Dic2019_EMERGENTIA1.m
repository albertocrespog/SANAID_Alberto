%% Defines mission properties
function [Data_ATM V_Performance] = Flight_Conditions_2020(h,V,V_max)

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
V_Performance.h = h;
V_Performance.Mach = Mach;
V_Performance.Temp = Temp;
V_Performance.rho = rho;
V_Performance.p = p;
V_Performance.a = a;
V_Performance.V = V;
V_Performance.V_max = V_max;
V_Performance.q_inf = 0.5*rho*V^2;