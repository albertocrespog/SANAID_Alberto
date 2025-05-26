function [rho_f,rho_fairing,rho_nose,rho_w,rho_HTP,rho_VTP,rho_tb]=CFIdensities

% Datos Céfiro I
S_w = 1.143;
S_h = 0.176;
S_v = 0.173;

% Volume Nose cone. approximated as cilinder(v3), truncated cone(v2), and a pyramid
% cone v1
r1 = 7.5/100;
h1 = 10/100;
r2 = 12.8/100;
h2 = 15/100;
r3 = 12.8/100;
h3 = 56.6/100;
% volume of cylinder
v3 = pi*r3^2*h3;
% Area
A_v3 = 2*pi*r3*h3;
% volume of truncated cylinder
v2 = (h2/3)*(pi*r1^2 + pi*r2^2 + sqrt(pi*(r1^2)*pi*(r2^2)));
% Area
s_2 = sqrt((r2-r1)^2 + h2^2); % generatriz
A_v2 = pi*(r1+r2)*s_2;
% volume of pyramid cone
v1 = (1/3)*pi*(r1^2)*h1;
% Area
s_1 = sqrt((r1)^2 + (h1^2)); % generatriz
A_v1 = pi*r1*s_1;
V_nose = v1 + v2 + v3;
A_nose = A_v1 + A_v2 + A_v3;

% Volume and area of center fuselage
% divided in a rectangular section with 
w_fuselage = 25/100; % width of fuselage
h_fuselage = (25/2)/100; % height of rectangular section
r_fuselage = (25/2)/100;
l_fuselage = 70/100;

v_fuselage_rec = w_fuselage*h_fuselage*l_fuselage;
A_fuselage_rec = w_fuselage*l_fuselage + 2*h_fuselage*l_fuselage;
v_fuselage_cyl = (pi*r_fuselage^2*l_fuselage)/2; % only half cylinder
A_fuselage_cyl = 2*pi*r_fuselage*l_fuselage/2; % only half cylinder

A_fuselage = A_fuselage_rec + A_fuselage_cyl;
v_fuselage = v_fuselage_rec + v_fuselage_cyl;

W_S_w = 5.022;
W_nose = 1.414;
W_fuse = 3.324;
W_S_hv = 1.523;
W_tb = 1.468;

rho_w = W_S_w/S_w;
rho_HTP = W_S_hv/(S_h+S_v);
rho_VTP = W_S_hv/(S_h+S_v);
rho_nose = W_nose/A_nose;
rho_f = W_fuse/A_fuselage;
% 
% S_w_e = 2.050131939508261;
% b_w = 5.225930044268314;
% b_w_e = 4.869930044268314;
% AR_w = 12.413793103448285;
% AR_w_e = 11.568142605375774;
% cR_w = 0.420977698010503;
% cT_w = 0.420977698010503;
% 
% cR_h = 0.259554525531993;
% cT_h = 0.259554525531993;
% b_h = 1.306482511067079;
% S_h = 0.339103448275862;
% AR_h = 5.033557046979866;
% 
% cR_v = 0.435494170355693;
% cT_v = 0.261296502213416;
% b_v = 0.478172599050551;
% S_v = 0.333186206896552;
% AR_v = 1.372500000000000;
% 
% Length Tail boom of metal
l_tb = 139/100;
% 
% % Fuselaje
% Surf_TOT = 2.700430097266073; % m^2
% Vol_TOT = 0.244461563000000; % m^3
% length_fus = 2.552000000000000; % m
% w_Area_b_max = 0.356000000000000;
% h_Area_b_max = 0.372000000000000;
% l_cajon = 1.750000000000000;
% Vol_cajon = l_cajon*w_Area_b_max*h_Area_b_max;
% 
% % Pesos
% W_S_w_e = 9.752;
% W_fus = 5.206;
% W_S_hv = 3.458;
% W_tb = 2.884;
% W_fairing = 4.149;
% 
% % Recalculo de las nuevas densidades
% rho_f = W_fus/Surf_TOT;
% rho_fairing = W_fairing/Surf_TOT;
% rho_w = W_S_w_e/S_w_e;
% rho_HTP = W_S_hv/(S_h+S_v);
% rho_VTP = W_S_hv/(S_h+S_v);
rho_tb = W_tb/(2*l_tb);

% Data fron Céfiro III for carbon fiber composite fairing 
Surf_TOT = 2.700430097266073; % m^2
W_fairing = 4.149;
W_fus = 5.206;
% Recalculo de las nuevas densidades
rho_fairing = W_fairing/Surf_TOT;

end
