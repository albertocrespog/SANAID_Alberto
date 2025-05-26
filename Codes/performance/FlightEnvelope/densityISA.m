function [rho,drho_dh] = densityISA(h)
%densityISA devuelve la densidad y su derivada con respecto a la altitud 
%calculadas de acuerdo al modelo ISA. Válida para troposfera

%Parámetros:
g       = 9.80665;   %[m/s^2] Gravedad
R_g     = 287.05287; %[m^2/(s^2K)] Constante del gas
alpha_T = -0.0065;   %[K/m] Variación de la temperatura con la altitud
Theta0  = 288.15;    %[K] Temperatura a nivel del mar
p0      = 101325;    %[Pa] Presión a nivel del mar

rho     = p0/(R_g*Theta0) * (1+alpha_T*h/Theta0).^(-g/(R_g*alpha_T)-1); 
drho_dh = p0/(R_g*Theta0) * (-g/(R_g*alpha_T)-1) * alpha_T/Theta0 *...
            (1+alpha_T*h/Theta0).^(-g/(R_g*alpha_T)-2); 