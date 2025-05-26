function ALT = densityISA2ALT(density)
%densityISA2ALT devuelve la altitud (según modelo ISA) teniendo commo
%entrada la densidad. Válida para troposfera

%Parámetros:
g       = 9.80665;   %[m/s^2] Gravedad
R_g     = 287.05287; %[m^2/(s^2K)] Constante del gas
alpha_T = -0.0065;   %[K/m] Variación de la temperatura con la altitud
Theta0  = 288.15;    %[K] Temeparatura a nivel del mar
p0      = 101325;    %[Pa] Presión a nivel del mar

ALT = Theta0/alpha_T*((density*R_g*(Theta0)/(p0)).^(1/(-g/(R_g*alpha_T)-1))-1); 