function ALT = densityISA2ALT(density)
%densityISA2ALT devuelve la altitud (seg�n modelo ISA) teniendo commo
%entrada la densidad. V�lida para troposfera

%Par�metros:
g       = 9.80665;   %[m/s^2] Gravedad
R_g     = 287.05287; %[m^2/(s^2K)] Constante del gas
alpha_T = -0.0065;   %[K/m] Variaci�n de la temperatura con la altitud
Theta0  = 288.15;    %[K] Temeparatura a nivel del mar
p0      = 101325;    %[Pa] Presi�n a nivel del mar

ALT = Theta0/alpha_T*((density*R_g*(Theta0)/(p0)).^(1/(-g/(R_g*alpha_T)-1))-1); 