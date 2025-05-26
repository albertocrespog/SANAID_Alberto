%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m_end = combustible_analisis(m0,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,tf,h_CR,h_f)
rho_SL = densityISA(0);   %[kg/m^3] Densidad a nivel del mar
g      = 9.80665;         %[m/s^2] Gravedad
misops = odeset('RelTol',1e-9,'AbsTol',1e-9); %Tolerancias para el integrador
misops_opt = optimset('Display','off','TolX',1e-9,'TolFun',1e-18); %Tolerancias para el optimizador. Ponemos una 
%TolFun exageradamente pequeña para que el criterio de parada sea siempre TolX

%Comparamos las velocidades mínima y de máxima autonomía a h_CR
rho = densityISA(h_CR);
V_stall = sqrt(m0*g./(0.5*rho*C_L_max_CR*S));
V_min = 1.3*V_stall;
%Utilizamos como primera estimación a la velocidad de máxima autonomía la 
%que se obtiene suponiendo que el consumo específico y el rendimiento 
%propulsivo no dependen de la velocidad de vuelo:
V_E0 = sqrt(m0*g/.5/rho/S)*sqrt(-(C_D1_CR/C_D0_CR/6) + sqrt((C_D1_CR/C_D0_CR/6)^2+C_D2_CR/C_D0_CR/3)); 
V_E = fmincon(@(V) 1e4*consumo_crucero(V,h_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,m0,g),V_E0,[],[],[],[],V_stall,2*V_stall,[],misops_opt);
V_EAS = max(V_min,V_E)*sqrt(rho/rho_SL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ascendemos hasta la altitud de crucero a velocidad EAS contante
h_A = 0;
h_B = h_CR;
m = m0;

[h,y] = ode45(@(h,y) subida(h,y,V_EAS,S,C_D0_CR,C_D1_CR,C_D2_CR,g,rho_SL),[h_A h_B],m,misops);

mF_subida = y(1) - y(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crucero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h   = h_CR;
rho = densityISA(h);
m   = m - mF_subida;

V = V_EAS *sqrt(rho_SL/rho);
if tf == 0;
    t = 0;
    y = m;
else
    [t,y] = ode45(@(t,y) crucero(t,y,V,S,C_D0_CR,C_D1_CR,C_D2_CR,g,rho,h),[0 tf],m,misops);
end
mF_crucero = y(1,1) - y(end,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descenso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Descendemos hasta 100 m
h_A = h_CR;
h_B = h_f;
m = m - mF_crucero;

[h,y] = ode45(@(h,y) descenso(h,y,V_EAS,S,C_D0_CR,C_D1_CR,C_D2_CR,g,rho_SL),[h_A h_B],m,misops);

mF_descenso = y(1) - y(end);

m_end = m - mF_descenso;