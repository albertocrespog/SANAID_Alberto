function [T_T0,C,T]=modprop_fan(h,delta_T,delta_T0,M,M_0,T_SL,c_SL,derivacion)
% Recibe las entradas en Sistema Internacional y devuelve las salidas en
% las mismas unidades

% 1 ft = 0.3048 m
% 1 lbf = 4.448221615255 N
% 1 lb/(lbf*hr) = 2.832546065 * 10^-5 kg/(N*s) = 0.1019716583 kg/(N*h)
gamma = 1.4;
%------------- Atmósfera internacional-------------------------------------
[Temp,rho,p,a]=atmos_inter_mio(h);
[Temp_SL,rho_SL,p_SL,a_SL]=atmos_inter_mio(0);
theta=(Temp/Temp_SL); 
delta_rho = rho/rho_SL;
%------------------EMPUJE--------------------------------------------------
T=delta_T*T_SL*((1+0.2*M^2)^(gamma./(gamma-1)))*delta_rho*(1-0.49*sqrt(M)); %[N]
%--------CORRECCIÓN DEL C_bhp CON LA POSICION DE PALANCA-------------------
a_1=3.559957437510763;
a_2=-10.739698199171459;
a_3= 11.989635150373475;
a_4=-5.869876557884609;
a_5=2.059994459180667;
c_SL=c_SL*(a_1*delta_T^4+a_2*delta_T^3+a_3*delta_T^2+a_4*delta_T+a_5); %[kg/(N*s)]
%---------CONSUMO DE COMBUSTIBLE POR UNIDAD DE PONTENCIA-------------------
switch derivacion
    case 1
        C = c_SL*(1+1.2*M)*sqrt(theta);
    case 2
        C = c_SL*(1+0.33*M)*sqrt(theta);
    case 3
        C = c_SL*(1+0.16875*M)*sqrt(theta);
end
%-------RELACIONES DE POTENCIA ENTRE SEGMENTOS----------------------------
T_T0=(delta_T/delta_T0)*(((1+0.2*M^2)/(1+0.2*M_0^2))^(gamma./(gamma-1)))*delta_rho*((1-0.49*sqrt(M))/(1-0.49*sqrt(M_0))); 
