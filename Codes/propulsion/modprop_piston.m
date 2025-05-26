function [T_T0,C,B,T]=modprop_piston(h,delta_P,delta_P0,M,M_0,B_SL,cp_P_SL)
% Recibe las entradas en Sistema Internacional y devuelve las salidas en
% las mismas unidades

% 1 ft = 0.3048 m
% 1 hp = 745.699872 W
% 1 lb/(hp*hr) = 0.6082773878418 * 10^-3 kg/(W*hr) = 1.68965941 * 10^-7 (kg/(W*s))
gamma = 1.4;

%------------- Atmósfera internacional-------------------------------------
[Temp,rho,p,a]=atmos_inter_mio(h);
[T_SL,rho_SL,p_SL,a_SL]=atmos_inter_mio(0);
theta=(Temp/T_SL); %Cociente de temperaturas 
delta_rho = rho/rho_SL;
%--------------Eficiencia propulsiva---------------------------------------
eta_inst=0.82;
if M<=0.1
    eta_p=eta_inst*M/0.1;
else
    eta_p=eta_inst;
end
if M_0<=0.1
    eta_p_0=eta_inst*M_0/0.1;
else
    eta_p_0=eta_inst;
end
%--------------POTENCIA----------------------------------------------------
B=delta_P*B_SL*((8.55*delta_rho - 1)/(7.55)); %[W]
%------------------EMPUJE--------------------------------------------------
V=a*M; %[m/s]
T=B*eta_p/V;%[N]
%---------CONSUMO DE COMBUSTIBLE POR UNIDAD DE PONTENCIA-------------------
C_bhp=cp_P_SL; %[kg/(W*s)]
%--------CORRECCIÓN DEL C_bhp CON LA POSICION DE PALANCA-------------------
a_1=3.559957437510763;
a_2=-10.739698199171459;
a_3= 11.989635150373475;
a_4=-5.869876557884609;
a_5=2.059994459180667;
C_bhp=C_bhp*(a_1*delta_T^4+a_2*delta_T^3+a_3*delta_T^2+a_4*delta_T+a_5);
C = (C_bhp * V/eta_p); %[kg/(N*s)]
%-------RELACIONES DE POTENCIA ENTRE SEGMENTOS----------------------------
P_3_P_0=(delta_P/delta_P0)*((8.55*delta_rho - 1)/(7.55)); 
T_T0=P_3_P_0*(a_SL*M_0/V)*(eta_p/eta_p_0);