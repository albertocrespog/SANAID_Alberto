%--------------------- Atmósfera internacional----------------------------%

function [Temp,rho,p,a]=atmos_inter(h)                                             %Altura en metros
%--------------------------- Constantes ----------------------------------%
g = 9.80665;                                                               %Aceleración de la gravedad(m/s^2)
R = 287.058;                                                               %Constante universal de los gases(J/K*kg)
gamma = 1.4;                                                               %Coeficiente de dilatación adiabática(-)
lambda = 6.5*10^-3;                                                        %Gradiente térmico(K/m)
Temp_SL = 288.15;                                                          %Temperatura a nivel del mar(K)
rho_SL = 1.225;                                                            %Densidad a nivel del mar(kg/m^3)
p_SL = rho_SL*R*Temp_SL;                                                   %Presión a nivel del mar(Pa)
a_SL = sqrt(gamma*R*Temp_SL);                                              %Velocidad del sonido a nivel del mar(m/s)

for i = 1:length(h),
%--------------------------- Troposfera ----------------------------------%
if h(i) <= 11000,
Temp(i) = Temp_SL - lambda * h(i);                                                       
rho(i) = rho_SL *(Temp(i)/Temp_SL)^(g/(R*lambda)-1);                                   
p(i) = rho(i)*R*Temp(i);
a(i) = sqrt(gamma*R*Temp(i));
end
%--------------------------- Tropopausa ----------------------------------%
if h(i) > 11000 && h(i) <= 20000,
    Temp(i) = 216.65;
    rho(i) = 0.3639 * exp((-g/(R*Temp(i)))*(h(i)-11000));
    p(i) = rho(i)*R*Temp(i);
    a(i) = sqrt(gamma*R*Temp(i));
end
end
