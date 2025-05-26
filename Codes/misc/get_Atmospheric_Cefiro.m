function Data_ATM = get_Atmospheric_Cefiro(h)

g_ATM=-9.81;                             %gravedad de la tierra 
R_ATM=287;                               %constante de los gases
T0_ATM=288;                                   %temperatura en Kelvin
  
if h<11000                           %densidad en la troposfera        
                                     %altura de la tropopausa
    alpha_ATM = -6.5e-3;                   %baja la temperatura 6,5º cada km                               
    T0_H      = T0_ATM + alpha_ATM*h;
    a         = sqrt(1.4*R_ATM*T0_H);
    rho       = 1.225*((T0_ATM + alpha_ATM*h)/T0_ATM)^(g_ATM/(R_ATM*alpha_ATM)-1);
                                     %número de Mach
else                                 %densidad en la estratosfera
    deltaestrato=exp(g_ATM*(h_ATM-11000)/(R_ATM*216.5));
    rho=0.3629*deltaestrato;
end

Data_ATM.rho = rho;
Data_ATM.a_speed = a;