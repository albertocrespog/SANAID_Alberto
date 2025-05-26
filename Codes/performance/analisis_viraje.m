function [fuel_viraje,tiempo_viraje,distancia_viraje,datos] = analisis_viraje(propul,aerodinamica,viraje,W_inicial,h_inicial,opcion,datos,i,Prop_data,data_electric)

% PROPUL
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR
% 4: CONSUMO ESPECIFICO
% 5: DERIVACION(TURBOFANES)
% 6: CAPACIDAD CALORIFICA DEL COMBUSTIBLE
% 7: DIAMETRO DEL MOTOR(TURBOHELICES)
% 8: EFICIENCIA DE LA HELICE (ETA_P) 
% 9: AVION CIVIL = 1/MILITAR = 2

% AERODINAMICA
% 1: SUPERFICIE
% 2: CD0
% 3: K1 = K
% 4: K2
% 5: CLMAX LIMPIO

% VIRAJE
% 1: TIEMPO FINAL       
% 2: MACH DE VUELO
% 3: PALANCA DE GASES
% 4: CL DE VIRAJE
% 5: ANGULO DE ALABEO
% 6: VELOCIDAD DE GUIÑADA
% 7: FACTOR DE CARGA
% 8: RADIO DE GIRO 

% AERODINAMICA
S = aerodinamica(1);
Cd0 = aerodinamica(2);
k1 = aerodinamica(3);
k2 = aerodinamica(4);
CLmax_limpio = aerodinamica(5);

% PROPULSION
tipo_motor = propul(1);
num_motores = propul(2);
T_SL = propul(3);
c_SL = propul(4);
civil_mil = propul(5);
eta_inst = propul(6);
derivacion = propul(7);

% ELECTRIC PROPULSION
e0 = data_electric(1)*3600; %from Wh to J
Dprop = data_electric(2); % in m
SF_prop = data_electric(3);
nps_max0 = data_electric(4);
nps_max = nps_max0*SF_prop;
eta_m = data_electric(7); % Will be replaced by excel
CT_Polyfit = Prop_data.CT_Polyfit;
CP_Polyfit = Prop_data.CP_Polyfit;
etamp_Polyfit = Prop_data.etamp_Polyfit;
tau = data_electric(8)/100;

% VIRAJE 
tiempo_final = viraje(1);
%Mach = viraje(2); Al final lo que se introduce es la velocidad
delta_T_gases = viraje(3);
C_l = viraje(4);
alabeo = viraje(5);
veloc_guiniada = viraje(6);
factor_carga = viraje(7);
radio_giro = viraje(8);


% ATMOSFERA ESTANDAR
[Temp_SL,rho_SL,p_SL,a_SL] = atmos_inter_mio(0);
[Temp,rho,p,a] = atmos_inter_mio(h_inicial);
Mach = viraje(2)/a;

tiempo_viraje = tiempo_final;
distancia_viraje = 0;

% ------------------------  CONSTANTES --------------------------
g = 9.80665;                                                                                     

W_ini = W_inicial;
h_ini = h_inicial;
x_ini = 0;

iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();   %ESTA LINEA PERMITE USAR CONDICIONES EN LAS FUNCIONES!!!

eta_p = eta_inst;

T_SL = T_SL * num_motores;

%-------------------------------------------------------------------------%
%ECUACIONES DE LA MECANICA DEL VUELO PARA EL VIRAJE

% V = CONSTANTE, H = CONSTANTE = DATO, VIRAJE HORIZONTAL --> GAMMA = 0

%           T(delta_T,V)         = D(V,CL) 
%           L(V,CL)*cos(mu)      = W
%           L(V,CL)*sen(mu)      = (W/g)V dXI/dt
%           dW/dt                = -C(V)*g*T(delta_T,V)

%-------------------------------------------------------------------------%

switch opcion
    
    case 1  % VIRAJE HORIZONTAL CON VELOCIDAD DADA Y PALANCA DE GASES
 
delta_T = delta_T_gases;
palanca = delta_T;
M = Mach;
V = M*a;

a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);

switch tipo_motor
    case 1
         T = delta_T*(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
        switch derivacion
            case 1
                C =  c_SL*correccion*(1+1.2*M)*sqrt(Temp/Temp_SL);
            case 2
                C =  c_SL*correccion*(1+0.33*M)*sqrt(Temp/Temp_SL);
            case 3
                C =  c_SL*correccion*(1+0.16875*M)*sqrt(Temp/Temp_SL);
        end
    case 2
        T = delta_T*(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
        C = c_SL*correccion*(1+1.44*M)*sqrt(Temp/Temp_SL)* V/eta_p;
    case 3
        T = delta_T*(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
        C = c_SL *correccion*V/eta_p;
end    

D = T;
CD = D/(0.5*rho*V^2*S);
CL = (k2 + sqrt(k2^2-4*k1*(Cd0-CD)))/(2*k1);
L = 0.5*rho*V^2*S*CL;

%Integrando la ecuacion de la variacion del peso
W   = @(t) W_ini - C*g*T.*t;

muu = @(t) acos(W(t)./L);  %2º ecuacion
N   = @(t) 1./cos(muu(t));

chi_punto = @(t) (L.*sin(muu(t))*g)./(W(t).*V);
R         = @(t)  V./chi_punto(t);

num = 10000;
TT = linspace(0,tiempo_final,num);

for k = 1:num
    peso(k) = W(TT(k));
    L_D(k) = CL/CD;
    empuje(k) = T;
    CL2(k) = CL;
    alabeo(k)    = muu(TT(k));
    n(k)           = N(TT(k));
    veloc_chi(k)   = chi_punto(TT(k));
    radio_giro(k)  = R(TT(k));
end

fuel_viraje = (W(0)-W(tiempo_final))/g;
vueltas = mean(veloc_chi)*tiempo_final/(2*pi);
W = W(tiempo_final);

datos(i).segmento.tiempo = linspace(0,tiempo_final,num);
datos(i).segmento.fuel = (peso(1) - peso)./g; %fuel consumido acumulado
datos(i).segmento.vueltas = vueltas; %escalar
datos(i).segmento.palanca = delta_T;
% datos(i).segmento.empuje = empuje;
datos(i).segmento.empuje = CL2;
datos(i).segmento.L_D = L_D;
datos(i).segmento.peso = peso;
datos(i).segmento.velocidad = V;
datos(i).segmento.alabeo = alabeo;
datos(i).segmento.factor_de_carga = n;
datos(i).segmento.velocidad_guiniada = veloc_chi;
datos(i).segmento.radio_giro = radio_giro;

%-------------------------------------------------------------------------%

    case 2 % VELOCIDAD DADA, CL DADO
        
CL = C_l;
M = Mach;
V = M*a;

CD = Cd0 -k2*CL + k1*CL^2;
D = 0.5*rho*V^2*S*CD;
L = 0.5*rho*V^2*S*CL;
T = D;

a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;

switch tipo_motor
    case 1
         delta_T = T/(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
         correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        switch derivacion
            case 1
                C =  c_SL*correccion*(1+1.2*M)*sqrt(Temp/Temp_SL);
            case 2
                C =  c_SL*correccion*(1+0.33*M)*sqrt(Temp/Temp_SL);
            case 3
                C =  c_SL*correccion*(1+0.16875*M)*sqrt(Temp/Temp_SL);
        end
    case 2
        delta_T = T/(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        C = c_SL*correccion*(1+1.44*M)*sqrt(Temp/Temp_SL)* V/eta_p;
    case 3
        delta_T = T/(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        C = c_SL *correccion*V/eta_p;
end    

palanca = delta_T;
%Integrando la ecuacion de la variacion del peso
W   = @(t) W_ini - C*g*T.*t;

muu = @(t) acos(W(t)./L);  %2º ecuacion
N   = @(t) 1./cos(muu(t));

chi_punto = @(t) (L.*sin(muu(t))*g)./(W(t).*V);
R         = @(t)  V./chi_punto(t);

num = 10000;
TT = linspace(0,tiempo_final,num);

for k = 1:num
    peso(k) = W(TT(k));
    L_D(k) = CL/CD;
    CL2(k) = CL;
    empuje(k) = T;
    alabeo(k)    = muu(TT(k));
    n(k)           = N(TT(k));
    veloc_chi(k)   = chi_punto(TT(k));
    radio_giro(k)  = R(TT(k));
end

fuel_viraje = (W(0)-W(tiempo_final))/g;
vueltas = mean(veloc_chi)*tiempo_final/(2*pi);
W = W(tiempo_final);

datos(i).segmento.tiempo = linspace(0,tiempo_final,num);
datos(i).segmento.fuel = (peso(1) - peso)./g; %fuel consumido acumulado
datos(i).segmento.vueltas = vueltas; %escalar
datos(i).segmento.palanca = delta_T;
% datos(i).segmento.empuje = empuje;
datos(i).segmento.empuje = CL2;
datos(i).segmento.L_D = L_D;
datos(i).segmento.peso = peso;
datos(i).segmento.velocidad = V;
datos(i).segmento.alabeo = alabeo;
datos(i).segmento.factor_de_carga = n;
datos(i).segmento.velocidad_guiniada = veloc_chi;
datos(i).segmento.radio_giro = radio_giro;


%-------------------------------------------------------------------------%


    case 3 % VELOCIDAD DADA, MU DADO

muu = alabeo;        
M = Mach;
V = M*a;

chi_punto = g*tan(muu)/V;

CL = @(W) W./(0.5*rho*V^2*S*cos(muu));
CD = @(W) Cd0 - k2.*CL(W) + k1.*CL(W).^2;
D  = @(W) 0.5*rho*V^2*S.*CD(W);
T  = @(W) D(W);

a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;

switch tipo_motor
    case 1
        delta_T = @(W) T(W)./(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        switch derivacion
            case 1
                C =  @(W) c_SL.*correccion(W).*(1+1.2*M).*sqrt(Temp/Temp_SL);
            case 2
                C =  @(W) c_SL.*correccion(W).*(1+0.33*M).*sqrt(Temp/Temp_SL);
            case 3
                C =  @(W) c_SL.*correccion(W).*(1+0.16875*M).*sqrt(Temp/Temp_SL);
        end
    case 2
        delta_T = @(W) T(W)./(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        C = @(W) c_SL.*correccion(W).*(1+1.44*M).*sqrt(Temp/Temp_SL).* V./eta_p;
    case 3
        delta_T = @(W) T(W)./(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        C = @(W) c_SL.*correccion(W).*V./eta_p;
end    

% Integramos  dW/dt = -C(W)gT(W)

num = 10000;
TT_v = linspace(0,tiempo_final,num);

dwdt = @(t,W) -C(W).*g.*T(W);
[TT,PESO] = ode45(dwdt,TT_v,W_ini);

R = V^2/(g*tan(muu));
N = 1/cos(muu);

for k = 1:num
    peso(k) = PESO(k);
    L_D(k) = CL(PESO(k))/CD(PESO(k));
    CL2(k) = CL(PESO(k));
    empuje(k) = T(PESO(k));
    palanca(k) = delta_T(PESO(k));
    n(k)           = N;
    veloc_chi(k)   = chi_punto;
    radio_giro(k)  = R;
end

fuel_viraje = (PESO(1)-PESO(end))/g;
W = PESO(end);
vueltas = chi_punto*tiempo_final/(2*pi);

datos(i).segmento.tiempo = TT_v;
datos(i).segmento.fuel = (peso(1) - peso)./g; %fuel consumido acumulado
datos(i).segmento.vueltas = vueltas; %escalar
datos(i).segmento.palanca = palanca;
%datos(i).segmento.empuje = empuje;
datos(i).segmento.empuje = CL2;
datos(i).segmento.L_D = L_D;
datos(i).segmento.peso = peso;
datos(i).segmento.velocidad = V;
datos(i).segmento.alabeo = alabeo;
datos(i).segmento.factor_de_carga = n;
datos(i).segmento.velocidad_guiniada = veloc_chi;
datos(i).segmento.radio_giro = radio_giro;

%-------------------------------------------------------------------------%

    case 4  %VELOCIDAD DADA, FACTOR DE CARGA DADO
        
N = factor_carga;        
M = Mach;
V = M*a;

muu = acos(1/N);

chi_punto = g*tan(muu)/V;

CL = @(W) W./(0.5*rho*V^2*S*cos(muu));
CD = @(W) Cd0 - k2.*CL(W) + k1.*CL(W).^2;
D  = @(W) 0.5*rho*V^2*S.*CD(W);
T  = @(W) D(W);

a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;

switch tipo_motor
    case 1
        delta_T = @(W) T(W)./(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        switch derivacion
            case 1
                C =  @(W) c_SL.*correccion(W).*(1+1.2*M).*sqrt(Temp/Temp_SL);
            case 2
                C =  @(W) c_SL.*correccion(W).*(1+0.33*M).*sqrt(Temp/Temp_SL);
            case 3
                C =  @(W) c_SL.*correccion(W).*(1+0.16875*M).*sqrt(Temp/Temp_SL);
        end
    case 2
        delta_T = @(W) T(W)./(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        C = @(W) c_SL.*correccion(W).*(1+1.44*M).*sqrt(Temp/Temp_SL).* V./eta_p;
    case 3
        delta_T = @(W) T(W)./(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        C = @(W) c_SL.*correccion(W).*V./eta_p;
end    

% Integramos  dW/dt = -C(W)gT(W)

num = 10000;
TT_v = linspace(0,tiempo_final,num);

dwdt = @(t,W) -C(W).*g.*T(W);
[TT,PESO] = ode45(dwdt,TT_v,W_ini);

R = V^2/(g*tan(muu));

for k = 1:num
    peso(k) = PESO(k);
    L_D(k) = CL(PESO(k))/CD(PESO(k));
    CL2(k) = CL(PESO(k));
    empuje(k) = T(PESO(k));
    palanca(k) = delta_T(PESO(k));
    n(k)           = N;
    alabeo(k)      = muu;
    veloc_chi(k)   = chi_punto;
    radio_giro(k)  = R;
end

fuel_viraje = (PESO(1)-PESO(end))/g;
W = PESO(end);
vueltas = chi_punto*tiempo_final/(2*pi);

datos(i).segmento.tiempo = TT_v;
datos(i).segmento.fuel = (peso(1) - peso)./g; %fuel consumido acumulado
datos(i).segmento.vueltas = vueltas; %escalar
datos(i).segmento.palanca = palanca;
% datos(i).segmento.empuje = empuje;
datos(i).segmento.empuje = CL2;
datos(i).segmento.L_D = L_D;
datos(i).segmento.peso = peso;
datos(i).segmento.velocidad = V;
datos(i).segmento.alabeo = alabeo;
datos(i).segmento.factor_de_carga = n;
datos(i).segmento.velocidad_guiniada = veloc_chi;
datos(i).segmento.radio_giro = radio_giro;

%-------------------------------------------------------------------------%

    case 5      %VELOCIDAD DADA, RADIO DE GIRO DADO
        
R = radio_giro;        
M = Mach;
V = M*a;

muu = atan(V^2/(R*g));

chi_punto = g*tan(muu)/V;

CL = @(W) W./(0.5*rho*V^2*S*cos(muu));
CD = @(W) Cd0 - k2.*CL(W) + k1.*CL(W).^2;
D  = @(W) 0.5*rho*V^2*S.*CD(W);
T  = @(W) D(W);

a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;

switch tipo_motor
    case 1
        delta_T = @(W) T(W)./(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        switch derivacion
            case 1
                C =  @(W) c_SL.*correccion(W).*(1+1.2*M).*sqrt(Temp/Temp_SL);
            case 2
                C =  @(W) c_SL.*correccion(W).*(1+0.33*M).*sqrt(Temp/Temp_SL);
            case 3
                C =  @(W) c_SL.*correccion(W).*(1+0.16875*M).*sqrt(Temp/Temp_SL);
        end
    case 2
        delta_T = @(W) T(W)./(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        C = @(W) c_SL.*correccion(W).*(1+1.44*M).*sqrt(Temp/Temp_SL).* V./eta_p;
    case 3
        delta_T = @(W) T(W)./(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        C = @(W) c_SL.*correccion(W).*V./eta_p;
end    

% Integramos  dW/dt = -C(W)gT(W)

num = 10000;
TT_v = linspace(0,tiempo_final,num);

dwdt = @(t,W) -C(W).*g.*T(W);
[TT,PESO] = ode45(dwdt,TT_v,W_ini);

N = 1/cos(muu);

for k = 1:num
    peso(k) = PESO(k);
    L_D(k) = CL(PESO(k))/CD(PESO(k));
    CL2(k) = CL(PESO(k));
    empuje(k) = T(PESO(k));
    palanca(k) = delta_T(PESO(k));
    alabeo(k)      = muu;
    n(k)           = N;
    veloc_chi(k)   = chi_punto;
    radio_giro(k)  = R;
end

fuel_viraje = (PESO(1)-PESO(end))/g;
W = PESO(end);
vueltas = chi_punto*tiempo_final/(2*pi);

datos(i).segmento.tiempo = TT_v;
datos(i).segmento.fuel = (peso(1) - peso)./g; %fuel consumido acumulado
datos(i).segmento.vueltas = vueltas; %escalar
datos(i).segmento.palanca = palanca;
% datos(i).segmento.empuje = empuje;
datos(i).segmento.empuje = CL2;
datos(i).segmento.L_D = L_D;
datos(i).segmento.peso = peso;
datos(i).segmento.velocidad = V;
datos(i).segmento.alabeo = alabeo;
datos(i).segmento.factor_de_carga = n;
datos(i).segmento.velocidad_guiniada = veloc_chi;
datos(i).segmento.radio_giro = radio_giro;

        
%-------------------------------------------------------------------------%

    case 6   %VELOCIDAD DADA, VELOCIDAD DE GUINIADA DADA
        
chi_punto = veloc_guiniada;        
M = Mach;
V = M*a;

muu = atan(V*chi_punto/g);

CL = @(W) W./(0.5*rho*V^2*S*cos(muu));
CD = @(W) Cd0 - k2.*CL(W) + k1.*CL(W).^2;
D  = @(W) 0.5*rho*V^2*S.*CD(W);
T  = @(W) D(W);

a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;

switch tipo_motor
    case 1
        delta_T = @(W) T(W)./(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        switch derivacion
            case 1
                C =  @(W) c_SL.*correccion(W).*(1+1.2*M).*sqrt(Temp/Temp_SL);
            case 2
                C =  @(W) c_SL.*correccion(W).*(1+0.33*M).*sqrt(Temp/Temp_SL);
            case 3
                C =  @(W) c_SL.*correccion(W).*(1+0.16875*M).*sqrt(Temp/Temp_SL);
        end
    case 2
        delta_T = @(W) T(W)./(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        C = @(W) c_SL.*correccion(W).*(1+1.44*M).*sqrt(Temp/Temp_SL).* V./eta_p;
    case 3
        delta_T = @(W) T(W)./(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
        correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
        C = @(W) c_SL.*correccion(W).*V./eta_p;
end    

% Integramos  dW/dt = -C(W)gT(W)

num = 10000;
TT_v = linspace(0,tiempo_final,num);

dwdt = @(t,W) -C(W).*g.*T(W);
[TT,PESO] = ode45(dwdt,TT_v,W_ini);

R = V^2/(g*tan(muu));
N = 1/cos(muu);

for k = 1:num
    peso(k) = PESO(k);
    L_D(k) = CL(PESO(k))/CD(PESO(k));
    CL2(k) = CL(PESO(k));
    empuje(k) = T(PESO(k));
    palanca(k) = delta_T(PESO(k));
    n(k)           = N;
    alabeo(k)      = muu;
    veloc_chi(k)   = chi_punto;
    radio_giro(k)  = R;
end

fuel_viraje = (PESO(1)-PESO(end))/g;
W = PESO(end);
vueltas = chi_punto*tiempo_final/(2*pi);

datos(i).segmento.tiempo = TT_v;
datos(i).segmento.fuel = (peso(1) - peso)./g; %fuel consumido acumulado
datos(i).segmento.vueltas = vueltas; %escalar
datos(i).segmento.palanca = palanca;
% datos(i).segmento.empuje = empuje;
datos(i).segmento.empuje = CL2;
datos(i).segmento.L_D = L_D;
datos(i).segmento.peso = peso;
datos(i).segmento.velocidad = V;
datos(i).segmento.alabeo = alabeo;
datos(i).segmento.factor_de_carga = n;
datos(i).segmento.velocidad_guiniada = veloc_chi;
datos(i).segmento.radio_giro = radio_giro;

        
%-------------------------------------------------------------------------%

    case 7 %PALANCA DE GASES DADA, HALLAR V TAL QUE N = N_MAX (EQUIVALENTE A MUU = MUU_MAX)
        
delta_T = delta_T_gases;
palanca = delta_T;

M_v = @(V) V./a;

a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);

switch tipo_motor
    case 1
        T_v = @(V) delta_T.*(T_SL*(1+0.2.*M_v(V).^2).^(1.4./0.4)*(1-0.49.*sqrt(M_v(V))).*(rho./rho_SL));
        switch derivacion
            case 1
                C_v = @(V) c_SL.*correccion.*(1+1.2.*M_v(V)).*sqrt(Temp/Temp_SL);
            case 2
                C_v = @(V) c_SL.*correccion.*(1+0.33.*M_v(V)).*sqrt(Temp/Temp_SL);
            case 3
                C_v = @(V) c_SL.*correccion.*(1+0.16875.*M_v(V)).*sqrt(Temp/Temp_SL);
        end
    case 2
        T_v = @(V) delta_T.*(T_SL.*eta_p./V.*(1+0.2.*M_v(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
        C_v = @(V) c_SL.*correccion.*(1+1.44.*M_v(V)).*sqrt(Temp/Temp_SL).* V./eta_p;
    case 3
        T_v = @(V) delta_T.*(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
        C_v = @(V) c_SL .*correccion.*V./eta_p;
end   

% El objetivo es despejar n = n(V) y derivar e igualar a 0 para obtener la 
% velocidad que hace el factor de carga máximo, iterando para cada peso

% D = 0.5*rho*V^2*S*Cd0 - k2*W*n + k1*n^2*W^2/(0.5*rho*V^2*S) =
% D = A_2 * V^2         - B_2*W*n + C_2*n^2*W^2/V^2

A_2 = 0.5*rho*S*Cd0;
B_2 = k2;
C_2 = k1/(0.5*rho*S);

N = @(V) (B_2.*V.^2 + V.*sqrt(B_2^2.*V.^2 + 4*C_2.*(T_v(V)-A_2.*V.^2)))./(2*C_2*W_ini);
syms vel positive
dNdv = diff(N(vel));

v = double(solve(dNdv));

% A partir de aqui estamos en el caso 1 
% La variacion de la velocidad con el factor de carga maxima es infima
% (del orden de 0.001 m/s)

V = v;
M = V/a;

switch tipo_motor
    case 1
         T = delta_T*(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
        switch derivacion
            case 1
                C =  c_SL*correccion.*(1+1.2*M)*sqrt(Temp/Temp_SL);
            case 2
                C =  c_SL*correccion.*(1+0.33*M)*sqrt(Temp/Temp_SL);
            case 3
                C =  c_SL*correccion.*(1+0.16875*M)*sqrt(Temp/Temp_SL);
        end
    case 2
        T = delta_T*(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
        C = c_SL*correccion.*(1+1.44*M)*sqrt(Temp/Temp_SL)* V/eta_p;
    case 3
        T = delta_T*(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
        C = c_SL *correccion.*V/eta_p;
end    

D = T;
CD = D/(0.5*rho*V^2*S);
CL = (k2 + sqrt(k2^2-4*k1*(Cd0-CD)))/(2*k1);
L = 0.5*rho*V^2*S*CL;

%Integrando la ecuacion de la variacion del peso
W   = @(t) W_ini - C*g*T.*t;

muu = @(t) acos(W(t)./L);  %2º ecuacion
N   = @(t) 1./cos(muu(t));

chi_punto = @(t) (L.*sin(muu(t))*g)./(W(t).*V);
R         = @(t)  V./chi_punto(t);

num = 1000;
TT = linspace(0,tiempo_final,num);

for k = 1:num
    peso(k) = W(TT(k));
    L_D(k) = CL/CD;
    CL2(k) = CL;
    empuje(k) = T;
    alabeo(k)    = muu(TT(k));
    n(k)           = N(TT(k));
    veloc_chi(k)   = chi_punto(TT(k));
    radio_giro(k)  = R(TT(k));
end

fuel_viraje = (W(0)-W(tiempo_final))/g;
vueltas = mean(veloc_chi)*tiempo_final/(2*pi);
W = W(tiempo_final);

datos(i).segmento.tiempo = TT;
datos(i).segmento.fuel = (peso(1) - peso)./g; %fuel consumido acumulado
datos(i).segmento.vueltas = vueltas; %escalar
datos(i).segmento.palanca = delta_T;
% datos(i).segmento.empuje = empuje;
datos(i).segmento.empuje = CL2;
datos(i).segmento.L_D = L_D;
datos(i).segmento.peso = peso;
datos(i).segmento.velocidad = v;
datos(i).segmento.alabeo = alabeo;
datos(i).segmento.factor_de_carga = n;
datos(i).segmento.velocidad_guiniada = veloc_chi;
datos(i).segmento.radio_giro = radio_giro;
    
%-------------------------------------------------------------------------%    
    
    case 8  %PALANCA DE GASES DADA, HALLAR V TAL QUE XI PUNTO = XI PUNTO_MAX

delta_T = delta_T_gases;

M = @(V) V./a;

a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);

switch tipo_motor
    case 1
        T = @(V) delta_T.*(T_SL*(1+0.2.*M(V).^2).^(1.4./0.4)*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
        switch derivacion
            case 1
                C = @(V) c_SL.*correccion.*(1+1.2.*M(V)).*sqrt(Temp/Temp_SL);
            case 2
                C = @(V) c_SL.*correccion.*(1+0.33.*M(V)).*sqrt(Temp/Temp_SL);
            case 3
                C = @(V) c_SL.*correccion.*(1+0.16875.*M(V)).*sqrt(Temp/Temp_SL);
        end
    case 2
        T = @(V) delta_T.*(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
        C = @(V) c_SL.*correccion.*(1+1.44.*M(V)).*sqrt(Temp/Temp_SL).* V./eta_p;
    case 3
        T = @(V) delta_T.*(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
        C = @(V) c_SL .*correccion.*V./eta_p;
end   

D = @(V) T(V);
CD = @(V) D(V)./(0.5.*rho.*V.^2*S);
L  = @(n,W) n.*W;
CL = @(n,V,W) L(n,W)./(0.5.*rho.*V.^2.*S);
muu = @(n) acos(1./n);
R = @(n,V) V.^2./(g.*sqrt(n.^2-1));

% El objetivo es despejar n = n(V) y derivar e igualar a 0 para obtener la 
% velocidad que hace el factor de carga máximo, iterando para cada peso

% D = 0.5*rho*V^2*S*Cd0 - k2*W*n + k1*n^2*W^2/(0.5*rho*V^2*S) =
% D = A_2 * V^2         - B_2*W*n + C_2*n^2*W^2/V^2

A_2 = 0.5*rho*S*Cd0;
B_2 = k2;
C_2 = k1/(0.5*rho*S);

num = 1000; %Numero de intervalos para discretizar el tiempo
num2 = 5;

TT = linspace(0,tiempo_final,num);
w(1) = W_ini;
k = 1;
epsilon = 10;
    N = @(V) (B_2.*V.^2 + V.*sqrt(B_2^2.*V.^2 + 4*C_2.*(T(V)-A_2.*V.^2)))./(2*C_2*w(k));
    chi_punto = @(V) g.*sqrt(N(V).^2-1)./V;    
    syms vel positive
    %dNdv   = diff(N(vel));
    dChidv = diff(chi_punto(vel));
    dChidv = matlabFunction(dChidv);
    d2Chidv = diff(dChidv(vel));
    %v(k) = double(solve(dChidv));
    [v(k),func,error_raph(k)] = Newton_Raphson(dChidv,0.001,10^-3,d2Chidv);
    n(k) = N(v(k));
    L_D(k) = CL(n(k),v(k),w(k))/CD(v(k));
    empuje(k) = T(v(k));
    CL2(k) = CL(n(k),v(k),w(k));
    palanca(k) = delta_T;
    alabeo(k)  = muu(n(k));
    veloc_chi(k)   = chi_punto(v(k));
    radio_giro(k)  = R(n(k),v(k));
    w(k+1) = w(k) -C(v(k))*g*T(v(k))*(TT(k+1)-TT(k));
    
for k = 2:(num-1),
    N = @(V) (B_2.*V.^2 + V.*sqrt(B_2^2.*V.^2 + 4*C_2.*(T(V)-A_2.*V.^2)))./(2*C_2*w(k));
    chi_punto = @(V) g.*sqrt(N(V).^2-1)./V;
    VV = linspace(v(k-1)-epsilon,v(k-1),num2); 
    for ii = 1:num2
        chi_prima(ii) = chi_punto(VV(ii));
    end
    [veloc_chi(k),indice] = max(chi_prima);
    v(k) = VV(indice);    
    n(k) = N(v(k));
    L_D(k) = CL(n(k),v(k),w(k))/CD(v(k));
    CL2(k) = CL(n(k),v(k),w(k));
    empuje(k) = T(v(k));
    palanca(k) = delta_T;
    alabeo(k)  = muu(n(k));
    radio_giro(k)  = R(n(k),v(k));
    w(k+1) = w(k) -C(v(k))*g*T(v(k))*(TT(k+1)-TT(k));
end

fuel_viraje = (w(1)-w(end))/g;
W = w(end);
vueltas = mean(veloc_chi)*tiempo_final/(2*pi);

datos(i).segmento.tiempo = TT;
datos(i).segmento.tiempo(end) = '';
datos(i).segmento.fuel = (w(1) - w)./g; %fuel consumido acumulado
datos(i).segmento.fuel(end) = '';
datos(i).segmento.vueltas = vueltas; %escalar
datos(i).segmento.palanca = palanca;
% datos(i).segmento.empuje = empuje;
datos(i).segmento.empuje = CL2;
datos(i).segmento.L_D = L_D;
datos(i).segmento.peso = w;
datos(i).segmento.peso(end) = '';
datos(i).segmento.velocidad = v;
datos(i).segmento.alabeo = alabeo;
datos(i).segmento.factor_de_carga = n;
datos(i).segmento.velocidad_guiniada = veloc_chi;
datos(i).segmento.radio_giro = radio_giro;        
    
%-------------------------------------------------------------------------%

case 9  %PALANCA DE GASES DADA, HALLAR V TAL QUE R=RMIN

delta_T = delta_T_gases;

M = @(V) V./a;

a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);

switch tipo_motor
    case 1
        T = @(V) delta_T.*(T_SL*(1+0.2.*M(V).^2).^(1.4./0.4)*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
        switch derivacion
            case 1
                C = @(V) c_SL.*correccion.*(1+1.2.*M(V)).*sqrt(Temp/Temp_SL);
            case 2
                C = @(V) c_SL.*correccion.*(1+0.33.*M(V)).*sqrt(Temp/Temp_SL);
            case 3
                C = @(V) c_SL.*correccion.*(1+0.16875.*M(V)).*sqrt(Temp/Temp_SL);
        end
    case 2
        T = @(V) delta_T.*(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
        C = @(V) c_SL.*correccion.*(1+1.44.*M(V)).*sqrt(Temp/Temp_SL).* V./eta_p;
    case 3
        T = @(V) delta_T.*(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
        C = @(V) c_SL .*correccion.*V./eta_p;
end   

D = @(V) T(V);
CD = @(V) D(V)./(0.5.*rho.*V.^2*S);
L  = @(n,W) n.*W;
CL = @(n,V,W) L(n,W)./(0.5.*rho.*V.^2.*S);
muu = @(n) acos(1./n);
chi_punto = @(n,V) g.*sqrt(n.^2-1)./V;

% El objetivo es despejar n = n(V) y derivar e igualar a 0 para obtener la 
% velocidad que hace el factor de carga máximo, iterando para cada peso

% D = 0.5*rho*V^2*S*Cd0 - k2*W*n + k1*n^2*W^2/(0.5*rho*V^2*S) =
% D = A_2 * V^2         - B_2*W*n + C_2*n^2*W^2/V^2

A_2 = 0.5*rho*S*Cd0;
B_2 = k2;
C_2 = k1/(0.5*rho*S);

num = 200; %Numero de intervalos para discretizar el tiempo
num2 = 5000;

TT = linspace(0,tiempo_final,num);
w(1) = W_ini;
k = 1;
epsilon = 1;
    N = @(V) (B_2.*V.^2 + V.*sqrt(B_2.^2.*V.^2 + 4*C_2.*(T(V)-A_2.*V.^2)))./(2*C_2*w(k));
%     N2 = @(V) (B_2.*V.^2 - V.*sqrt(B_2.^2.*V.^2 + 4*C_2.*(T(V)-A_2.*V.^2)))./(2*C_2*w(k));
    R = @(V) V.^2./(g.*sqrt(N(V).^2-1));
%     R2 = @(V) V.^2./(g.*sqrt(N2(V).^2-1));
    vvv = linspace(0,200,100000);
    radius = R(vvv);
%     radius2 = R2(vvv);
    longitud_radius = length(radius);
    ii = 1;
    
    while ii <= longitud_radius,
        if imag(radius(ii)) ~= 0,
            radius(ii) = '';
            vvv(ii) = '';
            longitud_radius = longitud_radius - 1;
        elseif real(radius(ii)) <= 0,
            radius(ii) = '';
            vvv(ii) = '';
            longitud_radius = longitud_radius - 1;
        elseif isnan(radius(ii)) == 1,
            radius(ii) = '';
            vvv(ii) = '';
            longitud_radius = longitud_radius - 1;
        else
            ii = ii + 1;
        end
        
    end   
        
%     delete(gcf);
%     figure
%     plot(vvv,real(R(vvv)),'b');
%     hold on;
%     plot(vvv,imag(R(vvv)),'r');

%    %dNdv   = diff(N(vel));
%    dRdv = diff(R(vel));
%    dRdv = matlabFunction(dRdv);
%    d2Rdv = diff(dRdv(vel));
%    %v(k) = double(solve(dRdv));
%    %v(k) = double(solve(N(vel)*vel*dNdv - 2*(N(vel)^2-1)));
%     [v(k),func,error_raph(k)] = Newton_Raphson(dRdv,1000,10^-4,d2Rdv);
    [radio_giro(k),indice] = min(radius);
%     [radio_giro2(k),indice2] = min(radius2);
    
    v(k) = vvv(indice);
    n(k) = N(v(k));
%     n2(k) = N2(v2(k));
    L_D(k) = CL(n(k),v(k),w(k))/CD(v(k));
    CL2(k) = CL(n(k),v(k),w(k));
    empuje(k) = T(v(k));
    palanca(k) = delta_T;
    alabeo(k)  = muu(n(k));
    veloc_chi(k)   = chi_punto(n(k),v(k));
    w(k+1) = w(k) -C(v(k))*g*T(v(k))*(TT(k+1)-TT(k));

for k = 2:(num-1),
    N = @(V) (B_2.*V.^2 + V.*sqrt(B_2^2.*V.^2 + 4*C_2.*(T(V)-A_2.*V.^2)))./(2*C_2*w(k));
    chi_punto = @(V) g.*sqrt(N(V).^2-1)./V;
    R = @(V) V.^2./(g.*sqrt(N(V).^2-1));
    VV = linspace(v(k-1)-epsilon,v(k-1)+epsilon,num2); 
    for ii = 1:num2
        radio_Giro(ii) = R(VV(ii));
    end
    
    longitud_radio = length(radio_Giro);
    ii = 1;
    while ii <= longitud_radio,
        if imag(radio_Giro(ii)) ~= 0,
            radio_Giro(ii) = '';
            VV(ii) = '';
            longitud_radio = longitud_radio - 1;
        elseif real(radio_Giro(ii)) <= 0,
            radio_Giro(ii) = '';
            VV(ii) = '';
            longitud_radio = longitud_radio - 1;
        elseif isnan(radio_Giro(ii)) == 1,
            radio_Giro(ii) = '';
            VV(ii) = '';
            longitud_radio = longitud_radio - 1;
        else
            ii = ii + 1;
        end
        
    end   
   
%     delete(gcf);
%     figure
%     plot(VV,real(radio_Giro),'b');
%     hold on;
%     plot(vvv,imag(R(vvv)),'r');
    
    [radio_giro(k),indice] = min(radio_Giro);
    v(k) = VV(indice);
    n(k) = N(v(k));
    L_D(k) = CL(n(k),v(k),w(k))/CD(v(k));
    CL2(k) = CL(n(k),v(k),w(k));
    empuje(k) = T(v(k));
    palanca(k) = delta_T;
    alabeo(k)  = muu(n(k));
    veloc_chi(k)   = chi_punto(v(k));
    w(k+1) = w(k) -C(v(k))*g*T(v(k))*(TT(k+1)-TT(k));
end

fuel_viraje = (w(1)-w(end))/g;
W = w(end);
vueltas = mean(veloc_chi)*tiempo_final/(2*pi);

datos(i).segmento.tiempo = TT;
datos(i).segmento.tiempo(end) = '';
datos(i).segmento.fuel = (w(1) - w)./g; %fuel consumido acumulado
datos(i).segmento.fuel(end) = '';
datos(i).segmento.vueltas = vueltas; %escalar
datos(i).segmento.palanca = palanca;
% datos(i).segmento.empuje = empuje;
datos(i).segmento.empuje = CL2;
datos(i).segmento.L_D = L_D;
datos(i).segmento.peso = w;
datos(i).segmento.peso(end) = '';
datos(i).segmento.velocidad = v;
datos(i).segmento.alabeo = alabeo;
datos(i).segmento.factor_de_carga = n;
datos(i).segmento.velocidad_guiniada = veloc_chi;
datos(i).segmento.radio_giro = radio_giro;        
%-------------------------------------------------------------------------%


end
datos(i).nombre = 'Viraje';
datos(i).segmento.altura = h_inicial;
datos(i).segmento.distancia = 0;

% datos(i).lista_variables = [{''};'Tiempo';'Fuel';'Peso';'Vueltas';'Velocidad';'Velocidad de guiñada';'Angulo de alabeo';...
%     'Factor de carga';'Radio de giro';'Palanca';'CL';'L/D';'Altura'];

global flag_palanca 
if max(palanca) > 1.2,
    flag_palanca = 1;
else
    flag_palanca = 0;
end

end
