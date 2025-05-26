function [fuel_total,tiempo_total,S_total,datos] = analisis_aterrizaje(aterrizaje_vec,propul,aerodinamica,aero_aterrizaje,W_inicial,datos,i,Prop_data,data_electric)

% PROPUL
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR
% 4: CONSUMO ESPECIFICO
% 5: AVION CIVIL =1/MILITAR = 2
% 6: EFICIENCIA DE LA HELICE (ETA_P)
% 7: DERIVACION(TURBOFANES)

% AERODINAMICA
% 1: SUPERFICIE
% 2: CD0
% 3: K1 = K
% 4: K2
% 5: CLMAX LIMPIO

% % --> %% AERODINAMICA_ATERRIZAJE
%       % 1: CD0 EN ATERRIZAJE (CONSIDERANDO TREN + FLAPS + SPOILERS)
%       % 2: CL EN ATERRIZAJE
%       % 3: K CON EFECTO SUELO (O ALTURA DEL ALA SOBRE EL SUELO Y
%       ENVERGADURA DEL AVION)
%       % 4: CLMAX SUCIO ATERRIZAJE (FLAPS = 30º)

% --> %% ATERRIZAJE
%       % 1: TEMPERATURA LOCAL
%       % 2: ALTURA LOCAL       
%       % 3: PRESION LOCAL
%       % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
%       % 5: PALANCA DE GAS PARA REVERSA
%       % 6: TIEMPO QUE SE TARDA EN ACTIVAR LOS FRENOS AL ATERRIZAR

aterrizaje(1) = aterrizaje_vec.temp_local;
aterrizaje(2) = aterrizaje_vec.h_inicial;
aterrizaje(3) = aterrizaje_vec.P_local;
aterrizaje(4) = aterrizaje_vec.mu_landing;
aterrizaje(5) = aterrizaje_vec.palanca;
aterrizaje(6) = aterrizaje_vec.t_brake;

R = 287.058;
g = 9.80665;
W = W_inicial;
%------------------- Recoleccion de datos ---------------------%

% ATERRIZAJE
temp_local = aterrizaje(1);
alt_local = aterrizaje(2);
pres_local = aterrizaje(3);
coef_fric = aterrizaje(4);
%delta_T = aterrizaje(5);
tiempo_frenos = aterrizaje(6);

% AERODINAMICA
S = aerodinamica(1);
Cd0 = aerodinamica(2);
k1 = aerodinamica(3);
k2 = aerodinamica(4);
CLmax_limpio = aerodinamica(5);

% AERODINAMICA_ATERRIZAJE
Cd0_land = aero_aterrizaje(1);
Cl_land = aero_aterrizaje(3);
k_land = aero_aterrizaje(2);
CLmax_land = aero_aterrizaje(4);

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

% ATMOSFERA ESTANDAR
[temp_SL,rho_SL,p_SL,a_SL] = atmos_inter_mio(0);
%-------------------------------------------------------------------%
% Calculo de la densidad local y los valores por defecto de la atmosfera
[temp_land, rho_land, p_land, a_land] = atmos_inter_mio(alt_local);
 
V_stall_land=sqrt(2*W/rho_land/S/CLmax_land);
a_land = sqrt(1.4*temp_local*R);
 
switch civil_mil
    case 1
        V_acer = 1.3 * V_stall_land;
        V_flare = 1.23 * V_stall_land;
        V_TD = 1.15 * V_stall_land;
    case 2
        V_acer = 1.2 * V_stall_land;
        V_flare = 1.15 * V_stall_land;
        V_TD = 1.1 * V_stall_land;
end

iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();   %ESTA LINEA PERMITE USAR CONDICIONES EN LAS FUNCIONES!!!

T_SL = T_SL * num_motores;

% -------------------- Acercamiento -------------------------%
                                                    
M = @(V) V./a_land;
eta_p = @(V) iif(M(V)<=0.1, @() eta_inst.*M(V)./0.1, ...
                 M(V)>0.1,@() eta_inst);

a_1=3.559957437510763;
a_2=-10.739698199171459;
a_3= 11.989635150373475;
a_4=-5.869876557884609;
a_5=2.059994459180667;

correccion = @(delta_T) (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);

switch tipo_motor
    case 1
        T = @(V,delta_T) delta_T.*T_SL.*(1+0.2.*M(V).^2).^(1.4/0.4).*(1-0.49.*sqrt(M(V))).*rho_land./rho_SL;
        switch derivacion
            case 1
                C = @(V,delta_T) c_SL.*correccion(delta_T).*(1+1.2.*M(V)).*sqrt(temp_local/temp_SL);
            case 2
                C = @(V,delta_T) c_SL.*correccion(delta_T).*(1+0.33.*M(V)).*sqrt(temp_local/temp_SL);
            case 3
                C = @(V,delta_T) c_SL.*correccion(delta_T).*(1+0.16875.*M(V)).*sqrt(temp_local/temp_SL);
        end
    case 2
        T  = @(V,delta_T) delta_T.*T_SL.*eta_p(V)./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*pres_local./p_SL;
        C = @(V,delta_T) c_SL.*correccion(delta_T).*(1+1.44.*M(V)).*sqrt(temp_local./temp_SL).* V./eta_p(V);
    case 3
        T = @(V,delta_T) delta_T.*T_SL.*eta_p(V)./V.*((8.55.*rho_land./rho_SL - 1)./7.55);
        C = @(V,delta_T) c_SL .*correccion(delta_T).*V./eta_p(V);
end

L_D_acer=1/(0.5*rho_land*V_acer^2*Cd0_land/(W/S)+...
    (W*k_land/(S*0.5*rho_land*V_acer^2)));

h_obst=15.24; %50 ft
           
% FLARE 

delta_T_flare = 1;

L_D_flare=1/(0.5*rho_land*V_flare^2*Cd0_land/(W/S)+...
    (W*k_land/(S*0.5*rho_land*V_flare^2)));

gamma_flare = abs(asin(T(V_flare,delta_T_flare)./W - 1/L_D_flare));
R_flare=V_flare^2/(0.2*g);
h_tr = R_flare*(1-cos(gamma_flare));

while h_tr > h_obst || gamma_flare > 0.087
delta_T_flare = delta_T_flare - 0.001;
L_D_flare=1/(0.5*rho_land*V_flare^2*Cd0_land/(W/S)+...
    (W*k_land/(S*0.5*rho_land*V_flare^2)));
gamma_flare = abs(asin(T(V_flare,delta_T_flare)./W - 1/L_D_flare));
h_tr = R_flare*(1-cos(gamma_flare));
end

if h_tr > h_obst
    h_tr = 10.866;
    R_flare = h_tr/(1-cos(gamma_flare));
    delta_T_flare = 0.05;
end

S_flare = R_flare * sin(gamma_flare);
tiempo_flare = S_flare/(V_flare * cos(gamma_flare));
fuel_flare = C(V_flare,delta_T_flare)*T(V_flare,delta_T_flare)*tiempo_flare;

gamma_acer= asin(T(V_acer,delta_T_flare)./W - 1/L_D_acer);
S_acer = abs((h_obst - h_tr)/tan(gamma_acer));
tiempo_acer = S_acer /(V_acer * cos(gamma_acer));
fuel_acer = C(V_acer,delta_T_flare)*T(V_acer,delta_T_flare)*tiempo_acer;

W = W - fuel_acer * g;

W = W  - fuel_flare*g;

% Freeroll
delta_T_free = 0.05;
tiempo_free = tiempo_frenos;
S_free = V_TD * tiempo_frenos;
fuel_free = C(V_TD,delta_T_free) * T(V_TD,delta_T_free) * tiempo_frenos;

W = W - fuel_free * g;
%---------------------------------- GROUND -------------------------------%
CD_land = Cd0_land + k_land * Cl_land^2;
CD_g = CD_land - coef_fric*Cl_land;

delta_T_at = -aterrizaje(5);

dvdt = @(V) g.*((T(V,delta_T_at)./W - coef_fric) - 0.5.*rho_land.*V.^2.*S.*CD_g./W);

dt = @(V) 1./dvdt(V);
dS = @(V) V./dvdt(V);
dW = @(V) C(V,0.05).*T(V,delta_T_at)./dvdt(V);

tiempo_ground = integral(dt,V_TD,0);
S_ground = integral(dS,V_TD,0);
fuel_ground = abs(integral(dW,V_TD,0));

W = W - fuel_ground * g;
%------------------------------- RESUMEN ---------------------------%
tiempo_total = tiempo_ground + tiempo_flare + tiempo_acer + tiempo_free;
S_total = S_ground + S_flare + S_acer + S_free;
fuel_total = fuel_ground + fuel_flare + fuel_acer + fuel_free;

% tiempo_aterrizaje = [tiempo_acer tiempo_flare tiempo_free tiempo_ground tiempo_total];
% S_aterrizaje = [S_acer  S_flare  S_free S_ground S_total];
% fuel_aterrizaje = [fuel_acer fuel_flare fuel_free fuel_ground fuel_total];

datos(i).segmento.acercamiento.tiempo = tiempo_acer;
datos(i).segmento.acercamiento.distancia = S_acer;
datos(i).segmento.acercamiento.fuel = fuel_acer;
datos(i).segmento.acercamiento.L_D = L_D_acer;
datos(i).segmento.acercamiento.altura = h_obst;
datos(i).segmento.acercamiento.gamma = gamma_acer;
datos(i).segmento.acercamiento.velocidad = V_acer;

datos(i).segmento.flare.tiempo = tiempo_flare;
datos(i).segmento.flare.distancia = S_flare;
datos(i).segmento.flare.fuel = fuel_flare;
datos(i).segmento.flare.L_D = L_D_flare;
datos(i).segmento.flare.velocidad = V_flare;
datos(i).segmento.flare.altura = h_tr;
datos(i).segmento.flare.gamma = gamma_flare;
datos(i).segmento.flare.radio = R_flare;

datos(i).segmento.free.tiempo = tiempo_free;
datos(i).segmento.free.distancia = S_free;
datos(i).segmento.free.fuel = fuel_free;
datos(i).segmento.free.velocidad = V_TD;

datos(i).segmento.ground.tiempo = tiempo_ground;
datos(i).segmento.ground.distancia = S_ground;
datos(i).segmento.ground.fuel = fuel_ground;

datos(i).segmento.total.tiempo = tiempo_total;
datos(i).segmento.total.distancia = S_total;
datos(i).segmento.total.fuel = fuel_total;

datos(i).nombre = 'Aterrizaje';
datos(i).lista_variables = [{''};'Acercamiento: tiempo';'Acercamiento: distancia';'Acercamiento: fuel';'Flare: tiempo';'Flare: distancia';'Flare: fuel';...
    'Free: tiempo';'Free: distancia';'Free: fuel';'Ground: tiempo';'Ground: distancia';'Ground: fuel';'Tiempo total';'Distancia total';'Fuel total'];




