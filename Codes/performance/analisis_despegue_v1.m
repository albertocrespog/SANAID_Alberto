function [fuel_total,tiempo_total,S_total,datos] = analisis_despegue(despegue,propul,aerodinamica,aero_despegue,W_inicial,i,datos)

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

% --> %% AERODINAMICA_DESPEGUE
%       % 1: CD0 EN DESPEGUE (CONSIDERANDO TREN + FLAPS)
%       % 2: CL EN DESPEGUE
%       % 3: K CON EFECTO SUELO (O ALTURA DEL ALA SOBRE EL SUELO Y
%       ENVERGADURA DEL AVION)
%       % 4: CLMAX SUCIO DESPEGUE (FLAPS = 20º)

% --> %% DESPEGUE
%       % 1: TEMPERATURA LOCAL
%       % 2: ALTURA LOCAL       
%       % 3: PRESION LOCAL
%       % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
%       % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
%       % 6: GAMMA DE SUBIDA MINIMO
%       % 7: PALANCA DE GASES PARA DESPEGUE
% --> %% DATOS DESPEGUE
%       % 1: PALANCA DE GASES PARA DESPEGUE
%       % 2: MACH DE DESPEGUE PARA VLO

R = 287.058;
g = 9.80665;
W = W_inicial;
%------------------- Recoleccion de datos ---------------------%

% DESPEGUE
temp_local = despegue(1);
alt_local = despegue(2);
pres_local = despegue(3);
coef_fric = despegue(4);
h_obst = (15*propul(5) + 20)*0.3048;

delta_T = despegue(7);

% AERODINAMICA
S = aerodinamica(1);
Cd0 = aerodinamica(2);
k1 = aerodinamica(3);
k2 = aerodinamica(4);
CLmax_limpio = aerodinamica(5);

% AERODINAMICA_DESPEGUE
Cd0_desp = aero_despegue(1);
Cl_desp = aero_despegue(3);
k_suelo = aero_despegue(2);
CLmax_desp = aero_despegue(4);

% PROPULSION
tipo_motor = propul(1);
num_motores = propul(2);
T_SL = propul(3);
c_SL = propul(4);
civil_mil = propul(5);
eta_inst = propul(6);
derivacion = propul(7);

% ATMOSFERA ESTANDAR
[temp_SL,rho_SL,p_SL,a_SL] = atmos_inter_mio(0);
%-------------------------------------------------------------------%
% Calculo de la densidad local y los valores por defecto de la atmosfera
[temp_desp, rho_desp, p_desp, a_desp] = atmos_inter_mio(alt_local);

%Calculo de M_0
V_stall_desp=sqrt(2*W_inicial/(S*CLmax_desp*rho_desp));
VLO=1.2*V_stall_desp;
a_desp = sqrt(1.4*R*temp_local);

M_0 = VLO/a_desp;

iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();   %ESTA LINEA PERMITE USAR CONDICIONES EN LAS FUNCIONES!!!

T_SL = T_SL * num_motores;

%Calculo tiempo de despegue
% ---------------------- GROUND ------------------------------------%
CD_desp = Cd0_desp + k_suelo * Cl_desp^2;
CD_g = CD_desp - coef_fric*Cl_desp;

%Codigo para integrar en velocidad el empuje, el consumo y la variacion
%de peso

M = @(V) V./a_desp;
eta_p = @(V) iif(M(V)<=0.1, @() eta_inst.*M(V)./0.1, ...
                 M(V)>0.1,@() eta_inst);

a_1=3.559957437510763;
a_2=-10.739698199171459;
a_3= 11.989635150373475;
a_4=-5.869876557884609;
a_5=2.059994459180667;
c_SL=c_SL*(a_1*delta_T^4+a_2*delta_T^3+a_3*delta_T^2+a_4*delta_T+a_5);

switch tipo_motor
    case 1
        T = @(V) delta_T.*T_SL.*(1+0.2.*M(V).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(V))).*(rho_desp./rho_SL);
        switch derivacion
            case 1
                C = @(V) c_SL.*(1+1.2.*M(V))*sqrt(temp_local/temp_SL);
            case 2
                C = @(V) c_SL.*(1+0.33.*M(V))*sqrt(temp_local/temp_SL);
            case 3
                C = @(V) c_SL.*(1+0.16875.*M(V))*sqrt(temp_local/temp_SL);
        end
    case 2
        T  = @(V) delta_T.*T_SL.*eta_p(V)./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*pres_local./p_SL;
        C = @(V) c_SL.*(1+1.44.*M(V)).*sqrt(temp_local/temp_SL).* V./eta_p(V);
    case 3
        T = @(V) delta_T.*T_SL*eta_p(V)./V.*((8.55*rho_desp/rho_SL - 1)/7.55);
        C = @(V) c_SL .*V/eta_p(V);
end

dvdt = @(V) g.*((T(V)./W - coef_fric) - 0.5.*rho_desp.*V.^2.*S.*CD_g./W);

dt = @(V) 1./dvdt(V);
dS = @(V) V./dvdt(V);
dW = @(V) C(V).*T(V)./dvdt(V);

tiempo_ground = integral(dt,0,VLO);
S_ground = integral(dS,0,VLO);
fuel_ground = integral(dW,0,VLO);

tiempo_ground = real(tiempo_ground);
S_ground = real(S_ground);
fuel_ground = real(fuel_ground);

W = W - fuel_ground*g;
% ---------------------- TRANSICION -------------------------%
V_tr = VLO;
M_tr = M_0;

R_tr=V_tr^2/(0.2*g);
D = (1/2) * rho_desp * V_tr^2 * S * Cd0_desp + k_suelo*W^2/((1/2) * rho_desp * V_tr^2 * S);
gamma_subida = real(asin((T(V_tr) - D)/W));
h_tr = R_tr*(1-cos(gamma_subida));

if h_tr > h_obst,    
h_tr = h_obst;
R_tr = h_obst/(1-cos(gamma_subida));
end
% indicador = 1;
% while h_tr > h_obst,
%     gamma_subida = real(asin((T(V_tr)/(indicador*1.2) - D)/W));
%     h_tr = R_tr*(1-cos(gamma_subida));
%     indicador = indicador + 0.2;
% end 

S_tr = R_tr * sin(gamma_subida);
tiempo_tr = S_tr/V_tr;
fuel_tr = C(V_tr)*T(V_tr)*tiempo_tr;

W = W - fuel_tr*g;
% -------------------------- CLIMB -------------------------------%
if (h_obst-h_tr) <= 0,
    S_c = 0;
    V_c = 0;
    tiempo_c = 0;
    fuel_c = 0;
    h_tr_muy_elevada = 'Altura de transicion superior a la altura del obstaculo';
else
    h_tr_muy_elevada = 'Altura de transicion correcta';
S_c = (h_obst - h_tr)/tan(gamma_subida);
V_c = V_tr;
tiempo_c = S_c/(V_c*cos(gamma_subida));
fuel_c = C(V_c)*T(V_c)*tiempo_c;
end

W = W - fuel_c*g;
%------------------------------- RESUMEN ---------------------------%
tiempo_total = tiempo_ground + tiempo_tr + tiempo_c;
S_total = S_ground + S_tr + S_c;
fuel_total = fuel_ground + fuel_tr + fuel_c;

% tiempo_despegue = [tiempo_ground tiempo_tr tiempo_c tiempo_total];
% S_despegue = [S_ground  S_tr  S_c S_total];
% fuel_despegue = [fuel_ground fuel_tr fuel_c fuel_total];

datos(i).despegue.ground.tiempo = tiempo_ground;
datos(i).despegue.ground.distancia = S_ground;
datos(i).despegue.ground.fuel = fuel_ground;

datos(i).despegue.transicion.tiempo = tiempo_tr;
datos(i).despegue.transicion.distancia = S_tr;
datos(i).despegue.transicion.fuel = fuel_tr;
datos(i).despegue.transicion.velocidad = V_tr;
datos(i).despegue.transicion.altura = h_tr;
datos(i).despegue.transicion.gamma = gamma_subida;

datos(i).despegue.climb.tiempo = tiempo_c;
datos(i).despegue.climb.distancia = S_c;
datos(i).despegue.climb.fuel = fuel_c;
datos(i).despegue.climb.velocidad = V_c;
datos(i).despegue.climb.altura = h_obst;
datos(i).despegue.climb.detalle = h_tr_muy_elevada;

datos(i).despegue.total.tiempo = tiempo_total;
datos(i).despegue.total.distancia = S_total;
datos(i).despegue.total.fuel = fuel_total;

datos(i).nombre = 'Despegue';
datos(i).lista_variables = [{''};'Ground: tiempo';'Ground: distancia';'Ground: fuel';'Transicion: tiempo';'Transicion: distancia';'Transicion: fuel';...
    'Climb: tiempo';'Climb; distancia';'Climb: fuel';'Tiempo total';'Distancia total';'Fuel total'];


