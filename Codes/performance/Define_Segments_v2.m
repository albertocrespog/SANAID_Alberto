function Segments = Define_Segments_v2(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE,OUTPUT_read_XLSX)

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
ft2m = conv_UNITS.ft2m;
m2ft = conv_UNITS.m2ft;

% climb_mode = FlightMODE.climb_mode;
% cruise_mode = FlightMODE.cruise_mode;
% turn_mode = FlightMODE.turn_mode;
% descent_mode = FlightMODE.descent_mode;

%% Taxy
temp_local_taxy = OUTPUT_read_XLSX.IPP_flags.temp_local_taxy;  % 1: TEMPERATURA LOCAL (Celsius)
h_inicial_taxy = OUTPUT_read_XLSX.IPP_flags.h_inicial_taxy; % (m)
P_local_taxy = OUTPUT_read_XLSX.IPP_flags.P_local_taxy; % (Pa)
delta_T_taxy = OUTPUT_read_XLSX.IPP_flags.delta_T_taxy; % 4: PALANCA DE RALENTI EN TAXI = 0.05
V_taxy = OUTPUT_read_XLSX.IPP_flags.V_taxy; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
t_taxy = OUTPUT_read_XLSX.IPP_flags.t_taxy; % 6: TIEMPO DE ESPERA EN TAXI (s)
% Store
taxy.T = temp_local_taxy;  % 1: TEMPERATURA LOCAL (Celsius)
taxy.h = h_inicial_taxy; % (m)
taxy.P = P_local_taxy; % (Pa)
taxy.dT = delta_T_taxy; % 4: PALANCA DE RALENTI EN TAXI = 0.05
taxy.V = V_taxy; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
taxy.t = t_taxy; % 6: TIEMPO DE ESPERA EN TAXI (s)

%% Takeoff
temp_local_TO = OUTPUT_read_XLSX.IPP_flags.temp_local_TO;  % 1: TEMPERATURA LOCAL (Celsius)
h_inicial_TO = OUTPUT_read_XLSX.IPP_flags.h_inicial_TO; % 2: ALTURA LOCAL (m)
P_local_TO = OUTPUT_read_XLSX.IPP_flags.P_local_TO; % 3: PRESION LOCAL (Pa)
mu_TO = OUTPUT_read_XLSX.IPP_flags.mu_TO; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
h_obstacle_TO = OUTPUT_read_XLSX.IPP_flags.h_obstacle_TO; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
gamma_climb_TO = OUTPUT_read_XLSX.IPP_flags.gamma_climb_TO; % 6: GAMMA DE SUBIDA MINIMO
delta_T_TO = OUTPUT_read_XLSX.IPP_flags.delta_T_TO; % 7: PALANCA DE GASES PARA DESPEGUE
% Store
despegue.T = temp_local_TO; % 1: TEMPERATURA LOCAL (Celsius)
despegue.h = h_inicial_TO; % 2: ALTURA LOCAL (m)
despegue.P = P_local_TO; % 3: PRESION LOCAL (Pa)
despegue.mu = mu_TO; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
despegue.ho = h_obstacle_TO; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
despegue.g = gamma_climb_TO; % 6: GAMMA DE SUBIDA MINIMO
despegue.dT = delta_T_TO; % 7: PALANCA DE GASES PARA DESPEGUE
%% Climb
h_inicial_cl = OUTPUT_read_XLSX.IPP_flags.h_inicial_cl; % 1: ALTURA INICIAL - [m]
h_final_cl = OUTPUT_read_XLSX.IPP_flags.h_final_cl; % 1: ALTURA FINAL - [m]
gamma_cl = OUTPUT_read_XLSX.IPP_flags.gamma_cl; % 2: GAMMA DE SUBIDA  - [-]
Mach_cl = OUTPUT_read_XLSX.IPP_flags.Mach_cl; % 3: MACH DE VUELO  - [-]
TAS_cl = OUTPUT_read_XLSX.IPP_flags.TAS_cl; % 4: VELOCIDAD TAS  - [m/s]
EAS_cl = OUTPUT_read_XLSX.IPP_flags.EAS_cl; % 5: VELOCIDAD EAS  - [m/s]
delta_T_cl = OUTPUT_read_XLSX.IPP_flags.delta_T_cl; % 6: PALANCA DE GASES  - [-]
V_ini_cl = OUTPUT_read_XLSX.IPP_flags.V_ini_cl; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
V_fin_cl = OUTPUT_read_XLSX.IPP_flags.V_fin_cl; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
% Stored
subida.hi = h_inicial_cl; % 1: ALTURA INICIAL - [m]
subida.hf = h_final_cl; % 1: ALTURA FINAL - [m]
subida.g = gamma_cl; % 2: GAMMA DE SUBIDA  - [-]
subida.M = Mach_cl; % 3: MACH DE VUELO  - [-]
subida.T = TAS_cl; % 4: VELOCIDAD TAS  - [m/s]
subida.E = EAS_cl; % 5: VELOCIDAD EAS  - [m/s]
subida.dT = delta_T_cl; % 6: PALANCA DE GASES  - [-]
subida.vi = V_ini_cl; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
subida.vf = V_fin_cl; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
%% VTOL Climb
% Stored
subidavt.hi = OUTPUT_read_XLSX.IPP_flags.h_inicial_vtcl; % 1: ALTURA INICIAL - [m]
subidavt.hf = OUTPUT_read_XLSX.IPP_flags.h_final_vtcl; % 2: ALTURA FINAL - [m]
subidavt.thover = OUTPUT_read_XLSX.IPP_flags.t_hover; % 3: TIEMPO EN HOVERING  - [s]
subidavt.dT = OUTPUT_read_XLSX.IPP_flags.delta_T_vtcl; % 4: PALANCA  - [-]
subidavt.vc = OUTPUT_read_XLSX.IPP_flags.vclimb_vtcl; % 5: VELOCIDAD CLIMB  - [m/s]
subidavt.mbathover = OUTPUT_read_XLSX.IPP_flags.mbat_vtcl; % 6: MASA BRUTA BATERÍAS  - [kg]
%% Cruise
h_inicial_cr = OUTPUT_read_XLSX.IPP_flags.h_inicial_cr;% - [m] % Altura inicial
dist_final_cr = OUTPUT_read_XLSX.IPP_flags.dist_final_cr;% - [m] % 1: DISTANCIA FINAL
V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
for i = 1:length(V_cr)
    [Data_ATM(i) Performance] = Flight_Conditions_2020_v1(h_inicial_cr(i),V_cr(i));
    Mach_cr(i) = V_cr(i)/Data_ATM(i).a;% - [-] % 2: MACH DE VUELO
end
CL_cr = sqrt(Aero_TH.CD0/Aero_TH.CD2);% - [-]% 3: CL DE CRUCERO


delta_T_cr = OUTPUT_read_XLSX.IPP_flags.delta_T_cr;% - []% 4: PALANCA DE GASES
V_ini_cr = OUTPUT_read_XLSX.IPP_flags.V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
V_fin_cr = OUTPUT_read_XLSX.IPP_flags.V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
fuel_cr = OUTPUT_read_XLSX.IPP_flags.fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
Cd0_cr = OUTPUT_read_XLSX.IPP_flags.Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
k1_cr = OUTPUT_read_XLSX.IPP_flags.k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
k2_cr = OUTPUT_read_XLSX.IPP_flags.k2_cr;% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
m_bat = OUTPUT_read_XLSX.IPP_flags.mbat;% - [kg] % 11: MASA EN BRUTO DE LAS BATERÍAS

% Stored
crucero.hi = h_inicial_cr;% - [m] % Altura inicial
crucero.d = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
crucero.vcr = V_cr; % Velocidad de crucero
crucero.M = Mach_cr;% - [-] % 2: MACH DE VUELO
crucero.CL = CL_cr;% - [-]% 3: CL DE CRUCERO
crucero.dT = delta_T_cr;% - []% 4: PALANCA DE GASES
crucero.vi = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
crucero.vf = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
crucero.f = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
crucero.cd0 = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
crucero.k1 = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
crucero.k2 = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
crucero.mbat = m_bat; % - [kg] % 11: MASA EN BRUTO DE LAS BATERÍAS
%% Load Deployment
carga_loadep = OUTPUT_read_XLSX.IPP_flags.carga_loadep;% - [kg]
loadep.carga = carga_loadep;
%% Turn
h_inicial_tr = OUTPUT_read_XLSX.IPP_flags.h_inicial_tr;% - [m]
t_final_tr = OUTPUT_read_XLSX.IPP_flags.t_final_tr;% 1: TIEMPO FINAL (seg)
V_turn = OUTPUT_read_XLSX.IPP_flags.V_turn; % turn velocity
% [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
% Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
for i = 1:length(V_turn)
    [Data_ATM(i) Performance] = Flight_Conditions_2020_v1(h_inicial_cr(i),V_turn(i));
    Mach_tr(i) = V_turn(i)/Data_ATM(i).a;% - [-] % 2: MACH DE VUELO
end
delta_T_tr = OUTPUT_read_XLSX.IPP_flags.delta_T_tr; % 3: PALANCA DE GASES
CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
phi_tr = OUTPUT_read_XLSX.IPP_flags.phi_tr;% 5: ANGULO DE ALABEO (rads)
V_psi = OUTPUT_read_XLSX.IPP_flags.V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
n_tr = OUTPUT_read_XLSX.IPP_flags.n_tr;% 7: FACTOR DE CARGA
R_tr = OUTPUT_read_XLSX.IPP_flags.R_tr;% 8: RADIO DE GIRO (m)
% Stored
viraje.hi = h_inicial_tr;% - [m]
viraje.tf = t_final_tr;% 1: TIEMPO FINAL (seg)
viraje.vtr = V_turn; % turn velocity
viraje.M = Mach_tr;% - [-] % 2: MACH DE VUELO
viraje.dT = delta_T_tr; % 3: PALANCA DE GASES
viraje.CL = CL_tr;% 4: CL DE VIRAJE
viraje.phi = phi_tr;% 5: ANGULO DE ALABEO (rads)
viraje.vpsi = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
viraje.n = n_tr;% 7: FACTOR DE CARGA
viraje.R = R_tr;% 8: RADIO DE GIRO (m)
%% Descent
h_inicial_d = OUTPUT_read_XLSX.IPP_flags.h_inicial_d;
h_final_d = OUTPUT_read_XLSX.IPP_flags.h_final_d;
gamma_d = OUTPUT_read_XLSX.IPP_flags.gamma_d;
Mach_d = OUTPUT_read_XLSX.IPP_flags.Mach_d;
EAS_d = OUTPUT_read_XLSX.IPP_flags.EAS_d;
TAS_d = OUTPUT_read_XLSX.IPP_flags.TAS_d;
delta_T_d = OUTPUT_read_XLSX.IPP_flags.delta_T_d;
V_ini_d = OUTPUT_read_XLSX.IPP_flags.V_ini_d;
V_fin_d = OUTPUT_read_XLSX.IPP_flags.V_fin_d;
% Stored
descenso.hi = h_inicial_d;
descenso.hf = h_final_d;
descenso.g = gamma_d;
descenso.M = Mach_d;
descenso.E = EAS_d;
descenso.T = TAS_d;
descenso.dT = delta_T_d;
descenso.vi = V_ini_d;
descenso.vf = V_fin_d;
%% VTOL Descent
% Stored
descensovt.hi = OUTPUT_read_XLSX.IPP_flags.h_inicial_vtd; % 1: ALTURA INICIAL - [m]
descensovt.hf = OUTPUT_read_XLSX.IPP_flags.h_final_vtd; % 2: ALTURA FINAL - [m]
descensovt.vd = OUTPUT_read_XLSX.IPP_flags.vdes_vtd; % 5: VELOCIDAD CLIMB  - [m/s]
%% Turn waiting area
% h_inicial_tr = OUTPUT_read_XLSX.IPP_flags.h_inicial_tr;% - [m]
% t_final_tr = OUTPUT_read_XLSX.IPP_flags.t_final_tr;% 1: TIEMPO FINAL (seg)
% V_turn = OUTPUT_read_XLSX.IPP_flags.V_turn; % turn velocity
% [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
% Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
% delta_T_tr = OUTPUT_read_XLSX.IPP_flags.delta_T_tr; % 3: PALANCA DE GASES
% CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
% phi_tr = OUTPUT_read_XLSX.IPP_flags.phi_tr;% 5: ANGULO DE ALABEO (rads)
% V_psi = OUTPUT_read_XLSX.IPP_flags.V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
% n_tr = OUTPUT_read_XLSX.IPP_flags.n_tr;% 7: FACTOR DE CARGA
% R_tr = OUTPUT_read_XLSX.IPP_flags.R_tr;% 8: RADIO DE GIRO (m)
% % Stored
% viraje_wt(1) = h_inicial_tr;% - [m]
% viraje_wt(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
% viraje_wt(3) = V_turn; % turn velocity
% viraje_wt(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
% viraje_wt(5) = delta_T_tr; % 3: PALANCA DE GASES
% viraje_wt(6) = CL_tr;% 4: CL DE VIRAJE
% viraje_wt(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
% viraje_wt(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
% viraje_wt(9) = n_tr;% 7: FACTOR DE CARGA
% viraje_wt(10) = R_tr;% 8: RADIO DE GIRO (m)
%% LANDING
temp_local_LND = OUTPUT_read_XLSX.IPP_flags.temp_local_TO; % 1: TEMPERATURA LOCAL (Celsius)
h_inicial_LND = OUTPUT_read_XLSX.IPP_flags.temp_local_TO; % 2: ALTURA LOCAL (m)
P_local_LND = OUTPUT_read_XLSX.IPP_flags.temp_local_TO; % 3: PRESION LOCAL (Pa)
mu_LND = OUTPUT_read_XLSX.IPP_flags.temp_local_TO; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
delta_T_LND = OUTPUT_read_XLSX.IPP_flags.temp_local_TO; % PALANCA DE REVERSA
t_brake = OUTPUT_read_XLSX.IPP_flags.temp_local_TO; % 'Tiempo en activar frenos' - [s]
% Stored
aterrizaje.T = temp_local_LND;  % 1: TEMPERATURA LOCAL (Celsius)
aterrizaje.h = h_inicial_LND; % 2: ALTURA LOCAL (m)
aterrizaje.P = P_local_LND; % 3: PRESION LOCAL (Pa)
aterrizaje.mu = mu_LND; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
aterrizaje.dT = delta_T_LND; % PALANCA DE REVERSA
aterrizaje.tb = t_brake; % 'Tiempo en activar frenos' - [s]

%% Dummy segment to account for 3 segments
dummy(1) = OUTPUT_read_XLSX.IPP_flags.dummy(1);% - [m] % Altura inicial
dummy(2) = OUTPUT_read_XLSX.IPP_flags.dummy(2);% - [m] % 1: DISTANCIA FINAL
dummy(3) = OUTPUT_read_XLSX.IPP_flags.dummy(3); % Velocidad de crucero
[Data_ATM Performance] = Flight_Conditions_2020_v1(dummy(1),dummy(3));
dummy(4) = dummy(3)/Data_ATM.a;% - [-] % 2: MACH DE VUELO

% Flight Modes
% MODES.cruise_mode = cruise_mode;
% MODES.climb_mode = climb_mode;
% MODES.turn_mode = turn_mode;
% MODES.descent_mode = descent_mode;

%% Stored
Segments{1}.taxy        = taxy;
Segments{1}.despegue    = despegue;
Segments{1}.subida      = subida;
Segments{1}.subidavt    = subidavt;
Segments{1}.crucero     = crucero;
Segments{1}.loadep      = loadep;
Segments{1}.viraje      = viraje;
% Segments{1}.viraje_wt = viraje_wt;
Segments{1}.descenso    = descenso;
Segments{1}.descensovt  = descensovt;
Segments{1}.aterrizaje  = aterrizaje;
Segments{1}.dummy       = dummy;
Segments{1}.MODES       = FlightMODE;

