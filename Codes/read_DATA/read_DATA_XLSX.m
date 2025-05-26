function OUTPUT_read_XLSX = read_DATA_XLSX(case_AC)

%% Pre-load variables stored in Excel
% Function that reads all the data for the different aircraft
get_add_path
filename = '../AIRCRAFT/SANAID_AIRCRAFT_DATA_v2.xlsx';
% [filename,path] = uigetfile('*.xlsx','../../AIRCRAFT','Select an Excel File');
% s = strcat(path,filename)
    
% Writes in the excel to assign the propper aircraft type
sheet = 3;
status = xlswrite(filename,case_AC,sheet,'D8');

%% MATLAB Flags
sheet = 1;
MATLAB_in = xlsread(filename,sheet,'D3');
CHECK_Efficiency = xlsread(filename,sheet,'D4');
Detailed_Profile = xlsread(filename,sheet,'D5');
save_MAT = xlsread(filename,sheet,'D6');

% Stores the flags
MATLAB_flags.MATLAB_in = MATLAB_in;
MATLAB_flags.CHECK_Efficiency = CHECK_Efficiency;
MATLAB_flags.Detailed_Profile = Detailed_Profile;
MATLAB_flags.save_MAT = save_MAT;
% Gathers all the flags
OUTPUT_read_XLSX.MATLAB_flags = MATLAB_flags;

%% Studies Flags
sheet = 2;
AERODYNAMIC_STUDY = xlsread(filename,sheet,'D3'); % Conducts Aerodynamic Studies
PROP_STUDY = xlsread(filename,sheet,'D4'); % Conducts Proppeller Studies
PERFORMANCE_STUDY = xlsread(filename,sheet,'D5'); % Conducts Performance Studies
STABILITY_STUDY = xlsread(filename,sheet,'D6'); % Conducts Stability Studies
STUDY_Weight_Fraction = xlsread(filename,sheet,'D7'); % Weight Fraction Study
ANALYSIS_PERFORMANCE_AP = xlsread(filename,sheet,'D8'); % Performance Analysis integrated with AP Codes
ANALYSIS_PERFORMANCE_AP_var = xlsread(filename,sheet,'D9'); % Performance Analysis integrated with AP Codes varying Conditions
variable_speed_AP = xlsread(filename,sheet,'D10'); % Variable Speed Performance Studies
variable_weight_AP = xlsread(filename,sheet,'D11'); % Variable Mass Performance Studies
STABILITY_STUDY_Trim = xlsread(filename,sheet,'D12'); % Trim conditions - one case
STABILITY_STUDY_Trim_varV_m0 = xlsread(filename,sheet,'D13'); % Trim conditions - variation V and mass
STABILITY_STUDY_Trim_var_XCG = xlsread(filename,sheet,'D14'); % Trim conditions - variation XCG
STABILITY_STUDY_Long_dyn = xlsread(filename,sheet,'D15'); % Longitudinal Stability Analysis
STABILITY_STUDY_LatDir_dyn = xlsread(filename,sheet,'D16'); % Lateral-Directional Stability Analysis
STABILITY_STUDY_Trim_lat = xlsread(filename,sheet,'D17'); % Lateral Directional Trim
STABILITY_STUDY_Turning = xlsread(filename,sheet,'D18'); % Turning Stability
STABILITY_STUDY_Stability_Derivatives_varV_m0 = xlsread(filename,sheet,'D19'); % Plots Stability Derivatives
STABILITY_STUDY_Stability_Analysis_varV_m0_long = xlsread(filename,sheet,'D20'); % Longitudinal Dynamic Response Plots
STABILITY_STUDY_Stability_Analysis_varV_m0_lat = xlsread(filename,sheet,'D21'); % Lateral-Directional Dynamic Response Plots
% Stores the flags
STUDY_flags.AERODYNAMIC_STUDY = AERODYNAMIC_STUDY; % Conducts Aerodynamic Studies
STUDY_flags.PROP_STUDY = PROP_STUDY; % Conducts Proppeller Studies
STUDY_flags.PERFORMANCE_STUDY = PERFORMANCE_STUDY; % Conducts Performance Studies
STUDY_flags.STABILITY_STUDY = STABILITY_STUDY; % Conducts Stability Studies
STUDY_flags.STUDY_Weight_Fraction = STUDY_Weight_Fraction;% Weight Fraction Study
STUDY_flags.ANALYSIS_PERFORMANCE_AP = ANALYSIS_PERFORMANCE_AP;% Performance Analysis integrated with AP Codes
STUDY_flags.ANALYSIS_PERFORMANCE_AP_var = ANALYSIS_PERFORMANCE_AP_var;% Performance Analysis integrated with AP Codes varying Conditions
STUDY_flags.variable_speed_AP = variable_speed_AP;% Variable Speed Performance Studies
STUDY_flags.variable_weight_AP = variable_weight_AP;% Variable Mass Performance Studies
STUDY_flags.STABILITY_STUDY_Trim = STABILITY_STUDY_Trim;% Trim conditions - one case
STUDY_flags.STABILITY_STUDY_Trim_varV_m0 = STABILITY_STUDY_Trim_varV_m0;% Trim conditions - variation V and mass
STUDY_flags.STABILITY_STUDY_Trim_var_XCG = STABILITY_STUDY_Trim_var_XCG;% Trim conditions - variation XCG
STUDY_flags.STABILITY_STUDY_Long_dyn = STABILITY_STUDY_Long_dyn;% Longitudinal Stability Analysis
STUDY_flags.STABILITY_STUDY_LatDir_dyn = STABILITY_STUDY_LatDir_dyn;% Lateral-Directional Stability Analysis
STUDY_flags.STABILITY_STUDY_Trim_lat = STABILITY_STUDY_Trim_lat;% Lateral Directional Trim
STUDY_flags.STABILITY_STUDY_Turning = STABILITY_STUDY_Turning;% Turning Stability
STUDY_flags.STABILITY_STUDY_Stability_Derivatives_varV_m0 = STABILITY_STUDY_Stability_Derivatives_varV_m0;% Plots Stability Derivatives
STUDY_flags.STABILITY_STUDY_Stability_Analysis_varV_m0_long = STABILITY_STUDY_Stability_Analysis_varV_m0_long;% Longitudinal Dynamic Response Plots
STUDY_flags.STABILITY_STUDY_Stability_Analysis_varV_m0_lat = STABILITY_STUDY_Stability_Analysis_varV_m0_lat;% Lateral-Directional Dynamic Response Plots
% Gathers all the flags
OUTPUT_read_XLSX.STUDY_flags = STUDY_flags;

%% Aircraft Data
sheet = 3;
SF = xlsread(filename,sheet,'D3'); % Scaling Factor
Weight_Estimation = xlsread(filename,sheet,'D4'); % Weight Estimation
AC_type = xlsread(filename,sheet,'D5'); % Aircaft Type
CASE_fuse = xlsread(filename,sheet,'D6');% Determines the fuselage than will be shown
ESCALADO = xlsread(filename,sheet,'D7'); % Determine the type of scaling
case_AC = xlsread(filename,sheet,'D8'); % Aircraft Model Analized
% Stores the flags
AC_Data_flags.SF = SF; %
AC_Data_flags.Weight_Estimation = Weight_Estimation; %
AC_Data_flags.AC_type = AC_type; % 
AC_Data_flags.CASE_fuse = CASE_fuse; % 
AC_Data_flags.ESCALADO = ESCALADO;% 
AC_Data_flags.case_AC = case_AC;% 
% Gathers all the flags
OUTPUT_read_XLSX.AC_Data_flags = AC_Data_flags;

%% Propulsive Data
sheet = 4;
propul(1) = xlsread(filename,sheet,'D3'); % Type of engine
propul(2) = xlsread(filename,sheet,'D4'); % Number of engines
propul(3) = xlsread(filename,sheet,'D5'); % EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine
propul(4) = xlsread(filename,sheet,'D6'); % CONSUMO ESPECIFICO: % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
propul(5) = xlsread(filename,sheet,'D7'); % DERIVACION(TURBOFANES) By-pass 
propul(6) = xlsread(filename,sheet,'D8'); % EFICIENCIA DE LA HELICE (ETA_P)
propul(7) = xlsread(filename,sheet,'D9'); % Normativa Propulsion
propul(8) = xlsread(filename,sheet,'D10'); % CAPACIDAD CALORIFICA DEL COMBUSTIBLE
propul(9) = xlsread(filename,sheet,'D11'); % DIAMETRO DEL MOTOR(TURBOHELICES)
type_battery = xlsread(filename,sheet,'D12'); % Type Battery Used for Propulsion
Engine_loc = xlsread(filename,sheet,'D13'); % Engine location
Engine_conf = xlsread(filename,sheet,'D14'); % Engine configuration
D_prop = xlsread(filename,sheet,'D15'); % Preliminary Prop Diameter (m)
type_prop = xlsread(filename,sheet,'D16'); % Type of propller
Prop_type = xlsread(filename,sheet,'D17'); % Propulsive Model- Number of Prop used of propller
l_eng = xlsread(filename,sheet,'D18'); % Engine geometry - length
d_eng = xlsread(filename,sheet,'D19'); % Engine geometry - diameter
l_nc = xlsread(filename,sheet,'D20'); % Nacelle geometry - length
d_nc = xlsread(filename,sheet,'D21'); % Nacelle geometry - diameter
model_prop = xlsread(filename,sheet,'D22'); %	Compares 3 prop models
prop_selec_APC = xlsread(filename,sheet,'D23'); %	Model 1 - APC data
prop_selec_WT1 = xlsread(filename,sheet,'D24'); %	Model 2 - Wind tunnel data for different props
prop_selec_WT2 = xlsread(filename,sheet,'D25'); %	Model 3 - Wind tunnel data for different angle of attack
eta_gear = xlsread(filename,sheet,'D26'); % Gearbox efficiency
eta_m = xlsread(filename,sheet,'D27'); % Motor efficiency
eta_esc = xlsread(filename,sheet,'D28'); % ESC efficiency
eta_dist = xlsread(filename,sheet,'D29'); % Electrical distribution efficiency
SE_LiFePO4 = xlsread(filename,sheet,'D30'); % Specific Energy of LiFePO4 batteries 90–160 Wh/kg (320–580 J/g or kJ/kg) 
SE_LiPo = xlsread(filename,sheet,'D31'); % Specific Energy of LiPo batteries (Wh/kg)
SE_FuelCells = xlsread(filename,sheet,'D32'); % Wh/kg Specific Energy for fuel cells
SP_motor = xlsread(filename,sheet,'D33'); % Motor Specific Power kW/kg - TIER 1
SP_ESC = xlsread(filename,sheet,'D34'); % ESC Specific Power kW/kg - TIER 1
m_prop_din = xlsread(filename,sheet,'D35'); % Mass of propeller as a function of diameter in inches

% Stores the flags
Propulsive_flags.propul = propul; %
Propulsive_flags.type_battery = type_battery; %
Propulsive_flags.Engine_loc = Engine_loc; %
Propulsive_flags.Engine_conf = Engine_conf; %
Propulsive_flags.D_prop = D_prop; %
Propulsive_flags.type_prop = type_prop; %
Propulsive_flags.Prop_type = Prop_type; %
Propulsive_flags.l_eng = l_eng; %
Propulsive_flags.d_eng = d_eng; %
Propulsive_flags.l_nc = l_nc; %
Propulsive_flags.d_nc = d_nc; %
Propulsive_flags.model_prop = model_prop; %
Propulsive_flags.prop_selec_APC = prop_selec_APC; %
Propulsive_flags.prop_selec_WT1 = prop_selec_WT1; %
Propulsive_flags.prop_selec_WT2 = prop_selec_WT2; %
Propulsive_flags.eta_gear = eta_gear; %
Propulsive_flags.eta_m = eta_m; %
Propulsive_flags.eta_esc = eta_esc; %
Propulsive_flags.eta_dist = eta_dist; %
Propulsive_flags.SE_LiFePO4 = SE_LiFePO4; %
Propulsive_flags.SE_LiPo = SE_LiPo; %
Propulsive_flags.SE_FuelCells = SE_FuelCells; %
Propulsive_flags.SP_motor = SP_motor; %
Propulsive_flags.SP_ESC = SP_ESC; %
Propulsive_flags.m_prop_din = m_prop_din; %

% Gathers all the flags
OUTPUT_read_XLSX.Propulsive_flags = Propulsive_flags;

%% Weights Data
sheet = 5;
W_eng = xlsread(filename,sheet,'D3'); % Mass of engine (kg)
m_prop = xlsread(filename,sheet,'D4'); % Mass of proppeller (kg)
m_ESC = xlsread(filename,sheet,'D5'); % Mass of ESC (kg)
n_servos = xlsread(filename,sheet,'D6'); % Number of servos 
m_servo = xlsread(filename,sheet,'D7'); % Mass per servo (kg)
m_wiring = xlsread(filename,sheet,'D8'); % Mass of electronic wiring (kg)
m_metal = xlsread(filename,sheet,'D9'); % Mass of misc (kg)
m_systems = xlsread(filename,sheet,'D10'); % Mass of Systems (kg)
m_landing_gear = xlsread(filename,sheet,'D11'); % Mass of landing gear (kg)
n_pax = xlsread(filename,sheet,'D12'); % Number of Passangers
n_crew = xlsread(filename,sheet,'D13'); % Number of Crew
m_pax = xlsread(filename,sheet,'D14'); % Mass per Passenger (kg)
m_crew = xlsread(filename,sheet,'D15'); % Mass per Crew (kg)
m_luggage = xlsread(filename,sheet,'D16'); % Mass of luggage (kg)
m_cargo = xlsread(filename,sheet,'D17'); % Mass of payload cargo (kg)
m_batteries = xlsread(filename,sheet,'D18'); % Mass of batteries (kg)
m_fuel = xlsread(filename,sheet,'D19'); % Mass of fuel (kg) 
flag_landing_gear = xlsread(filename,sheet,'D20'); % Determination of Landing Gear Estimation if available
f_f = xlsread(filename,sheet,'D21'); % Fudge Factor that increments the weight estimation
m_cell_A123 = xlsread(filename,sheet,'D22'); % Mass of one A123 cell
m_w1 = xlsread(filename,sheet,'D23'); % Mass of one wing
m_HTP = xlsread(filename,sheet,'D24'); % Mass of HTP
m_VTP = xlsread(filename,sheet,'D25'); % Mass of VTP
m_Can = xlsread(filename,sheet,'D26'); % Mass of Canard
m_Vee = xlsread(filename,sheet,'D27'); % Mass of Vee Tail
m_fus_fairing = xlsread(filename,sheet,'D28'); % Mass of fuselage and fairing
m_prop_fairing = xlsread(filename,sheet,'D29'); % Mass of propeller nacelle
I_xx = xlsread(filename,sheet,'D30'); % Moments of Inertia Ixx
I_yy = xlsread(filename,sheet,'D31'); % Moments of Inertia Ixx
I_zz = xlsread(filename,sheet,'D32'); % Moments of Inertia Ixx
I_xy = xlsread(filename,sheet,'D33'); % Moments of Inertia Ixx
I_xz = xlsread(filename,sheet,'D34'); % Moments of Inertia Ixx
I_yz = xlsread(filename,sheet,'D35'); % Moments of Inertia Ixx
flag_moments_inertia = xlsread(filename,sheet,'D36'); % Determination of Moments of innertia from literature	

% Stores the flags
Weights_flags.W_eng = W_eng; %
Weights_flags.m_prop = m_prop; %
Weights_flags.m_ESC = m_ESC; %
Weights_flags.n_servos = n_servos; %
Weights_flags.m_servo = m_servo; %
Weights_flags.m_wiring = m_wiring; %
Weights_flags.m_metal = m_metal; %
Weights_flags.m_systems = m_systems; %
Weights_flags.m_landing_gear = m_landing_gear; %
Weights_flags.n_pax = n_pax; %
Weights_flags.n_crew = n_crew; %
Weights_flags.m_pax = m_pax; %
Weights_flags.m_crew = m_crew; %
Weights_flags.m_luggage = m_luggage; %
Weights_flags.m_cargo = m_cargo; %
Weights_flags.m_batteries = m_batteries; %
Weights_flags.m_fuel = m_fuel; %
Weights_flags.flag_landing_gear = flag_landing_gear; %
Weights_flags.f_f = f_f; %
Weights_flags.m_cell_A123 = m_cell_A123; %
Weights_flags.m_w1 = m_w1; %
Weights_flags.m_HTP = m_HTP; %
Weights_flags.m_VTP = m_VTP; %
Weights_flags.m_Can = m_Can; %
Weights_flags.m_Vee = m_Vee; %
Weights_flags.m_fus_fairing = m_fus_fairing; %
Weights_flags.m_prop_fairing = m_prop_fairing; %
Weights_flags.I_xx = I_xx; %
Weights_flags.I_yy = I_yy; %
Weights_flags.I_zz = I_zz; %
Weights_flags.I_xy = I_xy; %
Weights_flags.I_xz = I_xz; %
Weights_flags.I_yz = I_yz; %
Weights_flags.flag_moments_inertia = flag_moments_inertia; %

% Gathers all the flags
OUTPUT_read_XLSX.Weights_flags = Weights_flags;


%% Fuselage Data
sheet = 6;
nSections = xlsread(filename,sheet,'D3'); % Number of sectins (x coordinate) (eliminated 1 for convergence)
nPoints = xlsread(filename,sheet,'D4'); % number of elements per section (y and z coordinates)
lecture_Geo = xlsread(filename,sheet,'D5'); % lineas a partir desde donde empieza a leer
STL_PLOT = xlsread(filename,sheet,'D6'); % Generates Fuselage mesh from CAD STL
[a,XFLR5_file] = xlsread(filename,sheet,'D7'); % Defines flies to be used for the Fuselage geometry - XFLR5
[a,STL_file] = xlsread(filename,sheet,'D8'); % Defines flies to be used for the Fuselage geometry - STL

% Stores the flags
Fuselage_flags.nSections = nSections; %
Fuselage_flags.nPoints = nPoints; %
Fuselage_flags.lecture_Geo = lecture_Geo; %
Fuselage_flags.STL_PLOT = STL_PLOT; %
Fuselage_flags.XFLR5_file = XFLR5_file; %
Fuselage_flags.STL_file = STL_file; %
% Gathers all the flags
OUTPUT_read_XLSX.Fuselage_flags = Fuselage_flags;

%% Performance Data Prelimina
sheet = 7;
h_climb = xlsread(filename,sheet,'D3'); % Preliminar Performance conditions: Climb altitude (m)
Endurance_v = xlsread(filename,sheet,'D4'); % Preliminar Performance conditions: Endurance in Hoovering (min)
h = xlsread(filename,sheet,'D5'); % Preliminar Performance conditions: Cruise initial altitude (m)
V = xlsread(filename,sheet,'D6'); % Preliminar Performance conditions: Cruise initial airspeed (m/s)
V_max = xlsread(filename,sheet,'D7'); % Preliminar Performance conditions: Cruise max speed
Range = xlsread(filename,sheet,'D8'); % Preliminar Performance conditions: Cruise Range (m)
Endurance = xlsread(filename,sheet,'D9'); % Preliminar Performance conditions: Cruise Endurance (min)
Flight_SF = xlsread(filename,sheet,'D10'); % Flight Safety Margin - normal 

% Stores the flags
Performance_pre_flags.h_climb = h_climb; %
Performance_pre_flags.Endurance_v = Endurance_v; %
Performance_pre_flags.h = h; %
Performance_pre_flags.V = V; %
Performance_pre_flags.V_max = V_max; %
Performance_pre_flags.Range = Range; %
Performance_pre_flags.Endurance = Endurance; %
Performance_pre_flags.Flight_SF = Flight_SF; %
% Gathers all the flags
OUTPUT_read_XLSX.Performance_pre_flags = Performance_pre_flags;

%% Aerodynamic Data
sheet = 8;
index_w1 = xlsread(filename,sheet,'D3'); % Selection of TXT that are used for the aerodynamic analysis (w1)
index_w2 = xlsread(filename,sheet,'D4'); % Selection of TXT that are used for the aerodynamic analysis (w2)
index_w3 = xlsread(filename,sheet,'D5'); % Selection of TXT that are used for the aerodynamic analysis (w3)
i_w1 = xlsread(filename,sheet,'D6'); % Selects the AoA of each surface relative to fuselage line (w1)
i_w2 = xlsread(filename,sheet,'D7'); % Selects the AoA of each surface relative to fuselage line (w2)
i_w3 = xlsread(filename,sheet,'D8'); % Selects the AoA of each surface relative to fuselage line (w3)
w1 = xlsread(filename,sheet,'D9');% Selection of elements for theoretical polar estimation - wing	w1
h = xlsread(filename,sheet,'D10');% Selection of elements for theoretical polar estimation - HTP	h
v = xlsread(filename,sheet,'D11');% Selection of elements for theoretical polar estimation - VTP	v
v2 = xlsread(filename,sheet,'D12');% Selection of elements for theoretical polar estimation - twin VTP	v2
vtail = xlsread(filename,sheet,'D13');% Selection of elements for theoretical polar estimation - Vtail	vtail
can = xlsread(filename,sheet,'D14');% Selection of elements for theoretical polar estimation - canard	can
fus = xlsread(filename,sheet,'D15');% Selection of elements for theoretical polar estimation - fuselage	fus
m_fus = xlsread(filename,sheet,'D16');% Selection of elements for theoretical polar estimation - multiple fuselage	m_fus
n_m_fus = xlsread(filename,sheet,'D17');% Selection of elements for theoretical polar estimation - number of multiple fuselage	n_m_fus
nac = xlsread(filename,sheet,'D18');% Selection of elements for theoretical polar estimation - nacelle	nac
landgear = xlsread(filename,sheet,'D19');% Selection of elements for theoretical polar estimation - landing gear	landgear

% Stores the flags
Aerodynamic_Data_flags.index_w1 = index_w1; %
Aerodynamic_Data_flags.index_w2 = index_w2; %
Aerodynamic_Data_flags.index_w3 = index_w3; %
Aerodynamic_Data_flags.i_w1 = i_w1; %
Aerodynamic_Data_flags.i_w2 = i_w2; %
Aerodynamic_Data_flags.i_w3 = i_w3; %
Aerodynamic_Data_flags.Conf.w1 = w1; %
Aerodynamic_Data_flags.Conf.h = h; %
Aerodynamic_Data_flags.Conf.v = v; %
Aerodynamic_Data_flags.Conf.v2 = v2; %
Aerodynamic_Data_flags.Conf.vtail = vtail; %
Aerodynamic_Data_flags.Conf.can = can; %
Aerodynamic_Data_flags.Conf.fus = fus; %
Aerodynamic_Data_flags.Conf.m_fus = m_fus; %
Aerodynamic_Data_flags.Conf.n_m_fus = n_m_fus; %
Aerodynamic_Data_flags.Conf.nac = nac; %
Aerodynamic_Data_flags.Conf.landgear = landgear; %

% Gathers all the flags
OUTPUT_read_XLSX.Aerodynamic_Data_flags = Aerodynamic_Data_flags;

%% Stability Data
sheet = 9;
SM = xlsread(filename,sheet,'D3'); % Static Margin
prop_wash_effect = xlsread(filename,sheet,'D4'); % Prop Wash Effect Considered on the dynamic pressure calculation
XCG_FF = xlsread(filename,sheet,'D5'); % elects XCG depending if it's from desired stability conditions in Forward Flight
StabilityModel = xlsread(filename,sheet,'D6'); % elects XCG depending if it's from desired stability conditions in Forward Flight
only_trim = xlsread(filename,sheet,'D7'); % elects XCG depending if it's from desired stability conditions in Forward Flight
V_low = xlsread(filename,sheet,'D8'); %	Trim Conditions Velocity Range - V low
V_high = xlsread(filename,sheet,'D9'); %	Trim Conditions Velocity Range - V high
N_V_VAR = xlsread(filename,sheet,'D10'); %	Number of iteration variation Speed
N_m_VAR = xlsread(filename,sheet,'D11'); %	Number of iteration variation Weight
beta = xlsread(filename,sheet,'D12'); %	Trim Lateral Stability Conditions -side slip angle (beta)
beta_i = xlsread(filename,sheet,'D13'); %	Trim Lateral Stability Conditions -Variable side slip angle (initial beta)
beta_f = xlsread(filename,sheet,'D14'); %	Trim Lateral Stability Conditions -Variable side slip angle (final beta)
N_Delta_beta = xlsread(filename,sheet,'D15'); %	Trim Lateral Stability Conditions -Variable side slip angle (Number of beta's)
phi = xlsread(filename,sheet,'D16'); %	Turning Condition Stability Conditions -bank angle (beta)
phi_i = xlsread(filename,sheet,'D17'); %	Turning Condition Stability Conditions -Variable bank angle (initial phi)
phi_f = xlsread(filename,sheet,'D18'); %	Turning Condition Stability Conditions -Variable bank angle (final phi)
N_Delta_phi = xlsread(filename,sheet,'D19'); %	Turning Condition Stability Conditions -Variable bank angle (Number of beta's)
n_viraje = xlsread(filename,sheet,'D20'); %	Loading factur during turning flight 
Munk_fuselage_constribution = xlsread(filename,sheet,'D21'); %	Include Munk's fuselage contribution todetermine Xac, & desired Xcg location for an SM  

% Stores the flags
Stability_flags.SM_des = SM/100; %
Stability_flags.prop_wash_effect = prop_wash_effect; %
Stability_flags.XCG_FF = XCG_FF; %
Stability_flags.StabilityModel = StabilityModel; %
Stability_flags.only_trim = only_trim; %
Stability_flags.V_low = V_low; %
Stability_flags.V_high = V_high; %
Stability_flags.N_V_VAR = N_V_VAR; %
Stability_flags.N_m_VAR = N_m_VAR; %
Stability_flags.beta = beta;
Stability_flags.beta_i = beta_i;
Stability_flags.beta_f = beta_f;
Stability_flags.N_Delta_beta = N_Delta_beta;
Stability_flags.phi = phi;
Stability_flags.phi_i = phi_i;
Stability_flags.phi_f = phi_f;
Stability_flags.N_Delta_phi = N_Delta_phi;
Stability_flags.n_viraje = n_viraje;
Stability_flags.Munk_fuselage_constribution = Munk_fuselage_constribution;

% Gathers all the flags
OUTPUT_read_XLSX.Stability_flags = Stability_flags;

%% Geometry Data
sheet = 10;

CASE = xlsread(filename,sheet,'D3'); % Model of Aircraft
Control_surface = xlsread(filename,sheet,'D4'); % Determines the ammount of control surface in the w2 aileron & elevator
AC_type = xlsread(filename,sheet,'D5'); % Aircraft type

% Stores the flags
Geometry_Data_flags.CASE = CASE; %
Geometry_Data_flags.Control_surface = Control_surface; %
Geometry_Data_flags.AC_type = AC_type; %

% Gathers all the flags
OUTPUT_read_XLSX.Geometry_Data_flags = Geometry_Data_flags;

%% Input Aircraft Geometry Data
sheet = 11;

x_offset_CAD = xlsread(filename,sheet,'D3'); % Distance from origin in CATIa to Fuselage (offset value) - x
z_offset_CAD = xlsread(filename,sheet,'D4'); % Distance from origin in CATIa to Fuselage (offset value) - z
w_fus = xlsread(filename,sheet,'D5'); % Fuselage geomtry (Units in m) - approx - width
h_fus = xlsread(filename,sheet,'D6'); % Fuselage geomtry (Units in m) - approx - heigth
l_fus = xlsread(filename,sheet,'D7'); % Fuselage geomtry (Units in m) - approx - length
d_fus = xlsread(filename,sheet,'D8'); % Fuselage geomtry (Units in m) - approx - diameter
%w1
y_loc_1R_y1_w1_CAD = xlsread(filename,sheet,'D9'); % y loc of wing (w1) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_w1_CAD = xlsread(filename,sheet,'D10'); %	y loc of wing (w1) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_w1_CAD = xlsread(filename,sheet,'D11'); %	z loc of wing (w1) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_w1_CAD = xlsread(filename,sheet,'D12'); %	x loc of wing (w1) root chord LE position (distance from CAD refference point)
Lambda_LE_w1_e = xlsread(filename,sheet,'D13'); %	Sweep of w1 (deg)
dihedral_w1_e = xlsread(filename,sheet,'D14'); %	Dihedral of w1 (deg)
%w2
y_loc_1R_y1_w2_CAD = xlsread(filename,sheet,'D15'); %	y loc of wing (w2) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_w2_CAD = xlsread(filename,sheet,'D16'); %	y loc of wing (w2) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_w2_CAD = xlsread(filename,sheet,'D17'); %	z loc of wing (w2) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_w2_CAD = xlsread(filename,sheet,'D18'); %	x loc of wing (w2) root chord LE position (distance from CAD refference point)
Lambda_LE_w2_e = xlsread(filename,sheet,'D19'); %	Sweep of w2 (deg)
dihedral_w2_e = xlsread(filename,sheet,'D20'); %	Dihedral of w2 (deg)
%can
y_loc_1R_y1_can_CAD = xlsread(filename,sheet,'D21'); %	y loc of wing (can) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_can_CAD = xlsread(filename,sheet,'D22'); %	y loc of wing (can) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_can_CAD = xlsread(filename,sheet,'D23'); %	z loc of wing (can) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_can_CAD = xlsread(filename,sheet,'D24'); %	x loc of wing (can) root chord LE position (distance from CAD refference point)
Lambda_LE_can_e = xlsread(filename,sheet,'D25'); %	Sweep of can (deg)
dihedral_can_e = xlsread(filename,sheet,'D26'); %	Dihedral of can (deg)
%HTP
y_loc_1R_y1_HTP_CAD = xlsread(filename,sheet,'D27'); %	y loc of wing (HTP) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_HTP_CAD = xlsread(filename,sheet,'D28'); %	y loc of wing (HTP) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_HTP_CAD = xlsread(filename,sheet,'D29'); %	z loc of wing (HTP) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_HTP_CAD = xlsread(filename,sheet,'D30'); %	x loc of wing (HTP) root chord LE position (distance from CAD refference point)
Lambda_LE_HTP_e = xlsread(filename,sheet,'D31'); %	Sweep of HTP (deg)
dihedral_HTP_e = xlsread(filename,sheet,'D32'); %	Dihedral of HTP (deg)
%VTP
y_loc_1R_y1_VTP_CAD = xlsread(filename,sheet,'D33'); %	y loc of wing (VTP) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_VTP_CAD = xlsread(filename,sheet,'D34'); %	y loc of wing (VTP) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_VTP_CAD = xlsread(filename,sheet,'D35'); %	z loc of wing (VTP) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_VTP_CAD = xlsread(filename,sheet,'D36'); %	x loc of wing (VTP) root chord LE position (distance from CAD refference point)
Lambda_LE_VTP_e = xlsread(filename,sheet,'D37'); %	Sweep of VTP (deg)
dihedral_VTP_e = xlsread(filename,sheet,'D38'); %	Dihedral of VTP (deg)
%VTAIL
y_loc_1R_y1_vee_CAD = xlsread(filename,sheet,'D39'); %	y loc of wing (V-tail) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_vee_CAD = xlsread(filename,sheet,'D40'); %	y loc of wing (V-tail) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_vee_CAD = xlsread(filename,sheet,'D41'); %	z loc of wing (V-tail) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_vee_CAD = xlsread(filename,sheet,'D42'); %	x loc of wing (V-tail) root chord LE position (distance from CAD refference point)
Lambda_LE_vee_e = xlsread(filename,sheet,'D43'); %	Sweep of V-tail (deg)
dihedral_vee_e = xlsread(filename,sheet,'D44'); %	Dihedral of V-tail (deg)
%Twin VTP
y_loc_1R_y1_2VTP_CAD = xlsread(filename,sheet,'D45'); %	y loc of Twin (VTP) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_2VTP_CAD = xlsread(filename,sheet,'D46'); %	y loc of Twin (VTP) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_2VTP_CAD = xlsread(filename,sheet,'D47'); %	z loc of wing (VTP) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_2VTP_CAD = xlsread(filename,sheet,'D48'); %	x loc of Twin (VTP) root chord LE position (distance from CAD refference point)
Lambda_LE_2VTP_e = xlsread(filename,sheet,'D49'); %	Sweep of Twin VTP (deg)
dihedral_2VTP_e = xlsread(filename,sheet,'D50'); %	Dihedral of Twin VTP (deg)
cR_w1 = xlsread(filename,sheet,'D51'); %	Root Chord w1
cT_w1 = xlsread(filename,sheet,'D52'); %	Tip Chord w1
cR_w2 = xlsread(filename,sheet,'D53'); %	Root Chord w2
cT_w2 = xlsread(filename,sheet,'D54'); %	Tip Chord w2
cR_can = xlsread(filename,sheet,'D55'); %	Root Chord canard
cT_can = xlsread(filename,sheet,'D56'); %	Tip Chord canard
cR_HTP = xlsread(filename,sheet,'D57'); %	Root Chord HTP
cT_HTP = xlsread(filename,sheet,'D58'); %	Tip Chord HTP
cR_VTP = xlsread(filename,sheet,'D59'); %	Root Chord VTP
cT_VTP = xlsread(filename,sheet,'D60'); %	Tip Chord VTP
cR_vee = xlsread(filename,sheet,'D61'); %	Root Chord vee
cT_vee = xlsread(filename,sheet,'D62'); %	Tip Chord vee
cR_2VTP = xlsread(filename,sheet,'D63'); %	Root Chord Twin VTP
cT_2VTP = xlsread(filename,sheet,'D64'); %	Tip Chord Twin VTP
% Perfecntage (%) of the control surface
K_y1_ail_w1 = xlsread(filename,sheet,'D65'); %	% of aileron from effective wing surface - inner position
K_y2_ail_w1 = xlsread(filename,sheet,'D66'); %	% of aileron from effective wing surface - outter position
K_y1_ele_w2 = xlsread(filename,sheet,'D67'); %	% of elevator from effective HTP surface - inner position
K_y2_ele_w2 = xlsread(filename,sheet,'D68'); %	% of elevator from effective HTP surface  - outter position
K_y1_elevon_w1 = xlsread(filename,sheet,'D69'); %	% of elevon from effective wing surface - inner position
K_y2_elevon_w1 = xlsread(filename,sheet,'D70'); %	% of elevon from effective wing surface - outter position
K_y1_flap_w1 = xlsread(filename,sheet,'D71'); %	% of flap from effective wing surface - inner position
K_y2_flap_w1 = xlsread(filename,sheet,'D72'); %	% of flap from effective wing surface - outter position
K_y1_rudder_VTP = xlsread(filename,sheet,'D73'); %	% of rudder from effective VTP surface - inner position
K_y2_rudder_VTP = xlsread(filename,sheet,'D74'); %	% of rudder from effective VTP surface  - outter position
K_y1_rudvtr_w2 = xlsread(filename,sheet,'D75'); %	% of ruddervator from effective V-tail surface - inner position
K_y2_rudvtr_w2 = xlsread(filename,sheet,'D76'); %	% of ruddervator from effective V-tail surface  - outter position
K_y1_can_w2 = xlsread(filename,sheet,'D77'); %	% of canard control surface from effective canard surface - inner position
K_y2_can_w2 = xlsread(filename,sheet,'D78'); %	% of canard control surface from effective canard surface - outter position
% Maximum and minimum control surface deflections
delta_ail_min = xlsread(filename,sheet,'D79'); % 	Minimum aileron deflection
delta_ail_max = xlsread(filename,sheet,'D80'); %	Maximum aileron deflection
delta_ele_min = xlsread(filename,sheet,'D81'); %	Minimum elevator deflection
delta_ele_max = xlsread(filename,sheet,'D82'); %	Maximum elevator deflection
delta_elevon_min = xlsread(filename,sheet,'D83'); %	Minimum elevon deflection
delta_elevon_max = xlsread(filename,sheet,'D84'); %	Maximum elevon deflection
delta_flap_min = xlsread(filename,sheet,'D85'); %	Minimum flap deflection
delta_flap_max = xlsread(filename,sheet,'D86'); %	Maximum flap deflection
delta_rudder_min = xlsread(filename,sheet,'D87'); %	Minimum rudder deflection
delta_rudder_max = xlsread(filename,sheet,'D88'); %	Maximum rudder deflection
delta_rudvtr_min = xlsread(filename,sheet,'D89'); %	Minimum ruddervator deflection
delta_rudvtr_max = xlsread(filename,sheet,'D90'); %	Maximum ruddervator deflection
delta_can_min = xlsread(filename,sheet,'D91'); %	Minimum canard deflection
delta_can_max = xlsread(filename,sheet,'D92'); %	Maximum canard deflection
% Control surface dimensions
cf_ail = xlsread(filename,sheet,'D93'); %	%  of control surface aileron (chrodwise)
t_c_ail = xlsread(filename,sheet,'D94'); %	Thinckness 2 chord ratio associated to the airfoil - aileron
cf_ele = xlsread(filename,sheet,'D95'); %	%  of control surface elevator (chrodwise)
t_c_ele = xlsread(filename,sheet,'D96'); %	Thinckness 2 chord ratio associated to the airfoil - elevator
cf_elevon = xlsread(filename,sheet,'D97'); %	%  of control surface elevon (chrodwise)
t_c_elevon = xlsread(filename,sheet,'D98'); %	Thinckness 2 chord ratio associated to the airfoil - elevon
cf_flap = xlsread(filename,sheet,'D99'); %	%  of control surface flap (chrodwise)
t_c_flap = xlsread(filename,sheet,'D100'); %	Thinckness 2 chord ratio associated to the airfoil - flap
cf_rudder = xlsread(filename,sheet,'D101'); %	%  of control surface rudder (chrodwise)
t_c_rudder = xlsread(filename,sheet,'D102'); %	Thinckness 2 chord ratio associated to the airfoil - rudder
cf_rudvtr = xlsread(filename,sheet,'D103'); %	%  of control surface ruddervator (chrodwise)
t_c_rudvtr = xlsread(filename,sheet,'D104'); %	Thinckness 2 chord ratio associated to the airfoil - ruddervator
cf_can = xlsread(filename,sheet,'D105'); %	%  of control surface canard (chrodwise)
t_c_can = xlsread(filename,sheet,'D106'); %	Thinckness 2 chord ratio associated to the airfoil - canard
% Initial center of gravity estimaton
x_XCG = xlsread(filename,sheet,'D107');
y_XCG = xlsread(filename,sheet,'D108');
z_XCG = xlsread(filename,sheet,'D109');

% Stores the flags
InputGeometry_Data_flags.x_offset_CAD = x_offset_CAD; %
InputGeometry_Data_flags.z_offset_CAD = z_offset_CAD; %
InputGeometry_Data_flags.w_fus = w_fus; %
InputGeometry_Data_flags.h_fus = h_fus; %
InputGeometry_Data_flags.l_fus = l_fus; %
InputGeometry_Data_flags.d_fus = d_fus; %
InputGeometry_Data_flags.d_fus = d_fus; %
InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD = y_loc_1R_y1_w1_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_w1_CAD = z_loc_1R_y1_w1_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; %
InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; %
InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; %
InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = y_loc_1R_y1_w2_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = y_loc_1R_y2_w2_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_w2_CAD = z_loc_1R_y1_w2_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD = x_loc_1R_y1_w2_CAD; %
InputGeometry_Data_flags.Lambda_LE_w2_e = Lambda_LE_w2_e; %
InputGeometry_Data_flags.dihedral_w2_e = dihedral_w2_e; %
InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = y_loc_1R_y1_can_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_can_CAD = y_loc_1R_y2_can_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_can_CAD = z_loc_1R_y1_can_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_can_CAD = x_loc_1R_y1_can_CAD; %
InputGeometry_Data_flags.Lambda_LE_can_e = Lambda_LE_can_e; %
InputGeometry_Data_flags.dihedral_can_e = dihedral_can_e; %
InputGeometry_Data_flags.y_loc_1R_y1_HTP_CAD = y_loc_1R_y1_HTP_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_HTP_CAD = y_loc_1R_y2_HTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_HTP_CAD = z_loc_1R_y1_HTP_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_HTP_CAD = x_loc_1R_y1_HTP_CAD; %
InputGeometry_Data_flags.Lambda_LE_HTP_e = Lambda_LE_HTP_e; %
InputGeometry_Data_flags.dihedral_HTP_e = dihedral_HTP_e; %
InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_VTP_CAD = y_loc_1R_y2_VTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; %
InputGeometry_Data_flags.Lambda_LE_VTP_e = Lambda_LE_VTP_e; %
InputGeometry_Data_flags.dihedral_VTP_e = dihedral_VTP_e; %
InputGeometry_Data_flags.y_loc_1R_y1_vee_CAD = y_loc_1R_y1_vee_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_vee_CAD = y_loc_1R_y2_vee_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_vee_CAD = z_loc_1R_y1_vee_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_vee_CAD = x_loc_1R_y1_vee_CAD; %
InputGeometry_Data_flags.Lambda_LE_vee_e = Lambda_LE_vee_e; %
InputGeometry_Data_flags.dihedral_vee_e = dihedral_vee_e; %
InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_2VTP_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_2VTP_CAD = y_loc_1R_y2_2VTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_2VTP_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_2VTP_CAD; %
InputGeometry_Data_flags.Lambda_LE_2VTP_e = Lambda_LE_2VTP_e; %
InputGeometry_Data_flags.dihedral_2VTP_e = dihedral_2VTP_e; %
InputGeometry_Data_flags.cR_w1 = cR_w1; %
InputGeometry_Data_flags.cT_w1 = cT_w1; %
InputGeometry_Data_flags.cR_w2 = cR_w2; %
InputGeometry_Data_flags.cT_w2 = cT_w2; %
InputGeometry_Data_flags.cR_can = cR_can; %
InputGeometry_Data_flags.cT_can = cT_can; %
InputGeometry_Data_flags.cR_HTP = cR_HTP; %
InputGeometry_Data_flags.cT_HTP = cT_HTP; %
InputGeometry_Data_flags.cR_VTP = cR_VTP; %
InputGeometry_Data_flags.cT_VTP = cT_VTP; %
InputGeometry_Data_flags.cR_vee = cR_vee; %
InputGeometry_Data_flags.cT_vee = cT_vee; %
InputGeometry_Data_flags.cR_2VTP = cR_2VTP; %
InputGeometry_Data_flags.cT_2VTP = cT_2VTP; %
InputGeometry_Data_flags.K_y1_ail_w1 = K_y1_ail_w1;
InputGeometry_Data_flags.K_y2_ail_w1 = K_y2_ail_w1;
InputGeometry_Data_flags.K_y1_ele_w2 = K_y1_ele_w2;
InputGeometry_Data_flags.K_y1_ele_w2 = K_y2_ele_w2;
InputGeometry_Data_flags.K_y1_elevon_w1 = K_y1_elevon_w1;
InputGeometry_Data_flags.K_y2_elevon_w1 = K_y2_elevon_w1;
InputGeometry_Data_flags.K_y1_flap_w1 = K_y1_flap_w1;
InputGeometry_Data_flags.K_y2_flap_w1 = K_y2_flap_w1;
InputGeometry_Data_flags.K_y1_rudder_VTP = K_y1_rudder_VTP;
InputGeometry_Data_flags.K_y2_rudder_VTP = K_y2_rudder_VTP;
InputGeometry_Data_flags.K_y1_rudvtr_w2 = K_y1_rudvtr_w2;
InputGeometry_Data_flags.K_y2_rudvtr_w2 = K_y2_rudvtr_w2;
InputGeometry_Data_flags.K_y1_can_w2 = K_y1_can_w2;
InputGeometry_Data_flags.K_y2_can_w2 = K_y2_can_w2;
InputGeometry_Data_flags.delta_ail_min = delta_ail_min; 
InputGeometry_Data_flags.delta_ail_max = delta_ail_max;
InputGeometry_Data_flags.delta_ele_min = delta_ele_min;
InputGeometry_Data_flags.delta_ele_max = delta_ele_max;
InputGeometry_Data_flags.delta_elevon_min = delta_elevon_min;
InputGeometry_Data_flags.delta_elevon_max = delta_elevon_max;
InputGeometry_Data_flags.delta_flap_min = delta_flap_min;
InputGeometry_Data_flags.delta_flap_max = delta_flap_max;
InputGeometry_Data_flags.delta_rudder_min = delta_rudder_min;
InputGeometry_Data_flags.delta_rudder_max = delta_rudder_max;
InputGeometry_Data_flags.delta_rudvtr_min = delta_rudvtr_min;
InputGeometry_Data_flags.delta_rudvtr_max = delta_rudvtr_max;
InputGeometry_Data_flags.delta_can_min = delta_can_min;
InputGeometry_Data_flags.delta_can_max = delta_can_max;
InputGeometry_Data_flags.cf_ail = cf_ail;
InputGeometry_Data_flags.t_c_ail = t_c_ail;
InputGeometry_Data_flags.cf_ele = cf_ele;
InputGeometry_Data_flags.t_c_ele = t_c_ele;
InputGeometry_Data_flags.cf_elevon = cf_elevon;
InputGeometry_Data_flags.t_c_elevon = t_c_elevon;
InputGeometry_Data_flags.cf_flap = cf_flap;
InputGeometry_Data_flags.t_c_flap = t_c_flap;
InputGeometry_Data_flags.cf_rudder = cf_rudder;
InputGeometry_Data_flags.t_c_rudder = t_c_rudder;
InputGeometry_Data_flags.cf_rudvtr = cf_rudvtr;
InputGeometry_Data_flags.t_c_rudvtr = t_c_rudvtr;
InputGeometry_Data_flags.cf_can = cf_can;
InputGeometry_Data_flags.t_c_can = t_c_can;
InputGeometry_Data_flags.x_XCG = x_XCG;
InputGeometry_Data_flags.y_XCG = y_XCG;
InputGeometry_Data_flags.z_XCG = z_XCG;

% Gathers all the flags
OUTPUT_read_XLSX.InputGeometry_Data_flags = InputGeometry_Data_flags;

%% PLOTS
sheet = 12;

plot(1) = xlsread(filename,sheet,'D3');% Prints plots - Aero
plot(2) = xlsread(filename,sheet,'D4');% Prints plots - Aero Polar
plot(3) = xlsread(filename,sheet,'D5');% Print Plots for Stimation ofXAC
plot(4) = xlsread(filename,sheet,'D6');% Print Plots for Propulsive Models
plot(5) = xlsread(filename,sheet,'D7');% Prints plots for 3D
plot(6) = xlsread(filename,sheet,'D8');% Prints plots of SM analysis
plot(7) = xlsread(filename,sheet,'D9');% Prints plots of longitudinal Trim
plot(8) = xlsread(filename,sheet,'D10');% Prints plots of longitudinal Trim with Variable mass & Variable Speed
plot(9) = xlsread(filename,sheet,'D11');% Prints plots of lateral Trim
plot(10) = xlsread(filename,sheet,'D12');% Prints plots of lateral Turning
plot(11) = xlsread(filename,sheet,'D13');% Prints PLOTS STABILITY DERIVATIVES FOR VAR MASS AND VELOCITY
plot(12) = xlsread(filename,sheet,'D14');% Prints PLOTS STABILITY ANALYSIS FOR VAR MASS AND VELOCITY
plot(13) = xlsread(filename,sheet,'D15');% Prints PLOTS PERFORMANCE STUDY
[a,prefix] = xlsread(filename,sheet,'D16');% Selection of String Characters that Define the AC selection

% Stores the flags
PLOT_flags.plot = plot; %
PLOT_flags.prefix = prefix; %

% Gathers all the flags
OUTPUT_read_XLSX.PLOT_flags = PLOT_flags;

%% Performance Inpute Preliminar
sheet = 13;
% Taxy
temp_local_taxy = xlsread(filename,sheet,'D3');% 	TAXY: TEMPERATURA LOCAL (Celsius)
h_inicial_taxy = xlsread(filename,sheet,'D4');% 	TAXY: INITIAL ALTITUDE
P_local_taxy = xlsread(filename,sheet,'D5');% 	TAXY: Local Pressure
delta_T_taxy = xlsread(filename,sheet,'D6');% 	TAXY: PALANCA DE RALENTI EN TAXI = 0.05
V_taxy = xlsread(filename,sheet,'D7');% 	TAXY:  VELOCIDAD A LA QUE HACE EL TAXI (m/s)
t_taxy = xlsread(filename,sheet,'D8');% 	TAXY: TIEMPO DE ESPERA EN TAXI (s)
% TakeOff
temp_local_TO = xlsread(filename,sheet,'D9');% 	TAKEOFF: TEMPERATURA LOCAL (Celsius)
h_inicial_TO = xlsread(filename,sheet,'D10');% 	TAKEOFF: ALTURA LOCAL (m)
P_local_TO = xlsread(filename,sheet,'D11');% 	TAKEOFF: PRESION LOCAL (Pa)
mu_TO = xlsread(filename,sheet,'D12');% 	TAKEOFF: COEFICIENTE DE FRICCION CON LA PISTA (MU)
h_obstacle_TO = xlsread(filename,sheet,'D13');% 	TAKEOFF: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
gamma_climb_TO = xlsread(filename,sheet,'D14');% 	TAKEOFF: GAMMA DE SUBIDA MINIMO
delta_T_TO = xlsread(filename,sheet,'D15');% 	TAKEOFF: PALANCA DE GASES PARA DESPEGUE
% Climb
h_inicial_cl = xlsread(filename,sheet,'D16');% 	CLIMB : ALTURA INICIAL - [m]
h_final_cl = xlsread(filename,sheet,'D17');% 	CLIMB : ALTURA FINAL - [m]
gamma_cl = xlsread(filename,sheet,'D18');% 	CLIMB : GAMMA DE SUBIDA - [-]
Mach_cl = xlsread(filename,sheet,'D19');% 	CLIMB : MACH DE VUELO  - [-]
TAS_cl = xlsread(filename,sheet,'D20');% 	CLIMB : VELOCIDAD TAS  - [m/s]
EAS_cl = xlsread(filename,sheet,'D21');% 	CLIMB : VELOCIDAD EAS  - [m/s]
delta_T_cl = xlsread(filename,sheet,'D22');%	CLIMB : PALANCA DE GASES  - [-]
V_ini_cl = xlsread(filename,sheet,'D23');% 	CLIMB : VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
V_fin_cl = xlsread(filename,sheet,'D24');% 	CLIMB : VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
% Cruise
h_inicial_cr = xlsread(filename,sheet,'D25');%	CRUISE: ALTURA INICIAL
dist_final_cr = xlsread(filename,sheet,'D26');%	CRUISE: DISTANCIA FINAL
V_cr = xlsread(filename,sheet,'D27');%	CRUISE: VELOCIDAD
delta_T_cr = xlsread(filename,sheet,'D28');%	CRUISE:  PALANCA DE GASES
V_ini_cr = xlsread(filename,sheet,'D29');%	CRUISE:  VELOCIDAD INICIAL
V_fin_cr = xlsread(filename,sheet,'D30');%	CRUISE: VELOCIDAD FINAL
fuel_cr = xlsread(filename,sheet,'D31');% - [kg] % 7: COMBUSTIBLE A QUEMAR
Cd0_cr = xlsread(filename,sheet,'D32');% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
k1_cr = xlsread(filename,sheet,'D33');% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
k2_cr = xlsread(filename,sheet,'D34');% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
% Turn
h_inicial_tr = xlsread(filename,sheet,'D35');%	TURN:  ALTURA INICIAL
t_final_tr = xlsread(filename,sheet,'D36');%	TURN: TIEMPO FINAL (seg)
V_turn  = xlsread(filename,sheet,'D37');% TURN: turn velocity
delta_T_tr = xlsread(filename,sheet,'D38');%	TURN: PALANCA DE GASES
phi_tr = xlsread(filename,sheet,'D39');% 	TURN: ANGULO DE ALABEO (rads)
V_psi = xlsread(filename,sheet,'D40');%	TURN: VELOCIDAD DE GUIÑADA (rads/seg)
n_tr = xlsread(filename,sheet,'D41');%	TURN: FACTOR DE CARGA
R_tr = xlsread(filename,sheet,'D42');%	TURN: RADIO DE GIRO (m)
% Descent
h_inicial_d = xlsread(filename,sheet,'D43');%	DESCENT: ALTURA INICIAL
h_final_d = xlsread(filename,sheet,'D44');%	DESCENT:  ALTURA FINAL
gamma_d = xlsread(filename,sheet,'D45');%	DESCENT: GAMMA
V_d = xlsread(filename,sheet,'D46');%	DESCENT: VELOCIDAD DE DESCENSI
EAS_d = xlsread(filename,sheet,'D47');%	DESCENT : VELOCIDAD TAS  - [m/s]
TAS_d = xlsread(filename,sheet,'D48');%	DESCENT : VELOCIDAD EAS  - [m/s]
delta_T_d = xlsread(filename,sheet,'D49');%	DESCENT : PALANCA DE GASES  - [-]
V_ini_d = xlsread(filename,sheet,'D50');%	DESCENT : VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
V_fin_d = xlsread(filename,sheet,'D51');%	DESCENT : VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
% Turn Wait
h_inicial_tr_wt = xlsread(filename,sheet,'D52');%	TURN WAIT:  ALTURA INICIAL
t_final_tr_wt = xlsread(filename,sheet,'D53');%	TURN WAIT: TIEMPO FINAL (seg)
V_turn_wt = xlsread(filename,sheet,'D54');%	TURN WAIT: turn velocity
delta_T_tr_wt = xlsread(filename,sheet,'D55');%	TURN WAIT: PALANCA DE GASES
phi_tr_wt = xlsread(filename,sheet,'D56');%	TURN WAIT: ANGULO DE ALABEO (rads)
V_psi_wt = xlsread(filename,sheet,'D57');%	TURN WAIT: VELOCIDAD DE GUIÑADA (rads/seg)
n_tr_wt = xlsread(filename,sheet,'D58');%	TURN WAIT: FACTOR DE CARGA
R_tr_wt = xlsread(filename,sheet,'D59');%	TURN WAIT: RADIO DE GIRO (m)
% lANDING
temp_local_LND = xlsread(filename,sheet,'D60');%	LANDING:  TEMPERATURA LOCAL (Celsius)
h_inicial_LND = xlsread(filename,sheet,'D61');%	LANDING: ALTURA LOCAL (m)
P_local_LND = xlsread(filename,sheet,'D62');%	LANDING:  PRESION LOCAL (Pa)
mu_LND = xlsread(filename,sheet,'D63');%	LANDING: COEFICIENTE DE FRICCION CON LA PISTA (MU)
delta_T_LND = xlsread(filename,sheet,'D64');%	LANDING:  PALANCA DE REVERSA
t_brake = xlsread(filename,sheet,'D65');%	LANDING: 'Tiempo en activar frenos' - [s]
% Dummy Segment
dummy(1) = xlsread(filename,sheet,'D66');%	DUMMY SEGMENT: ALTURA INICIAL (CRUISE)
dummy(2) = xlsread(filename,sheet,'D67');%	DUMMY SEGMENT: DISTANCIA FINAL
dummy(3) = xlsread(filename,sheet,'D68');%	DUMMY SEGMENT: VELOCIDAD (CRUISE)

% Stores the flags
% Taxy
IPP_flags.temp_local_taxy = temp_local_taxy;
IPP_flags.h_inicial_taxy = h_inicial_taxy;
IPP_flags.P_local_taxy = P_local_taxy;
IPP_flags.delta_T_taxy = delta_T_taxy;
IPP_flags.V_taxy = V_taxy;
IPP_flags.t_taxy = t_taxy;
% Takeoff
IPP_flags.temp_local_TO = temp_local_TO;
IPP_flags.h_inicial_TO = h_inicial_TO;
IPP_flags.P_local_TO = P_local_TO;
IPP_flags.mu_TO = mu_TO;
IPP_flags.h_obstacle_TO = h_obstacle_TO;
IPP_flags.gamma_climb_TO = gamma_climb_TO;
IPP_flags.delta_T_TO = delta_T_TO;
% Climb
IPP_flags.h_inicial_cl = h_inicial_cl; 
IPP_flags.h_final_cl = h_final_cl; 
IPP_flags.gamma_cl = gamma_cl; 
IPP_flags.Mach_cl = Mach_cl; 
IPP_flags.TAS_cl = TAS_cl; 
IPP_flags.EAS_cl = EAS_cl; 
IPP_flags.delta_T_cl = delta_T_cl;
IPP_flags.V_ini_cl = V_ini_cl; 
IPP_flags.V_fin_cl = V_fin_cl; 
% Cruise
IPP_flags.h_inicial_cr = h_inicial_cr;
IPP_flags.dist_final_cr = dist_final_cr;
IPP_flags.V_cr = V_cr;
IPP_flags.delta_T_cr = delta_T_cr;
IPP_flags.V_ini_cr = V_ini_cr;
IPP_flags.V_fin_cr = V_fin_cr;
IPP_flags.fuel_cr = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
IPP_flags.Cd0_cr = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
IPP_flags.k1_cr = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
IPP_flags.k2_cr = k2_cr;% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
% Turn
IPP_flags.h_inicial_tr = h_inicial_tr;
IPP_flags.t_final_tr = t_final_tr;
IPP_flags.V_turn = V_turn; 
IPP_flags.delta_T_tr = delta_T_tr;
IPP_flags.phi_tr = phi_tr; 
IPP_flags.V_psi = V_psi;
IPP_flags.n_tr = n_tr;
IPP_flags.R_tr = R_tr;
% Descent
IPP_flags.h_inicial_d = h_inicial_d;
IPP_flags.h_final_d = h_final_d;
IPP_flags.gamma_d = gamma_d;
IPP_flags.V_d = V_d;
IPP_flags.EAS_d = EAS_d;
IPP_flags.TAS_d = TAS_d;
IPP_flags.delta_T_d = delta_T_d;
IPP_flags.V_ini_d = V_ini_d;
IPP_flags.V_fin_d = V_fin_d;
% Turn Wait
IPP_flags.h_inicial_tr_wt = h_inicial_tr_wt;
IPP_flags.t_final_tr_wt = t_final_tr_wt;
IPP_flags.V_turn_wt = V_turn_wt;
IPP_flags.delta_T_tr_wt = delta_T_tr_wt;
IPP_flags.phi_tr_wt = phi_tr_wt;
IPP_flags.V_psi_wt = V_psi_wt;
IPP_flags.n_tr_wt = n_tr_wt;
IPP_flags.R_tr_wt = R_tr_wt;
% Landing
IPP_flags.temp_local_LND = temp_local_LND;
IPP_flags.h_inicial_LND = h_inicial_LND;
IPP_flags.P_local_LND = P_local_LND;
IPP_flags.mu_LND = mu_LND;
IPP_flags.delta_T_LND = delta_T_LND;
IPP_flags.t_brake = t_brake;
% Dummy Segment
IPP_flags.dummy = dummy;

% Gathers all the flags
OUTPUT_read_XLSX.IPP_flags = IPP_flags;

%% Performance Study Variblae conditions
sheet = 14;

type_missions_WF = xlsread(filename,sheet,'D3');%	Selection Type of Mission
num_missions_WF = xlsread(filename,sheet,'D4');%	Number of Segments
climb_mode = xlsread(filename,sheet,'D5');%	Selects Climb Mode Option
cruise_mode = xlsread(filename,sheet,'D6');%	Selects Cruise Mode Option
turn_mode = xlsread(filename,sheet,'D7');%	Selects Turn Mode Option
descent_mode = xlsread(filename,sheet,'D8');%	Selects Descent Mode Option

% Stores the flags
PerforMisionSelection_flags.type_missions_WF = type_missions_WF; 
PerforMisionSelection_flags.num_missions_WF = num_missions_WF;
PerforMisionSelection_flags.climb_mode = climb_mode;
PerforMisionSelection_flags.cruise_mode = cruise_mode;
PerforMisionSelection_flags.turn_mode = turn_mode;
PerforMisionSelection_flags.descent_mode = descent_mode;
% Gathers all the flags
OUTPUT_read_XLSX.PerforMisionSelection_flags = PerforMisionSelection_flags;


%% Performance Study Variblae conditions
sheet = 15;
V_low = xlsread(filename,sheet,'D3');% 	Performance Variable study - Low velocity  in m/s
V_high = xlsread(filename,sheet,'D4');% 	Performance Variable study - High velocity  in m/s
N_V_VAR_perf = xlsread(filename,sheet,'D5');% Performance Variable study for Velocity - Number of variable points
V_single = xlsread(filename,sheet,'D6');%	Performance Study - Single Speed
Wp_low = xlsread(filename,sheet,'D7');%	Performance Variable study - Low weight  in m/s
Wp_high = xlsread(filename,sheet,'D8');%	Performance Variable study - High weight  in m/s
N_Wp_VAR_perf = xlsread(filename,sheet,'D9');%	Performance Variable study for Weight - Number of variable points
W_single = xlsread(filename,sheet,'D10');%	Performance Study - Single Weight
Post_processing_PERFORMANCE = xlsread(filename,sheet,'D11');%	Conducts the Post processing PERFORMANCE
climb_mode = xlsread(filename,sheet,'D12');%	Selects Climb Mode Option
cruise_mode = xlsread(filename,sheet,'D13');%	Selects Cruise Mode Option
turn_mode = xlsread(filename,sheet,'D14');%	Selects Turn Mode Option
descent_mode = xlsread(filename,sheet,'D15');%	Selects Descent Mode Option

% Stores the flags
PerformanceStudy_flags.V_low = V_low; 
PerformanceStudy_flags.V_high = V_high;
PerformanceStudy_flags.N_V_VAR_perf = N_V_VAR_perf;
PerformanceStudy_flags.V_single = V_single;
PerformanceStudy_flags.Wp_low = Wp_low;
PerformanceStudy_flags.Wp_high = Wp_high;
PerformanceStudy_flags.N_Wp_VAR_perf = N_Wp_VAR_perf;
PerformanceStudy_flags.W_single = W_single;
PerformanceStudy_flags.Post_processing_PERFORMANCE = Post_processing_PERFORMANCE;
PerformanceStudy_flags.climb_mode = climb_mode;
PerformanceStudy_flags.cruise_mode = cruise_mode;
PerformanceStudy_flags.turn_mode = turn_mode;
PerformanceStudy_flags.descent_mode = descent_mode;

% Gathers all the flags
OUTPUT_read_XLSX.PerformanceStudy_flags = PerformanceStudy_flags;



