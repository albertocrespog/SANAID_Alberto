function OUTPUT_read_XLSX =  read_DATA_XLSX_advanced(case_AC)

%% Pre-load variables stored in Excel
% Function that reads all the data for the different aircraft
get_add_path
filename = '../AIRCRAFT/SANAID_AIRCRAFT_DATA_v6.xlsx';

excel_file = 'SANAID_AIRCRAFT_DATA_v6.xlsx';
    
% Writes in the excel to assign the propper aircraft type
sheet = 3;
status = xlswrite(filename,case_AC,sheet,'D8');

%% MATLAB Flags
% Sheet1 = readcell('SANAID_AIRCRAFT_DATA_v3.xlsx','Sheet','1MATLAB','Range','D3:D6');
Sheet1 = readcell(excel_file,'Sheet','1MATLAB','Range','D3:D8');

MATLAB_in = Sheet1{1};
CHECK_Efficiency = Sheet1{2};
Detailed_Profile = Sheet1{3};
save_MAT = Sheet1{4};
show_messages_screen = Sheet1{5};
write_data = Sheet1{6};

% Stores the flags
MATLAB_flags.MATLAB_in = MATLAB_in;
MATLAB_flags.CHECK_Efficiency = CHECK_Efficiency;
MATLAB_flags.Detailed_Profile = Detailed_Profile;
MATLAB_flags.save_MAT = save_MAT;
MATLAB_flags.show_messages_screen = show_messages_screen;
MATLAB_flags.write_data = write_data;

% Gathers all the flags
OUTPUT_read_XLSX.MATLAB_flags = MATLAB_flags;

%% Studies Flags
Sheet2 = readcell(excel_file,'Sheet','2Studies','Range','D3:D25');

AERODYNAMIC_STUDY = Sheet2{1}; % Conducts Aerodynamic Studies
PROP_STUDY = Sheet2{2}; % Conducts Proppeller Studies
PERFORMANCE_STUDY = Sheet2{3}; % Conducts Performance Studies
STABILITY_STUDY = Sheet2{4}; % Conducts Stability Studies
MISSIONS_STUDY = Sheet2{5}; % Conducts Mission Studies
STUDY_Weight_Fraction = Sheet2{6}; % Weight Fraction Study
ANALYSIS_PERFORMANCE_AP = Sheet2{7}; % Performance Analysis integrated with AP Codes
ANALYSIS_PERFORMANCE_AP_var = Sheet2{8}; % Performance Analysis integrated with AP Codes varying Conditions
variable_speed_AP = Sheet2{9}; % Variable Speed Performance Studies
variable_weight_AP = Sheet2{10}; % Variable Mass Performance Studies
STABILITY_STUDY_Trim = Sheet2{11}; % Trim conditions - one case
STABILITY_STUDY_Regular = Sheet2{12}; % Stability Analysis Regular Study
STABILITY_STUDY_Trim_varV_m0 = Sheet2{13}; % Trim conditions - variation V and mass
STABILITY_STUDY_Trim_var_XCG = Sheet2{14}; % Trim conditions - variation XCG
STABILITY_STUDY_Long_dyn = Sheet2{15}; % Longitudinal Stability Analysis
STABILITY_STUDY_LatDir_dyn = Sheet2{16}; % Lateral-Directional Stability Analysis
STABILITY_STUDY_Trim_lat = Sheet2{17}; % Lateral Directional Trim
% STABILITY_STUDY_Trim_lat_assymetric_drag = Sheet2{18}; % Turning Stability
% STABILITY_STUDY_Trim_lat_accelerations = Sheet2{18}; % Turning Stability
STABILITY_STUDY_Turning = Sheet2{18}; % Turning Stability
STABILITY_STUDY_Stability_Derivatives_varV_m0 = Sheet2{19}; % Plots Stability Derivatives
STABILITY_STUDY_Stability_Analysis_varV_m0_long = Sheet2{20}; % Longitudinal Dynamic Response Plots
STABILITY_STUDY_Stability_Analysis_varV_m0_lat = Sheet2{21}; % Lateral-Directional Dynamic Response Plots
STABILITY_STUDY_Stability_Analysis_dyna_impulse_long = Sheet2{22}; % Impulse Longitudinal Dynamic Response Plots
STABILITY_STUDY_Stability_Analysis_dyna_impulse_lat = Sheet2{23}; % Impulse Lateral-Directional Dynamic Response Plots

% Stores the flags
STUDY_flags.AERODYNAMIC_STUDY = AERODYNAMIC_STUDY; % Conducts Aerodynamic Studies
STUDY_flags.PROP_STUDY = PROP_STUDY; % Conducts Proppeller Studies
STUDY_flags.PERFORMANCE_STUDY = PERFORMANCE_STUDY; % Conducts Performance Studies
STUDY_flags.STABILITY_STUDY = STABILITY_STUDY; % Conducts Stability Studies
STUDY_flags.MISSIONS_STUDY = MISSIONS_STUDY; % Conducts Stability Studies
STUDY_flags.STUDY_Weight_Fraction = STUDY_Weight_Fraction;% Weight Fraction Study
STUDY_flags.ANALYSIS_PERFORMANCE_AP = ANALYSIS_PERFORMANCE_AP;% Performance Analysis integrated with AP Codes
STUDY_flags.ANALYSIS_PERFORMANCE_AP_var = ANALYSIS_PERFORMANCE_AP_var;% Performance Analysis integrated with AP Codes varying Conditions
STUDY_flags.variable_speed_AP = variable_speed_AP;% Variable Speed Performance Studies
STUDY_flags.variable_weight_AP = variable_weight_AP;% Variable Mass Performance Studies
STUDY_flags.STABILITY_STUDY_Trim = STABILITY_STUDY_Trim;% Trim conditions - one case
STUDY_flags.STABILITY_STUDY_Regular = STABILITY_STUDY_Regular;% Regular Stability Analysis
STUDY_flags.STABILITY_STUDY_Trim_varV_m0 = STABILITY_STUDY_Trim_varV_m0;% Trim conditions - variation V and mass
STUDY_flags.STABILITY_STUDY_Trim_var_XCG = STABILITY_STUDY_Trim_var_XCG;% Trim conditions - variation XCG
STUDY_flags.STABILITY_STUDY_Long_dyn = STABILITY_STUDY_Long_dyn;% Longitudinal Stability Analysis
STUDY_flags.STABILITY_STUDY_LatDir_dyn = STABILITY_STUDY_LatDir_dyn;% Lateral-Directional Stability Analysis
STUDY_flags.STABILITY_STUDY_Trim_lat = STABILITY_STUDY_Trim_lat;% Lateral Directional Trim
STUDY_flags.STABILITY_STUDY_Turning = STABILITY_STUDY_Turning;% Turning Stability
STUDY_flags.STABILITY_STUDY_Stability_Derivatives_varV_m0 = STABILITY_STUDY_Stability_Derivatives_varV_m0;% Plots Stability Derivatives
STUDY_flags.STABILITY_STUDY_Stability_Analysis_varV_m0_long = STABILITY_STUDY_Stability_Analysis_varV_m0_long;% Longitudinal Dynamic Response Plots
STUDY_flags.STABILITY_STUDY_Stability_Analysis_varV_m0_lat = STABILITY_STUDY_Stability_Analysis_varV_m0_lat;% Lateral-Directional Dynamic Response Plots
STUDY_flags.STABILITY_STUDY_Stability_Analysis_dyna_impulse_long = STABILITY_STUDY_Stability_Analysis_dyna_impulse_long;% Impulse Longitudinal Dynamic Response Plots
STUDY_flags.STABILITY_STUDY_Stability_Analysis_dyna_impulse_lat = STABILITY_STUDY_Stability_Analysis_dyna_impulse_lat;% Impulse Lateral-Directional Dynamic Response Plots

% Gathers all the flags
OUTPUT_read_XLSX.STUDY_flags = STUDY_flags;

%% Aircraft Data
Sheet3 = readcell(excel_file,'Sheet','3AC_Data','Range','D3:D9');

SF = Sheet3{1}; % Scaling Factor
Weight_Estimation = Sheet3{2}; % Weight Estimation
AC_type = Sheet3{3}; % Aircaft Type
CASE_fuse = Sheet3{4};% Determines the fuselage than will be shown
ESCALADO = Sheet3{5}; % Determine the type of scaling
case_AC = Sheet3{6}; % Aircraft Model Analized
twin_VTP = Sheet3{7}; % Defines twin VTP

% Stores the flags
AC_Data_flags.SF = SF; %
AC_Data_flags.Weight_Estimation = Weight_Estimation; %
AC_Data_flags.AC_type = AC_type; % 
AC_Data_flags.CASE_fuse = CASE_fuse; % 
AC_Data_flags.ESCALADO = ESCALADO;% 
AC_Data_flags.case_AC = case_AC;% 
AC_Data_flags.twin_VTP = twin_VTP;% 
% Gathers all the flags
OUTPUT_read_XLSX.AC_Data_flags = AC_Data_flags;

%% Propulsive Data
Sheet4 = readcell(excel_file,'Sheet','4Propulsion','Range','D3:D41');

propul(1) = Sheet4{1}; % Type of engine
propul(2) = Sheet4{2}; % Number of engines
propul(3) = Sheet4{3}; % EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine
propul(4) = Sheet4{4}; % CONSUMO ESPECIFICO: % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
propul(5) = Sheet4{5}; % DERIVACION(TURBOFANES) By-pass 
propul(6) = Sheet4{6}; % EFICIENCIA DE LA HELICE (ETA_P)
propul(7) = Sheet4{7}; % Normativa Propulsion
propul(8) = Sheet4{8}; % CAPACIDAD CALORIFICA DEL COMBUSTIBLE
propul(9) = Sheet4{9}; % DIAMETRO DEL MOTOR(TURBOHELICES)
type_battery = Sheet4{10}; % Type Battery Used for Propulsion
Engine_loc = Sheet4{11}; % Engine location
Engine_conf = Sheet4{12}; % Engine configuration
D_prop = Sheet4{13}; % Preliminary Prop Diameter (m)
prop_known = Sheet4{14}; % Flag that determines if the prop is known
prop_properties = Sheet4{15}; % model of prop:first 2 numbers diameter in inches last 2 numbers pitch
type_prop = Sheet4{16}; % Type of propller
Prop_type = Sheet4{17}; % Propulsive Model- Number of Prop used of propller
l_eng = Sheet4{18}; % Engine geometry - length
d_eng = Sheet4{19}; % Engine geometry - diameter
l_nc = Sheet4{20}; % Nacelle geometry - length
d_nc = Sheet4{21}; % Nacelle geometry - diameter
model_prop = Sheet4{22}; %	Compares 3 prop models
prop_selec_APC = Sheet4{23}; %	Model 1 - APC data
prop_selec_WT1 = Sheet4{24}; %	Model 2 - Wind tunnel data for different props
prop_selec_WT2 = Sheet4{25}; %	Model 3 - Wind tunnel data for different angle of attack
eta_gear = Sheet4{26}; % Gearbox efficiency
eta_m = Sheet4{27}; % Motor efficiency
eta_esc = Sheet4{28}; % ESC efficiency
eta_dist = Sheet4{29}; % Electrical distribution efficiency
SE_LiFePO4 = Sheet4{30}; % Specific Energy of LiFePO4 batteries 90–160 Wh/kg (320–580 J/g or kJ/kg) 
SE_LiPo = Sheet4{31}; % Specific Energy of LiPo batteries (Wh/kg)
SE_FuelCells = Sheet4{32}; % Wh/kg Specific Energy for fuel cells
SP_motor = Sheet4{33}; % Motor Specific Power kW/kg - TIER 1
SP_ESC = Sheet4{34}; % ESC Specific Power kW/kg - TIER 1
m_prop_din = Sheet4{35}; % Mass of propeller as a function of diameter in inches
compare_prop = Sheet4{36}; % Determinee Plots to compare different propulsion models
density_fuel = Sheet4{37}; % Fuel density
cost_fuel = Sheet4{38}; % Fuel Cost  cts/gal
CI = Sheet4{39}; % Cost Index Kg/s

% Stores the flags
Propulsive_flags.propul = propul; %
Propulsive_flags.type_battery = type_battery; %
Propulsive_flags.Engine_loc = Engine_loc; %
Propulsive_flags.Engine_conf = Engine_conf; %
Propulsive_flags.D_prop = D_prop; %
Propulsive_flags.prop_known = prop_known; %
Propulsive_flags.prop_properties = prop_properties; %
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
Propulsive_flags.compare_prop = compare_prop; %
Propulsive_flags.density_fuel = density_fuel; %
Propulsive_flags.cost_fuel = cost_fuel; %
Propulsive_flags.CI = CI; %

% Gathers all the flags
OUTPUT_read_XLSX.Propulsive_flags = Propulsive_flags;

%% Weights Data
Sheet5 = readcell(excel_file,'Sheet','5Weights','Range','D3:D46');

W_eng = Sheet5{1}; % Mass of engine (kg)
m_prop = Sheet5{2}; % Mass of proppeller (kg)
m_ESC = Sheet5{3}; % Mass of ESC (kg)
n_servos = Sheet5{4}; % Number of servos 
m_servo = Sheet5{5}; % Mass per servo (kg)
m_wiring = Sheet5{6}; % Mass of electronic wiring (kg)
m_metal = Sheet5{7}; % Mass of misc (kg)
m_systems = Sheet5{8}; % Mass of Systems (kg)
m_landing_gear = Sheet5{9}; % Mass of landing gear (kg)
n_pax = Sheet5{10}; % Number of Passangers
n_crew = Sheet5{11}; % Number of Crew
m_pax = Sheet5{12}; % Mass per Passenger (kg)
m_crew = Sheet5{13}; % Mass per Crew (kg)
m_luggage = Sheet5{14}; % Mass of luggage (kg)
m_cargo = Sheet5{15}; % Mass of payload cargo (kg)
m_batteries = Sheet5{16}; % Mass of batteries (kg)
m_fuel = Sheet5{17}; % Mass of fuel (kg) 
flag_landing_gear = Sheet5{18}; % Determination of Landing Gear Estimation if available
f_f = Sheet5{19}; % Fudge Factor that increments the weight estimation
m_cell_A123 = Sheet5{20}; % Mass of one A123 cell
m_w1 = Sheet5{21}; % Mass of one wing
m_HTP = Sheet5{22}; % Mass of HTP
m_VTP = Sheet5{23}; % Mass of VTP
m_Can = Sheet5{24}; % Mass of Canard
m_Vee = Sheet5{25}; % Mass of Vee Tail
m_fus_fairing = Sheet5{26}; % Mass of fuselage and fairing
m_prop_fairing = Sheet5{27}; % Mass of propeller nacelle
I_xx = Sheet5{28}; % Moments of Inertia Ixx
I_yy = Sheet5{29}; % Moments of Inertia Ixx
I_zz = Sheet5{30}; % Moments of Inertia Ixx
I_xy = Sheet5{31}; % Moments of Inertia Ixx
I_xz = Sheet5{32}; % Moments of Inertia Ixx
I_yz = Sheet5{33}; % Moments of Inertia Ixx
flag_moments_inertia = Sheet5{34}; % Determination of Moments of innertia from literature	
flag_total_weights = Sheet5{35}; % Determination if using the true weights
MTOW_true = Sheet5{36}; % True MTOW - Maximum Takeke Off Weight
ME_true = Sheet5{37}; % True ME - Empty Weight
MCREW_true = Sheet5{38}; % True MCREW - Crew Weight
MLW_true = Sheet5{39}; % True MLW - Máxcimum landing Weight (85% MTOW)
MF_true = Sheet5{40}; % True MTOW - Maximum Fuel Weight
MPL_true = Sheet5{41}; % True MTOW - Maximum Paylod Weigh
fraction_MF = Sheet5{42}; % Fraction of fuel
fraction_PL = Sheet5{43}; % Fraction of payload
get_shift_XCG_variation = Sheet5{44}; % Variation of XCG study	

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
Weights_flags.flag_total_weights = flag_total_weights; % Determination if using the true weights
Weights_flags.MTOW_true = MTOW_true; % Truen MTOW - Maximum Takeke Off Weight
Weights_flags.ME_true = ME_true; % Truen ME - Empty Weight
Weights_flags.MCREW_true = MCREW_true; % Truen MCREW - Crew Weight
Weights_flags.MLW_true = MLW_true; % Truen MLW - Máxcimum landing Weight (85% MTOW)
Weights_flags.MF_true = MF_true; % Truen MTOW - Maximum Fuel Weight
Weights_flags.MPL_true = MPL_true; % Truen MTOW - Maximum Paylod Weight
Weights_flags.fraction_MF = fraction_MF; % Fraction of fuel
Weights_flags.fraction_PL = fraction_PL; % Fraction of payload
Weights_flags.get_shift_XCG_variation = get_shift_XCG_variation; % 

% Gathers all the flags
OUTPUT_read_XLSX.Weights_flags = Weights_flags;


%% Fuselage Data
Sheet6 = readcell(excel_file,'Sheet','6CAD_data','Range','D3:D19');

nSections = Sheet6{1}; % Number of sectins (x coordinate) (eliminated 1 for convergence)
nPoints = Sheet6{2}; % number of elements per section (y and z coordinates)
lecture_Geo = Sheet6{3}; % lineas a partir desde donde empieza a leer
STL_PLOT = Sheet6{4}; % Generates Fuselage mesh from CAD STL
XFLR5_file = Sheet6{5}; % Defines flies to be used for the Fuselage geometry - XFLR5
STL_ac = Sheet6{6}; % Defines flies to be used for the aircraft geometry - STL
STL_fus = Sheet6{7}; % Defines flies to be used for the Fuselage geometry - STL
STL_wing = Sheet6{8}; % Defines flies to be used for the wing geometry - STL
STL_canard = Sheet6{9}; % Defines flies to be used for the canard geometry - STL
STL_HTP = Sheet6{10}; % Defines flies to be used for the Fuselage geometry - STL
STL_VTP = Sheet6{11}; % Defines flies to be used for the HTP geometry - STL
STL_Vee = Sheet6{12}; % Defines flies to be used for the VTP geometry - STL
STL_engine = Sheet6{13}; % Defines flies to be used for the engine geometry - STL
STL_nacelle = Sheet6{14}; % Defines flies to be used for the nacelle geometry - STL
SF_CAD = Sheet6{15}; % Defines Scaling Factor for STL files for units conversion (example mm 2 m)	
AC_STL_Compare = Sheet6{16}; % % Plots both Generated Geometry and Aircraft STL comparisso	AC_STL_Compare
PLOT_Slices_fuselage_STL = Sheet6{17}; % Plotsfuselage slices to check if they are rotated accordingly	

% Stores the flags
Fuselage_flags.nSections = nSections; %
Fuselage_flags.nPoints = nPoints; %
Fuselage_flags.lecture_Geo = lecture_Geo; %
Fuselage_flags.STL_PLOT = STL_PLOT; %
Fuselage_flags.XFLR5_file = XFLR5_file; %
Fuselage_flags.STL_ac = STL_ac; %
Fuselage_flags.STL_fus = STL_fus; %
Fuselage_flags.STL_wing = STL_wing; %
Fuselage_flags.STL_canard = STL_canard; %
Fuselage_flags.STL_HTP = STL_HTP; %
Fuselage_flags.STL_VTP = STL_VTP; %
Fuselage_flags.STL_Vee = STL_Vee; %
Fuselage_flags.STL_engine = STL_engine; %
Fuselage_flags.STL_nacelle = STL_nacelle; %
Fuselage_flags.SF_CAD = SF_CAD; %
Fuselage_flags.AC_STL_Compare = AC_STL_Compare; %
Fuselage_flags.PLOT_Slices_fuselage_STL = PLOT_Slices_fuselage_STL; %

% Gathers all the flags
OUTPUT_read_XLSX.Fuselage_flags = Fuselage_flags;

%% Performance Data Prelimina
Sheet7 = readcell(excel_file,'Sheet','7Performance','Range','D3:D13');

h_climb = Sheet7{1}; % Preliminar Performance conditions: Climb altitude (m)
Endurance_v = Sheet7{2}; % Preliminar Performance conditions: Endurance in Hoovering (min)
h = Sheet7{3}; % Preliminar Performance conditions: Cruise initial altitude (m)
V = Sheet7{4}; % Preliminar Performance conditions: Cruise initial airspeed (m/s)
V_max = Sheet7{5}; % Preliminar Performance conditions: Cruise max speed
Range = Sheet7{6}; % Preliminar Performance conditions: Cruise Range (m)
Endurance = Sheet7{7}; % Preliminar Performance conditions: Cruise Endurance (min)
Flight_SF = Sheet7{8}; % Flight Safety Margin - normal 
Flight_cruise = Sheet7{9}; % Flight Condition - Cruise 	 
Flight_takeoff = Sheet7{10}; % Flight Condition - TakeOff	 
Flight_climb = Sheet7{11}; % Flight Condition - Climb

% Stores the flags
Performance_pre_flags.h_climb = h_climb; %
Performance_pre_flags.Endurance_v = Endurance_v; %
Performance_pre_flags.h = h; %
Performance_pre_flags.V = V; %
Performance_pre_flags.V_max = V_max; %
Performance_pre_flags.Range = Range; %
Performance_pre_flags.Endurance = Endurance; %
Performance_pre_flags.Flight_SF = Flight_SF; %
Performance_pre_flags.Flight_cruise = Flight_cruise; %
Performance_pre_flags.Flight_takeoff = Flight_takeoff; %
Performance_pre_flags.Flight_climb = Flight_climb; %

% Gathers all the flags
OUTPUT_read_XLSX.Performance_pre_flags = Performance_pre_flags;

%% Aerodynamic Data
Sheet8 = readcell(excel_file,'Sheet','8Aero','Range','D3:D49');

index_w1 = Sheet8{1}; % Selection of TXT that are used for the aerodynamic analysis (w1) - LLT
index_w2 = Sheet8{2}; % Selection of TXT that are used for the aerodynamic analysis (w2) - LLT
index_w3 = Sheet8{3}; % Selection of TXT that are used for the aerodynamic analysis (w3) - LLT
index_VTP = Sheet8{4}; % Selection of TXT that are used for the aerodynamic analysis (vtp) - LLT
index_w1_VLM = Sheet8{5}; % Selection of TXT that are used for the aerodynamic analysis (w1) - VLM
index_w2_VLM = Sheet8{6}; % Selection of TXT that are used for the aerodynamic analysis (w2) - VLM
index_w3_VLM = Sheet8{7}; % Selection of TXT that are used for the aerodynamic analysis (w3) - VLM
index_VTP_VLM = Sheet8{8}; % Selection of TXT that are used for the aerodynamic analysis (VTP) - VLM
i_w1 = Sheet8{9}; % Selects the AoA of each surface relative to fuselage line (w1)
i_w2 = Sheet8{10}; % Selects the AoA of each surface relative to fuselage line (w2)
i_w3 = Sheet8{11}; % Selects the AoA of each surface relative to fuselage line (w3)
i_VTP = Sheet8{12}; % Selects the AoA of each surface relative to fuselage line (w3)
w1 = Sheet8{13};% Selection of elements for theoretical polar estimation - wing	w1
h = Sheet8{14};% Selection of elements for theoretical polar estimation - HTP	h
v = Sheet8{15};% Selection of elements for theoretical polar estimation - VTP	v
v2 = Sheet8{16};% Selection of elements for theoretical polar estimation - twin VTP	v2
vtail = Sheet8{17};% Selection of elements for theoretical polar estimation - Vtail	vtail
can = Sheet8{18};% Selection of elements for theoretical polar estimation - canard	can
fus = Sheet8{19};% Selection of elements for theoretical polar estimation - fuselage	fus
m_fus = Sheet8{20};% Selection of elements for theoretical polar estimation - multiple fuselage	m_fus
n_m_fus = Sheet8{21};% Selection of elements for theoretical polar estimation - number of multiple fuselage	n_m_fus
nac = Sheet8{22};% Selection of elements for theoretical polar estimation - nacelle	nac
landgear = Sheet8{23};% Selection of elements for theoretical polar estimation - landing gear	landgear
tailboom = Sheet8{24};% Selection of elements for theoretical polar estimation - tailboom
m_tailboom = Sheet8{25};% Selection of elements for theoretical polar estimation - multiple tailboom
n_m_tailboom = Sheet8{26};% Selection of elements for theoretical polar estimation - number of multiple tailboom
missile = Sheet8{27};% Selection of elements for theoretical polar estimation - misile
n_missile = Sheet8{28};% Selection of elements for theoretical polar estimation - misile
pod = Sheet8{29};% Selection of elements for theoretical polar estimation - pod
flap = Sheet8{30};% Selection of elements for theoretical polar estimation - flap
compare_plot_aero = Sheet8{31};% Determines the aero plots to compare
read_XFLR5 = Sheet8{32};% Reads Aero data from XFLR5	
read_FLOW5 = Sheet8{33}; % Reads Aero data from FLOW5	
airfoil_w1 = Sheet8{34}; % Airfoil for wing (primary)	
airfoil_w2 = Sheet8{35}; % Airfoil for wing (secondary)	
airfoil_c1 = Sheet8{36}; % Airfoil for canard (primary)	
airfoil_c2 = Sheet8{37}; % Airfoil for canard (secondary)	
airfoil_HTP1 = Sheet8{38}; % Airfoil for HTP (primary)	
airfoil_HTP2 = Sheet8{39}; % Airfoil for HTP (secondary)
airfoil_VTP1 = Sheet8{40}; % Airfoil for VTP (primary)	
airfoil_VTP2 = Sheet8{41}; % Airfoil for VTP (secondary)	
airfoil_Vee1 = Sheet8{42}; % Airfoil for VeeTail (primary)	
airfoil_Vee2 = Sheet8{43}; % Airfoil for VeeTail (secondary)	
polar_model = Sheet8{44}; % Use Global Weights	flag_total_weight	
fuse_aero_FLOW_and_CBM = Sheet8{45}; % Use Global Weights	flag_total_weight	
Flap_type = Sheet8{46}; % High Lift Devices  - Leading Edge
LED_type = Sheet8{47}; %High Lift Devices  - Trailing  Edge (Flaps)	

% Stores the flags
Aerodynamic_Data_flags.index_w1 = index_w1; %
Aerodynamic_Data_flags.index_w2 = index_w2; %
Aerodynamic_Data_flags.index_w3 = index_w3; %
Aerodynamic_Data_flags.index_VTP = index_VTP; %
Aerodynamic_Data_flags.index_w1_VLM = index_w1_VLM; %
Aerodynamic_Data_flags.index_w2_VLM = index_w2_VLM; %
Aerodynamic_Data_flags.index_w3_VLM = index_w3_VLM; %
Aerodynamic_Data_flags.index_VTP_VLM = index_VTP_VLM; %
Aerodynamic_Data_flags.i_w1 = i_w1; %
Aerodynamic_Data_flags.i_w2 = i_w2; %
Aerodynamic_Data_flags.i_w3 = i_w3; %
Aerodynamic_Data_flags.i_VTP = i_VTP; %
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
Aerodynamic_Data_flags.Conf.tailboom = tailboom; %
Aerodynamic_Data_flags.Conf.m_tailboom = m_tailboom; %
Aerodynamic_Data_flags.Conf.n_m_tailboom = n_m_tailboom; %
Aerodynamic_Data_flags.Conf.missile = missile; %
Aerodynamic_Data_flags.Conf.n_missile = n_missile; %
Aerodynamic_Data_flags.Conf.pod = pod; %
Aerodynamic_Data_flags.Conf.flap = flap; %
Aerodynamic_Data_flags.compare_plot_aero = compare_plot_aero; %
Aerodynamic_Data_flags.read_XFLR5 = read_XFLR5; %
Aerodynamic_Data_flags.read_FLOW5 = read_FLOW5; %
Aerodynamic_Data_flags.airfoil_w1 = airfoil_w1; % Airfoil for wing (primary)	
Aerodynamic_Data_flags.airfoil_w2 = airfoil_w2; % Airfoil for wing (secondary)	
Aerodynamic_Data_flags.airfoil_c1 = airfoil_c1; % Airfoil for canard (primary)	
Aerodynamic_Data_flags.airfoil_c2 = airfoil_c2; % Airfoil for canard (secondary)	
Aerodynamic_Data_flags.airfoil_HTP1 = airfoil_HTP1; % Airfoil for HTP (primary)	
Aerodynamic_Data_flags.airfoil_HTP2 = airfoil_HTP2; % Airfoil for HTP (secondary)
Aerodynamic_Data_flags.airfoil_VTP1 = airfoil_VTP1; % Airfoil for VTP (primary)	
Aerodynamic_Data_flags.airfoil_VTP2 = airfoil_VTP2; % Airfoil for VTP (secondary)	
Aerodynamic_Data_flags.airfoil_Vee1 = airfoil_Vee1; % Airfoil for VeeTail (primary)	
Aerodynamic_Data_flags.airfoil_Vee2 = airfoil_Vee2; % 
% Airfoil for VeeTail (secondary)	
Aerodynamic_Data_flags.polar_model = polar_model; %
Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %% High Lift Devices  - Leading Edge
Aerodynamic_Data_flags.Flap_type = Flap_type; %% High Lift Devices  - Leading Edge
Aerodynamic_Data_flags.LED_type = LED_type; %%High Lift Devices  - Trailing  Edge (Flaps)	

% Gathers all the flags
OUTPUT_read_XLSX.Aerodynamic_Data_flags = Aerodynamic_Data_flags;

%% Stability Data
Sheet9 = readcell(excel_file,'Sheet','9Stability','Range','D3:D35');

SM = Sheet9{1}; % Static Margin
prop_wash_effect = Sheet9{2}; % Prop Wash Effect Considered on the dynamic pressure calculation
XCG_FF = Sheet9{3}; % elects XCG depending if it's from desired stability conditions in Forward Flight
StabilityModel = Sheet9{4}; % elects XCG depending if it's from desired stability conditions in Forward Flight
only_trim = Sheet9{5}; % elects XCG depending if it's from desired stability conditions in Forward Flight
V_low = Sheet9{6}; %	Trim Conditions Velocity Range - V low
V_high = Sheet9{7}; %	Trim Conditions Velocity Range - V high
N_V_VAR = Sheet9{8}; %	Number of iteration variation Speed
N_m_VAR = Sheet9{9}; %	Number of iteration variation Weight
beta = Sheet9{10}; %	Trim Lateral Stability Conditions -side slip angle (beta)
beta_i = Sheet9{11}; %	Trim Lateral Stability Conditions -Variable side slip angle (initial beta)
beta_f = Sheet9{12}; %	Trim Lateral Stability Conditions -Variable side slip angle (final beta)
N_Delta_beta = Sheet9{13}; %	Trim Lateral Stability Conditions -Variable side slip angle (Number of beta's)
phi = Sheet9{14}; %	Turning Condition Stability Conditions -bank angle (beta)
phi_i = Sheet9{15}; %	Turning Condition Stability Conditions -Variable bank angle (initial phi)
phi_f = Sheet9{16}; %	Turning Condition Stability Conditions -Variable bank angle (final phi)
N_Delta_phi = Sheet9{17}; %	Turning Condition Stability Conditions -Variable bank angle (Number of beta's)
n_viraje = Sheet9{18}; %	Loading factur during turning flight 
Munk_fuselage_constribution = Sheet9{19}; %	Include Munk's fuselage contribution todetermine Xac, & desired Xcg location for an SM  
tf1_long = Sheet9{20}; %	Time period ploting Phugoid Mode
tf2_long = Sheet9{21}; %	Time period ploting Short Period Mode
tf1_lat = Sheet9{22}; %	Time period ploting Dutch Roll
Du = Sheet9{23}; %	Perturbation in forward speed velocity (percentage of trim velocity)
Dalpha = Sheet9{24}; %	Perturbation in angle of attack
Dq = Sheet9{25}; %	Perturbation in pitch rate
Dtheta = Sheet9{26}; %	Perturbation in pitch angle
Dbeta = Sheet9{27}; %	Perturbation in side slip angle
Dp = Sheet9{28}; %	Perturbation in roll rate
Dr = Sheet9{29}; %	Perturbation in yaw rate
Dphi = Sheet9{30}; % Perturbation in bank angle
n_min = Sheet9{31}; % Minimum Load Factor
n_max = Sheet9{32}; % Maximum Load Factor
N_n_VAR = Sheet9{33}; %	Number of iterations for variation of load factor

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
Stability_flags.tf1_long = tf1_long;
Stability_flags.tf2_long = tf2_long;
Stability_flags.tf1_lat = tf1_lat;
Stability_flags.Du = Du;
Stability_flags.Dalpha = Dalpha;
Stability_flags.Dq = Dq;
Stability_flags.Dtheta = Dtheta;
Stability_flags.Dbeta = Dbeta;
Stability_flags.Dp = Dp;
Stability_flags.Dr = Dr;
Stability_flags.Dphi = Dphi;
Stability_flags.n_min = n_min;
Stability_flags.n_max = n_max;
Stability_flags.N_n_VAR = N_n_VAR;

% Gathers all the flags
OUTPUT_read_XLSX.Stability_flags = Stability_flags;

%% Geometry Data
Sheet10 = readcell(excel_file,'Sheet','10Configuration','Range','D3:D5');

CASE = Sheet10{1}; % Model of Aircraft
Control_surface = Sheet10{2}; % Determines the ammount of control surface in the w2 aileron & elevator
AC_type = Sheet10{3}; % Aircraft type

% Stores the flags
Geometry_Data_flags.CASE = CASE; %
Geometry_Data_flags.Control_surface = Control_surface; %
Geometry_Data_flags.AC_type = AC_type; %

% Gathers all the flags
OUTPUT_read_XLSX.Geometry_Data_flags = Geometry_Data_flags;

%% Input Aircraft Geometry Data
Sheet11 = readcell(excel_file,'Sheet','11InputGeometry','Range','D3:D137');

x_offset_CAD = Sheet11{1}; % Distance from origin in CATIa to Fuselage (offset value) - x
z_offset_CAD = Sheet11{2}; % Distance from origin in CATIa to Fuselage (offset value) - z
w_fus = Sheet11{3}; % Fuselage geomtry (Units in m) - approx - width
h_fus = Sheet11{4}; % Fuselage geomtry (Units in m) - approx - heigth
l_fus = Sheet11{5}; % Fuselage geomtry (Units in m) - approx - length
d_fus = Sheet11{6}; % Fuselage geomtry (Units in m) - approx - diameter
%w1
y_loc_1R_y1_w1_CAD = Sheet11{7}; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_w1_CAD = Sheet11{8}; %	y loc of wing (w1) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_w1_CAD = Sheet11{9}; %	z loc of wing (w1) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_w1_CAD = Sheet11{10}; %	x loc of wing (w1) root chord LE position (distance from CAD refference point)
Lambda_LE_w1_e = Sheet11{11}; %	Sweep of w1 (deg)
dihedral_w1_e = Sheet11{12}; %	Dihedral of w1 (deg)
%w2
y_loc_1R_y1_w2_CAD = Sheet11{13}; %	y loc of wing (w2) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_w2_CAD = Sheet11{14}; %	y loc of wing (w2) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_w2_CAD = Sheet11{15}; %	z loc of wing (w2) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_w2_CAD = Sheet11{16}; %	x loc of wing (w2) root chord LE position (distance from CAD refference point)
Lambda_LE_w2_e = Sheet11{17}; %	Sweep of w2 (deg)
dihedral_w2_e = Sheet11{18}; %	Dihedral of w2 (deg)
%can
y_loc_1R_y1_can_CAD = Sheet11{19}; %	y loc of wing (can) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_can_CAD = Sheet11{20}; %	y loc of wing (can) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_can_CAD = Sheet11{21}; %	z loc of wing (can) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_can_CAD = Sheet11{22}; %	x loc of wing (can) root chord LE position (distance from CAD refference point)
Lambda_LE_can_e = Sheet11{23}; %	Sweep of can (deg)
dihedral_can_e = Sheet11{24}; %	Dihedral of can (deg)
%HTP
y_loc_1R_y1_HTP_CAD = Sheet11{25}; %	y loc of wing (HTP) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_HTP_CAD = Sheet11{26}; %	y loc of wing (HTP) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_HTP_CAD = Sheet11{27}; %	z loc of wing (HTP) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_HTP_CAD = Sheet11{28}; %	x loc of wing (HTP) root chord LE position (distance from CAD refference point)
Lambda_LE_HTP_e = Sheet11{29}; % Sweep of HTP (deg)
dihedral_HTP_e = Sheet11{30}; %	Dihedral of HTP (deg)
%VTP
y_loc_1R_y1_VTP_CAD = Sheet11{31}; %	y loc of wing (VTP) root chord LE position (distance from CAD refference point)
z_loc_1R_y1_VTP_CAD = Sheet11{32}; %	y loc of wing (VTP) tip chord LE position (distance from CAD refference point)
z_loc_1R_y2_VTP_CAD = Sheet11{33}; %	z loc of wing (VTP) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_VTP_CAD = Sheet11{34}; %	x loc of wing (VTP) root chord LE position (distance from CAD refference point)
Lambda_LE_VTP_e = Sheet11{35}; %	Sweep of VTP (deg)
dihedral_VTP_e = Sheet11{36}; %	Dihedral of VTP (deg)
%VTAIL
y_loc_1R_y1_vee_CAD = Sheet11{37}; %	y loc of wing (V-tail) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_vee_CAD = Sheet11{38}; %	y loc of wing (V-tail) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_vee_CAD = Sheet11{39}; %	z loc of wing (V-tail) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_vee_CAD = Sheet11{40}; %	x loc of wing (V-tail) root chord LE position (distance from CAD refference point)
Lambda_LE_vee_e = Sheet11{42}; %	Sweep of V-tail (deg)
dihedral_vee_e = Sheet11{43}; %	Dihedral of V-tail (deg)
%Twin VTP
y_loc_1R_y1_2VTP_CAD = Sheet11{43}; %	y loc of Twin (VTP) root chord LE position (distance from CAD refference point)
z_loc_1R_y1_2VTP_CAD = Sheet11{44}; %	y loc of Twin (VTP) tip chord LE position (distance from CAD refference point)
z_loc_1R_y2_2VTP_CAD = Sheet11{45}; %	z loc of wing (VTP) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_2VTP_CAD = Sheet11{46}; %	x loc of Twin (VTP) root chord LE position (distance from CAD refference point)
Lambda_LE_2VTP_e = Sheet11{47}; %	Sweep of Twin VTP (deg)
dihedral_2VTP_e = Sheet11{48}; %	Dihedral of Twin VTP (deg)
cR_w1 = Sheet11{49}; %	Root Chord w1
cT_w1 = Sheet11{50}; %	Tip Chord w1
cR_w2 = Sheet11{51}; %	Root Chord w2
cT_w2 = Sheet11{52}; %	Tip Chord w2
cR_can = Sheet11{53}; %	Root Chord canard
cT_can = Sheet11{54}; %	Tip Chord canard
cR_HTP = Sheet11{55}; %	Root Chord HTP
cT_HTP = Sheet11{56}; %	Tip Chord HTP
cR_VTP = Sheet11{57}; %	Root Chord VTP
cT_VTP = Sheet11{58}; %	Tip Chord VTP
cR_vee = Sheet11{59}; %	Root Chord vee
cT_vee = Sheet11{60}; %	Tip Chord vee
cR_2VTP = Sheet11{61}; %	Root Chord Twin VTP
cT_2VTP = Sheet11{62}; %	Tip Chord Twin VTP
% Perfecntage (%) of the control surface
K_y1_ail_w1 = Sheet11{63}; %	% of aileron from effective wing surface - inner position
K_y2_ail_w1 = Sheet11{64}; %	% of aileron from effective wing surface - outter position
K_y1_ele_w2 = Sheet11{65}; %	% of elevator from effective HTP surface - inner position
K_y2_ele_w2 = Sheet11{66}; %	% of elevator from effective HTP surface  - outter position
K_y1_elevon_w1 = Sheet11{67}; %	% of elevon from effective wing surface - inner position
K_y2_elevon_w1 = Sheet11{68}; %	% of elevon from effective wing surface - outter position
K_y1_flap_w1 = Sheet11{69}; %	% of flap from effective wing surface - inner position
K_y2_flap_w1 = Sheet11{70}; %	% of flap from effective wing surface - outter position
K_y1_rudder_VTP = Sheet11{71}; %	% of rudder from effective VTP surface - inner position
K_y2_rudder_VTP = Sheet11{72}; %	% of rudder from effective VTP surface  - outter position
K_y1_rudvtr_w2 = Sheet11{73}; %	% of ruddervator from effective V-tail surface - inner position
K_y2_rudvtr_w2 = Sheet11{74}; %	% of ruddervator from effective V-tail surface  - outter position
K_y1_canard_can = Sheet11{75}; %	% of canard control surface from effective canard surface - inner position
K_y2_canard_can = Sheet11{76}; %	% of canard control surface from effective canard surface - outter position
% Maximum and minimum control surface deflections
delta_ail_min = Sheet11{77}; % 	Minimum aileron deflection
delta_ail_max = Sheet11{78}; %	Maximum aileron deflection
delta_ele_min = Sheet11{79}; %	Minimum elevator deflection
delta_ele_max = Sheet11{80}; %	Maximum elevator deflection
delta_elevon_min = Sheet11{81}; %	Minimum elevon deflection
delta_elevon_max = Sheet11{82}; %	Maximum elevon deflection
delta_flap_min = Sheet11{83}; %	Minimum flap deflection
delta_flap_max = Sheet11{84}; %	Maximum flap deflection
delta_rudder_min = Sheet11{85}; %	Minimum rudder deflection
delta_rudder_max = Sheet11{86}; %	Maximum rudder deflection
delta_rudvtr_min = Sheet11{87}; %	Minimum ruddervator deflection
delta_rudvtr_max = Sheet11{88}; %	Maximum ruddervator deflection
delta_can_min = Sheet11{89}; %	Minimum canard deflection
delta_can_max = Sheet11{90}; %	Maximum canard deflection
% Control surface dimensions
cf_ail = Sheet11{91}; %	%  of control surface aileron (chrodwise)
t_c_ail = Sheet11{92}; %	Thinckness 2 chord ratio associated to the airfoil - aileron
cf_ele = Sheet11{93}; %	%  of control surface elevator (chrodwise)
t_c_ele = Sheet11{94}; %	Thinckness 2 chord ratio associated to the airfoil - elevator
cf_elevon = Sheet11{95}; %	%  of control surface elevon (chrodwise)
t_c_elevon = Sheet11{96}; %	Thinckness 2 chord ratio associated to the airfoil - elevon
cf_flap = Sheet11{97}; %	%  of control surface flap (chrodwise)
t_c_flap = Sheet11{98}; %	Thinckness 2 chord ratio associated to the airfoil - flap
cf_rudder = Sheet11{99}; %	%  of control surface rudder (chrodwise)
t_c_rudder = Sheet11{100}; %	Thinckness 2 chord ratio associated to the airfoil - rudder
cf_rudvtr = Sheet11{101}; %	%  of control surface ruddervator (chrodwise)
t_c_rudvtr = Sheet11{102}; %	Thinckness 2 chord ratio associated to the airfoil - ruddervator
cf_canard = Sheet11{103}; %	%  of control surface canard (chrodwise)
t_c_canard = Sheet11{104}; %	Thinckness 2 chord ratio associated to the airfoil - canard
% Initial center of gravity estimaton
x_XCG = Sheet11{105};
y_XCG = Sheet11{106};
z_XCG = Sheet11{107};
% Location of Engine type 1 (symetry applied)
x_eng_ybar1 = Sheet11{108};
y_eng_ybar1 = Sheet11{109};
z_eng_ybar1 = Sheet11{110};
% Location of Engine type 1 (symetry applied)
x_eng_ybar2 = Sheet11{111};
y_eng_ybar2 = Sheet11{112};
z_eng_ybar2 = Sheet11{113};
% Location of the center of Leading Edge Nacelle of propulsion system 1
x_nac_ybar1 = Sheet11{114}; % Location of the center of LE (z) Nacelle of propulsion system 1	x_nac_ybar1
y_nac_ybar1 = Sheet11{115}; % Location of the center of LE (y) Nacelle of propulsion system 1	y_nac_ybar1
z_nac_ybar1 = Sheet11{116}; % Location of the center of LE (z) Nacelle of propulsion system 1	z_nac_ybar1
% Location of the center of Leading Edge Nacelle of propulsion system 2
x_nac_ybar2 = Sheet11{117}; % Location of the center of LE (z) Nacelle of propulsion system 2	x_nac_ybar1
y_nac_ybar2 = Sheet11{118}; % Location of the center of LE (y) Nacelle of propulsion system 2	y_nac_ybar1
z_nac_ybar2 = Sheet11{119}; % Location of the center of LE (z) Nacelle of propulsion system 2	z_nac_ybar1
% Wing Offset: Determines if the wing geometrý includes center section offset	
wing_offset_w1 = Sheet11{120}; % Wing Offset: Determines if the wing geometrý includes center section offset
canard_offset_can = Sheet11{121};% Canard Offset: Determines if the canard geometrý includes center section offset
zoffset_fuselage = Sheet11{122};% Fuselage Offset: Determines if the canard geometrý includes center section offset
l_tailboom = Sheet11{123};% Length of Tailboom
w_tailboom = Sheet11{124};% Width of Tailboom
h_tailboom = Sheet11{125};% Height of Tailboom
l_pod = Sheet11{126};% Length of nacelle
w_pod = Sheet11{127};% Width of nacelle
h_pod = Sheet11{128};% Height of nacelle
% Available control surfaces
d_ail = Sheet11{129}; %Definition of available control surface - aileron
d_ele = Sheet11{130}; %Definition of available control surface - elevator
d_elevon = Sheet11{131}; %Definition of available control surface - elevon
d_flap = Sheet11{132}; %Definition of available control surface - flap
d_rudder = Sheet11{133}; %Definition of available control surface - rudder
d_rudvtr = Sheet11{134}; %Definition of available control surface - ruddervator
d_can = Sheet11{135}; %Definition of available control surface - canard

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
InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; %
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
InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_2VTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_2VTP_CAD; %
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
InputGeometry_Data_flags.K_y2_ele_w2 = K_y2_ele_w2;
InputGeometry_Data_flags.K_y1_elevon_w1 = K_y1_elevon_w1;
InputGeometry_Data_flags.K_y2_elevon_w1 = K_y2_elevon_w1;
InputGeometry_Data_flags.K_y1_flap_w1 = K_y1_flap_w1;
InputGeometry_Data_flags.K_y2_flap_w1 = K_y2_flap_w1;
InputGeometry_Data_flags.K_y1_rudder_VTP = K_y1_rudder_VTP;
InputGeometry_Data_flags.K_y2_rudder_VTP = K_y2_rudder_VTP;
InputGeometry_Data_flags.K_y1_rudvtr_w2 = K_y1_rudvtr_w2;
InputGeometry_Data_flags.K_y2_rudvtr_w2 = K_y2_rudvtr_w2;
InputGeometry_Data_flags.K_y1_canard_can = K_y1_canard_can;
InputGeometry_Data_flags.K_y2_canard_can = K_y2_canard_can;
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
InputGeometry_Data_flags.cf_canard = cf_canard;
InputGeometry_Data_flags.t_c_canard = t_c_canard;
InputGeometry_Data_flags.x_XCG = x_XCG;
InputGeometry_Data_flags.y_XCG = y_XCG;
InputGeometry_Data_flags.z_XCG = z_XCG;
% Location of Engine type 1 (symetry applied)
InputGeometry_Data_flags.x_eng_ybar1 = x_eng_ybar1;
InputGeometry_Data_flags.y_eng_ybar1 = y_eng_ybar1;
InputGeometry_Data_flags.z_eng_ybar1 = z_eng_ybar1;
% Location of Engine type 1 (symetry applied)
InputGeometry_Data_flags.x_eng_ybar2 = x_eng_ybar2;
InputGeometry_Data_flags.y_eng_ybar2 = y_eng_ybar2;
InputGeometry_Data_flags.z_eng_ybar2 = z_eng_ybar2;
% Location of the center of Leading Edge Nacelle of propulsion system 1
InputGeometry_Data_flags.x_nac_ybar1 = x_nac_ybar1;
InputGeometry_Data_flags.y_nac_ybar1 = y_nac_ybar1;
InputGeometry_Data_flags.z_nac_ybar1 = z_nac_ybar1;
% Location of the center of Leading Edge Nacelle of propulsion system 1
InputGeometry_Data_flags.x_nac_ybar2 = x_nac_ybar2;
InputGeometry_Data_flags.y_nac_ybar2 = y_nac_ybar2;
InputGeometry_Data_flags.z_nac_ybar2 = z_nac_ybar2;
InputGeometry_Data_flags.wing_offset_w1 = wing_offset_w1;
InputGeometry_Data_flags.canard_offset_can = canard_offset_can;
InputGeometry_Data_flags.zoffset_fuselage = zoffset_fuselage;
% Tailboom geometry
InputGeometry_Data_flags.l_tailboom = l_tailboom;
InputGeometry_Data_flags.w_tailboom = w_tailboom;
InputGeometry_Data_flags.h_tailboom = h_tailboom;
% Nacelle Geometry
InputGeometry_Data_flags.l_pod = l_pod;
InputGeometry_Data_flags.w_pod = w_pod;
InputGeometry_Data_flags.h_nac = h_pod;
% Available control surfaces
InputGeometry_Data_flags.d_ail = d_ail; %Definition of available control surface - aileron
InputGeometry_Data_flags.d_ele = d_ele; %Definition of available control surface - elevator
InputGeometry_Data_flags.d_elevon = d_elevon; %Definition of available control surface - elevon
InputGeometry_Data_flags.d_flap = d_flap; %Definition of available control surface - flap
InputGeometry_Data_flags.d_rudder = d_rudder; %Definition of available control surface - rudder
InputGeometry_Data_flags.d_rudvtr = d_rudvtr; %Definition of available control surface - ruddervator
InputGeometry_Data_flags.d_can = d_can; %Definition of available control surface - canard

% Gathers all the flags
OUTPUT_read_XLSX.InputGeometry_Data_flags = InputGeometry_Data_flags;

%% PLOTS
Sheet12 = readcell(excel_file,'Sheet','12Plots','Range','D3:D25');

plot(1) = Sheet12{1};% Prints plots - Aero
plot(2) = Sheet12{2};% Prints plots - Aero Polar
plot(3) = Sheet12{3};% Print Plots for Stimation ofXAC
plot(4) = Sheet12{4};% Print Plots for Propulsive Models
plot(5) = Sheet12{5};% Prints plots for 3D
plot(6) = Sheet12{6};% Prints plots of SM analysis
plot(7) = Sheet12{7};% Prints plots of longitudinal Trim
plot(8) = Sheet12{8};% Prints plots of longitudinal Trim with Variable mass & Variable Speed
plot(9) = Sheet12{9};% Prints plots of lateral Trim
plot(10) = Sheet12{10};% Prints plots of lateral Turning
plot(11) = Sheet12{11};% Prints PLOTS STABILITY DERIVATIVES FOR VAR MASS AND VELOCITY
plot(12) = Sheet12{12};% Prints PLOTS STABILITY ANALYSIS FOR VAR MASS AND VELOCITY
plot(13) = Sheet12{13};% Prints PLOTS DYNAMIC STABILITY ANALYSIS AFTER IMPULSE
plot(14) = Sheet12{14};% Prints PLOTS PERFORMANCE STUDY
plot(15) = Sheet12{15};% Prints plots of longitudinal Trim with V-n diagram n & Variable Speed
prefix = Sheet12{16};% Selection of String Characters that Define the AC selection
plot_individuals = Sheet12{17};% Shows plots individually or collecting similar plots	
fname = Sheet12{18};% Selection of String Characters that Define the AC selection storing location	
LS = Sheet12{19}; % Line size	LS
FS = Sheet12{20}; % Text Font size	FS
LFS = Sheet12{21}; % Legend Font Size	LFS
Video_3D = Sheet12{22}; % Saves video 3D	Video_3D
SAVE_FIGS = Sheet12{23}; % saves the plots: fig, jpg and pdf	SAVE_FIGS

% Stores the flags
PLOT_flags.plot = plot; %
PLOT_flags.prefix = prefix; %
PLOT_flags.plot_individuals = plot_individuals; %
PLOT_flags.fname = fname; %
PLOT_flags.LS = LS; % Line size	LS
PLOT_flags.FS = FS; % Text Font size	FS
PLOT_flags.LFS = LFS; % Legend Font Size	LFS
PLOT_flags.Video_3D = Video_3D; % Saves video 3D	Video_3D
PLOT_flags.SAVE_FIGS = SAVE_FIGS; % saves the plots: fig, jpg and pdf	SAVE_FIGS

% Gathers all the flags
OUTPUT_read_XLSX.PLOT_flags = PLOT_flags;

%% Performance Inpute Preliminar
Sheet13 = readcell(excel_file,'Sheet','13InputPerforInitial','Range','D3:D68');

% Taxy
temp_local_taxy = Sheet13{1};% 	TAXY: TEMPERATURA LOCAL (Celsius)
h_inicial_taxy = Sheet13{2};% 	TAXY: INITIAL ALTITUDE
P_local_taxy = Sheet13{3};% 	TAXY: Local Pressure
delta_T_taxy = Sheet13{4};% 	TAXY: PALANCA DE RALENTI EN TAXI = 0.05
V_taxy = Sheet13{5};% 	TAXY:  VELOCIDAD A LA QUE HACE EL TAXI (m/s)
t_taxy = Sheet13{6};% 	TAXY: TIEMPO DE ESPERA EN TAXI (s)
% TakeOff
temp_local_TO = Sheet13{7};% 	TAKEOFF: TEMPERATURA LOCAL (Celsius)
h_inicial_TO = Sheet13{8};% 	TAKEOFF: ALTURA LOCAL (m)
P_local_TO = Sheet13{9};% 	TAKEOFF: PRESION LOCAL (Pa)
mu_TO = Sheet13{10};% 	TAKEOFF: COEFICIENTE DE FRICCION CON LA PISTA (MU)
h_obstacle_TO = Sheet13{11};% 	TAKEOFF: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
gamma_climb_TO = Sheet13{12};% 	TAKEOFF: GAMMA DE SUBIDA MINIMO
delta_T_TO = Sheet13{13};% 	TAKEOFF: PALANCA DE GASES PARA DESPEGUE
% Climb
h_inicial_cl = Sheet13{14};% 	CLIMB : ALTURA INICIAL - [m]
h_final_cl = Sheet13{15};% 	CLIMB : ALTURA FINAL - [m]
gamma_cl = Sheet13{16};% 	CLIMB : GAMMA DE SUBIDA - [-]
Mach_cl = Sheet13{17};% 	CLIMB : MACH DE VUELO  - [-]
TAS_cl = Sheet13{18};% 	CLIMB : VELOCIDAD TAS  - [m/s]
EAS_cl = Sheet13{19};% 	CLIMB : VELOCIDAD EAS  - [m/s]
delta_T_cl = Sheet13{20};%	CLIMB : PALANCA DE GASES  - [-]
V_ini_cl = Sheet13{21};% 	CLIMB : VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
V_fin_cl = Sheet13{22};% 	CLIMB : VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
% Cruise
h_inicial_cr = Sheet13{23};%	CRUISE: ALTURA INICIAL
dist_final_cr = Sheet13{24};%	CRUISE: DISTANCIA FINAL
V_cr = Sheet13{25};%	CRUISE: VELOCIDAD
delta_T_cr = Sheet13{26};%	CRUISE:  PALANCA DE GASES
V_ini_cr = Sheet13{27};%	CRUISE:  VELOCIDAD INICIAL
V_fin_cr = Sheet13{28};%	CRUISE: VELOCIDAD FINAL
fuel_cr = Sheet13{29};% - [kg] % 7: COMBUSTIBLE A QUEMAR
Cd0_cr = Sheet13{30};% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
k1_cr = Sheet13{31};% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
k2_cr = Sheet13{32};% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
% Turn
h_inicial_tr = Sheet13{33};%	TURN:  ALTURA INICIAL
t_final_tr = Sheet13{34};%	TURN: TIEMPO FINAL (seg)
V_turn  = Sheet13{35};% TURN: turn velocity
delta_T_tr = Sheet13{36};%	TURN: PALANCA DE GASES
phi_tr = Sheet13{37};% 	TURN: ANGULO DE ALABEO (rads)
V_psi = Sheet13{38};%	TURN: VELOCIDAD DE GUIÑADA (rads/seg)
n_tr = Sheet13{39};%	TURN: FACTOR DE CARGA
R_tr = Sheet13{40};%	TURN: RADIO DE GIRO (m)
% Descent
h_inicial_d = Sheet13{41};%	DESCENT: ALTURA INICIAL
h_final_d = Sheet13{42};%	DESCENT:  ALTURA FINAL
gamma_d = Sheet13{43};%	DESCENT: GAMMA
V_d = Sheet13{44};%	DESCENT: VELOCIDAD DE DESCENSI
EAS_d = Sheet13{45};%	DESCENT : VELOCIDAD TAS  - [m/s]
TAS_d = Sheet13{46};%	DESCENT : VELOCIDAD EAS  - [m/s]
delta_T_d = Sheet13{47};%	DESCENT : PALANCA DE GASES  - [-]
V_ini_d = Sheet13{48};%	DESCENT : VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
V_fin_d = Sheet13{49};%	DESCENT : VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
% Turn Wait
h_inicial_tr_wt = Sheet13{50};%	TURN WAIT:  ALTURA INICIAL
t_final_tr_wt = Sheet13{51};%	TURN WAIT: TIEMPO FINAL (seg)
V_turn_wt = Sheet13{52};%	TURN WAIT: turn velocity
delta_T_tr_wt = Sheet13{53};%	TURN WAIT: PALANCA DE GASES
phi_tr_wt = Sheet13{54};%	TURN WAIT: ANGULO DE ALABEO (rads)
V_psi_wt = Sheet13{55};%	TURN WAIT: VELOCIDAD DE GUIÑADA (rads/seg)
n_tr_wt = Sheet13{56};%	TURN WAIT: FACTOR DE CARGA
R_tr_wt = Sheet13{57};%	TURN WAIT: RADIO DE GIRO (m)
% lANDING
temp_local_LND = Sheet13{58};%	LANDING:  TEMPERATURA LOCAL (Celsius)
h_inicial_LND = Sheet13{59};%	LANDING: ALTURA LOCAL (m)
P_local_LND = Sheet13{60};%	LANDING:  PRESION LOCAL (Pa)
mu_LND = Sheet13{61};%	LANDING: COEFICIENTE DE FRICCION CON LA PISTA (MU)
delta_T_LND = Sheet13{62};%	LANDING:  PALANCA DE REVERSA
t_brake = Sheet13{63};%	LANDING: 'Tiempo en activar frenos' - [s]
% Dummy Segment
dummy(1) = Sheet13{64};%	DUMMY SEGMENT: ALTURA INICIAL (CRUISE)
dummy(2) = Sheet13{65};%	DUMMY SEGMENT: DISTANCIA FINAL
dummy(3) = Sheet13{66};%	DUMMY SEGMENT: VELOCIDAD (CRUISE)

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
Sheet14 = readcell(excel_file,'Sheet','14PerfoMision_Selection','Range','D3:D8');

type_missions_WF = Sheet14{1};%	Selection Type of Mission
num_missions_WF = Sheet14{2};%	Number of Segments
climb_mode = Sheet14{3};%	Selects Climb Mode Option
cruise_mode = Sheet14{4};%	Selects Cruise Mode Option
turn_mode = Sheet14{5};%	Selects Turn Mode Option
descent_mode = Sheet14{6};%	Selects Descent Mode Option

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
Sheet15 = readcell(excel_file,'Sheet','15PerformanceVariable','Range','D3:D15');

V_low = Sheet15{1};% 	Performance Variable study - Low velocity  in m/s
V_high = Sheet15{2};% 	Performance Variable study - High velocity  in m/s
N_V_VAR_perf = Sheet15{3};% Performance Variable study for Velocity - Number of variable points
V_single = Sheet15{4};%	Performance Study - Single Speed
Wp_low = Sheet15{5};%	Performance Variable study - Low weight  in m/s
Wp_high = Sheet15{6};%	Performance Variable study - High weight  in m/s
N_Wp_VAR_perf = Sheet15{7};%	Performance Variable study for Weight - Number of variable points
W_single = Sheet15{8};%	Performance Study - Single Weight
Post_processing_PERFORMANCE = Sheet15{9};%	Conducts the Post processing PERFORMANCE
climb_mode = Sheet15{10};%	Selects Climb Mode Option
cruise_mode = Sheet15{11};%	Selects Cruise Mode Option
turn_mode = Sheet15{12};%	Selects Turn Mode Option
descent_mode = Sheet15{13};%	Selects Descent Mode Option

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
