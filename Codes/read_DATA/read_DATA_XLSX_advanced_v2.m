function OUTPUT_read_XLSX =  read_DATA_XLSX_advanced_v2(case_AC)

%% Pre-load variables stored in Excel
% Function that reads all the data for the different aircraft
get_add_path
filename = '../AIRCRAFT/SANAID_AIRCRAFT_DATA_v9.xlsx';

excel_file = 'SANAID_AIRCRAFT_DATA_v9.xlsx';
    
% Writes in the excel to assign the propper aircraft type
sheet = 3;
status = xlswrite(filename,case_AC,sheet,'D8');

%% MATLAB Flags
% Sheet1 = readcell('SANAID_AIRCRAFT_DATA_v3.xlsx','Sheet','1MATLAB','Range','D3:D6');
Sheet1 = readcell(excel_file,'Sheet','1MATLAB','Range','D3:D8');
el=1;

MATLAB_in = Sheet1{el};el=el+1;
CHECK_Efficiency = Sheet1{el};el=el+1;
Detailed_Profile = Sheet1{el};el=el+1;
save_MAT = Sheet1{el};el=el+1;
show_messages_screen = Sheet1{el};el=el+1;
write_data = Sheet1{el};el=el+1;

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
Sheet2 = readcell(excel_file,'Sheet','2Studies','Range','D3:D28');
el=1;

AERODYNAMIC_STUDY = Sheet2{el};el=el+1; % Conducts Aerodynamic Studies
PROP_STUDY = Sheet2{el};el=el+1; % Conducts Proppeller Studies
PERFORMANCE_STUDY = Sheet2{el};el=el+1; % Conducts Performance Studies
STABILITY_STUDY = Sheet2{el};el=el+1; % Conducts Stability Studies
MISSIONS_STUDY = Sheet2{el};el=el+1; % Conducts Mission Studies
STUDY_Weight_Fraction = Sheet2{el};el=el+1; % Weight Fraction Study
ANALYSIS_PERFORMANCE_AP = Sheet2{el};el=el+1; % Performance Analysis integrated with AP Codes
ANALYSIS_PERFORMANCE_AP_var = Sheet2{el};el=el+1; % Performance Analysis integrated with AP Codes varying Conditions
variable_speed_AP = Sheet2{el};el=el+1; % Variable Speed Performance Studies
variable_weight_AP = Sheet2{el};el=el+1; % Variable Mass Performance Studies
STABILITY_STUDY_Trim = Sheet2{el};el=el+1; % Trim conditions - one case
STABILITY_STUDY_Regular = Sheet2{el};el=el+1; % Stability Analysis Regular Study
STABILITY_STUDY_Trim_varV_m0 = Sheet2{el};el=el+1; % Trim conditions - variation V and mass
STABILITY_STUDY_Trim_var_XCG = Sheet2{el};el=el+1; % Trim conditions - variation XCG
STABILITY_STUDY_Long_dyn = Sheet2{el};el=el+1; % Longitudinal Stability Analysis
STABILITY_STUDY_LatDir_dyn = Sheet2{el};el=el+1; % Lateral-Directional Stability Analysis
STABILITY_STUDY_Trim_lat = Sheet2{el};el=el+1; % Lateral Directional Trim
STABILITY_STUDY_Trim_lat_asymmetries = Sheet2{el};el=el+1; % Lateral Directional Trim - Assymetries
STABILITY_STUDY_Trim_lat_accelerations = Sheet2{el};el=el+1; % Lateral Directional Trim - Accelertations
STABILITY_STUDY_Trim_lat_Trim_TAB = Sheet2{el};el=el+1; % Lateral Directional Trim - Trim Tab
STABILITY_STUDY_Turning = Sheet2{el};el=el+1; % Turning Stability
STABILITY_STUDY_Stability_Derivatives_varV_m0 = Sheet2{el};el=el+1; % Plots Stability Derivatives
STABILITY_STUDY_Stability_Analysis_varV_m0_long = Sheet2{el};el=el+1; % Longitudinal Dynamic Response Plots
STABILITY_STUDY_Stability_Analysis_varV_m0_lat = Sheet2{el};el=el+1; % Lateral-Directional Dynamic Response Plots
STABILITY_STUDY_Stability_Analysis_dyna_impulse_long = Sheet2{el};el=el+1; % Impulse Longitudinal Dynamic Response Plots
STABILITY_STUDY_Stability_Analysis_dyna_impulse_lat = Sheet2{el};el=el+1; % Impulse Lateral-Directional Dynamic Response Plots

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
STUDY_flags.STABILITY_STUDY_Trim_lat_asymmetries = STABILITY_STUDY_Trim_lat_asymmetries;% Lateral Directional Trim - Assymetries 
STUDY_flags.STABILITY_STUDY_Trim_lat_accelerations = STABILITY_STUDY_Trim_lat_accelerations;% Lateral Directional Trim - Accelertations
STUDY_flags.STABILITY_STUDY_Trim_lat_Trim_TAB = STABILITY_STUDY_Trim_lat_Trim_TAB;% Lateral Directional Trim - Trim Tab
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
el=1;

SF = Sheet3{el};el=el+1; % Scaling Factor
Weight_Estimation = Sheet3{el};el=el+1; % Weight Estimation
AC_type = Sheet3{el};el=el+1; % Aircaft Type
CASE_fuse = Sheet3{el};el=el+1;% Determines the fuselage than will be shown
ESCALADO = Sheet3{el};el=el+1; % Determine the type of scaling
case_AC = Sheet3{el};el=el+1; % Aircraft Model Analized
twin_VTP = Sheet3{el};el=el+1; % Defines twin VTP

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
Sheet4 = readcell(excel_file,'Sheet','4Propulsion','Range','D3:D46');
el=1;

propul(1) = Sheet4{el};el=el+1; % Type of engine
propul(2) = Sheet4{el};el=el+1; % Number of engines
propul(3) = Sheet4{el};el=el+1; % EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine
propul(4) = Sheet4{el};el=el+1; % CONSUMO ESPECIFICO: % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
propul(5) = Sheet4{el};el=el+1; % DERIVACION(TURBOFANES) By-pass 
propul(6) = Sheet4{el};el=el+1; % EFICIENCIA DE LA HELICE (ETA_P)
propul(7) = Sheet4{el};el=el+1; % Normativa Propulsion
propul(8) = Sheet4{el};el=el+1; % CAPACIDAD CALORIFICA DEL COMBUSTIBLE
propul(9) = Sheet4{el};el=el+1; % DIAMETRO DEL MOTOR(TURBOHELICES)
type_battery = Sheet4{el};el=el+1; % Type Battery Used for Propulsion
Engine_loc = Sheet4{el};el=el+1; % Engine location
Engine_conf = Sheet4{el};el=el+1; % Engine configuration
D_prop = Sheet4{el};el=el+1; % Preliminary Prop Diameter (m)
RPM_min = Sheet4{el};el=el+1; % RPM limits (if known) - LOW LIMIT	
RPM_max = Sheet4{el};el=el+1; % RPM limits (if known) - HIGH LIMIT	
beta_pitch = Sheet4{el};el=el+1; % Pitch angle	
beta_variable_pitch = Sheet4{el};el=el+1; % Use of tables of family of variable pitch propellers
b_prop = Sheet4{el};el=el+1; % Number of blades	to account for correction for more than 2 blades
prop_known = Sheet4{el};el=el+1; % Flag that determines if the prop is known
prop_properties = Sheet4{el};el=el+1; % model of prop:first 2 numbers diameter in inches last 2 numbers pitch
type_prop = Sheet4{el};el=el+1; % Type of propller
Prop_type = Sheet4{el};el=el+1; % Propulsive Model- Number of Prop used of propller
l_eng = Sheet4{el};el=el+1; % Engine geometry - length
d_eng = Sheet4{el};el=el+1; % Engine geometry - diameter
l_nc = Sheet4{el};el=el+1; % Nacelle geometry - length
d_nc = Sheet4{el};el=el+1; % Nacelle geometry - diameter
model_prop = Sheet4{el};el=el+1; %	Compares 3 prop models
prop_selec_APC = Sheet4{el};el=el+1; %	Model 1 - APC data
prop_selec_WT1 = Sheet4{el};el=el+1; %	Model 2 - Wind tunnel data for different props
prop_selec_WT2 = Sheet4{el};el=el+1; %	Model 3 - Wind tunnel data for different angle of attack
eta_gear = Sheet4{el};el=el+1; % Gearbox efficiency
eta_m = Sheet4{el};el=el+1; % Motor efficiency
eta_esc = Sheet4{el};el=el+1; % ESC efficiency
eta_dist = Sheet4{el};el=el+1; % Electrical distribution efficiency
SE_LiFePO4 = Sheet4{el};el=el+1; % Specific Energy of LiFePO4 batteries 90–160 Wh/kg (320–580 J/g or kJ/kg) 
SE_LiPo = Sheet4{el};el=el+1; % Specific Energy of LiPo batteries (Wh/kg)
SE_FuelCells = Sheet4{el};el=el+1; % Wh/kg Specific Energy for fuel cells
SP_motor = Sheet4{el};el=el+1; % Motor Specific Power kW/kg - TIER 1
SP_ESC = Sheet4{el};el=el+1; % ESC Specific Power kW/kg - TIER 1
m_prop_din = Sheet4{el};el=el+1; % Mass of propeller as a function of diameter in inches
compare_prop = Sheet4{el};el=el+1; % Determinee Plots to compare different propulsion models
density_fuel = Sheet4{el};el=el+1; % Fuel density
cost_fuel = Sheet4{el};el=el+1; % Fuel Cost  cts/gal
CI = Sheet4{el};el=el+1; % Cost Index Kg/s

% Stores the flags
Propulsive_flags.propul = propul; %
Propulsive_flags.type_battery = type_battery; %
Propulsive_flags.Engine_loc = Engine_loc; %
Propulsive_flags.Engine_conf = Engine_conf; %
Propulsive_flags.D_prop = D_prop; %
Propulsive_flags.RPM_min = RPM_min; %
Propulsive_flags.RPM_max = RPM_max; %
Propulsive_flags.beta_pitch = beta_pitch; %
Propulsive_flags.beta_variable_pitch = beta_variable_pitch; %
Propulsive_flags.b_prop = b_prop; %
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
el=1;

W_eng = Sheet5{el};el=el+1; % Mass of engine (kg)
m_prop = Sheet5{el};el=el+1; % Mass of proppeller (kg)
m_ESC = Sheet5{el};el=el+1; % Mass of ESC (kg)
n_servos = Sheet5{el};el=el+1; % Number of servos 
m_servo = Sheet5{el};el=el+1; % Mass per servo (kg)
m_wiring = Sheet5{el};el=el+1; % Mass of electronic wiring (kg)
m_metal = Sheet5{el};el=el+1; % Mass of misc (kg)
m_systems = Sheet5{el};el=el+1; % Mass of Systems (kg)
m_landing_gear = Sheet5{el};el=el+1; % Mass of landing gear (kg)
n_pax = Sheet5{el};el=el+1; % Number of Passangers
n_crew = Sheet5{el};el=el+1; % Number of Crew
m_pax = Sheet5{el};el=el+1; % Mass per Passenger (kg)
m_crew = Sheet5{el};el=el+1; % Mass per Crew (kg)
m_luggage = Sheet5{el};el=el+1; % Mass of luggage (kg)
m_cargo = Sheet5{el};el=el+1; % Mass of payload cargo (kg)
m_batteries = Sheet5{el};el=el+1; % Mass of batteries (kg)
m_fuel = Sheet5{el};el=el+1; % Mass of fuel (kg) 
flag_landing_gear = Sheet5{el};el=el+1; % Determination of Landing Gear Estimation if available
f_f = Sheet5{el};el=el+1; % Fudge Factor that increments the weight estimation
m_cell_A123 = Sheet5{el};el=el+1; % Mass of one A123 cell
m_w1 = Sheet5{el};el=el+1; % Mass of one wing
m_HTP = Sheet5{el};el=el+1; % Mass of HTP
m_VTP = Sheet5{el};el=el+1; % Mass of VTP
m_Can = Sheet5{el};el=el+1; % Mass of Canard
m_Vee = Sheet5{el};el=el+1; % Mass of Vee Tail
m_Vee2 = Sheet5{el};el=el+1; % Mass of Vee Tail
m_fus_fairing = Sheet5{el};el=el+1; % Mass of fuselage and fairing
m_prop_fairing = Sheet5{el};el=el+1; % Mass of propeller nacelle
I_xx = Sheet5{el};el=el+1; % Moments of Inertia Ixx
I_yy = Sheet5{el};el=el+1; % Moments of Inertia Ixx
I_zz = Sheet5{el};el=el+1; % Moments of Inertia Ixx
I_xy = Sheet5{el};el=el+1; % Moments of Inertia Ixx
I_xz = Sheet5{el};el=el+1; % Moments of Inertia Ixx
I_yz = Sheet5{el};el=el+1; % Moments of Inertia Ixx
flag_moments_inertia = Sheet5{el};el=el+1; % Determination of Moments of innertia from literature	
flag_total_weights = Sheet5{el};el=el+1; % Determination if using the true weights
MTOW_true = Sheet5{el};el=el+1; % True MTOW - Maximum Takeke Off Weight
ME_true = Sheet5{el};el=el+1; % True ME - Empty Weight
MCREW_true = Sheet5{el};el=el+1; % True MCREW - Crew Weight
MLW_true = Sheet5{el};el=el+1; % True MLW - Máxcimum landing Weight (85% MTOW)
MF_true = Sheet5{el};el=el+1; % True MTOW - Maximum Fuel Weight
MPL_true = Sheet5{el};el=el+1; % True MTOW - Maximum Paylod Weigh
fraction_MF = Sheet5{el};el=el+1; % Fraction of fuel
fraction_PL = Sheet5{el};el=el+1; % Fraction of payload
get_shift_XCG_variation = Sheet5{el};el=el+1; % Variation of XCG study	

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
Weights_flags.m_Vee2 = m_Vee2; %
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
Sheet6 = readcell(excel_file,'Sheet','6CAD_data','Range','D3:D22');
el=1;

nSections = Sheet6{el};el=el+1; % Number of sectins (x coordinate) (eliminated 1 for convergence)
nPoints = Sheet6{el};el=el+1; % number of elements per section (y and z coordinates)
lecture_Geo = Sheet6{el};el=el+1; % lineas a partir desde donde empieza a leer
STL_PLOT = Sheet6{el};el=el+1; % Generates Fuselage mesh from CAD STL
XFLR5_file = Sheet6{el};el=el+1; % Defines flies to be used for the Fuselage geometry - XFLR5
STL_ac = Sheet6{el};el=el+1; % Defines flies to be used for the aircraft geometry - STL
STL_fus = Sheet6{el};el=el+1; % Defines flies to be used for the Fuselage geometry - STL
STL_wing = Sheet6{el};el=el+1; % Defines flies to be used for the wing geometry - STL
STL_canard = Sheet6{el};el=el+1; % Defines flies to be used for the canard geometry - STL
STL_HTP = Sheet6{el};el=el+1; % Defines flies to be used for the Fuselage geometry - STL
STL_VTP = Sheet6{el};el=el+1; % Defines flies to be used for the HTP geometry - STL
STL_Vee = Sheet6{el};el=el+1; % Defines flies to be used for the VTP geometry - STL
STL_Vee2 = Sheet6{el};el=el+1; % Defines flies to be used for the VTP geometry - STL
STL_engine = Sheet6{el};el=el+1; % Defines flies to be used for the engine geometry - STL
STL_nacelle = Sheet6{el};el=el+1; % Defines flies to be used for the nacelle geometry - STL
SF_CAD = Sheet6{el};el=el+1; % Defines Scaling Factor for STL files for units conversion (example mm 2 m)	
AC_STL_Compare = Sheet6{el};el=el+1; % % Plots both Generated Geometry and Aircraft STL comparisso	AC_STL_Compare
PLOT_Slices_fuselage_STL = Sheet6{el};el=el+1; % Plotsfuselage slices to check if they are rotated accordingly	
x_CAD_OFSET_Plot = Sheet6{el};el=el+1; % CAD Offset for plotting		
CAD_Kink_Wing = Sheet6{el};el=el+1; % Plots Wing with the Kink	
STL_reading_test = Sheet6{el};el=el+1; % New STL Method - Ion test only for Cessna 208	

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
Fuselage_flags.STL_Vee2 = STL_Vee2; %
Fuselage_flags.STL_engine = STL_engine; %
Fuselage_flags.STL_nacelle = STL_nacelle; %
Fuselage_flags.SF_CAD = SF_CAD; %
Fuselage_flags.AC_STL_Compare = AC_STL_Compare; %
Fuselage_flags.PLOT_Slices_fuselage_STL = PLOT_Slices_fuselage_STL; %
Fuselage_flags.x_CAD_OFSET_Plot = x_CAD_OFSET_Plot; %
Fuselage_flags.CAD_Kink_Wing = CAD_Kink_Wing; %	
Fuselage_flags.STL_reading_test = STL_reading_test; %	

% Gathers all the flags
OUTPUT_read_XLSX.Fuselage_flags = Fuselage_flags;

%% Performance Data Prelimina
Sheet7 = readcell(excel_file,'Sheet','7Performance','Range','D3:D15');
el=1;

h_climb = Sheet7{el};el=el+1; % Preliminar Performance conditions: Climb altitude (m)
Endurance_v = Sheet7{el};el=el+1; % Preliminar Performance conditions: Endurance in Hoovering (min)
h = Sheet7{el};el=el+1; % Preliminar Performance conditions: Cruise initial altitude (m)
V = Sheet7{el};el=el+1; % Preliminar Performance conditions: Cruise initial airspeed (m/s)
V_max = Sheet7{el};el=el+1; % Preliminar Performance conditions: Cruise max speed
Range = Sheet7{el};el=el+1; % Preliminar Performance conditions: Cruise Range (m)
Endurance = Sheet7{el};el=el+1; % Preliminar Performance conditions: Cruise Endurance (min)
Flight_SF = Sheet7{el};el=el+1; % Flight Safety Margin - normal 
Flight_cruise = Sheet7{el};el=el+1; % Flight Condition - Cruise 	 
Flight_takeoff = Sheet7{el};el=el+1; % Flight Condition - TakeOff	 
Flight_climb = Sheet7{el};el=el+1; % Flight Condition - Climb
N_V_VAR_perfo = Sheet7{el};el=el+1; % Number of iteration variation Speed
N_m_VAR_perfo = Sheet7{el};el=el+1; % Number of iteration variation Weight

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
Performance_pre_flags.N_V_VAR_perfo = N_V_VAR_perfo; %
Performance_pre_flags.N_m_VAR_perfo = N_m_VAR_perfo; %


% Gathers all the flags
OUTPUT_read_XLSX.Performance_pre_flags = Performance_pre_flags;

%% Aerodynamic Data
Sheet8A = readcell(excel_file,'Sheet','8AeroA','Range','D3:D37');

el=1;

index_w1 = Sheet8A{el};el=el+1; % Selection of TXT that are used for the aerodynamic analysis (w1) - LLT
index_w2 = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (w2) - LLT
index_vee = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (vee) - LLT
index_vee2 = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (vee) - LLT
index_w3 = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (w3) - LLT
index_VTP = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (vtp) - LLT
index_vee_cy = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (vtail-sideslip) - LLT
index_vee2_cy = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (vtail-sideslip) - LLT
index_w1_VLM = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (w1) - VLM
index_w2_VLM = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (w2) - VLM
index_vee_VLM = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (w2) - VLM
index_vee2_VLM = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (w2) - VLM
index_w3_VLM = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (w3) - VLM
index_VTP_VLM = Sheet8A{el}; el=el+1; % Selection of TXT that are used for the aerodynamic analysis (VTP) - VLM
i_w1 = Sheet8A{el}; el=el+1; % Selects the AoA of each surface relative to fuselage line (w1)
i_w2 = Sheet8A{el}; el=el+1; % Selects the AoA of each surface relative to fuselage line (w2)
i_vee = Sheet8A{el}; el=el+1; % Selects the AoA of each surface relative to fuselage line (w2)
i_vee2 = Sheet8A{el}; el=el+1; % Selects the AoA of each surface relative to fuselage line (w2)
i_w3 = Sheet8A{el}; el=el+1; % Selects the AoA of each surface relative to fuselage line (w3)
i_VTP = Sheet8A{el}; el=el+1; % Selects the AoA of each surface relative to fuselage line (w3)
compare_plot_aero = Sheet8A{el}; el=el+1;% Determines the aero plots to compare
read_XFLR5 = Sheet8A{el}; el=el+1; % Reads Aero data from XFLR5	
read_FLOW5 = Sheet8A{el}; el=el+1; % Reads Aero data from FLOW5	
airfoil_w1 = Sheet8A{el}; el=el+1; % Airfoil for wing (primary)	
airfoil_w2 = Sheet8A{el}; el=el+1; % Airfoil for wing (secondary)	
airfoil_c1 = Sheet8A{el}; el=el+1; % Airfoil for canard (primary)	
airfoil_c2 = Sheet8A{el}; el=el+1; % Airfoil for canard (secondary)	
airfoil_HTP1 = Sheet8A{el}; el=el+1; % Airfoil for HTP (primary)	
airfoil_HTP2 = Sheet8A{el}; el=el+1; % Airfoil for HTP (secondary)
airfoil_VTP1 = Sheet8A{el}; el=el+1; % Airfoil for VTP (primary)	
airfoil_VTP2 = Sheet8A{el}; el=el+1; % Airfoil for VTP (secondary)	
airfoil_Vee1 = Sheet8A{el}; el=el+1; % Airfoil for VeeTail (primary)	
airfoil_Vee2 = Sheet8A{el}; el=el+1; % Airfoil for VeeTail (secondary)	
polar_model = Sheet8A{el}; el=el+1; % Use Global Weights	flag_total_weight	
fuse_aero_FLOW_and_CBM = Sheet8A{el}; % Use Global Weights	flag_total_weight	

% Stores the flags
Aerodynamic_Data_flags.index_w1 = index_w1; %
Aerodynamic_Data_flags.index_w2 = index_w2; %
Aerodynamic_Data_flags.index_vee = index_vee; %
Aerodynamic_Data_flags.index_vee2 = index_vee2; %
Aerodynamic_Data_flags.index_w3 = index_w3; %
Aerodynamic_Data_flags.index_VTP = index_VTP; %
Aerodynamic_Data_flags.index_vee_cy = index_vee_cy; %
Aerodynamic_Data_flags.index_vee2_cy = index_vee2_cy; %
Aerodynamic_Data_flags.index_w1_VLM = index_w1_VLM; %
Aerodynamic_Data_flags.index_w2_VLM = index_w2_VLM; %
Aerodynamic_Data_flags.index_vee_VLM = index_vee_VLM; %
Aerodynamic_Data_flags.index_vee2_VLM = index_vee2_VLM; %
Aerodynamic_Data_flags.index_w3_VLM = index_w3_VLM; %
Aerodynamic_Data_flags.index_VTP_VLM = index_VTP_VLM; %
Aerodynamic_Data_flags.i_w1 = i_w1; %
Aerodynamic_Data_flags.i_w2 = i_w2; %
Aerodynamic_Data_flags.i_vee = i_vee; %
Aerodynamic_Data_flags.i_vee2 = i_vee2; %
Aerodynamic_Data_flags.i_w3 = i_w3; %
Aerodynamic_Data_flags.i_VTP = i_VTP; %
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

%% Aerodynamic Data
Sheet8B = readcell(excel_file,'Sheet','8AeroB','Range','D3:D22');
el=1;

w1 = Sheet8B{el}; el=el+1; % Selection of elements for theoretical polar estimation - wing	w1
h = Sheet8B{el}; el=el+1; % Selection of elements for theoretical polar estimation - HTP	h
v = Sheet8B{el}; el=el+1; % Selection of elements for theoretical polar estimation - VTP	v
v2 = Sheet8B{el}; el=el+1; % Selection of elements for theoretical polar estimation - twin VTP	v2
vtail = Sheet8B{el}; el=el+1; % Selection of elements for theoretical polar estimation - Vtail	vtail
vtail2 = Sheet8B{el}; el=el+1; % Selection of elements for theoretical polar estimation - Vtail	vtail
can = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - canard	can
fus = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - fuselage	fus
m_fus = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - multiple fuselage	m_fus
n_m_fus = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - number of multiple fuselage	n_m_fus
nac = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - nacelle	nac
landgear = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - landing gear	landgear
tailboom = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - tailboom
m_tailboom = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - multiple tailboom
n_m_tailboom = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - number of multiple tailboom
missile = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - misile
n_missile = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - misile
pod = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - pod
flap = Sheet8B{el}; el=el+1;% Selection of elements for theoretical polar estimation - flap
Flap_type = Sheet8B{el}; el=el+1; % High Lift Devices  - Leading Edge
LED_type = Sheet8B{el}; %High Lift Devices  - Trailing  Edge (Flaps)	

% Stores the flags
Aerodynamic_Data_flags.Conf.w1 = w1; %
Aerodynamic_Data_flags.Conf.h = h; %
Aerodynamic_Data_flags.Conf.v = v; %
Aerodynamic_Data_flags.Conf.v2 = v2; %
Aerodynamic_Data_flags.Conf.vtail = vtail; %
Aerodynamic_Data_flags.Conf.vtail2 = vtail2; %
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
Aerodynamic_Data_flags.Flap_type = Flap_type; %% High Lift Devices  - Leading Edge
Aerodynamic_Data_flags.LED_type = LED_type; %%High Lift Devices  - Trailing  Edge (Flaps)	

% Gathers all the flags
OUTPUT_read_XLSX.Aerodynamic_Data_flags = Aerodynamic_Data_flags;

%% Stability Data
Sheet9 = readcell(excel_file,'Sheet','9Stability','Range','D3:D36');
el=1;

SM = Sheet9{el}; el=el+1; % Static Margin
prop_wash_effect = Sheet9{el}; el=el+1; % Prop Wash Effect Considered on the dynamic pressure calculation
XCG_FF = Sheet9{el}; el=el+1; % elects XCG depending if it's from desired stability conditions in Forward Flight
StabilityModel = Sheet9{el}; el=el+1; % elects XCG depending if it's from desired stability conditions in Forward Flight
only_trim = Sheet9{el}; el=el+1; % elects XCG depending if it's from desired stability conditions in Forward Flight
V_low = Sheet9{el}; el=el+1; %	Trim Conditions Velocity Range - V low
V_high = Sheet9{el}; el=el+1; %	Trim Conditions Velocity Range - V high
N_V_VAR = Sheet9{el}; el=el+1; %	Number of iteration variation Speed
N_m_VAR = Sheet9{el}; el=el+1; %	Number of iteration variation Weight
beta = Sheet9{el}; el=el+1; %	Trim Lateral Stability Conditions -side slip angle (beta)
beta_i = Sheet9{el}; el=el+1; %	Trim Lateral Stability Conditions -Variable side slip angle (initial beta)
beta_f = Sheet9{el}; el=el+1; %	Trim Lateral Stability Conditions -Variable side slip angle (final beta)
N_Delta_beta = Sheet9{el}; el=el+1; %	Trim Lateral Stability Conditions -Variable side slip angle (Number of beta's)
phi = Sheet9{el}; el=el+1; %	Turning Condition Stability Conditions -bank angle (beta)
phi_i = Sheet9{el}; el=el+1; %	Turning Condition Stability Conditions -Variable bank angle (initial phi)
phi_f = Sheet9{el}; el=el+1; %	Turning Condition Stability Conditions -Variable bank angle (final phi)
N_Delta_phi = Sheet9{el}; el=el+1; %	Turning Condition Stability Conditions -Variable bank angle (Number of beta's)
n_viraje = Sheet9{el}; el=el+1; %	Loading factur during turning flight 
Munk_fuselage_constribution = Sheet9{el}; el=el+1; %	Include Munk's fuselage contribution todetermine Xac, & desired Xcg location for an SM  
flagwingspan2bodydiam = Sheet9{el}; el=el+1; %	Takes into account thew wing span to body diameter interference 
tf1_long = Sheet9{el}; el=el+1; %	Time period ploting Phugoid Mode
tf2_long = Sheet9{el}; el=el+1; %	Time period ploting Short Period Mode
tf1_lat = Sheet9{el}; el=el+1; %	Time period ploting Dutch Roll
Du = Sheet9{el}; el=el+1; %	Perturbation in forward speed velocity (percentage of trim velocity)
Dalpha = Sheet9{el}; el=el+1; %	Perturbation in angle of attack
Dq = Sheet9{el}; el=el+1; %	Perturbation in pitch rate
Dtheta = Sheet9{el}; el=el+1; %	Perturbation in pitch angle
Dbeta = Sheet9{el}; el=el+1; %	Perturbation in side slip angle
Dp = Sheet9{el}; el=el+1; %	Perturbation in roll rate
Dr = Sheet9{el}; el=el+1; %	Perturbation in yaw rate
Dphi = Sheet9{el}; el=el+1; % Perturbation in bank angle
n_min = Sheet9{el}; el=el+1; % Minimum Load Factor
n_max = Sheet9{el}; el=el+1; % Maximum Load Factor
N_n_VAR = Sheet9{el}; el=el+1; %	Number of iterations for variation of load factor

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
Stability_flags.flagwingspan2bodydiam = flagwingspan2bodydiam;
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
el=1;

CASE = Sheet10{el}; el=el+1; % Model of Aircraft
Control_surface = Sheet10{el}; el=el+1; % Determines the ammount of control surface in the w2 aileron & elevator
AC_type = Sheet10{el}; el=el+1; % Aircraft type

% Stores the flags
Geometry_Data_flags.CASE = CASE; %
Geometry_Data_flags.Control_surface = Control_surface; %
Geometry_Data_flags.AC_type = AC_type; %

% Gathers all the flags
OUTPUT_read_XLSX.Geometry_Data_flags = Geometry_Data_flags;

%% Input Aircraft Geometry Data
Sheet11 = readcell(excel_file,'Sheet','11InputGeometry','Range','D3:D209');
el=1;
x_offset_CAD = Sheet11{el}; el=el+1; % Distance from origin in CATIa to Fuselage (offset value) - x
z_offset_CAD = Sheet11{el}; el=el+1; % Distance from origin in CATIa to Fuselage (offset value) - z
w_fus = Sheet11{el}; el=el+1; % Fuselage geomtry (Units in m) - approx - width
h_fus = Sheet11{el}; el=el+1; % Fuselage geomtry (Units in m) - approx - heigth
l_fus = Sheet11{el}; el=el+1; % Fuselage geomtry (Units in m) - approx - length
d_fus = Sheet11{el}; el=el+1; % Fuselage geomtry (Units in m) - approx - diameter
%w1
y_loc_1R_y1_w1_CAD = Sheet11{el}; el=el+1; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
% Kink wings
y_loc_1R_yB1_w1_CAD = Sheet11{el}; el=el+1; % y loc of wing (w1) kink1 chord LE position (distance from CAD refference point)	
y_loc_1R_yB2_w1_CAD = Sheet11{el}; el=el+1; % y loc of wing (w1) kink2 chord LE position (distance from CAD refference point)	
y_loc_1R_y2_w1_CAD = Sheet11{el}; el=el+1; %	y loc of wing (w1) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_w1_CAD = Sheet11{el}; el=el+1; %	z loc of wing (w1) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_w1_CAD = Sheet11{el}; el=el+1; %	x loc of wing (w1) root chord LE position (distance from CAD refference point)
Lambda_LE_w1_e = Sheet11{el}; el=el+1; %	Sweep of w1 (deg)
Lambda_LE_w1_k1_e  = Sheet11{el}; el=el+1; % Sweep of kink1 (deg)	
Lambda_LE_w1_k2_e  = Sheet11{el}; el=el+1; % Sweep of kink2 (deg)	
dihedral_w1_e = Sheet11{el}; el=el+1; %	Dihedral of w1 (deg)
dihedral_w1_k1_e = Sheet11{el}; el=el+1; % Dihedral of kink1 (deg)	
dihedral_w1_k2_e = Sheet11{el}; el=el+1; % Dihedral of kink2 (deg)	
%w2
y_loc_1R_y1_w2_CAD = Sheet11{el}; el=el+1; %	y loc of wing (w2) root chord LE position (distance from CAD refference point)
% Kink surface
y_loc_1R_yB1_w2_CAD = Sheet11{el}; el=el+1; % y loc of wing (w2) kink1 chord LE position (distance from CAD refference point)	
y_loc_1R_yB2_w2_CAD = Sheet11{el}; el=el+1; % y loc of wing (w2) kink2 chord LE position (distance from CAD refference point)	
y_loc_1R_y2_w2_CAD = Sheet11{el}; el=el+1; %	y loc of wing (w2) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_w2_CAD = Sheet11{el}; el=el+1; %	z loc of wing (w2) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_w2_CAD = Sheet11{el}; el=el+1; %	x loc of wing (w2) root chord LE position (distance from CAD refference point)
Lambda_LE_w2_e = Sheet11{el}; el=el+1; %	Sweep of w2 (deg)
Lambda_LE_w2_k1_e  = Sheet11{el}; el=el+1; % Sweep of kink1 (deg)	
Lambda_LE_w2_k2_e  = Sheet11{el}; el=el+1; % Sweep of kink2 (deg)	
dihedral_w2_e = Sheet11{el}; el=el+1; %	Dihedral of w2 (deg)
dihedral_w2_k1_e = Sheet11{el}; el=el+1; % Dihedral of kink1 (deg)	
dihedral_w2_k2_e = Sheet11{el}; el=el+1; % Dihedral of kink2 (deg)	
%can
y_loc_1R_y1_can_CAD = Sheet11{el}; el=el+1; %	y loc of wing (can) root chord LE position (distance from CAD refference point)
% Kink surface
y_loc_1R_yB1_can_CAD = Sheet11{el}; el=el+1; % y loc of can kink1 chord LE position (distance from CAD refference point)	
y_loc_1R_yB2_can_CAD = Sheet11{el}; el=el+1; % y loc of can kink2 chord LE position (distance from CAD refference point)	
y_loc_1R_y2_can_CAD = Sheet11{el}; el=el+1; %	y loc of wing (can) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_can_CAD = Sheet11{el}; el=el+1; %	z loc of wing (can) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_can_CAD = Sheet11{el}; el=el+1; %	x loc of wing (can) root chord LE position (distance from CAD refference point)
Lambda_LE_can_e = Sheet11{el}; el=el+1; %	Sweep of can (deg)
Lambda_LE_can_k1_e  = Sheet11{el}; el=el+1; % Sweep of kink1 (deg)	
Lambda_LE_can_k2_e  = Sheet11{el}; el=el+1; % Sweep of kink2 (deg)	
dihedral_can_e = Sheet11{el}; el=el+1; %	Dihedral of can (deg)
dihedral_can_k1_e = Sheet11{el}; el=el+1; % Dihedral of kink1 (deg)	
dihedral_can_k2_e = Sheet11{el}; el=el+1; % Dihedral of kink2 (deg)	
%HTP
y_loc_1R_y1_HTP_CAD = Sheet11{el}; el=el+1; %	y loc of wing (HTP) root chord LE position (distance from CAD refference point)
% Kink surface
y_loc_1R_yB1_HTP_CAD = Sheet11{el}; el=el+1; % y loc of HTP kink1 chord LE position (distance from CAD refference point)	
y_loc_1R_yB2_HTP_CAD = Sheet11{el}; el=el+1; % y loc of HTP kink2 chord LE position (distance from CAD refference point)	
y_loc_1R_y2_HTP_CAD = Sheet11{el}; el=el+1; %	y loc of wing (HTP) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_HTP_CAD = Sheet11{el}; el=el+1; %	z loc of wing (HTP) root chord LE position (distance from CAD refference point)
 x_loc_1R_y1_HTP_CAD = Sheet11{el}; el=el+1; %	x loc of wing (HTP) root chord LE position (distance from CAD refference point)
Lambda_LE_HTP_e = Sheet11{el}; el=el+1; % Sweep of HTP (deg)
Lambda_LE_HTP_k1_e  = Sheet11{el}; el=el+1; % Sweep of kink1 (deg)	
Lambda_LE_HTP_k2_e  = Sheet11{el}; el=el+1; % Sweep of kink2 (deg)	
dihedral_HTP_e = Sheet11{el}; el=el+1; %	Dihedral of HTP (deg)
dihedral_HTP_k1_e = Sheet11{el}; el=el+1; % Dihedral of kink1 (deg)	
dihedral_HTP_k2_e = Sheet11{el}; el=el+1; % Dihedral of kink2 (deg)	
%VTP
y_loc_1R_y1_VTP_CAD = Sheet11{el}; el=el+1; %	y loc of wing (VTP) root chord LE position (distance from CAD refference point)
z_loc_1R_y1_VTP_CAD = Sheet11{el}; el=el+1; %	y loc of wing (VTP) tip chord LE position (distance from CAD refference point)
z_loc_1R_y2_VTP_CAD = Sheet11{el}; el=el+1; %	z loc of wing (VTP) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_VTP_CAD = Sheet11{el}; el=el+1; %	x loc of wing (VTP) root chord LE position (distance from CAD refference point)
Lambda_LE_VTP_e = Sheet11{el}; el=el+1; %	Sweep of VTP (deg)
dihedral_VTP_e = Sheet11{el}; el=el+1; %	Dihedral of VTP (deg)
%VTAIL
y_loc_1R_y1_vee_CAD = Sheet11{el}; el=el+1; %	y loc of wing (V-tail) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_vee_CAD = Sheet11{el}; el=el+1; %	y loc of wing (V-tail) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_vee_CAD = Sheet11{el}; el=el+1; %	z loc of wing (V-tail) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_vee_CAD = Sheet11{el}; el=el+1; %	x loc of wing (V-tail) root chord LE position (distance from CAD refference point)
Lambda_LE_vee_e = Sheet11{el}; el=el+1; %	Sweep of V-tail (deg)
dihedral_vee_e = Sheet11{el}; el=el+1; %	Dihedral of V-tail (deg)
%VTAIL2
y_loc_1R_y1_vee2_CAD = Sheet11{el}; el=el+1; %	y loc of wing (V-tail) root chord LE position (distance from CAD refference point)
y_loc_1R_y2_vee2_CAD = Sheet11{el}; el=el+1; %	y loc of wing (V-tail) tip chord LE position (distance from CAD refference point)
z_loc_1R_y1_vee2_CAD = Sheet11{el}; el=el+1; %	z loc of wing (V-tail) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_vee2_CAD = Sheet11{el}; el=el+1; %	x loc of wing (V-tail) root chord LE position (distance from CAD refference point)
Lambda_LE_vee2_e = Sheet11{el}; el=el+1; %	Sweep of V-tail (deg)
dihedral_vee2_e = Sheet11{el}; el=el+1; %	Dihedral of V-tail (deg)
%Twin VTP
y_loc_1R_y1_2VTP_CAD = Sheet11{el}; el=el+1; %	y loc of Twin (VTP) root chord LE position (distance from CAD refference point)
z_loc_1R_y1_2VTP_CAD = Sheet11{el}; el=el+1; %	y loc of Twin (VTP) tip chord LE position (distance from CAD refference point)
z_loc_1R_y2_2VTP_CAD = Sheet11{el}; el=el+1; %	z loc of wing (VTP) root chord LE position (distance from CAD refference point)
x_loc_1R_y1_2VTP_CAD = Sheet11{el}; el=el+1; %	x loc of Twin (VTP) root chord LE position (distance from CAD refference point)
Lambda_LE_2VTP_e = Sheet11{el}; el=el+1; %	Sweep of Twin VTP (deg)
dihedral_2VTP_e = Sheet11{el}; el=el+1; %	Dihedral of Twin VTP (deg)
% Chord geometry
cR_w1 = Sheet11{el}; el=el+1; %	Root Chord w1
cB_k1_w1 = Sheet11{el}; el=el+1; %	Kink1 Chord w1	
cB_k2_w1 = Sheet11{el}; el=el+1; %	Kink2 Chord w1	
cT_w1 = Sheet11{el}; el=el+1; %	Tip Chord w1
cR_w2 = Sheet11{el}; el=el+1; %	Root Chord w2
cB_k1_w2 = Sheet11{el}; el=el+1; %	Kink1 Chord w2	
cB_k2_w2 = Sheet11{el}; el=el+1; %	Kink2 Chord w2	
cT_w2 = Sheet11{el}; el=el+1; %	Tip Chord w2
cR_can = Sheet11{el}; el=el+1; %	Root Chord canard
cB_k1_can = Sheet11{el}; el=el+1; %	Kink1 Chord canard	
cB_k2_can = Sheet11{el}; el=el+1; %	Kink2 Chord canard	
cT_can = Sheet11{el}; el=el+1; %	Tip Chord canard
cR_HTP = Sheet11{el}; el=el+1; %	Root Chord HTP
cB_k1_HTP = Sheet11{el}; el=el+1; %	Kink1 Chord HTP
cB_k2_HTP = Sheet11{el}; el=el+1; %	Kink2 Chord HTP	
cT_HTP = Sheet11{el}; el=el+1; %	Tip Chord HTP
cR_VTP = Sheet11{el}; el=el+1; %	Root Chord VTP
cT_VTP = Sheet11{el}; el=el+1; %	Tip Chord VTP
cR_vee = Sheet11{el}; el=el+1; %	Root Chord vee
cT_vee = Sheet11{el}; el=el+1; %	Tip Chord vee
cR_vee2 = Sheet11{el}; el=el+1; %	Root Chord vee
cT_vee2 = Sheet11{el}; el=el+1; %	Tip Chord vee
cR_2VTP = Sheet11{el}; el=el+1; %	Root Chord Twin VTP
cT_2VTP = Sheet11{el}; el=el+1; %	Tip Chord Twin VTP
% Perfecntage (%) of the control surface
K_y1_ail_w1 = Sheet11{el}; el=el+1; %	% of aileron from effective wing surface - inner position
K_y2_ail_w1 = Sheet11{el}; el=el+1; %	% of aileron from effective wing surface - outter position
K_y1_ele_w2 = Sheet11{el}; el=el+1; %	% of elevator from effective HTP surface - inner position
K_y2_ele_w2 = Sheet11{el}; el=el+1; %	% of elevator from effective HTP surface  - outter position
K_y1_elevon_w1 = Sheet11{el}; el=el+1; %	% of elevon from effective wing surface - inner position
K_y2_elevon_w1 = Sheet11{el}; el=el+1; %	% of elevon from effective wing surface - outter position
K_y1_flap_w1 = Sheet11{el}; el=el+1; %	% of flap from effective wing surface - inner position
K_y2_flap_w1 = Sheet11{el}; el=el+1; %	% of flap from effective wing surface - outter position
K_y1_rudder_VTP = Sheet11{el}; el=el+1; %	% of rudder from effective VTP surface - inner position
K_y2_rudder_VTP = Sheet11{el}; el=el+1; %	% of rudder from effective VTP surface  - outter position
K_y1_rudvtr_w2 = Sheet11{el}; el=el+1; %	% of ruddervator from effective V-tail surface - inner position
K_y2_rudvtr_w2 = Sheet11{el}; el=el+1; %	% of ruddervator from effective V-tail surface  - outter position
K_y1_rudvtr2_w2 = Sheet11{el}; el=el+1; %	% of ruddervator from effective V-tail surface - inner position
K_y2_rudvtr2_w2 = Sheet11{el}; el=el+1; %	% of ruddervator from effective V-tail surface  - outter position
K_y1_canard_can = Sheet11{el}; el=el+1; %	% of canard control surface from effective canard surface - inner position
K_y2_canard_can = Sheet11{el}; el=el+1; %	% of canard control surface from effective canard surface - outter position
% Trim Tabs
K_y1_TT_ail_w1 = Sheet11{el}; el=el+1;% of Trim Tab aileron from effective wing surface - inner position	
K_y2_TT_ail_w1 = Sheet11{el}; el=el+1;% of Trim Tab  aileron from effective wing surface - outter position	
K_y1_TT_ele_w2 = Sheet11{el}; el=el+1;% of Trim Tab  elevator from effective HTP surface - inner position	
K_y2_TT_ele_w2 = Sheet11{el}; el=el+1;% of Trim Tab  elevator from effective HTP surface  - outter position	
K_y1_TT_rudder_VTP = Sheet11{el}; el=el+1;% of Trim Tab rudder from effective VTP surface - inner position	
K_y2_TT_rudder_VTP = Sheet11{el}; el=el+1;% of Trim Tab rudder from effective VTP surface  - outter position	
K_y1_TT_sp_w1 = Sheet11{el}; el=el+1;% of Trim Tab spoiler from effective wing surface - inner position	
K_y2_TT_sp_w1 = Sheet11{el}; el=el+1;% of Trim Tab spoiler from effective wing surface  - outter position	
% Maximum and minimum control surface deflections
delta_ail_min = Sheet11{el}; el=el+1; % 	Minimum aileron deflection
delta_ail_max = Sheet11{el}; el=el+1; %	Maximum aileron deflection
delta_ele_min = Sheet11{el}; el=el+1; %	Minimum elevator deflection
delta_ele_max = Sheet11{el}; el=el+1; %	Maximum elevator deflection
delta_elevon_min = Sheet11{el}; el=el+1; %	Minimum elevon deflection
delta_elevon_max = Sheet11{el}; el=el+1; %	Maximum elevon deflection
delta_flap_min = Sheet11{el}; el=el+1; %	Minimum flap deflection
delta_flap_max = Sheet11{el}; el=el+1; %	Maximum flap deflection
delta_rudder_min = Sheet11{el}; el=el+1; %	Minimum rudder deflection
delta_rudder_max = Sheet11{el}; el=el+1; %	Maximum rudder deflection
delta_rudvtr_min = Sheet11{el}; el=el+1; %	Minimum ruddervator deflection
delta_rudvtr_max = Sheet11{el}; el=el+1; %	Maximum ruddervator deflection
delta_can_min = Sheet11{el}; el=el+1; %	Minimum canard deflection
delta_can_max = Sheet11{el}; el=el+1; %	Maximum canard deflection
% Max defelction Trim Tabs
delta_TT_ail_min = Sheet11{el}; el=el+1; % Minimum Trim Tab aileron deflection
delta_TT_ail_max = Sheet11{el}; el=el+1; %	Maximum Trim Tab aileron deflection
delta_TT_ele_min = Sheet11{el}; el=el+1; %	Minimum Trim Tab elevator deflection
delta_TT_ele_max = Sheet11{el}; el=el+1; %	Maximum Trim Tab elevator deflection
delta_TT_rudder_min = Sheet11{el}; el=el+1; %	Minimum Trim Tab rudder deflection
delta_TT_rudder_max = Sheet11{el}; el=el+1; %	Maximum Trim Tab rudder deflection
delta_TT_sp_min = Sheet11{el}; el=el+1; %	Minimum Trim Tab spolier deflection
delta_TT_sp_max = Sheet11{el}; el=el+1; %	Maximum Trim Tab spoiler deflection

% Control surface dimensions
cf_ail = Sheet11{el}; el=el+1; %	%  of control surface aileron (chrodwise)
t_c_ail = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - aileron
cf_ele = Sheet11{el}; el=el+1; %	%  of control surface elevator (chrodwise)
t_c_ele = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - elevator
cf_elevon = Sheet11{el}; el=el+1; %	%  of control surface elevon (chrodwise)
t_c_elevon = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - elevon
cf_flap = Sheet11{el}; el=el+1; %	%  of control surface flap (chrodwise)
t_c_flap = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - flap
cf_rudder = Sheet11{el}; el=el+1; %	%  of control surface rudder (chrodwise)
t_c_rudder = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - rudder
cf_rudvtr = Sheet11{el}; el=el+1; %	%  of control surface ruddervator (chrodwise)
t_c_rudvtr = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - ruddervator
cf_rudvtr2 = Sheet11{el}; el=el+1; %	%  of control surface ruddervator (chrodwise)
t_c_rudvtr2 = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - ruddervator
cf_canard = Sheet11{el}; el=el+1; %	%  of control surface canard (chrodwise)
t_c_canard = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - canard
% Geometry of Trim Tabs control surfaces
cf_TT_ail = Sheet11{el}; el=el+1; %	%  of control surface Trim Tab  aileron (chrodwise)
t_c_TT_ail = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - Trim Tab  aileron
cf_TT_ele = Sheet11{el}; el=el+1; %	%  of control surface Trim Tab  elevator (chrodwise)
t_c_TT_ele = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - Trim Tab  elevator
cf_TT_rudder = Sheet11{el}; el=el+1; %	%  of control surface Trim Tab  rudder (chrodwise)
t_c_TT_rudder = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - Trim Tab  rudder
cf_TT_sp = Sheet11{el}; el=el+1; %	%  of control surface Trim Tab  spoiler (chrodwise)
t_c_TT_sp = Sheet11{el}; el=el+1; %	Thinckness 2 chord ratio associated to the airfoil - Trim Tab  spoiler
% Initial center of gravity estimaton
x_XCG = Sheet11{el}; el=el+1;
y_XCG = Sheet11{el}; el=el+1;
z_XCG = Sheet11{el}; el=el+1;
% Location of Engine type 1 (symetry applied)
x_eng_ybar1 = Sheet11{el}; el=el+1;
y_eng_ybar1 = Sheet11{el}; el=el+1;
z_eng_ybar1 = Sheet11{el}; el=el+1;
% Location of Engine type 2 (symetry applied)
x_eng_ybar2 = Sheet11{el}; el=el+1;
y_eng_ybar2 = Sheet11{el}; el=el+1;
z_eng_ybar2 = Sheet11{el}; el=el+1;
% % Location of Engine type 1 (symetry applied)
% x_eng_ybar3 = Sheet11{el}; el=el+1;
% y_eng_ybar3 = Sheet11{el}; el=el+1;
% z_eng_ybar3 = Sheet11{el}; el=el+1;
% % Location of Engine type 2 (symetry applied)
% x_eng_ybar4 = Sheet11{el}; el=el+1;
% y_eng_ybar4 = Sheet11{el}; el=el+1;
% z_eng_ybar4 = Sheet11{el}; el=el+1;
% Location of the center of Leading Edge Nacelle of propulsion system 1
x_nac_ybar1 = Sheet11{el}; el=el+1; % Location of the center of LE (z) Nacelle of propulsion system 1	x_nac_ybar1
y_nac_ybar1 = Sheet11{el}; el=el+1; % Location of the center of LE (y) Nacelle of propulsion system 1	y_nac_ybar1
z_nac_ybar1 = Sheet11{el}; el=el+1; % Location of the center of LE (z) Nacelle of propulsion system 1	z_nac_ybar1
% Location of the center of Leading Edge Nacelle of propulsion system 2
x_nac_ybar2 = Sheet11{el}; el=el+1; % Location of the center of LE (z) Nacelle of propulsion system 2	x_nac_ybar1
y_nac_ybar2 = Sheet11{el}; el=el+1; % Location of the center of LE (y) Nacelle of propulsion system 2	y_nac_ybar1
z_nac_ybar2 = Sheet11{el}; el=el+1; % Location of the center of LE (z) Nacelle of propulsion system 2	z_nac_ybar1
% % Location of the center of Leading Edge Nacelle of propulsion system 1
% x_nac_ybar3 = Sheet11{el}; el=el+1; % Location of the center of LE (z) Nacelle of propulsion system 1	x_nac_ybar1
% y_nac_ybar3 = Sheet11{el}; el=el+1; % Location of the center of LE (y) Nacelle of propulsion system 1	y_nac_ybar1
% z_nac_ybar3 = Sheet11{el}; el=el+1; % Location of the center of LE (z) Nacelle of propulsion system 1	z_nac_ybar1
% % Location of the center of Leading Edge Nacelle of propulsion system 2
% x_nac_ybar4 = Sheet11{el}; el=el+1; % Location of the center of LE (z) Nacelle of propulsion system 2	x_nac_ybar1
% y_nac_ybar4 = Sheet11{el}; el=el+1; % Location of the center of LE (y) Nacelle of propulsion system 2	y_nac_ybar1
% z_nac_ybar4 = Sheet11{el}; el=el+1; % Location of the center of LE (z) Nacelle of propulsion system 2	z_nac_ybar1
% Wing Offset: Determines if the wing geometrý includes center section offset	
wing_offset_w1 = Sheet11{el}; el=el+1; % Wing Offset: Determines if the wing geometrý includes center section offset
canard_offset_can = Sheet11{el}; el=el+1;% Canard Offset: Determines if the canard geometrý includes center section offset
zoffset_fuselage = Sheet11{el}; el=el+1;% Fuselage Offset: Determines if the canard geometrý includes center section offset
l_tailboom = Sheet11{el}; el=el+1;% Length of Tailboom
w_tailboom = Sheet11{el}; el=el+1;% Width of Tailboom
h_tailboom = Sheet11{el}; el=el+1;% Height of Tailboom
l_pod = Sheet11{el}; el=el+1;% Length of nacelle
w_pod = Sheet11{el}; el=el+1;% Width of nacelle
h_pod = Sheet11{el}; el=el+1;% Height of nacelle
% Available control surfaces
d_ail = Sheet11{el}; el=el+1; %Definition of available control surface - aileron
d_ele = Sheet11{el}; el=el+1; %Definition of available control surface - elevator
d_elevon = Sheet11{el}; el=el+1; %Definition of available control surface - elevon
d_flap = Sheet11{el}; el=el+1; %Definition of available control surface - flap
d_rudder = Sheet11{el}; el=el+1; %Definition of available control surface - rudder
d_rudvtr = Sheet11{el}; el=el+1; %Definition of available control surface - ruddervator
d_can = Sheet11{el}; el=el+1; %Definition of available control surface - canard
% Available Trim Tabs
d_TT_ail = Sheet11{el}; el=el+1; %Definition of available control surface - Trim Tab aileron
d_TT_ele = Sheet11{el}; el=el+1; %Definition of available control surface - Trim Tab elevator
d_TT_rudder = Sheet11{el}; el=el+1; %Definition of available control surface - Trim Tab rudder
d_TT_sp = Sheet11{el}; el=el+1; %Definition of available control surface - Trim Tab Spolier

% Stores the flags
InputGeometry_Data_flags.x_offset_CAD = x_offset_CAD; %
InputGeometry_Data_flags.z_offset_CAD = z_offset_CAD; %
InputGeometry_Data_flags.w_fus = w_fus; %
InputGeometry_Data_flags.h_fus = h_fus; %
InputGeometry_Data_flags.l_fus = l_fus; %
InputGeometry_Data_flags.d_fus = d_fus; %
InputGeometry_Data_flags.d_fus = d_fus; %
% wing
InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD = y_loc_1R_y1_w1_CAD; %
% Kink wings
InputGeometry_Data_flags.y_loc_1R_yB1_w1_CAD = y_loc_1R_yB1_w1_CAD; %
InputGeometry_Data_flags.y_loc_1R_yB2_w1_CAD = y_loc_1R_yB2_w1_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_w1_CAD = z_loc_1R_y1_w1_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; %
InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; %
InputGeometry_Data_flags.Lambda_LE_w1_k1_e = Lambda_LE_w1_k1_e; %
InputGeometry_Data_flags.Lambda_LE_w1_k2_e = Lambda_LE_w1_k2_e; %
InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; %
InputGeometry_Data_flags.dihedral_w1_k1_e = dihedral_w1_k1_e; %
InputGeometry_Data_flags.dihedral_w1_k2_e = dihedral_w1_k2_e; %
% wing 2
InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = y_loc_1R_y1_w2_CAD; %
% Kink wings
InputGeometry_Data_flags.y_loc_1R_yB1_w2_CAD = y_loc_1R_yB1_w2_CAD; %
InputGeometry_Data_flags.y_loc_1R_yB2_w2_CAD = y_loc_1R_yB2_w2_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = y_loc_1R_y2_w2_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_w2_CAD = z_loc_1R_y1_w2_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD = x_loc_1R_y1_w2_CAD; %
InputGeometry_Data_flags.Lambda_LE_w2_e = Lambda_LE_w2_e; %
InputGeometry_Data_flags.Lambda_LE_w2_k1_e = Lambda_LE_w2_k1_e; %
InputGeometry_Data_flags.Lambda_LE_w2_k2_e = Lambda_LE_w2_k2_e; %
InputGeometry_Data_flags.dihedral_w2_e = dihedral_w2_e; %
InputGeometry_Data_flags.dihedral_w2_k1_e = dihedral_w2_k1_e; %
InputGeometry_Data_flags.dihedral_w2_k2_e = dihedral_w2_k2_e; %

% Can
InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = y_loc_1R_y1_can_CAD; %
% Kink wings
InputGeometry_Data_flags.y_loc_1R_yB1_can_CAD = y_loc_1R_yB1_can_CAD; %
InputGeometry_Data_flags.y_loc_1R_yB2_can_CAD = y_loc_1R_yB2_can_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_can_CAD = y_loc_1R_y2_can_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_can_CAD = z_loc_1R_y1_can_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_can_CAD = x_loc_1R_y1_can_CAD; %
InputGeometry_Data_flags.Lambda_LE_can_e = Lambda_LE_can_e; %
InputGeometry_Data_flags.Lambda_LE_can_k1_e = Lambda_LE_can_k1_e; %
InputGeometry_Data_flags.Lambda_LE_can_k2_e = Lambda_LE_can_k2_e; %
InputGeometry_Data_flags.dihedral_can_e = dihedral_can_e; %
InputGeometry_Data_flags.dihedral_can_k1_e = dihedral_can_k1_e; %
InputGeometry_Data_flags.dihedral_can_k2_e = dihedral_can_k2_e; %

% HTP
InputGeometry_Data_flags.y_loc_1R_y1_HTP_CAD = y_loc_1R_y1_HTP_CAD; %
% Kink wings
InputGeometry_Data_flags.y_loc_1R_yB1_HTP_CAD = y_loc_1R_yB1_HTP_CAD; %
InputGeometry_Data_flags.y_loc_1R_yB2_HTP_CAD = y_loc_1R_yB2_HTP_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_HTP_CAD = y_loc_1R_y2_HTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_HTP_CAD = z_loc_1R_y1_HTP_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_HTP_CAD = x_loc_1R_y1_HTP_CAD; %
InputGeometry_Data_flags.Lambda_LE_HTP_e = Lambda_LE_HTP_e; %
InputGeometry_Data_flags.Lambda_LE_HTP_k1_e = Lambda_LE_HTP_k1_e; %
InputGeometry_Data_flags.Lambda_LE_HTP_k2_e = Lambda_LE_HTP_k2_e; %
InputGeometry_Data_flags.dihedral_HTP_e = dihedral_HTP_e; %
InputGeometry_Data_flags.dihedral_HTP_k1_e = dihedral_HTP_k1_e; %
InputGeometry_Data_flags.dihedral_HTP_k2_e = dihedral_HTP_k2_e; %

% VTP
InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; %
InputGeometry_Data_flags.Lambda_LE_VTP_e = Lambda_LE_VTP_e; %
InputGeometry_Data_flags.dihedral_VTP_e = dihedral_VTP_e; %
% Vee Tail 
InputGeometry_Data_flags.y_loc_1R_y1_vee_CAD = y_loc_1R_y1_vee_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_vee_CAD = y_loc_1R_y2_vee_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_vee_CAD = z_loc_1R_y1_vee_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_vee_CAD = x_loc_1R_y1_vee_CAD; %
InputGeometry_Data_flags.Lambda_LE_vee_e = Lambda_LE_vee_e; %
InputGeometry_Data_flags.dihedral_vee_e = dihedral_vee_e; %
% Vee Tail 2 
InputGeometry_Data_flags.y_loc_1R_y1_vee2_CAD = y_loc_1R_y1_vee2_CAD; %
InputGeometry_Data_flags.y_loc_1R_y2_vee2_CAD = y_loc_1R_y2_vee2_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_vee2_CAD = z_loc_1R_y1_vee2_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_vee2_CAD = x_loc_1R_y1_vee2_CAD; %
InputGeometry_Data_flags.Lambda_LE_vee2_e = Lambda_LE_vee2_e; %
InputGeometry_Data_flags.dihedral_vee2_e = dihedral_vee2_e; %

% Second VTP
InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_2VTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_2VTP_CAD; %
InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_2VTP_CAD; %
InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_2VTP_CAD; %
InputGeometry_Data_flags.Lambda_LE_2VTP_e = Lambda_LE_2VTP_e; %
InputGeometry_Data_flags.dihedral_2VTP_e = dihedral_2VTP_e; %
% Chord geometry
%wing 1
InputGeometry_Data_flags.cR_w1 = cR_w1; %
InputGeometry_Data_flags.cB_k1_w1 = cB_k1_w1; %
InputGeometry_Data_flags.cB_k2_w1 = cB_k2_w1; %
InputGeometry_Data_flags.cT_w1 = cT_w1; %
% Wing2
InputGeometry_Data_flags.cR_w2 = cR_w2; %
InputGeometry_Data_flags.cB_k1_w2 = cB_k1_w2; %
InputGeometry_Data_flags.cB_k2_w2 = cB_k2_w2; %
InputGeometry_Data_flags.cT_w2 = cT_w2; %
% canard
InputGeometry_Data_flags.cR_can = cR_can; %
InputGeometry_Data_flags.cB_k1_can = cB_k1_can; %
InputGeometry_Data_flags.cB_k2_can = cB_k2_can; %
InputGeometry_Data_flags.cT_can = cT_can; %
% HTP
InputGeometry_Data_flags.cR_HTP = cR_HTP; %
InputGeometry_Data_flags.cB_k1_HTP = cB_k1_HTP; %
InputGeometry_Data_flags.cB_k2_HTP = cB_k2_HTP; %
InputGeometry_Data_flags.cT_HTP = cT_HTP; %
% VTP
InputGeometry_Data_flags.cR_VTP = cR_VTP; %
InputGeometry_Data_flags.cT_VTP = cT_VTP; %
%Vee Tail
InputGeometry_Data_flags.cR_vee = cR_vee; %
InputGeometry_Data_flags.cT_vee = cT_vee; %
%Vee Tail2
InputGeometry_Data_flags.cR_vee2 = cR_vee2; %
InputGeometry_Data_flags.cT_vee2 = cT_vee2; %
% Twin VTP
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
InputGeometry_Data_flags.K_y1_rudvtr2_w2 = K_y1_rudvtr2_w2;
InputGeometry_Data_flags.K_y2_rudvtr2_w2 = K_y2_rudvtr2_w2;
InputGeometry_Data_flags.K_y1_canard_can = K_y1_canard_can;
InputGeometry_Data_flags.K_y2_canard_can = K_y2_canard_can;
% Trim Tabs
InputGeometry_Data_flags.K_y1_TT_ail_w1 = K_y1_TT_ail_w1;
InputGeometry_Data_flags.K_y2_TT_ail_w1 = K_y2_TT_ail_w1;
InputGeometry_Data_flags.K_y1_TT_ele_w2 = K_y1_TT_ele_w2;
InputGeometry_Data_flags.K_y2_TT_ele_w2 = K_y2_TT_ele_w2;
InputGeometry_Data_flags.K_y1_TT_rudder_VTP = K_y1_TT_rudder_VTP;
InputGeometry_Data_flags.K_y2_TT_rudder_VTP = K_y2_TT_rudder_VTP;
InputGeometry_Data_flags.K_y1_TT_sp_w1 = K_y1_TT_sp_w1;
InputGeometry_Data_flags.K_y2_TT_sp_w1 = K_y2_TT_sp_w1;
% Max Deflections
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
% TRaim Tabs Max and min deflections
InputGeometry_Data_flags.delta_TT_ail_min = delta_TT_ail_min; 
InputGeometry_Data_flags.delta_TT_ail_max = delta_TT_ail_max;
InputGeometry_Data_flags.delta_TT_ele_min = delta_TT_ele_min;
InputGeometry_Data_flags.delta_TT_ele_max = delta_TT_ele_max;
InputGeometry_Data_flags.delta_TT_rudder_min = delta_TT_rudder_min;
InputGeometry_Data_flags.delta_TT_rudder_max = delta_TT_rudder_max;
InputGeometry_Data_flags.delta_TT_sp_min = delta_TT_sp_min;
InputGeometry_Data_flags.delta_TT_sp_max = delta_TT_sp_max;
% Geometry of control surfaces
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
InputGeometry_Data_flags.cf_rudvtr2 = cf_rudvtr2;
InputGeometry_Data_flags.t_c_rudvtr2 = t_c_rudvtr2;
InputGeometry_Data_flags.cf_canard = cf_canard;
InputGeometry_Data_flags.t_c_canard = t_c_canard;
% Trim Tabs
InputGeometry_Data_flags.cf_TT_ail = cf_TT_ail;
InputGeometry_Data_flags.t_c_TT_ail = t_c_TT_ail;
InputGeometry_Data_flags.cf_TT_ele = cf_TT_ele;
InputGeometry_Data_flags.t_c_TT_ele = t_c_TT_ele;
InputGeometry_Data_flags.cf_TT_rudder = cf_TT_rudder;
InputGeometry_Data_flags.t_c_TT_rudder = t_c_TT_rudder;
InputGeometry_Data_flags.cf_TT_sp = cf_TT_sp;
InputGeometry_Data_flags.t_c_TT_sp = t_c_TT_sp;
InputGeometry_Data_flags.x_XCG = x_XCG;
InputGeometry_Data_flags.y_XCG = y_XCG;
InputGeometry_Data_flags.z_XCG = z_XCG;
% Location of Engine type 1 (symetry applied)
InputGeometry_Data_flags.x_eng_ybar1 = x_eng_ybar1;
InputGeometry_Data_flags.y_eng_ybar1 = y_eng_ybar1;
InputGeometry_Data_flags.z_eng_ybar1 = z_eng_ybar1;
% Location of Engine type 2 (symetry applied)
InputGeometry_Data_flags.x_eng_ybar2 = x_eng_ybar2;
InputGeometry_Data_flags.y_eng_ybar2 = y_eng_ybar2;
InputGeometry_Data_flags.z_eng_ybar2 = z_eng_ybar2;
% % Location of Engine type 3 (symetry applied)
% InputGeometry_Data_flags.x_eng_ybar3 = x_eng_ybar3;
% InputGeometry_Data_flags.y_eng_ybar3 = y_eng_ybar3;
% InputGeometry_Data_flags.z_eng_ybar3 = z_eng_ybar3;
% % Location of Engine type 4 (symetry applied)
% InputGeometry_Data_flags.x_eng_ybar4 = x_eng_ybar4;
% InputGeometry_Data_flags.y_eng_ybar4 = y_eng_ybar4;
% InputGeometry_Data_flags.z_eng_ybar4 = z_eng_ybar4;
% Location of the center of Leading Edge Nacelle of propulsion system 1
InputGeometry_Data_flags.x_nac_ybar1 = x_nac_ybar1;
InputGeometry_Data_flags.y_nac_ybar1 = y_nac_ybar1;
InputGeometry_Data_flags.z_nac_ybar1 = z_nac_ybar1;
% Location of the center of Leading Edge Nacelle of propulsion system 2
InputGeometry_Data_flags.x_nac_ybar2 = x_nac_ybar2;
InputGeometry_Data_flags.y_nac_ybar2 = y_nac_ybar2;
InputGeometry_Data_flags.z_nac_ybar2 = z_nac_ybar2;
% % Location of the center of Leading Edge Nacelle of propulsion system 3
% InputGeometry_Data_flags.x_nac_ybar3 = x_nac_ybar3;
% InputGeometry_Data_flags.y_nac_ybar3 = y_nac_ybar3;
% InputGeometry_Data_flags.z_nac_ybar3 = z_nac_ybar3;
% % Location of the center of Leading Edge Nacelle of propulsion system 4
% InputGeometry_Data_flags.x_nac_ybar4 = x_nac_ybar4;
% InputGeometry_Data_flags.y_nac_ybar4 = y_nac_ybar4;
% InputGeometry_Data_flags.z_nac_ybar4 = z_nac_ybar4;
% Flags
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
% Trim Tabs
% Available control surfaces
InputGeometry_Data_flags.d_TT_ail = d_TT_ail; %Definition of available control surface - aileron
InputGeometry_Data_flags.d_TT_ele = d_TT_ele; %Definition of available control surface - elevator
InputGeometry_Data_flags.d_TT_rudder = d_TT_rudder; %Definition of available control surface - rudder
InputGeometry_Data_flags.d_TT_sp = d_TT_sp; %Definition of available control surface - rudder



% Gathers all the flags
OUTPUT_read_XLSX.InputGeometry_Data_flags = InputGeometry_Data_flags;

%% PLOTS
Sheet12 = readcell(excel_file,'Sheet','12Plots','Range','D3:D29');
el=1;

plot(el) = Sheet12{el}; el=el+1;% Prints plots - Aero
plot(el) = Sheet12{el}; el=el+1;% Prints plots - Aero Polar
plot(el) = Sheet12{el}; el=el+1;% Print Plots for Stimation ofXAC
plot(el) = Sheet12{el}; el=el+1;% Print Plots for Propulsive Models
plot(el) = Sheet12{el}; el=el+1;% Prints plots for 3D
plot(el) = Sheet12{el}; el=el+1;% Prints plots of SM analysis
plot(el) = Sheet12{el}; el=el+1;% Prints plots of longitudinal Trim
plot(el) = Sheet12{el}; el=el+1;% Prints plots of longitudinal Trim with Variable mass & Variable Speed
plot(el) = Sheet12{el}; el=el+1;% Prints plots of lateral Trim
plot(el) = Sheet12{el}; el=el+1;% Prints plots of lateral Trim - Asymmetries
plot(el) = Sheet12{el}; el=el+1;% Prints plots of lateral Trim - Accelerations
plot(el) = Sheet12{el}; el=el+1;% Prints plots of lateral Trim - Trim Tabs
plot(el) = Sheet12{el}; el=el+1;% Prints plots of lateral Turning
plot(el) = Sheet12{el}; el=el+1;% Prints PLOTS STABILITY DERIVATIVES FOR VAR MASS AND VELOCITY
plot(el) = Sheet12{el}; el=el+1;% Prints PLOTS STABILITY ANALYSIS FOR VAR MASS AND VELOCITY
plot(el) = Sheet12{el}; el=el+1;% Prints PLOTS DYNAMIC STABILITY ANALYSIS AFTER IMPULSE
plot(el) = Sheet12{el}; el=el+1;% Prints PLOTS PERFORMANCE STUDY
plot(el) = Sheet12{el}; el=el+1;% Prints plots of longitudinal Trim with V-n diagram n & Variable Speed
plot(el) = Sheet12{el}; el=el+1; % Prints plots of Performance for Variable V and mass
prefix = Sheet12{el}; el=el+1;% Selection of String Characters that Define the AC selection
plot_individuals = Sheet12{el}; el=el+1;% Shows plots individually or collecting similar plots	
fname = Sheet12{el}; el=el+1;% Selection of String Characters that Define the AC selection storing location	
LS = Sheet12{el}; el=el+1; % Line size	LS
FS = Sheet12{el}; el=el+1; % Text Font size	FS
LFS = Sheet12{el}; el=el+1; % Legend Font Size	LFS
Video_3D = Sheet12{el}; el=el+1; % Saves video 3D	Video_3D
SAVE_FIGS = Sheet12{el}; el=el+1; % saves the plots: fig, jpg and pdf	SAVE_FIGS

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
el=1;

% Taxy
temp_local_taxy = Sheet13{el}; el=el+1;% 	TAXY: TEMPERATURA LOCAL (Celsius)
h_inicial_taxy = Sheet13{el}; el=el+1;% 	TAXY: INITIAL ALTITUDE
P_local_taxy = Sheet13{el}; el=el+1;% 	TAXY: Local Pressure
delta_T_taxy = Sheet13{el}; el=el+1;% 	TAXY: PALANCA DE RALENTI EN TAXI = 0.05
V_taxy = Sheet13{el}; el=el+1;% 	TAXY:  VELOCIDAD A LA QUE HACE EL TAXI (m/s)
t_taxy = Sheet13{el}; el=el+1;% 	TAXY: TIEMPO DE ESPERA EN TAXI (s)
% TakeOff
temp_local_TO = Sheet13{el}; el=el+1;% 	TAKEOFF: TEMPERATURA LOCAL (Celsius)
h_inicial_TO = Sheet13{el}; el=el+1;% 	TAKEOFF: ALTURA LOCAL (m)
P_local_TO = Sheet13{el}; el=el+1;% 	TAKEOFF: PRESION LOCAL (Pa)
mu_TO = Sheet13{el}; el=el+1;% 	TAKEOFF: COEFICIENTE DE FRICCION CON LA PISTA (MU)
h_obstacle_TO = Sheet13{el}; el=el+1;% 	TAKEOFF: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
gamma_climb_TO = Sheet13{el}; el=el+1;% 	TAKEOFF: GAMMA DE SUBIDA MINIMO
delta_T_TO = Sheet13{el}; el=el+1;% 	TAKEOFF: PALANCA DE GASES PARA DESPEGUE
% Climb
h_inicial_cl = Sheet13{el}; el=el+1;% 	CLIMB : ALTURA INICIAL - [m]
h_final_cl = Sheet13{el}; el=el+1;% 	CLIMB : ALTURA FINAL - [m]
gamma_cl = Sheet13{el}; el=el+1;% 	CLIMB : GAMMA DE SUBIDA - [-]
Mach_cl = Sheet13{el}; el=el+1;% 	CLIMB : MACH DE VUELO  - [-]
TAS_cl = Sheet13{el}; el=el+1;% 	CLIMB : VELOCIDAD TAS  - [m/s]
EAS_cl = Sheet13{el}; el=el+1;% 	CLIMB : VELOCIDAD EAS  - [m/s]
delta_T_cl = Sheet13{el}; el=el+1;%	CLIMB : PALANCA DE GASES  - [-]
V_ini_cl = Sheet13{el}; el=el+1;% 	CLIMB : VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
V_fin_cl = Sheet13{el}; el=el+1;% 	CLIMB : VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
% Cruise
h_inicial_cr = Sheet13{el}; el=el+1;%	CRUISE: ALTURA INICIAL
dist_final_cr = Sheet13{el}; el=el+1;%	CRUISE: DISTANCIA FINAL
V_cr = Sheet13{el}; el=el+1;%	CRUISE: VELOCIDAD
delta_T_cr = Sheet13{el}; el=el+1;%	CRUISE:  PALANCA DE GASES
V_ini_cr = Sheet13{el}; el=el+1;%	CRUISE:  VELOCIDAD INICIAL
V_fin_cr = Sheet13{el}; el=el+1;%	CRUISE: VELOCIDAD FINAL
fuel_cr = Sheet13{el}; el=el+1;% - [kg] % 7: COMBUSTIBLE A QUEMAR
Cd0_cr = Sheet13{el}; el=el+1;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
k1_cr = Sheet13{el}; el=el+1;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
k2_cr = Sheet13{el}; el=el+1;% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
% Turn
h_inicial_tr = Sheet13{el}; el=el+1;%	TURN:  ALTURA INICIAL
t_final_tr = Sheet13{el}; el=el+1;%	TURN: TIEMPO FINAL (seg)
V_turn  = Sheet13{el}; el=el+1;% TURN: turn velocity
delta_T_tr = Sheet13{el}; el=el+1;%	TURN: PALANCA DE GASES
phi_tr = Sheet13{el}; el=el+1;% 	TURN: ANGULO DE ALABEO (rads)
V_psi = Sheet13{el}; el=el+1;%	TURN: VELOCIDAD DE GUIÑADA (rads/seg)
n_tr = Sheet13{el}; el=el+1;%	TURN: FACTOR DE CARGA
R_tr = Sheet13{el}; el=el+1;%	TURN: RADIO DE GIRO (m)
% Descent
h_inicial_d = Sheet13{el}; el=el+1;%	DESCENT: ALTURA INICIAL
h_final_d = Sheet13{el}; el=el+1;%	DESCENT:  ALTURA FINAL
gamma_d = Sheet13{el}; el=el+1;%	DESCENT: GAMMA
V_d = Sheet13{el}; el=el+1;%	DESCENT: VELOCIDAD DE DESCENSI
EAS_d = Sheet13{el}; el=el+1;%	DESCENT : VELOCIDAD TAS  - [m/s]
TAS_d = Sheet13{el}; el=el+1;%	DESCENT : VELOCIDAD EAS  - [m/s]
delta_T_d = Sheet13{el}; el=el+1;%	DESCENT : PALANCA DE GASES  - [-]
V_ini_d = Sheet13{el}; el=el+1;%	DESCENT : VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
V_fin_d = Sheet13{el}; el=el+1;%	DESCENT : VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
% Turn Wait
h_inicial_tr_wt = Sheet13{el}; el=el+1;%	TURN WAIT:  ALTURA INICIAL
t_final_tr_wt = Sheet13{el}; el=el+1;%	TURN WAIT: TIEMPO FINAL (seg)
V_turn_wt = Sheet13{el}; el=el+1;%	TURN WAIT: turn velocity
delta_T_tr_wt = Sheet13{el}; el=el+1;%	TURN WAIT: PALANCA DE GASES
phi_tr_wt = Sheet13{el}; el=el+1;%	TURN WAIT: ANGULO DE ALABEO (rads)
V_psi_wt = Sheet13{el}; el=el+1;%	TURN WAIT: VELOCIDAD DE GUIÑADA (rads/seg)
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
