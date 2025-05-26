%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANAID - Stability ANAlysis Interactive Design Tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
path(pathdef)
clc

% Call function that Defines Generic Path so that files can be organized in folders
get_add_path

% Load the data from SANAID for the Original Aircraft Configuration and
% Scaled to 25%
% load('Geo_tier.mat')

% Units conversion
conv_UNITS = conversion_UNITS;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
g = conv_UNITS.g;

%% Prompts User for reading process
[case_AC OUTPUT_read_XLSX] = Initial_prompt_advanced;

%% MATLAB Compatibility
% Flag that determines MATLAB incompatibility (For example functions that do not work in MATLAB 2017 )
MATLAB_in = OUTPUT_read_XLSX.MATLAB_flags.MATLAB_in;
CHECK_Efficiency = OUTPUT_read_XLSX.MATLAB_flags.CHECK_Efficiency;
Detailed_Profile = OUTPUT_read_XLSX.MATLAB_flags.Detailed_Profile;

%% Check the efficiency of the coding
% CHECK_Efficiency = 0;
% Detailed_Profile = 0;
% if CHECK_Efficiency == 1
%     if Detailed_Profile
%         profile on -history
%     else
%         profile on
%     end
% end

%% Plots graphics
%% Defines the plots that will be plotted, and in which order
% Flags for plotting the figures, defiens the order and the number of plots
% Case 1 PRINT_PLOTS_XFLR5 = 0; % Prints plots - Aero
% Case 2 PRINT_PLOTS_XFLR5_POLAR = 0; % Prints plots - Aero
% Case 3 PRINT_PLOTS_XAC = 0; % print Plots for Stimation ofXAC
% Case 4 PRINT_PLOTS_CT_Model = 0; % print Plots for Propulsive Models
% Case 5 PRINT_PLOTS_3D = 0; % Prints plots for 3D
% Case 6 PRINT_PLOTS_STABILITY_SM = 0; % prints plots of SM analysis
% Case 7 PRINT_PLOTS_TRIM_LON = 0; % prints plots of longitudinal Trim
% Case 8 PRINT_PLOTS_TRIM_LON_VAR = 0; % prints plots of longitudinal Trim with Variable mass & Variable Speed
% Case 9 PRINT_PLOTS_TRIM_LAT = 0; % prints plots of lateral Trim
% Case 10 PRINT_PLOTS_TURNING_LAT = 0; % prints plots of lateral Turning
% Case 11 PLOTS STABILITY DERIVATIVES FOR VAR MASS AND VELOCITY
% Case 12 PLOTS STABILITY ANALYSIS FOR VAR MASS AND VELOCITY
% Case 13 PERFORMANCE_STUDY
% assigns from the Excel the Figures that will be plotted
plot = OUTPUT_read_XLSX.PLOT_flags.plot;
for i=1:length(plot)
    if plot(i) == 1
        PLOTS(i) = i;
    end
end
%% PLOTS Can be overwritten by defining PLOTS
% PLOTS = [5];

%% Defines options for the plots
[mark_Type COLOR_scheme, Plot_Options] = Generate_Plot_Options_init;
Plot_Options.MATLAB_in = MATLAB_in;
% Scaling Factor
SF = OUTPUT_read_XLSX.AC_Data_flags.SF;

% Defines the flag to determine method to Estimate weights
% CASE1 = Weight_Estimation from composite - Cefiro III
% CASE2 = Weight_Estimation from wood - Cefiro I
% CASE3 = Factores lineales
% CASE4 = Exact Extimation
Weight_Estimation = OUTPUT_read_XLSX.AC_Data_flags.Weight_Estimation;

%% Estimation of prop diameter, just preliminary for geometric conditions
% Stores the flags
% Propulsive_flags.type_battery = type_battery; %
% Propulsive_flags.Engine_loc = Engine_loc; %
% Propulsive_flags.Engine_conf = Engine_conf; %
D_prop = OUTPUT_read_XLSX.Propulsive_flags.D_prop; %

%% Aircraft type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
AC_type = OUTPUT_read_XLSX.AC_Data_flags.AC_type;

% [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_EMERGENTIA(SF,conv_UNITS,AC_type,case_AC)
[AC_CONFIGURATION] = Generation_AC_configuration(SF,conv_UNITS,OUTPUT_read_XLSX,case_AC);
AC_CONFIGURATION.type_battery = OUTPUT_read_XLSX.Propulsive_flags.type_battery;

% Geo_input_tier = Generation_Input_Geometric_Data_EMERGENTIA(conv_UNITS,AC_CONFIGURATION,SF);
Geo_input_tier = Generation_Input_Geometric_Data_v2(conv_UNITS,AC_CONFIGURATION,SF,OUTPUT_read_XLSX);

% Determines the fuselage than will be shown
% CASE_fuse = 1 - Just STL of fuselage
% CASE_fuse = 2 - STL of the entire aircraft
% CASE_fuse = 3 - Fuselage from XFLR5
CASE_fuse = OUTPUT_read_XLSX.AC_Data_flags.CASE_fuse;
% Determine the type of scaling
% ESCALADO = 0; % no scaling Ratio =1
% ESCALADO = 1; % scaling acording to desired for each 3 axix
% ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
ESCALADO = OUTPUT_read_XLSX.AC_Data_flags.ESCALADO;
% [XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_EMERGENTIA;
[XFLR5_DATA,XFLR5_file,STL_file] = Extraction_fuselage_data(OUTPUT_read_XLSX);

Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX,case_AC);
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    save Geo_tier.mat Geo_tier
end
        
% type_battery used
% case 1 LiFePO4
% case 2 LiPo
% case 3 FuelCells
type_battery = OUTPUT_read_XLSX.Propulsive_flags.type_battery; %
AC_CONFIGURATION.type_battery = type_battery;
alpha=0;
beta =0;

%% identifies the aerodynamic surfaces being used
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Nac = AC_CONFIGURATION.Nac;

% Defines S ref
S_ref = Geo_tier.S_ref;

%% Generation of Geometry
% [Body_Geo,meshData] = Read_Fuselage_Data(Geo_tier,OUTPUT_read_XLSX,XFLR5_DATA,XFLR5_file,STL_file); % Defines Fuselage DATA
% [Body_Geo,meshData] = Generation_Fuselage_Data(Geo_tier,OUTPUT_read_XLSX,XFLR5_DATA,XFLR5_file,STL_file,conversion_UNITS); % Defines Fuselage DATA
[Body_Geo,meshData] = Generation_Fuselage_Data_old(Geo_tier,XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file,SF);% Defines Propulsion DATA

%% Defines Estimation of Weights according to Cefiro III densities
Weight_tier = Generation_Weight_Data(Geo_tier,Body_Geo,AC_CONFIGURATION,conv_UNITS,Weight_Estimation,OUTPUT_read_XLSX);
m_TOW = Weight_tier.m_TOW;
% Save MAT Structure
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    save Weight_tier.mat Weight_tier
end

% Propulsion data for combustion Engine
Propulsion = Propulsion_Data_CE(AC_CONFIGURATION);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stores information for the plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [mark_Type] = Generates_plot_Info; % Generates the style of the lines and the Legend

%% Data for analysis of XFLR5 DATA
% Determines the plots that want to show
% Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for W_1
% Defines the number of Aerodynamic Cases to be analyzed: See "read_aero_files_Aug2018.m" to select them
% The VECTOR in order to analyze the results of XFLR consist of 3 grous of
% vector each for: w1, w2, vtp, and full airplane
% compare = 1 elements: w1,w2, vtp... alone
% compare = 2 elements: w1,w2, vtp... alone
% compare = 3 elements: w1,w2, vtp... alone
compare = 1;
VECTOR_XFLR5.compare = compare;

%% Generates the file for Aerodynamic Data
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        switch compare
            case 1 % compare = 1 elements: w1,w2, vtp... alone
                % W1 cases to be eplotted
                VECTOR_XFLR5.v1 = [14];
            case 2 % compare = 2 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [15,16];
                VECTOR_XFLR5.v2 = [20,21,22,23,24,25];
            case 3 % compare = 3 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [15];
                VECTOR_XFLR5.v2 = [16];
                VECTOR_XFLR5.v3 = [18];
        end
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        switch compare
            case 1 % compare = 1 elements: w1,w2, vtp... alone
                % W1 cases to be eplotted
                VECTOR_XFLR5.v1 = [14];
            case 2 % compare = 2 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [15,16];
                VECTOR_XFLR5.v2 = [20,21,22,23,24,25];
            case 3 % compare = 3 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [15];
                VECTOR_XFLR5.v2 = [16];
                VECTOR_XFLR5.v3 = [18];
        end
    case 3 % case_AC = 3 - PEPIÑO XXL
        switch compare
            case 1 % compare = 1 elements: w1,w2, vtp... alone
                % W1 cases to be eplotted
                VECTOR_XFLR5.v1 = [1,2];
            case 2 % compare = 2 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [1];
                VECTOR_XFLR5.v2 = [2];
            case 3 % compare = 3 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [];
                VECTOR_XFLR5.v2 = [];
                VECTOR_XFLR5.v3 = [];
        end
    case 4 % case_AC = 4 COMERCIAL
        switch compare
            case 1 % compare = 1 elements: w1,w2, vtp... alone
                % W1 cases to be eplotted
                VECTOR_XFLR5.v1 = [1,2];
            case 2 % compare = 2 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [1];
                VECTOR_XFLR5.v2 = [2];
            case 3 % compare = 3 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [];
                VECTOR_XFLR5.v2 = [];
                VECTOR_XFLR5.v3 = [];
        end
    case 5 % case_AC = 5 WIG
        switch compare
            case 1 % compare = 1 elements: w1,w2, vtp... alone
                % W1 cases to be eplotted
                VECTOR_XFLR5.v1 = [1,2];
            case 2 % compare = 2 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [1];
                VECTOR_XFLR5.v2 = [2];
            case 3 % compare = 3 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [1];
                VECTOR_XFLR5.v2 = [2];
                VECTOR_XFLR5.v3 = [3];
        end
    case 6 % case_AC = 1 - CERVERA
        switch compare
            case 1 % compare = 1 elements: w1,w2, vtp... alone
                % W1 cases to be eplotted
                VECTOR_XFLR5.v1 = [14];
            case 2 % compare = 2 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [15,16];
                VECTOR_XFLR5.v2 = [20,21,22,23,24,25];
            case 3 % compare = 3 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [15];
                VECTOR_XFLR5.v2 = [16];
                VECTOR_XFLR5.v3 = [18];
        end
    case 12 % case_AC = 1 - ALO
        switch compare
            case 1 % compare = 1 elements: w1,w2, vtp... alone
                % W1 cases to be eplotted
                VECTOR_XFLR5.v1 = [1];
            case 2 % compare = 2 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [1];
                VECTOR_XFLR5.v2 = [2];
            case 3 % compare = 3 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [1];
                VECTOR_XFLR5.v2 = [2];
                VECTOR_XFLR5.v3 = [1];
        end
    case 13 % case_AC = 13 - ALO Fuel CEll
        switch compare
            case 1 % compare = 1 elements: w1,w2, vtp... alone
                % W1 cases to be eplotted
                VECTOR_XFLR5.v1 = [1];
            case 2 % compare = 2 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [1];
                VECTOR_XFLR5.v2 = [2];
            case 3 % compare = 3 elements: w1,w2, vtp... alone
                VECTOR_XFLR5.v1 = [1];
                VECTOR_XFLR5.v2 = [2];
                VECTOR_XFLR5.v3 = [1];
        end
end

% Desired Static Margin
SM_des = OUTPUT_read_XLSX.Stability_flags.SM_des;

%% Generates the file for Aerodynamic Data
index_w1 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1;
index_w2 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2;
index_w3 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w3;
i_w1 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1;
i_w2 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2;
i_w3 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w3;
Conf = OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf;
        
if W1 == 1
    %% Stores data
    alpha_selected_w1 = i_w1; % (degs)
    Design_criteria.index_w1 = index_w1;
    Design_criteria.alpha_selected_w1 = alpha_selected_w1;
    % Design choices for the incidence of the different surfaces
    Design_criteria.i_w1 = i_w1*D2R; % incidence of Front Wing
end

if HTP == 1 || Vee == 1
    %% Stores data
    alpha_selected_w2 = i_w2; % (degs)
    Design_criteria.index_w2 = index_w2;
    Design_criteria.alpha_selected_w2 = alpha_selected_w2;
    % Design choices for the incidence of the different surfaces
    Design_criteria.i_w2 = i_w2*D2R; % incidence of Rear Wing
end

if VTP == 1
    Design_criteria.index_VTP = index_VTP;
    % Design choices for the incidence of the different surfaces
    Design_criteria.i_VTP = 0*D2R; % incidence of Rear Wing
end

%% Flight Safe Margin to calculate aerodynamic properties
Flight_SF = OUTPUT_read_XLSX.Performance_pre_flags.Flight_SF;
Design_criteria.Flight_SF = Flight_SF;
%% Revisar si pasar al Excel

%% Propulsion Generation
alpha_f = 0*D2R;
beta_f = 0*D2R;

%% Initial Estimate of Design altitude and velocity
%% Needs to introduce preliminary results for initial estimate

%% Generates the file for Aerodynamic Data
% Indentifies preliminary performance results
Performance_preliminar = OUTPUT_read_XLSX.Performance_pre_flags;

% Atmospheric conditions
[Temp_init,rho_init,p_init,a_init]=atmos_inter_mio(Performance_preliminar.h);
Performance_preliminar.Temp = Temp_init;
Performance_preliminar.rho = rho_init;
Performance_preliminar.p = p_init;
Performance_preliminar.a = a_init;
Mach_init = Performance_preliminar.V/a_init;
Performance_preliminar.Mach = Mach_init;
q_inf_init = 0.5*rho_init*(Performance_preliminar.V)^2;
Performance_preliminar.q_inf = q_inf_init;

%% Defines the Reading files
if OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY == 1
    %% Generates the file for Aerodynamic Data
    switch case_AC
        case 1 % case_AC = 1 - EMERGENTIA 1:1
            [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 2 % case_AC = 2 - EMERGENTIA 1:2
            [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 3 % case_AC = 3 - PEPIÑO XXL
            [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_PEPINOXXL(Performance_preliminar);
        case 4 % case_AC = 4 - COMERCIAL
            [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 5 % case_AC = 5 - WIGL
            [casos,prefix,mark_legend,X_OC] = read_aero_files_Feb2020_v3_WIG(Performance_preliminar);
        case 6 % case_AC = 6 - CERVERA
            [casos prefix mark_legend] = read_aero_files_CERVERA(Performance_preliminar);
        case 7 % Existing Aircraft = 7
            [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 8 % case_AC = 8 - TAMIZ
            [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
            [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
            [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 11 % case_AC = 11 - EMERGENTIA Manufacturing
            [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 12 % case_AC = 12 - ALO
            [casos,prefix,mark_legend,X_OC] = read_aero_files_ALO_v1(Performance_preliminar);
        case 13 % case_AC = 13 - ALO Fuel Cell
            [casos,prefix,mark_legend,X_OC] = read_aero_files_ALO_v1(Performance_preliminar);
    end
        prefix = OUTPUT_read_XLSX.PLOT_flags.prefix;

    % Read Aerodynamic Data
    [DATA_Ae] = Read_data_Aero(casos);
    
    % Fusion Aerodynamic Propert1ies
    [Aero_TH] = Polar_Literatura_v1(Performance_preliminar,Geo_tier,conv_UNITS,Body_Geo,AC_CONFIGURATION,Conf);
    % Save MAT Structure
    if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
        save Aero_TH.mat Aero_TH
    end

    %% Initializes figures
    Fig = 0;
    
    % Premiliminary Aerodynamic Design
%     [Aero,DATA_PL,Performance] = Generate_Aero_Data_2020_v1(DATA_Ae,Design_criteria,Performance_preliminar,Geo_tier,Weight_tier,conv_UNITS,AC_CONFIGURATION);
    [Aero,DATA_PL,Performance] = Generate_Aero_Data_2021_v1(DATA_Ae,Design_criteria,Performance_preliminar,Geo_tier,Weight_tier,conv_UNITS,AC_CONFIGURATION, Body_Geo);
pause
    % [Aero,DATA_PL,Fig] = Aerodynamic_Design_2020_v1(Geo_tier,...
    %     Weight_tier,conv_UNITS,Design_criteria,Performance,CASOS,degrees_XAC,Fig,XFLR5_DATA,...
    %     MAC_Estimation,DATA_Ae,X_OC,Plot_Options,VECTOR_XFLR5);
    % Save MAT Structure
    if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
        save Aero.mat Aero
    end
end

%% Data for analysis of Prop DATA
% Determines the plots that want to show
% Compares Propeller properties for different sources
% The VECTOR in order to analyze the results of Props consist of 2 grous of
% vector each prop model
% compare_prop = 1 compares APC Models alone
% compare_prop = 2 compares APC with wind tunnel models
% VECTOR_Prop.compare_props = compare_props;

%% Propulsive Model Wind Tunnel Model 1
% Number of Prop used
% 1 - APC 20x8
% 2 - APC 22x10
% 3 - APC 22x12
% 4 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
% 5 - APC 21x13W
% 6 - APC 21x14

%% Propulsive Model Wind Tunnel Model 2 (Rai models with alpha)
% Number of Prop used
% 1 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/

%% Generates the file for Aerodynamic Data
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 3 % case_AC = 3 - PEPIÑO XXL
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 4 % case_AC = 4 COMERCIAL
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 5 % case_AC = 5 WIG
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 6 % case_AC = 1 - CERVERA
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 7 % case_AC = 7 - Existing Aircraft
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 8 % case_AC = 8 - TAMIZ
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 11 % case_AC = 11 - EMERGENTIA Manufactured
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 12 % case_AC = 12 - ALO 
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 13 % case_AC = 13 - ALO Fuel Cell
        switch compare
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
end

%% Generates the file for Prop Data to be used in the codes
% actualizes the Prop geometry
% Select:
% Propulsion model: SEE read_prop_files_May2020.m FOR DETAILS OF PROP MODEL!!
Flag.APC_model = 1;
Flag.Wind_tunnel_model = 1;
Flag.compare_prop_models = 1;

%% Propulsive Model
%% Compares 3 prop models
% - Model 1 - APC data
% - Model 2 - Wind tunnel data for different props
% 1 - APC 20x8
% 2 - APC 22x10
% 3 - APC 22x12
% 4 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
% 5 - APC 21x13W
% 6 - APC 21x14
% - Model 3 - Wind tunnel data for different angle of attack for APC 22x12W
model_prop = OUTPUT_read_XLSX.Propulsive_flags.model_prop; %
prop_selec_APC = OUTPUT_read_XLSX.Propulsive_flags.prop_selec_APC; %
prop_selec_WT1 = OUTPUT_read_XLSX.Propulsive_flags.prop_selec_WT1; %
prop_selec_WT2 = OUTPUT_read_XLSX.Propulsive_flags.prop_selec_WT2; %
Prop_selection.model_prop = model_prop;
% Selects the prop for
% - Model 1 - APC data
% - Model 2 - Wind tunnel data for different props
% Stores info
Prop_selection.prop_selec_APC = prop_selec_APC;
Prop_selection.prop_selec_WT1 = prop_selec_WT1;
Prop_selection.prop_selec_WT2 = prop_selec_WT2;
[Prop_data] = Generation_Propulsion_Data(SF,D_prop,AC_CONFIGURATION,Prop_selection);
% Defines Propulsion DATA
% [Fig] = plot_prop_APC(Data_P,Plot_Options,Fig,prefix,VECTOR_Prop)

if OUTPUT_read_XLSX.STUDY_flags.PROP_STUDY == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Determine Minim Prop Diameter %%%%%%%%%%%%%%%
    %%
    %% Lets user select if wants to conduct study for the propellers
    answer = questdlg('Would you like to conduct PROPELLER''s STUDY?', ...
        'Prop Limits', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            %         disp([answer ' coming right up.'])
            Study_Prop_Limits = 1;
        case 'No'
            %         disp([answer ' coming right up.'])
            Study_Prop_Limits = 0;
    end
    
    if Study_Prop_Limits == 1
        % Determination
        switch case_AC
            case 1 % case_AC = 1 - EMERGENTIA 1:1
                pp = 0.90; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 1; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1;
            case 2 % case_AC = 1 - EMERGENTIA 1:2
                pp = 0.80; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 0.01; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1;
            case 3 % case_AC = 3 - PEPIÑO XXL
                pp = 0.95; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 0;
                Study_horizontal = 0;
            case 6 % case_AC = 6 - CERVERA
                pp = 1; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 1; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1    ;
            case 12 % case_AC = 12 - ALO
                pp = 1; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 1; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1    ;
            case 13 % case_AC = 13 - ALO Fuel Cell
                pp = 1; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 1; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1    ;
        end
        % uses fsolve to find associated rpms
        solve_w_fero = 0;
        
        %% Conducts study of sensitivity analysis of the Performance
        %% Lets user select if wants to conduct study for the propellers
        answer = questdlg('Would you like to conduct Propeller Performance?', ...
            'Prop Limits', ...
            'Yes','No','No');
        % Handle response
        switch answer
            case 'Yes'
                % Selects if 3D Plots are shown
                answer1 = questdlg('Would you like to conduct horizontal flight limits?', ...
                    'Horizontal limits', ...
                    'Yes','No','No');
                % Handle response
                switch answer1
                    case 'Yes'
                        Study_vertical = 1;
                    case 'No'
                        Study_vertical = 0;
                end
                
                % Selects if Contour Plots are shown
                answer2 = questdlg('Would you like to conduct vertical flight limits?', ...
                    '3D plots', ...
                    'Yes','No','No');
                % Handle response
                switch answer2
                    case 'Yes'
                        Study_horizontal = 1;
                    case 'No'
                        Study_horizontal = 0;
                end
                flight_minits_study = 1;
            case 'No'
                %         disp([answer ' coming right up.'])
                flight_minits_study = 0;
        end
        
        if flight_minits_study == 1
            % Flag to determine if vertical flight & horizontal flight
            % prop limits are calculated
            if Study_vertical == 1;
                % Selects to plot Get Vertical Flight Limits that determines the
                % minimum prop diameter for the different maneuvers
                PLOT_Get_Vertical_Flight_Limits = 1;
                [Fig] = Get_Vertical_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Performance,Geo_tier,...
                    AC_CONFIGURATION,N_Vv,Vv_vect,Fig,PLOT_Get_Vertical_Flight_Limits,solve_w_fero,Plot_Options);
                Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
                disp(Warning)
                pause
            end
            
            % Flag to determine if vertical flight & horizontal flight
            % prop limits are calculated
            if Study_horizontal == 1
                % Selects to plot Get Vertical Flight Limits that determines the
                % minimum prop diameter for the different maneuvers
                PLOT_Get_Horizontal_Flight_Limits = 1;
                [Fig] = Get_Horizontal_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Performance,Geo_tier,...
                    AC_CONFIGURATION,N_Vh,Vh_vect,Fig,PLOT_Get_Horizontal_Flight_Limits,solve_w_fero,Plot_Options);
                Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
                disp(Warning)
                pause
            end
        end
    end
    
    %% Conducts study of sensitivity analysis of the Performance
    %% Lets user select if wants to conduct study for the propellers
    answer = questdlg('Would you like to conduct Sensitivity Study To Determine the Propeller Diameter?', ...
        'Prop Limits', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            
            % Selects if 3D Plots are shown
            answer1 = questdlg('Would you like to show Hover Performence vs. D_{prop}?', ...
                'Hover plots', ...
                'Yes','No','No');
            % Handle response
            switch answer1
                case 'Yes'
                    Hover_PLOTS = 1;
                case 'No'
                    Hover_PLOTS = 0;
            end
            
            % Selects if 3D Plots are shown
            answer2 = questdlg('Would you like to show 3D plots?', ...
                '3D plots', ...
                'Yes','No','No');
            % Handle response
            switch answer2
                case 'Yes'
                    plots_3d_sensitivity = 1;
                case 'No'
                    plots_3d_sensitivity = 0;
            end
            
            % Selects if Contour Plots are shown
            answer3 = questdlg('Would you like to show contour plots?', ...
                '3D plots', ...
                'Yes','No','No');
            % Handle response
            switch answer3
                case 'Yes'
                    plots_contour_sensitivity = 1;
                case 'No'
                    plots_contour_sensitivity = 0;
            end
            
            % Selects the minimum Energy plots
            answer4 = questdlg('Would you like to show plots searching for minimum Energy?', ...
                'Minimum Energy', ...
                'Yes','No','No');
            % Handle response
            switch answer4
                case 'Yes'
                    special_PLOTS = 1; % represents the plots searching for minimum Energy
                case 'No'
                    special_PLOTS = 0; % represents the plots searching for minimum Energy
            end
            
            % Selects the minimum Energy plots
            answer5 = questdlg('Would you like to show plots for Performance Study vs battery mass variation for hover flight?', ...
                'Mass Variation', ...
                'Yes','No','No');
            % Handle response
            switch answer5
                case 'Yes'
                    variation_mass = 1; % represents the plots searching for minimum Energy
                case 'No'
                    variation_mass = 0; % represents the plots searching for minimum Energy
            end
            
            % Selects if Intersection curves between Vertical and
            % horizontal Flight are shown
            
            answer3 =questdlg ('Would you like to perform the Diameter optimization?',...
                '3D plots',...
                'Yes','No','No');
            %Handle response
            switch answer3
                case 'Yes'
                    diameter_optim=1;
                case 'No'
                    diameter_optim=0;
            end
            
            PERFORMANCE_sensitivity = 1;
            % Saves selection of plotting options
            PLOTS_Prop_Performance.Hover_PLOTS = Hover_PLOTS;
            PLOTS_Prop_Performance.PERFORMANCE_sensitivity = PERFORMANCE_sensitivity;
            PLOTS_Prop_Performance.plots_contour_sensitivity = plots_contour_sensitivity;
            PLOTS_Prop_Performance.plots_3d_sensitivity = plots_3d_sensitivity;
            PLOTS_Prop_Performance.special_PLOTS = special_PLOTS;
            PLOTS_Prop_Performance.variation_mass = variation_mass;
            
        case 'No'
            %         disp([answer ' coming right up.'])
            PERFORMANCE_sensitivity = 0;
            % Saves selection of plotting options
            PLOTS_Prop_Performance.Hover_PLOTS = 0;
            PLOTS_Prop_Performance.PERFORMANCE_sensitivity = 0;
            PLOTS_Prop_Performance.plots_contour_sensitivity = 0;
            PLOTS_Prop_Performance.plots_3d_sensitivity = 0;
            PLOTS_Prop_Performance.special_PLOTS = 0;
            
    end
    
    PLOTS_SENSITIVITY_VERTICAL = 1;
    PLOTS_SENSITIVITY_HORIZONTAL = 1;
    % Uses Fzero which takes longer time
    solve_w_fero = 0;
    
    if PERFORMANCE_sensitivity == 1
        % Selects the Engine-Prop configuration according to the Limit Study
        N_prop = 500; % number of prop elements to analyze between the limits
        pp_D_prop_min = 0.75;
        pp_D_prop_max = 1.5;
        N_contour_lines = 25; % number of contour lines
        
        switch case_AC
            case 1 % case_AC = 1 - EMERGENTIA 1:1
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 2 % case_AC = 1 - EMERGENTIA 1:2
                D_prop = 22*2.54/100;
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 3 % case_AC = 3 - PEPIÑO XXL
                D_prop = 32*2.54/100;
                % Actualizes Prop Data
                SF_prop = 0.85;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                % Flag to determine if vertical flight
                sensitivity_vertical = 0;
                sensitivity_horizontal = 1;
            case 6 % case_AC = 6 - CERVERA
                D_prop = 30*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = pp; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 9 % case_AC = 1 - EMERGENTIA Wind Tunnel w2
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel w1 (sweep)
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 11 % case_AC = 11 - EMERGENTIA Manufactured
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 12 % case_AC = 12 - ALO
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = pp; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 13 % case_AC = 13 - ALO Fuel Cell
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = pp; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
        end
        
        % Loop that conducts sensitivity study for vertical flight
        if sensitivity_vertical == 1;
            % Sensitivity analysis for the Vertical Performance
            [Xv,Yv,Zv_study,Fig] = Sensitivity_Study_Vertical_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,...
                Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_VERTICAL,Performance,solve_w_fero,...
                N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines,Vv_min,PLOTS_Prop_Performance,Prop_selection,prefix,Plot_Options);
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            query_properties = 1;
            while query_properties == 1
                % Lets user select if wants to conduct study
                answer_prop = questdlg('Would you like to calculate properties for a given prop diameter and vertical speed?', ...
                    'Prop Limits', ...
                    'Yes','No','Yes');
                % Handle response
                switch answer_prop
                    case 'Yes'
                        prompt = {'Enter the prop diameter (inches):','Enter vertical VTOl speed (m/s):'};
                        dlgtitle = 'Performance properties';
                        dims = [1 35];
                        definput = {'30','5'};
                        answer_sv = inputdlg(prompt,dlgtitle,dims,definput);
                        
                        D_propin_q = str2num(answer_sv{1});
                        Vv_q = str2num(answer_sv{2});
                        
                        z_Tv = interp2(Xv,Yv,Zv_study.Z_T,D_propin_q,Vv_q,'cubic');
                        z_Pev = interp2(Xv,Yv,Zv_study.Z_Pe,D_propin_q,Vv_q,'cubic');
                        %                         z_Pe_CHv = interp2(Xv,Yv,Zv_study.Z_Pe_CH,D_propin_q,Vv_q,'cubic');
                        z_Energyv = interp2(Xv,Yv,Zv_study.Z_Energy,D_propin_q,Vv_q,'cubic');
                        z_Energy_CHv = interp2(Xv,Yv,Zv_study.Z_Energy_CH,D_propin_q,Vv_q,'cubic');
                        z_Qv = interp2(Xv,Yv,Zv_study.Z_Q,D_propin_q,Vv_q,'cubic');
                        z_deltav = interp2(Xv,Yv,Zv_study.Z_delta,D_propin_q,Vv_q,'cubic');
                        z_RPMv = interp2(Xv,Yv,Zv_study.Z_RPM,D_propin_q,Vv_q,'cubic');
                        z_ethav = interp2(Xv,Yv,Zv_study.Z_etha,D_propin_q,Vv_q,'cubic');
                        z_Battery_massv = interp2(Xv,Yv,Zv_study.Z_Battery_mass,D_propin_q,Vv_q,'cubic');
                        z_Battery_mass_CHv = interp2(Xv,Yv,Zv_study.Z_Battery_mass_CH,D_propin_q,Vv_q,'cubic');
                        z_Pe_CHv = 0;
                        
                        uiwait(msgbox({['Total Thrust (climb) = ',num2str(z_Tv),' (N)']...
                            ['Total Electric Power (climb) = ',num2str(z_Pev),' (W)']...
                            ['Total Electric Power (climb & hover) = ',num2str(z_Pe_CHv),' (W)']...
                            ['Total Energy (climb) = ',num2str(z_Energyv),' (kW-h)']...
                            ['Total Energy (climb & hover) = ',num2str(z_Energy_CHv),' (kW-h)']...
                            ['Torque = ',num2str(z_Qv),' (Nm)']...
                            ['Throttle = ',num2str(z_deltav),' (%)']...
                            ['RPM = ',num2str(z_RPMv),' (RPM)']...
                            ['\eta = ',num2str(z_ethav)]...
                            ['m_{bat} (climb) = ',num2str(z_Battery_massv),' (kg)']...
                            ['m_{bat} (climb & ) hover= ',num2str(z_Battery_mass_CHv),' (kg)']}));
                        query_properties = 1;
                    case 'No'
                        query_properties = 0;
                end
            end
        end
        
        % Loop that conducts sensitivity study for horizontal flight
        if sensitivity_horizontal == 1
            [Xh,Yh,Zh_study,Fig] = Sensitivity_Study_Horizontal_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,...
                Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_VERTICAL,Performance,solve_w_fero,...
                N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines,PLOTS_Prop_Performance,Prop_selection,prefix,Plot_Options);
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            query_properties = 1;
            Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
            disp(Warning)
            pause

            while query_properties == 1
                % Lets user select if wants to conduct study
                answer_prop = questdlg('Would you like to calculate properties for a given prop diameter and horizontal speed?', ...
                    'Prop Limits', ...
                    'Yes','No','Yes');
                % Handle response
                switch answer_prop
                    case 'Yes'
                        prompt = {'Enter the prop diameter (inches):','Enter horizontal speed (m/s):'};
                        dlgtitle = 'Performance properties';
                        dims = [1 35];
                        definput = {'30','30'};
                        answer_sv = inputdlg(prompt,dlgtitle,dims,definput);
                        
                        D_propin_q = str2num(answer_sv{1});
                        Vv_q = str2num(answer_sv{2});
                        
                        z_Th = interp2(Xh,Yh,Zh_study.Z_T,D_propin_q,Vv_q,'cubic');
                        z_Peh = interp2(Xh,Yh,Zh_study.Z_Pe,D_propin_q,Vv_q,'cubic');
                        z_Energyh = interp2(Xh,Yh,Zh_study.Z_Energy,D_propin_q,Vv_q,'cubic');
                        z_Qh = interp2(Xh,Yh,Zh_study.Z_Q,D_propin_q,Vv_q,'cubic');
                        z_deltah = interp2(Xh,Yh,Zh_study.Z_delta,D_propin_q,Vv_q,'cubic');
                        z_RPMh = interp2(Xh,Yh,Zh_study.Z_RPM,D_propin_q,Vv_q,'cubic');
                        z_ethah = interp2(Xh,Yh,Zh_study.Z_etha,D_propin_q,Vv_q,'cubic');
                        z_Battery_massh = interp2(Xh,Yh,Zh_study.Z_Battery_mass,D_propin_q,Vv_q,'cubic');
                        
                        uiwait(msgbox({['Total Thrust = ',num2str(z_Th),' (N)']...
                            ['Total Electric Power = ',num2str(z_Peh),' (W)']...
                            ['Total Energy = ',num2str(z_Energyh),' (kW-h)']...
                            ['Torque = ',num2str(z_Qh),' (Nm)']...
                            ['Throttle = ',num2str(z_deltah),' (%)']...
                            ['RPM = ',num2str(z_RPMh),' (RPM)']...
                            ['\eta = ',num2str(z_ethah)]...
                            ['m_{bat} (climb) = ',num2str(z_Battery_massh),' (kg)']}));
                        query_properties = 1;
                    case 'No'
                        query_properties = 0;
                end
            end
        end
        
        if sensitivity_horizontal == 1 && sensitivity_vertical == 1 && diameter_optim== 1
            
            [Vvopt_mc,Vhopt_mc,Dopt_mc] = diameteroptimization(Xh,Yv,Yh,Zh_study.Z_T,Zh_study.Z_Pe,Zh_study.Z_Energy,Zh_study.Z_etha,...
                Zv_study.Z_T,Zv_study.Z_Pe,Zv_study.Z_Energy,Zv_study.Z_etha);
            
        end
        
    end
end

if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
    %% Enter the number of mission segments
    %% Generates the file for Aerodynamic Data
    type_missions_WF = OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF;
    num_missions_WF = OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF;
    climb_mode = OUTPUT_read_XLSX.PerforMisionSelection_flags.climb_mode;
    cruise_mode = OUTPUT_read_XLSX.PerforMisionSelection_flags.cruise_mode;
    turn_mode = OUTPUT_read_XLSX.PerforMisionSelection_flags.turn_mode;
    descent_mode = OUTPUT_read_XLSX.PerforMisionSelection_flags.descent_mode;
    
    FlightMODE_WF.climb_mode = climb_mode;
    FlightMODE_WF.cruise_mode = cruise_mode;
    FlightMODE_WF.turn_mode = turn_mode;
    FlightMODE_WF.descent_mode = descent_mode;
    
    Segments = Define_Segments_v2(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF,OUTPUT_read_XLSX);
    
    % Weight_tier
    if OUTPUT_read_XLSX.STUDY_flags.STUDY_Weight_Fraction == 1
        [seg_WF] = Generation_Mission_Segments_v3(conv_UNITS,num_missions_WF,type_missions_WF,Performance,case_AC,Segments,FlightMODE_WF);
        Weights = Weight_Fraction(seg_WF,Propulsion,conv_UNITS,Aero_TH,type_missions_WF,Weight_tier,Geo_tier,AC_CONFIGURATION);
    end
    
    %% Performance
    if OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var ==1

        %% Generates the file for Aerodynamic Data
        V_low = OUTPUT_read_XLSX.PerformanceStudy_flags.V_low;
        V_high = OUTPUT_read_XLSX.PerformanceStudy_flags.V_high;
        N_V_VAR_perf = OUTPUT_read_XLSX.PerformanceStudy_flags.N_V_VAR_perf;
        V_single = OUTPUT_read_XLSX.PerformanceStudy_flags.V_single;
        Wp_low = OUTPUT_read_XLSX.PerformanceStudy_flags.Wp_low;
        Wp_high = OUTPUT_read_XLSX.PerformanceStudy_flags.Wp_high;
        N_Wp_VAR_perf = OUTPUT_read_XLSX.PerformanceStudy_flags.N_Wp_VAR_perf;
        W_single = OUTPUT_read_XLSX.PerformanceStudy_flags.W_single;
        Post_processing_PERFORMANCE = OUTPUT_read_XLSX.PerformanceStudy_flags.Post_processing_PERFORMANCE;

%         switch case_AC
            %% PROPULSION DATA
            % 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON
            % 2: NUMERO DE MOTORES
            % 3: EMPUJE/POTENCIA A NIVEL DEL MAR
            % 4: CONSUMO ESPECIFICO
            % 5: AVION CIVIL =1/MILITAR = 2
            % 6: EFICIENCIA DE LA HELICE (ETA_P)
            % 7: DERIVACION(TURBOFANES)
            
                % Enter type of mission segments betwee brackets being
                type_missions(1) = OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF;
                type_missions(2) = 13;
                type_missions(3) = 13;
                num_missions = length(type_missions);
                %% User needs to define the inputs for each segment
                FlightMODE_var.climb_mode(1) = OUTPUT_read_XLSX.PerforMisionSelection_flags.climb_mode;
                FlightMODE_var.climb_mode(2) = OUTPUT_read_XLSX.PerforMisionSelection_flags.climb_mode;
                FlightMODE_var.climb_mode(3) = OUTPUT_read_XLSX.PerforMisionSelection_flags.climb_mode;
                FlightMODE_var.cruise_mode(1) = OUTPUT_read_XLSX.PerforMisionSelection_flags.cruise_mode;
                FlightMODE_var.cruise_mode(2) = OUTPUT_read_XLSX.PerforMisionSelection_flags.cruise_mode;
                FlightMODE_var.cruise_mode(3) = OUTPUT_read_XLSX.PerforMisionSelection_flags.cruise_mode;
                FlightMODE_var.turn_mode(1) = OUTPUT_read_XLSX.PerforMisionSelection_flags.turn_mode;
                FlightMODE_var.turn_mode(2) = OUTPUT_read_XLSX.PerforMisionSelection_flags.turn_mode;
                FlightMODE_var.turn_mode(3) = OUTPUT_read_XLSX.PerforMisionSelection_flags.turn_mode;
                FlightMODE_var.descent_mode(1) = OUTPUT_read_XLSX.PerforMisionSelection_flags.descent_mode;
                FlightMODE_var.descent_mode(2) = OUTPUT_read_XLSX.PerforMisionSelection_flags.descent_mode;
                FlightMODE_var.descent_mode(3) = OUTPUT_read_XLSX.PerforMisionSelection_flags.descent_mode;
                
                % propulsive selected data
                propul = OUTPUT_read_XLSX.Propulsive_flags.propul;
                prop_data(1) = propul(1); % Type of engine
                prop_data(2) = propul(2); % Number of engines
                prop_data(3) = propul(3); % Thrust (lbf) or Power (shp) per engine
                prop_data(4) = propul(4); % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
                prop_data(5) = propul(7); % Normativa
                prop_data(6) = propul(6); % Prop efficiency
                prop_data(7) = propul(5); % By-pass
                Prop_sel = prop_data;
                
                % Range of low and high speed to be analized
                V_low = OUTPUT_read_XLSX.PerformanceStudy_flags.V_low; 
                V_high = OUTPUT_read_XLSX.PerformanceStudy_flags.V_high;
                N_V_VAR_perf = OUTPUT_read_XLSX.PerformanceStudy_flags.N_V_VAR_perf;
                % Value of single speed
                V_single = OUTPUT_read_XLSX.PerformanceStudy_flags.V_single;
                
                % Range of low and high payloads to be analized
                Wp_low = OUTPUT_read_XLSX.PerformanceStudy_flags.Wp_low;
                Wp_high = OUTPUT_read_XLSX.PerformanceStudy_flags.Wp_high;
                N_Wp_VAR_perf = OUTPUT_read_XLSX.PerformanceStudy_flags.N_Wp_VAR_perf;
                % Value of single payload
                W_single = OUTPUT_read_XLSX.PerformanceStudy_flags.W_single;
                Post_processing_PERFORMANCE = OUTPUT_read_XLSX.PerformanceStudy_flags.Post_processing_PERFORMANCE;

        % Creates the variable or single condictions
        if OUTPUT_read_XLSX.STUDY_flags.variable_speed_AP == 1
            V_VAR_perf = linspace(V_low,V_high,N_V_VAR_perf);
        else
            N_V_VAR_perf = 1; % number of iterations
            V_VAR_perf = V_single;
        end
        Plot_Options.V_VAR_perf = V_VAR_perf;
        Plot_Options.N_V_VAR_perf = N_V_VAR_perf;
        
        % Creates the variable or single condictions
        if OUTPUT_read_XLSX.STUDY_flags.variable_weight_AP == 1
            Wp_VAR_perf = linspace(Wp_low,Wp_high,N_Wp_VAR_perf);
        else
            N_Wp_VAR_perf = 1; % number of iterations
            Wp_VAR_perf = W_single;
        end
        Plot_Options.Wp_VAR_perf = Wp_VAR_perf;
        Plot_Options.N_Wp_VAR_perf = N_Wp_VAR_perf;
        
        % defines which series plot as a surface plot
        N_Wp_plot = 1;
        if N_Wp_plot > N_Wp_VAR_perf
            Warning = 'Warning, the Number of iterations is smaller that the number of figures - Paused, Press a key to continue';
            disp(Warning)
            pause
            N_Wp_plot = N_Wp_VAR_perf;
        end
        Plot_Options.N_Wp_plot = N_Wp_plot;
        
        % Conducts Performance studies
        for i = 1:N_V_VAR_perf
            Segments{1}.crucero(3) = V_VAR_perf(i);
            for j = 1:N_Wp_VAR_perf
                
                handles = caracteristicas_avanz_v3(Aero_TH,Prop_data,Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Aero,OUTPUT_read_XLSX);
                handles.pesos(2) = Wp_VAR_perf(j);
                PROP_VAR = prop_data;
                
                empuje = PROP_VAR(3);
                consumo = PROP_VAR(4);
                
                if PROP_VAR(1) == 1,
                    handles.propul(3) = empuje * 4.448221615255; %[N]
                    handles.propul(4) = consumo * 2.832546065 * 10^-5; %[kg/(N*s)]
                else
                    handles.propul(3) = empuje * 745.699872; %[W]
                    handles.propul(4) = consumo * 1.68965941 * 10^-7; %[kg/(W*s)]
                end
                
                % Flight Conditions For Cruise
                h_inicial_cr = Segments{1}.crucero(1);% - [m] % Altura inicial
                dist_final_cr = Segments{1}.crucero(2);% - [m] % 1: DISTANCIA FINAL
                V_cr = V_VAR_perf(i);
                Performance_prelimin = Performance;
                [Data_ATM,Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_cr,Performance_prelimin);
                Mach_cr = V_cr/Data_ATM.a;% - [-] % 2: MACH DE VUELO
                
                Segments = Define_Segments_v2(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF,OUTPUT_read_XLSX);

                Segments{1}.crucero(1) = h_inicial_cr;
                Segments{1}.crucero(2) = dist_final_cr;
                Segments{1}.crucero(4) = Mach_cr;
                CL_opt = sqrt(Aero_TH.CD0/Aero_TH.CD2);
                Segments{1}.crucero(5) = CL_opt;
                
                % Generates the different segments
                [seg] = Generation_Mission_Segments_v3(conv_UNITS,num_missions,type_missions,Performance,case_AC,Segments,FlightMODE_var);
                
                [Weights_AP,Total_Datos,datos] = procesar_mision_v3(seg,num_missions,handles,Segments,FlightMODE_var,Geo_tier,Prop_data,AC_CONFIGURATION);
                Weights_AP_var{i,j} = Weights_AP;
                fuel_total_var{i,j} = Total_Datos.fuel_total;
                tiempo_total_var{i,j} = Total_Datos.tiempo_total;
                distancia_total_var{i,j} = Total_Datos.distancia_total;
                W_var{i,j} = Total_Datos.W;
                datos_var{i,j}.crucero = datos(1).crucero;
            end
            i
            j
        end
        % Cálculo de Datos de Performance
        CI = 0; % kg/s
        density_fuel =  0.809; % kg/l
        cost_fuel = 161.69;%  cts/gal - 7 Feb 2020
        
        if Post_processing_PERFORMANCE == 1
            l2gal = conv_UNITS.l2gal;
            m2km = conv_UNITS.m2km;
            kg2Tm = conv_UNITS.kg2Tm;
            % Generation of variable to plot
            for i=1:N_V_VAR_perf
                for j=1:N_Wp_VAR_perf
                    
                    m_f_PLOT(i,j) = Weights_AP_var{i,j}.m_f;
                    tiempo_total_PLOT(i,j) = tiempo_total_var{i,j};
                    distancia_total_PLOT(i,j) = distancia_total_var{i,j};
                    Cost_fuel(i,j) = (m_f_PLOT(i,j)/density_fuel)*l2gal*cost_fuel;
                    DOC_PLOT(i,j) = (tiempo_total_PLOT(i,j)*CI + Cost_fuel(i,j));
                    CAPM_PLOT(i,j) = DOC_PLOT(i,j)/((distancia_total_PLOT(i,j)*m2km)*(Wp_VAR_perf(j)*kg2Tm));
                end
            end
            
            for j=1:length(datos_var{i}.crucero.tiempo)
                distancia_PLOT(j) = datos_var{i}.crucero.distancia(j);
            end
            
            % Generation of variable to plot for Cruise
            for i=1:N_V_VAR_perf
                for j=1:length(datos_var{i,N_Wp_plot}.crucero.L_D)
                    L_D_PLOT(i,j) = datos_var{i,N_Wp_plot}.crucero.L_D(j);
                    tiempo_PLOT(i,j) = datos_var{i,N_Wp_plot}.crucero.tiempo(j);
                    CL_PLOT(i,j)  = datos_var{i,N_Wp_plot}.crucero.CL(j);
                    palanca_PLOT(i,j) = datos_var{i,N_Wp_plot}.crucero.palanca(j);
                end
            end
            
            Plots_performance.m_f_PLOT = m_f_PLOT;
            Plots_performance.tiempo_total_PLOT = tiempo_total_PLOT;
            Plots_performance.distancia_total_PLOT = distancia_total_PLOT;
            Plots_performance.Cost_fuel = Cost_fuel;
            Plots_performance.DOC_PLOT = DOC_PLOT;
            Plots_performance.CAPM_PLOT = CAPM_PLOT;
            Plots_performance.distancia_PLOT = distancia_PLOT;
            Plots_performance.L_D_PLOT = L_D_PLOT;
            Plots_performance.tiempo_PLOT = tiempo_PLOT;
            Plots_performance.CL_PLOT = CL_PLOT;
            Plots_performance.palanca_PLOT = palanca_PLOT;
            Plots_performance.CI = CI;
        else
            % Dummy Variable
            Plots_performance = 0;
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CONFIGURACIÓN PUSHER. MODELO DE HELICÓPTERO%%%%%%%%%%%%%%%
Posicion_Palanca = 1;
% Will be used in future version of propulsion
conditions.alpha_f = alpha_f;
conditions.beta_f = beta_f;
conditions.h = Performance.h;
conditions.V = Performance.V;
conditions.study_var_xcg = 0;
conditions.x_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
conditions.x_XCG = 0.64500;

% Selects XCG depending if it's from desired stability conditions in Forward Flight (XCG_FF = 1) or
% AXIAL Flight XCG_FF = 0)
XCG_FF =1;
XCG_FF = OUTPUT_read_XLSX.Stability_flags.XCG_FF;

if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
    %% Dynamic Stability Formulation
    % Stability_Pamadi = 0 - Fran's Formulation
    % Stability_Pamadi = 1 - Pamadi
    % Stability_Pamadi = 2 - Roskam
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim == 1
        conditions.m_TOW = m_TOW;
        only_trim = 0;
        only_trim = OUTPUT_read_XLSX.Stability_flags.only_trim;
        StabilityModel =1;
        StabilityModel = OUTPUT_read_XLSX.Stability_flags;
%         [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts]  = Calculo_Stability_Derivatives_Abril2020...
%             (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
%             Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim);
        [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts] = Calculo_Stability_Derivatives_Noviembre2021_v2...
            (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
            Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim);
        %% XCG range
        x_XCG_fwd = TRIM_RESULTS.x_XCG_fwd;
        x_XCG_rwd = TRIM_RESULTS.x_XCG_rwd;
        N_x_XCG_VAR = 100;
        x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);
        Plot_Options.x_XCG_VAR = x_XCG_VAR;
        if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
            save Stability_Derivatives.mat TRIM_RESULTS Trim_ITER Stab_Der Stab_Der_parts
        end
            
    end
    
        StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
    Variable_Study = 0;
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Stab_Dyn_Long = longitudinal_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
            Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions);
    end
    
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Stab_Dyn_LatDir = lateral_directional_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
            Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions);

    end
       
    %% Speed range
%     V_low = Performance.V_min;
%     V_high = Performance.V_max;
    V_low = OUTPUT_read_XLSX.Stability_flags.V_low;
    V_high = OUTPUT_read_XLSX.Stability_flags.V_high;
    N_V_VAR = OUTPUT_read_XLSX.Stability_flags.N_V_VAR;
    V_VAR = linspace(V_low,V_high,N_V_VAR);
    Plot_Options.V_VAR = V_VAR;
    
    %% Weight range
    % m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_systems + Weight_tier.m_batteries/2 + Weight_tier.m_fuel/2);
    m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy);
    m_mid = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy)/2;
    m_high = Weight_tier.m_TOW;
    N_m_VAR = 20;
    N_m_VAR = OUTPUT_read_XLSX.Stability_flags.N_m_VAR;
    m_VAR = linspace(m_low,m_high,N_m_VAR);
    Plot_Options.W_VAR = m_VAR;
    
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
        % reduces the calculations since only for trim conditions
        only_trim = OUTPUT_read_XLSX.Stability_flags.only_trim;
        %% Study of variation of Velocity and mass with XCG selected for trim conditions at desired Static Margin
        for i=1:N_V_VAR
            
            conditions.x_XCG = TRIM_RESULTS.x_XCG_des;
            conditions.V = V_VAR(i);
            conditions.x_XCG_var = TRIM_RESULTS.x_XCG_des;
            for j=1:N_m_VAR
                conditions.m_TOW = m_VAR(j);
%                 [TRIM_RESULTS_calc,Trim_ITER_calc,Stab_Der_calc,Stab_Der_parts_calc] = Calculo_Stability_Derivatives_Abril2020...
%                     (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
%                     Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim);
% KIKO's
                [TRIM_RESULTS_calc,Trim_ITER_calc,Stab_Der_calc,Stab_Der_parts_calc] = Calculo_Stability_Derivatives_Noviembre2021_v2...
                    (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
                    Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim);
                TRIM_RESULTS_var_V_XCG{i,j} = TRIM_RESULTS_calc;
                Trim_ITER_var_V_XCG{i,j} = Trim_ITER_calc;
                Stab_Der_var_V_XCG{i,j} = Stab_Der_calc;
                Stab_Der_parts_V_XCG{i,j} = Stab_Der_parts_calc;
                % Variable longitudinal analysis of the poles
                % Estimation through Pamadi if Stability_Pamadi = 1;
                % Estimation through ROSKAM if Stability_Pamadi = 2;
                StabilityModel = 2;
                StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
                Variable_Study = 1;
                Stab_Dyn_Long_var_V_XCG{i,j} = longitudinal_analysis(Performance,Stab_Der_calc,conv_UNITS,StabilityModel,...
                    Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions);
                Stab_Dyn_LatDir_var_V_XCG{i,j} = lateral_directional_analysis(Performance,Stab_Der_calc,conv_UNITS,StabilityModel,...
                    Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions);
            end
        end
        
        %Saves the data to a mat file
        switch case_AC
            case 1
                save Study_var_m_varV_EMERGENTIA.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 2
                save Study_var_m_varV_EMERGENTIA_SC.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 3
                save Study_var_m_varV_PEPINOXXL.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 4
                save Study_var_m_varV_COMERCIAL.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 5
                save Study_var_m_varV_WIG.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 6
                save Study_var_m_varV_CERVERA.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 7
                save Study_var_m_varV_EXISTING_AC.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 8
                save Study_var_m_varV_TAMIZ.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 9 % EMERGENTIA Wind Tunnel w2 (No sweep)
                save Study_var_m_varV_EMERGENTIA_WT2.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 10 % EMERGENTIA Wind Tunnel w1 (No sweep)
                save Study_var_m_varV_EMERGENTIA_WT1.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 11 % EMERGENTIA Manufacturing
                save Study_var_m_varV_EMERGENTIA_MFT.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 12 % ALO
                save Study_var_m_varV_ALO.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
            case 13 % ALO FC
                save Study_var_m_varV_ALO_FC.mat TRIM_RESULTS_var_V_XCG Trim_ITER_var_V_XCG Stab_Der_var_V_XCG Stab_Der_parts_V_XCG
        end
    end
    
    
    %% Study of variation of Velocity with XCG selected for trim conditions at desired Static Margin
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_var_XCG
        conditions.study_var_xcg = 1;
        for i=1:N_x_XCG_VAR
%             x_XCG_var = TRIM_RESULTS.x_XCG_des;
            conditions.V = Performance.V;
            conditions.V = 30;
            conditions.x_XCG = x_XCG_VAR(i);
%             conditions.x_XCG_var = TRIM_RESULTS.x_XCG_des;
            conditions.m_TOW = m_TOW;
%             [TRIM_RESULTS_calc,Trim_ITER_calc,Stab_Der_calc,Stab_Der_parts_calc] = Calculo_Stability_Derivatives_Abril2020...
%                 (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
%                 Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim);
            [TRIM_RESULTS_calc,Trim_ITER_calc,Stab_Der_calc,Stab_Der_parts_calc] = Calculo_Stability_Derivatives_Noviembre2021_v2...
                (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
                Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim);

            TRIM_RESULTS_var_XCG{i} = TRIM_RESULTS_calc;
            Trim_ITER_var_XCG{i} = Trim_ITER_calc;
            Stab_Der_var_XCG{i} = Stab_Der_calc;
            Stab_Der_parts_var_XCG{i} = Stab_Der_parts_calc;
            
        end
        % Saves the data to a mat file
        switch case_AC
            case 1
                save Study_varXCG_EMERGENTIA.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 2
                save Study_varXCG_EMERGENTIA_SC.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 3
                save Study_varXCG_PEPINOXXL.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 4
                save Study_varXCG_COMERCIAL.ma  t TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 5
                save Study_varXCG_WIG.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 6
                save Study_varXCG_CERVERA.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 7
                save Study_varXCG_EXISTING_AC.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 8
                save Study_varXCG_TAMIZ.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 9 % EMERGENTIA Wind Tunnel w2 (No sweep)
                save Study_varXCG_EMERGENTIA_WT2.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 10 % EMERGENTIA Wind Tunnel w1 (No sweep)
                save Study_varXCG_EMERGENTIA_WT1.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 11 % EMERGENTIA Manufacturing
                save Study_varXCG_EMERGENTIA_MFT.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 12 % ALO
                save Study_varXCG_ALO.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
            case 13 % ALO Fuel Cell
                save Study_varXCG_ALO_FC.mat TRIM_RESULTS_var_XCG Trim_ITER_var_XCG Stab_Der_var_XCG Stab_Der_var_XCG
        end
    end
    
    % Estimation through Pamadi if Stability_Pamadi = 1;
    % Estimation through ROSKAM if Stability_Pamadi = 2;
    StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
    Variable_Study = 0;
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Stab_Dyn_Long = longitudinal_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
            Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions);
    end
    
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Stab_Dyn_LatDir = lateral_directional_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
            Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions);

    end
    
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat == 1
        conditions.study_var_xcg = 0;
        % %%%%%%%%%%%%%%%%%%%%TRIM LATERAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Cosntant beta
        beta = OUTPUT_read_XLSX.Stability_flags.beta;
        % variable beta
        beta_i = OUTPUT_read_XLSX.Stability_flags.beta_i;
        beta_f = OUTPUT_read_XLSX.Stability_flags.beta_f;
        N_Delta_beta = OUTPUT_read_XLSX.Stability_flags.N_Delta_beta;
        beta_vec = linspace(beta_i,beta_f,N_Delta_beta);
        % Storing study conditions
        conditions_TRIM_lat.beta = beta;
        conditions_TRIM_lat.beta_vec = beta_vec;
        [Trim_ITER_LAT] = Calculo_Trim_ITER_LAT_v3(Performance,Geo_tier,conv_UNITS,Trim_ITER,...
            Stab_Der,Weight_tier,conditions_TRIM_lat);
    end
    
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Turning == 1
        % %%%%%%%%%%%%%%%%%%%% Turnng conditioms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Constant phi
        phi = OUTPUT_read_XLSX.Stability_flags.phi;
        % variable phi
        phi_i = OUTPUT_read_XLSX.Stability_flags.phi_i;
        phi_f = OUTPUT_read_XLSX.Stability_flags.phi_f;
        N_Delta_phi = OUTPUT_read_XLSX.Stability_flags.N_Delta_phi;
        n_viraje = OUTPUT_read_XLSX.Stability_flags.n_viraje;
        phi_vec = linspace(phi_i,phi_f,N_Delta_phi);
        % Storing study conditions
        conditions_TRIM_turning.phi = phi;
        conditions_TRIM_turning.phi_vec = phi_vec;
        conditions_TRIM_turning.n_viraje = n_viraje;
        [Trim_ITER_LAT_Viraje] = Calculo_Trim_ITER_LAT_Viraje_v3(Performance,Geo_tier,conv_UNITS,...
            Trim_ITER,Stab_Der,Weight_tier,conditions_TRIM_turning);
    end
    
else
    Plot_Options.dummy = 1;
end


%% Plots graphics
%% Defines the plots that will be plotted, and in which order
% Flags for plotting the figures, defiens the order and the number of plots
% Case 1 PRINT_PLOTS_XFLR5 = 0; % Prints plots - Aero
% Case 2 PRINT_PLOTS_XFLR5_POLAR = 0; % Prints plots - Aero
% Case 3 PRINT_PLOTS_XAC = 0; % print Plots for Stimation ofXAC
% Case 4 PRINT_PLOTS_CT_Model = 0; % print Plots for Propulsive Models
% Case 5 PRINT_PLOTS_3D = 0; % Prints plots for 3D
% Case 6 PRINT_PLOTS_STABILITY_SM = 0; % prints plots of SM analysis
% Case 7 PRINT_PLOTS_TRIM_LON = 0; % prints plots of longitudinal Trim
% Case 8 PRINT_PLOTS_TRIM_LON_VAR = 0; % prints plots of longitudinal Trim with Variable mass & Variable Speed
% Case 9 PRINT_PLOTS_TRIM_LAT = 0; % prints plots of lateral Trim
% Case 10 PRINT_PLOTS_TURNING_LAT = 0; % prints plots of lateral Turning
% Case 11 PLOTS STABILITY DERIVATIVES FOR VAR MASS AND VELOCITY
% Case 12 PLOTS STABILITY ANALYSIS FOR VAR MASS AND VELOCITY
% Case 13 PERFORMANCE_STUDY

%% Defines options for the plots
[mark_Type COLOR_scheme, Plot_Options] = Generate_Plot_Options(prefix,mark_legend,VECTOR_XFLR5,Plot_Options);

for i=1:length(PLOTS)
    switch PLOTS(i)
        case 1 % compare = 1 elements: w1,w2, vtp... alone
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY == 1
                [Fig] = plot_XFLR5_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix);
            else
                Warning = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 2 % compare = 1 elements: w1,w2, vtp... alone
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY == 1
                [Fig] = plot_XFLR5_Polar_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix);
            else
                Warning = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 3
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY == 1
%                 [XAC,Fig] = plot_Generate_XAC(DATA_Ae,Aero,DATA_PL,CASOS,degrees_XAC,Fig,...
%                     XFLR5_DATA,Performance,X_OC,Plot_Options,VECTOR_XFLR5);
                
            else
                Warning = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
%         case 4
%             if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
%                 if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_var_XCG == 1
%                     [Fig] = Generates_Plots_PropulsionModels(Prop_data,Data_Trim,V,alpha,D,Plot_Options,Fig);
%                 else
%                     Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Trim_var_XCG == 1  - Code in PAUSE';
%                     disp(Warning)
%                     pause
%                 end
%             else
%                 Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
%                 disp(Warning)
%                 pause
%             end
            
        case 5
%             [Fig] = plot_GEOMETRY_2018(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC,Performance);
            
        case 6
            %             % Logic ensures that the studies have been conducted
            %             if STABILITY_STUDY == 1
            %                 if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim == 1
            %                     [Fig] = Generates_Plots_Long_Stability(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var,Trim_ITER_var,Geo_tier,Plot_Options,Fig);
            %                 else
            %                     Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Trim == 1  - Code in PAUSE';
            %                     disp(Warning)
            %                     pause
            %                 end
            %             else
            %                 Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
            %                 disp(Warning)
            %                 pause
            %             end
            %
        case 7
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_var_XCG == 1
                    [Fig] = Generates_Plots_Longitudinal_Trim(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var_XCG,Trim_ITER_var_XCG,Geo_tier,Plot_Options,Fig);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Trim_var_XCG == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 8
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
                    [Fig] = Generates_Plots_Longitudinal_Trim_VAR(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var_V_XCG,Trim_ITER_var_V_XCG,Stab_Der_var_V_XCG,Geo_tier,Plot_Options,Fig);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Trim_varV_m0 == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 9
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat == 1
                    [Fig] = Generates_Plots_Lateral_Trim(TRIM_RESULTS,Trim_ITER,Trim_ITER_LAT,Geo_tier,...
                        Plot_Options,conv_UNITS,conditions_TRIM_lat,Fig);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Trim_lat == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 10
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Turning == 1
                    [Fig] = Generates_Plots_Lateral_Turning(TRIM_RESULTS,Trim_ITER,Trim_ITER_LAT_Viraje,Geo_tier,...
                        Plot_Options,conv_UNITS,conditions_TRIM_turning,Fig);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Turning == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 11
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Stability_Derivatives_varV_m0 == 1
                    [Fig] = Generates_Plots_Derivatives_VAR(TRIM_RESULTS_var_XCG,Trim_ITER_var_XCG,Stab_Der_var_V_XCG,Geo_tier,Plot_Options,Fig);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Stability_Derivatives_varV_m0 == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
        case 12
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
                    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Stability_Analysis_varV_m0_long == 1
                        [Fig] = Generates_Plots_StabilityAnalysis_long_VAR(Geo_tier,Plot_Options,Stab_Dyn_Long_var_V_XCG,Fig);
                    end
                    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Stability_Analysis_varV_m0_lat == 1
                        [Fig] = Generates_Plots_StabilityAnalysis_lat_VAR(Geo_tier,Plot_Options,Stab_Dyn_LatDir_var_V_XCG,Fig);
                    end 
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Stability_Analysis_varV_m0 == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
        case 13
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var == 1
                    [Fig] = Generates_Plots_Performance_v1(Geo_tier,Plot_Options,conv_UNITS,Fig,...
                        Segments, handles,Weights_AP_var,fuel_total_var,tiempo_total_var,distancia_total_var,W_var,datos_var,case_AC,Plots_performance,Post_processing_PERFORMANCE);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that ANALYSIS_PERFORMANCE_AP_var == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that PERFORMANCE_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
    end
end

%% Checks Efficiency of code
if CHECK_Efficiency == 1
    if Detailed_Profile ==  1
        p = profile('info');
        for n = 1:size(p.FunctionHistory,2)
            if p.FunctionHistory(1,n)==0
                str = 'entering function: ';
            else
                str = ' exiting function: ';
            end
            disp([str p.FunctionTable(p.FunctionHistory(2,n)).FunctionName]);
        end
    else
        profile viewer
        profsave(profile('info'),'profile_results')
    end
end