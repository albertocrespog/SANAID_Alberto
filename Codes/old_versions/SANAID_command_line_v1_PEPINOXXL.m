%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANAID - Stability ANAlysis Interactive Design Tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
% Defines Generic Path so that files can be organized in folders
% change for the different versions
addpath(genpath('../../MATLAB/src'))
addpath(genpath('../../MATLAB/XFLR5'))
addpath(genpath('../../MATLAB/CAD'))
addpath(genpath('../../MATLAB/Prop_data'))
addpath(genpath('../../MATLAB/Fuselage'))

% Frag that determines MATLAB incompatibility (For example functions that do not work in MATLAB 2017 )
% Case MATLAB_in = 1; Modified codes for compatibility
% Case MATLAB_in = 0; Normal codes
MATLAB_in = 0;
% Units conversion
conv_UNITS = conversion_UNITS;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
g = conv_UNITS.g;

%% Check the efficiency of the coding
CHECK_Efficiency = 0;
Detailed_Profile = 0;
if CHECK_Efficiency == 1
    if Detailed_Profile
        profile on -history
    else
        profile on
    end
end

%%
% Defines de aircraft to be analyzed
% case_AC = 1 - EMERGENTIA 1:1
% case_AC = 2 - EMERGENTIA 1:2
% case_AC = 3 - PEPIÑOXXL
% case_AC = 4 - Comercial AC
% case_AC = 5 - WIG
% case_AC = 6 - CERVERA
% case_AC = 7 - Existing AC - B747
case_AC = 3;

%% Studies
% Conducts Performance Study
AERODYNAMIC_STUDY = 1;
% Propulsion Analysis
PROP_STUDY = 0;
% Performance Analysis
PERFORMANCE_STUDY = 0;
%Stability Analysis
STABILITY_STUDY = 1;

% % Conducts study of the limits of the props
% Study_Prop_Limits = 1;
% % Conducts study of sensitivity analysis of the Performance
% PERFORMANCE_sensitivity = 1;

% Weight Fraction Study
STUDY_Weight_Fraction = 0;
% Analizes poerformance integrated with Academic Performance Codes
ANALYSIS_PERFORMANCE_AP = 1;
% Analizes poerformance integrated with Academic Performance Codes varying
% conditions
ANALYSIS_PERFORMANCE_AP_var = 1;
% Determines if variable study is conducted, if different thant 1 then uses
% just one value
variable_speed_AP = 1;
variable_weight_AP = 1;

%%
% Trim conditions - one case
STABILITY_STUDY_Trim = 1;
% Trim conditions - variation V and mass
STABILITY_STUDY_Trim_varV_m0 = 1;
% Trim conditions - variation XCG
STABILITY_STUDY_Trim_var_XCG = 1;
% Longitudinal Stability Analysis
STABILITY_STUDY_Long_dyn = 1;
% Lateral-Directional Stability Analysis
STABILITY_STUDY_LatDir_dyn = 1;
% Lateral Directional Trim
STABILITY_STUDY_Trim_lat = 1;
% Turning Stability
STABILITY_STUDY_Turning = 1;

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
% Case 11 PERFORMANCE_STUDY
PLOTS = [5,6,8,9];

%% Defines options for the plots
[mark_Type COLOR_scheme, Plot_Options] = Generate_Plot_Options_init;
Plot_Options.MATLAB_in = MATLAB_in;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Geometry %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PERFORMANCE ANALYSIS %%%%%%%%%%%%%%%
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        % Scaling Factor
        SF = 1.0;
        % Defines the flag to determine method to Estimate weights
        % CASE1 = Weight_Estimation from composite - Cefiro III
        % CASE2 = Weight_Estimation from wood - Cefiro I
        % CASE3 = Factores lineales
        Weight_Estimation = 1;
        %% Estimation of prop diameter, just preliminary for geometric
        % conditions
        D_prop=28*2.54/100;
        %% Aircraft type
        % AC_type = 1 - flying wing
        % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        % AC_type = 4 - 2 surface: wing + V-tail
        % AC_type = 5 - 3 surface: cannard + wing + V-tail
        AC_type = 4;
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_EMERGENTIA(SF,conv_UNITS,AC_type,case_AC);
        Geo_input_tier = Generation_Input_Geometric_Data_EMERGENTIA(conv_UNITS,AC_CONFIGURATION,SF);
        % Determines the fuselage than will be shown
        % CASE_fuse = 1 - Just STL of fuselage
        % CASE_fuse = 2 - STL of the entire aircraft
        % CASE_fuse = 3 - Fuselage from XFLR5
        CASE_fuse = 3;
        % Determine the type of scaling
        % ESCALADO = 0; % no scaling Ratio =1
        % ESCALADO = 1; % scaling acording to desired for each 3 axix
        % ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
        ESCALADO = 0;
        [XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_EMERGENTIA;
        Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data,case_AC);
        
        % type_battery used
        % case 1 LiFePO4
        % case 2 LiPo
        % case 3 FuelCells
        type_battery = 2;
        AC_CONFIGURATION.type_battery = type_battery;
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        % Scaling Factor
        SF = 0.50;
        % Defines the flag to determine method to Estimate weights
        % CASE1 = Weight_Estimation from composite - Cefiro III
        % CASE2 = Weight_Estimation from wood - Cefiro I
        Weight_Estimation = 1;
        %% Estimation of prop diameter, just preliminary for geometric
        D_prop=22*2.54/100;
        %% Aircraft type
        % AC_type = 1 - flying wing
        % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        % AC_type = 4 - 2 surface: wing + V-tail
        % AC_type = 5 - 3 surface: cannard + wing + V-tail
        AC_type = 4;
        [XCG_data, AC_CONFIGURATION] = Generation_AC_configuration_EMERGENTIA(SF,conv_UNITS,AC_type,case_AC);
        Geo_input_tier = Generation_Input_Geometric_Data_EMERGENTIA(conv_UNITS,AC_CONFIGURATION,SF);
        % Determines the fuselage than will be shown
        % CASE_fuse = 1 - Just STL of fuselage
        % CASE_fuse = 2 - STL of the entire aircraft
        % CASE_fuse = 3 - Fuselage from XFLR5
        CASE_fuse = 1;
        % Determine the type of scaling
        % ESCALADO = 0; % no scaling Ratio =1
        % ESCALADO = 1; % scaling acording to desired for each 3 axix
        % ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
        ESCALADO = 2;
        [XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_EMERGENTIA;
        Geo_tier = Generation_Geometric_Data_v2(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data,case_AC);
        
        % type_battery used
        % case 1 LiFePO4
        % case 2 LiPo
        % case 3 FuelCells
        type_battery = 2;
        AC_CONFIGURATION.type_battery = type_battery;
    case 3 % case_AC = 3 - PEPIÑO XXL
        % Scaling Factor
        SF = 1;
        % Defines the flag to determine method to Estimate weights
        % CASE1 = Weight_Estimation from composite - Cefiro III
        % CASE2 = Weight_Estimation from wood - Cefiro I
        % CASE3 = Factores lineales
        Weight_Estimation = 1;
        %% Estimation of prop diameter, just preliminary for geometric
        % conditions
        D_prop=32*2.54/100;
        %% Aircraft type
        % AC_type = 1 - flying wing
        % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        % AC_type = 4 - 2 surface: wing + V-tail
        % AC_type = 5 - 3 surface: cannard + wing + V-tail
        AC_type = 4;
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_PEPINOXXL(SF,conv_UNITS,AC_type,case_AC);
        Geo_input_tier = Generation_Input_Geometric_Data_PEPINOXXL(conv_UNITS,AC_CONFIGURATION,SF);
        % Determines the fuselage than will be shown
        % CASE_fuse = 1 - Just STL of fuselage
        % CASE_fuse = 2 - STL of the entire aircraft
        % CASE_fuse = 3 - Fuselage from XFLR5
        CASE_fuse = 3;
        % Determine the type of scaling
        % ESCALADO = 0; % no scaling Ratio =1
        % ESCALADO = 1; % scaling acording to desired for each 3 axix
        % ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
        ESCALADO = 0;
        [XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_PEPINOXXL;
        Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data,case_AC);
        
        % type_battery used
        % case 1 LiFePO4
        % case 2 LiPo
        % case 3 FuelCells
        type_battery = 2;
        AC_CONFIGURATION.type_battery = type_battery;
    case 4 % Commercial example
        SF = sqrt(30);
        % Defines the flag to determine method to Estimate weights
        % CASE1 = Weight_Estimation from composite - Cefiro III
        % CASE2 = Weight_Estimation from wood - Cefiro I
        % CASE3 = Factores lineales
        Weight_Estimation = 3;
        %% Estimation of prop diameter, just preliminary for geometric
        % conditions
        D_prop=60*2.54/100;
        %% Aircraft type
        % AC_type = 1 - flying wing
        % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        % AC_type = 4 - 2 surface: wing + V-tail
        % AC_type = 5 - 3 surface: cannard + wing + V-tail
        AC_type = 2;
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_Comercial(SF,conv_UNITS,AC_type,case_AC);
        Geo_input_tier = Generation_Input_Geometric_Data_EMERGENTIA(conv_UNITS,AC_CONFIGURATION,SF);
        %--------------------------- CALCULATES GEOMETRY DATA ---------------------------------
        % Determines the fuselage than will be shown
        % CASE_fuse = 1 - Just STL of fuselage
        % CASE_fuse = 2 - STL of the entire aircraft
        % CASE_fuse = 3 - Fuselage from XFLR5
        CASE_fuse = 1;
        % Determine the type of scaling
        % ESCALADO = 0; % no scaling Ratio =1
        % ESCALADO = 1; % scaling acording to desired for each 3 axix
        % ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
        ESCALADO = 1;
        [XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_Comercial;
        Geo_tier = Generation_Geometric_Data_v2(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data,case_AC);
        
        % type_battery used
        % case 1 LiFePO4
        % case 2 LiPo
        % case 3 FuelCells
        type_battery = 2;
        AC_CONFIGURATION.type_battery = type_battery;
    case 5 % WIG Aircraft
        SF = 0.42; % Para dimensiones de
        % Defines the flag to determine method to Estimate weights
        % CASE1 = Weight_Estimation from composite - Cefiro III
        % CASE2 = Weight_Estimation from wood - Cefiro I
        % CASE3 = Factores lineales
        Weight_Estimation = 3;
        %% Estimation of prop diameter, just preliminary for geometric
        % conditions
        D_prop=60*2.54/100;
        %% Aircraft type
        % AC_type = 1 - flying wing
        % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        % AC_type = 4 - 2 surface: wing + V-tail
        % AC_type = 5 - 3 surface: cannard + wing + V-tail
        AC_type = 2;
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_WIG(SF,conv_UNITS,AC_type,case_AC);
        Geo_input_tier = Generation_Input_Geometric_Data_WIG(conv_UNITS,AC_CONFIGURATION,SF);
        %--------------------------- CALCULATES GEOMETRY DATA ---------------------------------
        % Determines the fuselage than will be shown
        % CASE_fuse = 1 - Just STL of fuselage
        % CASE_fuse = 2 - STL of the entire aircraft
        % CASE_fuse = 3 - Fuselage from XFLR5
        CASE_fuse = 2;
        % Determine the type of scaling
        % ESCALADO = 0; % no scaling Ratio =1
        % ESCALADO = 1; % scaling acording to desired for each 3 axix
        % ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
        ESCALADO = 1;
        STL_PLOT = 1; % Generates Fuselage mesh from CAD STL is flag STL_PLOT = 1;
        XFLR5_DATA.STL_PLOT = STL_PLOT;
        [XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_WIG(XFLR5_DATA);
        Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data,case_AC);
        
        % type_battery used
        % case 1 LiFePO4
        % case 2 LiPo
        % case 3 FuelCells
        type_battery = 2;
        AC_CONFIGURATION.type_battery = type_battery;
    case 6 % case_AC = 6 CERVERA
        % Scaling Factor
        SF = 1.0;
        % Defines the flag to determine method to Estimate weights
        % CASE1 = Weight_Estimation from composite - Cefiro III
        % CASE2 = Weight_Estimation from wood - Cefiro I
        % CASE3 = Factores lineales
        Weight_Estimation = 2;
        %% Estimation of prop diameter, just preliminary for geometric
        % conditions
        D_prop=28*2.54/100;
        %% Aircraft type
        % AC_type = 1 - flying wing
        % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        % AC_type = 4 - 2 surface: wing + V-tail
        % AC_type = 5 - 3 surface: cannard + wing + V-tail
        AC_type = 4;
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_PEPINOXXL(SF,conv_UNITS,AC_type,case_AC);
        Geo_input_tier = Generation_Input_Geometric_Data_PEPINCOXXL(conv_UNITS,AC_CONFIGURATION,SF);
        % Determines the fuselage than will be shown
        % CASE_fuse = 1 - Just STL of fuselage
        % CASE_fuse = 2 - STL of the entire aircraft
        % CASE_fuse = 3 - Fuselage from XFLR5
        CASE_fuse = 3;
        % Determine the type of scaling
        % ESCALADO = 0; % no scaling Ratio =1
        % ESCALADO = 1; % scaling acording to desired for each 3 axix
        % ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
        ESCALADO = 0;
        [XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_PEPINOXXL;
        Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data,case_AC);
        
        % type_battery used
        % case 1 LiFePO4
        % case 2 LiPo
        % case 3 FuelCells
        type_battery = 3;
        AC_CONFIGURATION.type_battery = type_battery;
end

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
[Body_Geo,meshData] = Generation_Fuselage_Data(Geo_tier,XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file,SF); % Defines Propulsion DATA

%% Defines Estimation of Weights according to Cefiro III densities
f_f = 1.5; % fudge factor to increment the estimation associated to linear fraction estimations
Weight_tier = Generation_Weight_Data(Geo_tier,Body_Geo,AC_CONFIGURATION,conv_UNITS,Weight_Estimation,SF,f_f);
m_TOW = Weight_tier.m_TOW;

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
end

% Valor deseado del Punto Neutro
SM_des = 0.25;

%% Generates the file for Aerodynamic Data
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        %% Selection of TXT that are used for the aerodynamic analysis
        % Select txt files associated with each aerodynamic study
        % Make sure to check with case asssigned in read_aero_files_XXX.m
        index_w1 = 14; % Front Wing W1 @ 0.045 m from the wing root chord LE
        index_w2 = 23; % Rear Wing W2
        % Allows the user to select the min AoA used to determine Curve Lift Slope
        i_w1 = 4;
        i_w2 = -3;
        
        %% Defines elements to analyze Polar Composite Build-Up-Method
        Conf.w1 = AC_CONFIGURATION.W1; % wing
        Conf.h = AC_CONFIGURATION.HTP; % HTP
        if AC_CONFIGURATION.VTP == 1
            if AC_CONFIGURATION.twin_VTP == 1
                Conf.v = 0;% VTP
                Conf.v2 = 1;% double vertical
            else
                Conf.v = 1;% VTP
                Conf.v2 = 0;% double vertical
            end
        else
            Conf.v = 0;% VTP
            Conf.v2 = 0;% double vertical
        end
        Conf.vtail = AC_CONFIGURATION.Vee;% V-tail
        Conf.can   = AC_CONFIGURATION.Can;% Canard
        Conf.fus = 1;% fuselage
        Conf.m_fus = 0;% multiple fuselage
        if Conf.m_fus == 1
            Conf.n_m_fus = 0;% number of multiple
        else
            Conf.n_m_fus = 0;% number of multiple
        end
        Conf.nac = 1;% nacelles
        Conf.landgear = 0;
        
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        %% Selection of TXT that are used for the aerodynamic analysis
        % Select txt files associated with each aerodynamic study
        % Make sure to check with case asssigned in read_aero_files_XXX.m
        index_w1 = 14; % Front Wing W1 @ 0.045 m from the wing root chord LE
        index_w2 = 23; % Rear Wing W2
        % Allows the user to select the min AoA used to determine Curve Lift Slope
        i_w1 = 4;
        i_w2 = -6;
        
        %% Defines elements to analyze Polar Composite Build-Up-Method
        Conf.w1 = AC_CONFIGURATION.W1; % wing
        Conf.h = AC_CONFIGURATION.HTP; % HTP
        if AC_CONFIGURATION.VTP == 1
            if AC_CONFIGURATION.twin_VTP == 1
                Conf.v = 0;% VTP
                Conf.v2 = 1;% double vertical
            else
                Conf.v = 1;% VTP
                Conf.v2 = 0;% double vertical
            end
        end
        Conf.vtail = AC_CONFIGURATION.Vee;% V-tail
        Conf.can   = AC_CONFIGURATION.Can;% Canard
        Conf.fus = 1;% fuselage
        Conf.m_fus = 0;% multiple fuselage
        if Conf.m_fus == 1
            Conf.n_m_fus = 0;% number of multiple
        else
            Conf.n_m_fus = 0;% number of multiple
        end
        Conf.nac = 1;% nacelles
        Conf.landgear = 0;
        
    case 3 % case_AC = 3 - PEPIÑO XXL
        %% Selection of TXT that are used for the aerodynamic analysis
        % Select txt files associated with each aerodynamic study
        % Make sure to check with case asssigned in read_aero_files_XXX.m
        index_w1 = 1; % Front Wing W1 @ 0.045 m from the wing root chord LE
        index_w2 = 2; % Rear Wing W2
        % Allows the user to select the min AoA used to determine Curve Lift Slope
        i_w1 = 0;
        i_w2 = -1;
        
        %% Defines elements to analyze Polar Composite Build-Up-Method
        Conf.w1 = AC_CONFIGURATION.W1; % wing
        Conf.h = AC_CONFIGURATION.HTP; % HTP
        if AC_CONFIGURATION.VTP == 1
            if AC_CONFIGURATION.twin_VTP == 1
                Conf.v = 0;% VTP
                Conf.v2 = 1;% double vertical
            else
                Conf.v = 1;% VTP
                Conf.v2 = 0;% double vertical
            end
        else
            Conf.v = 0;% VTP
            Conf.v2 = 0;% double vertical
        end
        Conf.vtail = AC_CONFIGURATION.Vee;% V-tail
        Conf.can   = AC_CONFIGURATION.Can;% Canard
        Conf.fus = 1;% fuselage
        Conf.m_fus = 0;% multiple fuselage
        if Conf.m_fus == 1
            Conf.n_m_fus = 0;% number of multiple
        else
            Conf.n_m_fus = 0;% number of multiple
        end
        Conf.nac = 1;% nacelles
        Conf.landgear = 0;
        
    case 4 % case_AC = 4 - COMERCIAL
        %% Selection of TXT that are used for the aerodynamic analysis
        % Select txt files associated with each aerodynamic study
        % Make sure to check with case asssigned in read_aero_files_XXX.m
        index_w1 = 14; % Front Wing W1 @ 0.045 m from the wing root chord LE
        index_w2 = 23; % Rear Wing W2
        % Allows the user to select the min AoA used to determine Curve Lift Slope
        i_w1 = 4;
        i_w2 = -3;
        
        %% Defines elements to analyze Polar Composite Build-Up-Method
        Conf.w1 = AC_CONFIGURATION.W1; % wing
        Conf.h = AC_CONFIGURATION.HTP; % HTP
        if AC_CONFIGURATION.VTP == 1
            if AC_CONFIGURATION.twin_VTP == 1
                Conf.v = 0;% VTP
                Conf.v2 = 1;% double vertical
            else
                Conf.v = 1;% VTP
                Conf.v2 = 0;% double vertical
            end
        end
        Conf.vtail = AC_CONFIGURATION.Vee;% V-tail
        Conf.can   = AC_CONFIGURATION.Can;% Canard
        Conf.fus = 1;% fuselage
        Conf.m_fus = 0;% multiple fuselage
        if Conf.m_fus == 1
            Conf.n_m_fus = 0;% number of multiple
        else
            Conf.n_m_fus = 0;% number of multiple
        end
        Conf.nac = 1;% nacelles
        Conf.landgear = 0;
        
    case 5 % case_AC = 5 - WIG
        %% Selection of TXT that are used for the aerodynamic analysis
        % Select txt files associated with each aerodynamic study
        % Make sure to check with case asssigned in read_aero_files_XXX.m
        index_w1 = 1; % Front Wing W1 @ 0.045 m from the wing root chord LE
        index_w2 = 2; % Rear Wing W2
        index_VTP = 3; % Rear Wing W2
        % Allows the user to select the min AoA used to determine Curve Lift Slope
        i_w1 = 1;
        i_w2 = -1;
        
        %% Defines elements to analyze Polar Composite Build-Up-Method
        Conf.w1 = AC_CONFIGURATION.W1; % wing
        Conf.h = AC_CONFIGURATION.HTP; % HTP
        if AC_CONFIGURATION.VTP == 1
            if AC_CONFIGURATION.twin_VTP == 1
                Conf.v = 0;% VTP
                Conf.v2 = 1;% double vertical
            else
                Conf.v = 1;% VTP
                Conf.v2 = 0;% double vertical
            end
        end
        Conf.vtail = AC_CONFIGURATION.Vee;% V-tail
        Conf.can   = AC_CONFIGURATION.Can;% Canard
        Conf.fus = 1;% fuselage
        Conf.m_fus = 0;% multiple fuselage
        Conf.m_fus = 1;% multiple fuselage
        if Conf.m_fus == 1
            Conf.n_m_fus = 0;% number of multiple
        else
            Conf.n_m_fus = 0;% number of multiple
        end
        Conf.nac = 1;% nacelles
        Conf.landgear = 0;
    case 6 % case_AC = 6 CERVERA
        %% Selection of TXT that are used for the aerodynamic analysis
        % Select txt files associated with each aerodynamic study
        % Make sure to check with case asssigned in read_aero_files_XXX.m
        index_w1 = 14; % Front Wing W1 @ 0.045 m from the wing root chord LE
        index_w2 = 23; % Rear Wing W2
        % Allows the user to select the min AoA used to determine Curve Lift Slope
        i_w1 = 4;
        i_w2 = -3;
        
        %% Defines elements to analyze Polar Composite Build-Up-Method
        Conf.w1 = AC_CONFIGURATION.W1; % wing
        Conf.h = AC_CONFIGURATION.HTP; % HTP
        if AC_CONFIGURATION.VTP == 1
            if AC_CONFIGURATION.twin_VTP == 1
                Conf.v = 0;% VTP
                Conf.v2 = 1;% double vertical
            else
                Conf.v = 1;% VTP
                Conf.v2 = 0;% double vertical
            end
        else
            Conf.v = 0;% VTP
            Conf.v2 = 0;% double vertical
        end
        Conf.vtail = AC_CONFIGURATION.Vee;% V-tail
        Conf.can   = AC_CONFIGURATION.Can;% Canard
        Conf.fus = 1;% fuselage
        Conf.m_fus = 0;% multiple fuselage
        if Conf.m_fus == 1
            Conf.n_m_fus = 0;% number of multiple
        else
            Conf.n_m_fus = 0;% number of multiple
        end
        Conf.nac = 1;% nacelles
        Conf.landgear = 0;
end

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

% Flight Safe Margin to calculate aerodynamic properties
Flight_SF = 1.2; % Stall Safe Margin Conditions
Design_criteria.Flight_SF = Flight_SF;
%% Propulsion Generation
% D = Prop_data.D_prop;
alpha_f = 0*D2R;
beta_f = 0*D2R;

%% Initial Estimate of Design altitude and velocity
%% Needs to introduce preliminary results for initial estimate

%% Generates the file for Aerodynamic Data
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        % VTOL Flight
        h_climb = 100; % (m)
        Endurance_v = 0.5; % (min)
        % Cruise Flight
        h_init = 200; % (m)
        V_init = 30; % speed
        V_max_init = 1.25*V_init;
        Range = 22*1000; % m
        Endurance = 60; % (min)
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        % VTOL Flight
        h_climb = 100; % (m)
        Endurance_v = 0.5; % (min)
        % Cruise Flight
        h_init = 200; % (m)
        V_init = 30; % speed
        V_max_init = 1.25*V_init;
        Range = 22*1000; % m
        Endurance = 60; % (min)
    case 3 % case_AC = 3 - PEPIÑO XXL
        Endurance_v = 0.0; % (min)
        % Cruise Flight
        h_init = 3000; % (m)
        h_climb = 3000; % (m)
        V_init = 80; % speed
        V_max_init = 1.25*V_init;
        Range = 50*1000; % m
        Endurance = 4*60; % (min)
        h_climb = 3000; % (m)
    case 4 % case_AC = 4 - COMERCIAL
        h_init = 3000; % (m)
        V_init = 90; % speed
        V_max_init = 1.25*V_init;
        Range = 22*1000; % m
        h_climb = 10000; % (m)
    case 5 % case_AC = 5 - WIG
        Endurance_v = 0; % (min)
        % Cruise Flight
        Endurance = 4*60; % (min)
        h_init = 5; % (m)
        V_init = 70; % speed
        V_max_init = 1.25*V_init;
        Range = 2000*1000; % m
        h_climb = 5; % (m)
    case 6 % case_AC = 6 - CERVERA
        % VTOL Flight
        h_climb = 100; % (m)
        Endurance_v = 20; % (min)
        % Cruise Flight
        h_init = 100; % (m)
        V_init = 30; % speed
        V_max_init = 1.25*V_init;
        Range = 22*1000; % m
        Endurance = 60; % (min)
end

% Indentifies preliminary performance results
Performance_preliminar.h_climb = h_climb;
Performance_preliminar.Endurance_v = Endurance_v;
Performance_preliminar.h = h_init;
Performance_preliminar.V = V_init;
Performance_preliminar.V_max = V_max_init;
Performance_preliminar.Range = Range;
Performance_preliminar.Endurance = Endurance;

%v Atmospheric conditions
[Temp_init,rho_init,p_init,a_init]=atmos_inter_mio(h_init);
Performance_preliminar.Temp = Temp_init;
Performance_preliminar.rho = rho_init;
Performance_preliminar.p = p_init;
Performance_preliminar.a = a_init;
Mach_init = V_init/a_init;
Performance_preliminar.Mach = Mach_init;
q_inf_init = 0.5*rho_init*V_init^2;
Performance_preliminar.q_inf = q_inf_init;

% Defines the Prefix for plotting options
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        prefix = strcat('CERVERA_');
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        prefix = strcat('CERVERA_');
    case 3 % case_AC = 3 - PEPIÑO XXL
        prefix = strcat('PEPINO_XXXL_');;
    case 4 % case_AC = 4 - COMERCIAL
        prefix = strcat('COMERCIAL_');
    case 5 % case_AC = 5 - WIGL
        prefix = strcat('WIG_');;
    case 6 % case_AC = 6 - CERVERA
        prefix = strcat('CERVERA_');
end

if AERODYNAMIC_STUDY == 1
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
    end
    %% Read Aerodynamic Data
    [DATA_Ae] = Read_data_Aero(casos);
    
    %% Fusion Aerodynamic Propert1ies
    [Aero_TH] = Polar_Literatura_v1(Performance_preliminar,Geo_tier,conv_UNITS,Body_Geo,AC_CONFIGURATION,Conf);
    
    %% Initializes figures
    Fig = 0;
    
    % Premiliminary Aerodynamic Design
    [Aero,DATA_PL,Performance] = Generate_Aero_Data_2020_v1(DATA_Ae,Design_criteria,Performance_preliminar,Geo_tier,Weight_tier,conv_UNITS,AC_CONFIGURATION);
    % [Aero,DATA_PL,Fig] = Aerodynamic_Design_2020_v1(Geo_tier,...
    %     Weight_tier,conv_UNITS,Design_criteria,Performance,CASOS,degrees_XAC,Fig,XFLR5_DATA,...
    %     MAC_Estimation,DATA_Ae,X_OC,Plot_Options,VECTOR_XFLR5);
end

%% Data for analysis of Prop DATA
% Determines the plots that want to show
% Compares Propeller properties for different sources
% The VECTOR in order to analyze the results of Props consist of 2 grous of
% vector each prop model
% compare_prop = 1 compares APC Models alone
% compare_prop = 2 compares APC with wind tunnel models
compare_props = 1;
VECTOR_Prop.compare_props = compare_props;

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
end

%% Generates the file for Prop Data to be used in the codes
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        %% Selection of TXT that are used for the prop analysis
        % Select txt files associated with each prop study
        % Make sure to check with case asssigned in read_prop_files_XXX.m
        index_prop = 1;
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        %% Selection of TXT that are used for the prop analysis
        % Select txt files associated with each prop study
        % Make sure to check with case asssigned in read_prop_files_XXX.m
        index_prop = 1;
    case 3 % case_AC = 3 - PEPIÑO XXL
        %% Selection of TXT that are used for the prop analysis
        % Select txt files associated with each prop study
        % Make sure to check with case asssigned in read_prop_files_XXX.m
        index_prop = 1;
    case 4 % case_AC = 4 - COMERCIAL
        %% Selection of TXT that are used for the prop analysis
        % Select txt files associated with each prop study
        % Make sure to check with case asssigned in read_prop_files_XXX.m
        index_prop = 1;
    case 5 % case_AC = 5 - WIG
        %% Selection of TXT that are used for the prop analysis
        % Select txt files associated with each prop study
        % Make sure to check with case asssigned in read_prop_files_XXX.m
        index_prop = 1;
    case 6 % case_AC = 6 CERVERA
        %% Selection of TXT that are used for the prop analysis
        % Select txt files associated with each prop study
        % Make sure to check with case asssigned in read_prop_files_XXX.m
        index_prop = 1;
end

% actualizes the Prop geometry
% Select:
% Propulsion model: SEE read_prop_files_May2020.m FOR DETAILS OF PROP MODEL!!
Flag.APC_model = 1;
Flag.Wind_tunnel_model = 1;
Flag.compare_prop_models = 1;
% [casos_prop_APC prefix_prop mark_legend_prop]= read_prop_files_May2020;

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
Prop_selection.model = 2;
% Selects the prop for
% - Model 1 - APC data
% - Model 2 - Wind tunnel data for different props
prop_selec_APC = 1;
prop_selec_WT1 = 4;
% Stores info
Prop_selection.prop_selec_APC = prop_selec_APC;
Prop_selection.prop_selec_WT1 = prop_selec_WT1;
[Prop_data] = Generation_Propulsion_Data(SF,D_prop,AC_CONFIGURATION,Prop_selection);
% Defines Propulsion DATA
% [Fig] = plot_prop_APC(Data_P,Plot_Options,Fig,prefix,VECTOR_Prop)

% actualizes the Prop geometry
% Prop_data = Generation_Propulsion_Data(SF,D_prop,AC_CONFIGURATION); % Defines Propulsion DATA

if PROP_STUDY== 1
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


if PERFORMANCE_STUDY == 1
    %% Enter the number of mission segments
    % - 1 Taxy
    % - 2 TakeOff
    % - 3 = Climb
    % - 4 = VTOL Climb
    % - 5 = Cruise
    % - 6 = Load Deployment
    % - 7 = Turn
    % - 8 = Descent
    % - 9 = Descent (VTOL)
    % - 10 = Alternative Airport climb to 3000ft
    % - 11 = Turn loitter 45 min
    % - 12 = Landing
    % - 13 = Dummy to complete the 3 segment requirement for AP
    
    %% Climb options
    % case 2 % 'Subida dados M y gamma';
    % case 3 % 'Subida dados EAS y gamma';
    % case 4 % 'Subida dados TAS y gamma';
    % case 5 % 'Subida dados M y palanca';
    % case 6 % 'Subida dados EAS y palanca';
    % case 7 % 'Subida dados TAS y palanca';
    % case 8 % 'Subida dados V inicial,final y gamma';
    % case 9 % 'Subida steppest climb';
    % case 10 % 'Subida fastest climb';
    % case 11 % 'Subida dados V inicial,final y palanca'
    %% Cruise options
    % case 2 % 'Crucero dado M y distancia'
    % case 3 % ;'Crucero dado CL y distancia'
    % case 4 % ;'Crucero dados V inicial,final y palanca'
    % case 5 % 'Crucero con polar en funcion de M';
    % case 6 % 'Crucero de max alcance dado peso final'
    % case 7 % 'Crucero de max autonomia dado peso final'
    %% Turn options
    % case 2 % 'Viraje horizontal dado V y palanca'
    % case 3 % 'Viraje horizontal dado V y CL'
    % case 4 % 'Viraje horizontal dado V y balance'
    % case 5 % 'Viraje horizontal dado V y n'
    % case 6 % 'V.H dado V y radio de giro '
    % case 7 % 'V.H dado V y velocidad de guiñada';...
    % case 8 % 'V.H dado palanca y a factor de carga max'
    % case 9 % 'V.H dado palanca y a v de guiñada max'
    % case 10 % 'V.H dado palanca y a radio min'
    %% Descent options
    % case 2 % 'Descenso dados M y gamma'
    % case 3 % 'Descenso dados EAS y gamma';
    % case 4 % 'Descenso dados TAS y gamma'
    % case 5 % 'Descenso dados M y palanca'
    % case 6 % 'Descenso dados EAS y palanca'
    % case 7 % 'Descenso dados TAS y palanca'
    % case 8 % 'Descenso dados V inicial,final y gamma'
    % case 9 % 'Descenso a minimo gamma'
    % case 10 % 'Descenso "slowest sink"'
    % case 11 % 'Descenso dados V inicial,final y palanca'
    %% Generates the file for Aerodynamic Data
    switch case_AC
        case 1 % case_AC = 1 - EMERGENTIA 1:1
            % Enter type of mission segments between brackets being
            type_missions_WF = [5]; % mission for the weight fraction method
            num_missions_WF = length(type_missions_WF);
            %% User needs to define the inputs for each segment
            FlightMODE_WF.climb_mode = 2;
            FlightMODE_WF.turn_mode = 5;
            FlightMODE_WF.descent_mode = 2;
            FlightMODE_WF.cruise_mode = 2;
            Segments = Define_Segments(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF);
        case 2 % case_AC = 2 - EMERGENTIA 1:2
            % Enter type of mission segments between brackets being
            type_missions_WF = [5]; % mission for the weight fraction method
            num_missions_WF = length(type_missions_WF);
            %% User needs to define the inputs for each segment
            FlightMODE_WF.climb_mode = 2;
            FlightMODE_WF.turn_mode = 5;
            FlightMODE_WF.descent_mode = 2;
            FlightMODE_WF.cruise_mode = 2;
            Segments = Define_Segments(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF);
        case 3 % case_AC = 3 - PEPIÑO XXL
            % Enter type of mission segments between brackets being
            type_missions_WF = [5]; % mission for the weight fraction method
            num_missions_WF = length(type_missions_WF);
            %% User needs to define the inputs for each segment
            FlightMODE_WF.climb_mode = 2;
            FlightMODE_WF.turn_mode = 5;
            FlightMODE_WF.descent_mode = 2;
            FlightMODE_WF.cruise_mode = 2;
            Segments = Define_Segments(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF);
        case 4 % case_AC = 4 - COMERCIAL
            % Enter type of mission segments between brackets being
            type_missions_WF = [5]; % mission for the weight fraction method
            num_missions_WF = length(type_missions_WF);
            %% User needs to define the inputs for each segment
            FlightMODE_WF.climb_mode = 2;
            FlightMODE_WF.turn_mode = 5;
            FlightMODE_WF.descent_mode = 2;
            FlightMODE_WF.cruise_mode = 2;
            Segments = Define_Segments(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF);
        case 5 % case_AC = 5 - WIG
            % Enter type of mission segments between brackets being
            type_missions_WF = [5]; % mission for the weight fraction method
            num_missions_WF = length(type_missions_WF);
            %% User needs to define the inputs for each segment
            FlightMODE_WF.climb_mode = 2;
            FlightMODE_WF.turn_mode = 5;
            FlightMODE_WF.descent_mode = 2;
            FlightMODE_WF.cruise_mode = 2;
            Segments = Define_Segments(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF);
        case 6 % case_AC = 6 - CERVERA
            % Enter type of mission segments between brackets being
            type_missions_WF = [5]; % mission for the weight fraction method
            num_missions_WF = length(type_missions_WF);
            %% User needs to define the inputs for each segment
            FlightMODE_WF.climb_mode = 2;
            FlightMODE_WF.turn_mode = 5;
            FlightMODE_WF.descent_mode = 2;
            FlightMODE_WF.cruise_mode = 2;
            Segments = Define_Segments(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF);
    end
    
    %     Weight_tier
    if STUDY_Weight_Fraction == 1
        [seg_WF] = Generation_Mission_Segments_v3(conv_UNITS,num_missions_WF,type_missions_WF,Performance,case_AC,Segments,FlightMODE_WF);
        Weights = Weight_Fraction(seg_WF,Propulsion,conv_UNITS,Aero_TH,type_missions_WF,Weight_tier,Geo_tier,AC_CONFIGURATION);
    end
    %     % Performance Analysis
    % %     type_missions = [1 2 3 5 7 8 10 11 12];
    %     type_missions = [5 13 13];
    %     num_missions = length(type_missions);
    %     %% Define the modes of mission according to thenumber of segments, making all segments eqaul in mode of flight
    %     FlightMODE.climb_mode = 2*ones(1,num_missions);
    %     FlightMODE.turn_mode = 5*ones(1,num_missions);
    %     FlightMODE.descent_mode = 2*ones(1,num_missions);
    %     FlightMODE.cruise_mode = 3*ones(1,num_missions);
    %     [seg] = Generation_Mission_Segments_v3(conv_UNITS,num_missions,type_missions,Performance,case_AC,Segments,FlightMODE);
    %
    %     %% Performance
    %     if ANALYSIS_PERFORMANCE_AP ==1
    %         output_caracteristics = caracteristicas_avanz_v1(Aero_TH,Prop_data,Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Aero);
    %         [Weights_AP,fuel_total,tiempo_total,distancia_total,W,datos] = procesar_mision_v2(seg,num_missions,output_caracteristics,Segments);
    %     end
    
    %% Performance
    if ANALYSIS_PERFORMANCE_AP_var ==1
        
        %% Generates the file for Aerodynamic Data
        switch case_AC
            %% PROPULSION DATA
            % 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON
            % 2: NUMERO DE MOTORES
            % 3: EMPUJE/POTENCIA A NIVEL DEL MAR
            % 4: CONSUMO ESPECIFICO
            % 5: AVION CIVIL =1/MILITAR = 2
            % 6: EFICIENCIA DE LA HELICE (ETA_P)
            % 7: DERIVACION(TURBOFANES)
            
            case 1 % case_AC = 1 - EMERGENTIA 1:1
                % Enter type of mission segments betwee brackets being
                type_missions = [5 13 13]; % mission for the weight fraction method
                num_missions = length(type_missions);
                %% User needs to define the inputs for each segment
                FlightMODE_var.climb_mode = [2 2 2];
                FlightMODE_var.turn_mode = [5 5 5];
                FlightMODE_var.descent_mode = [2 2 2];
                FlightMODE_var.cruise_mode = [2 2 2];
                
                prop_data(1) = 2; % Type of engine
                prop_data(2) = AC_CONFIGURATION.n_eng; % Number of engines
                prop_data(3) = 30; % Thrust (lbf) or Power (shp) per engine
                prop_data(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
                prop_data(5) = 2; % Normativa
                prop_data(6) = 0.82; % Prop efficiency
                prop_data(7) = 1; % By-pass
                Prop_sel = prop_data;
                
                % Range of low and high speed to be analized
                V_low = 250/3.6; % low value in m/s
                V_high = 450/3.6; % high value in m/s
                N_V_VAR_perf = 5; % number of iterations
                % Value of single speed
                V_single = 250/3.6; % singe value in m/s
                
                % Range of low and high payloads to be analized
                Wp_low = 10*1000; % low value in kg
                Wp_high = 30*1000; % high value in kg
                N_Wp_VAR_perf = 2;
                % Value of single payload
                W_single = 23*1000; % single value in kg
                Post_processing_PERFORMANCE = 0;

            case 2 % case_AC = 2 - EMERGENTIA 1:2
                % Enter type of mission segments betwee brackets being
                type_missions = [5 13 13]; % mission for the weight fraction method
                num_missions = length(type_missions);
                %% User needs to define the inputs for each segment
                FlightMODE_var.climb_mode = [2 2 2];
                FlightMODE_var.turn_mode = [5 5 5];
                FlightMODE_var.descent_mode = [2 2 2];
                FlightMODE_var.cruise_mode = [2 2 2];
                
                prop_data(1) = 2; % Type of engine
                prop_data(2) = AC_CONFIGURATION.n_eng; % Number of engines
                prop_data(3) = 30; % Thrust (lbf) or Power (shp) per engine
                prop_data(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
                prop_data(5) = 2; % Normativa
                prop_data(6) = 0.82; % Prop efficiency
                prop_data(7) = 1; % By-pass
                Prop_sel = prop_data;
                
                % Range of low and high speed to be analized
                V_low = 250/3.6; % low value in m/s
                V_high = 450/3.6; % high value in m/s
                N_V_VAR_perf = 5; % number of iterations
                % Value of single speed
                V_single = 250/3.6; % single value in m/s
                
                % Range of low and high payloads to be analized
                Wp_low = 10*1000; % low value in kg
                Wp_high = 30*1000; % high value in kg
                N_Wp_VAR_perf = 2;
                % Value of single payload
                W_single = 23*1000; % single value in kg
                Post_processing_PERFORMANCE = 0;

            case 3 % case_AC = 3 - PEPIÑO XXL
                % Enter type of mission segments betwee brackets being
                type_missions = [5 13 13]; % mission for the weight fraction method
                num_missions = length(type_missions);
                %% User needs to define the inputs for each segment
                FlightMODE_var.climb_mode = [2 2 2];
                FlightMODE_var.turn_mode = [5 5 5];
                FlightMODE_var.descent_mode = [2 2 2];
                FlightMODE_var.cruise_mode = [2 2 2];
                
                prop_data(1) = 2; % Type of engine
                prop_data(2) = AC_CONFIGURATION.n_eng; % Number of engines
                prop_data(3) = 30; % Thrust (lbf) or Power (shp) per engine
                prop_data(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
                prop_data(5) = 2; % Normativa
                prop_data(6) = 0.82; % Prop efficiency
                prop_data(7) = 1; % By-pass
                Prop_sel = prop_data;
                
                % Range of low and high speed to be analized
                V_low = 50; % low value in m/s
                V_high = 120; % high value in m/s
                N_V_VAR_perf = 4; % number of iterations
                % Value of single speed
                V_single = 80; % single value in m/s
                
                % Range of low and high payloads to be analized
                Wp_low = 10; % low value in kg
                Wp_high = 30; % high value in kg
                N_Wp_VAR_perf = 1;
                % Value of single payload
                W_single = 10; % single value in kg
                Post_processing_PERFORMANCE = 0;

            case 4 % case_AC = 4 - COMERCIAL
                % Enter type of mission segments betwee brackets being
                type_missions = [5 13 13]; % mission for the weight fraction method
                num_missions = length(type_missions);
                %% User needs to define the inputs for each segment
                FlightMODE_var.climb_mode = [2 2 2];
                FlightMODE_var.turn_mode = [5 5 5];
                FlightMODE_var.descent_mode = [2 2 2];
                FlightMODE_var.cruise_mode = [2 2 2];
                
                prop_data(1) = 2; % Type of engine
                prop_data(2) = AC_CONFIGURATION.n_eng; % Number of engines
                prop_data(3) = 30; % Thrust (lbf) or Power (shp) per engine
                prop_data(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
                prop_data(5) = 2; % Normativa
                prop_data(6) = 0.82; % Prop efficiency
                prop_data(7) = 1; % By-pass
                Prop_sel = prop_data;
                
                % Range of low and high speed to be analized
                V_low = 250/3.6; % low value in m/s
                V_high = 450/3.6; % high value in m/s
                N_V_VAR_perf = 5; % number of iterations
                % Value of single speed
                V_single = 250/3.6; % single value in m/s
                
                % Range of low and high payloads to be analized
                Wp_low = 10*1000; % low value in kg
                Wp_high = 30*1000; % high value in kg
                N_Wp_VAR_perf = 2;
                % Value of single payload
                W_single = 23*1000; % single value in kg
                Post_processing_PERFORMANCE = 0;

            case 5 % case_AC = 5 - WIG
                % Enter type of mission segments betwee brackets being
                type_missions = [5 13 13]; % mission for the weight fraction method
                num_missions = length(type_missions);
                %% User needs to define the inputs for each segment
                FlightMODE_var.climb_mode = [2 2 2];
                FlightMODE_var.turn_mode = [5 5 5];
                FlightMODE_var.descent_mode = [2 2 2];
                FlightMODE_var.cruise_mode = [2 2 2];
                
                prop_data(1) = 2; % Type of engine
                prop_data(2) = AC_CONFIGURATION.n_eng; % Number of engines
                prop_data(3) = 8500; % Thrust (lbf) or Power (shp) per engine
                prop_data(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
                prop_data(5) = 2; % Normativa
                prop_data(6) = 0.82; % Prop efficiency
                prop_data(7) = 1; % By-pass
                Prop_sel = prop_data;
                
                % Range of low and high speed to be analized
                V_low = 250/3.6; % low value in m/s
                V_high = 450/3.6; % high value in m/s
                N_V_VAR_perf = 3; % number of iterations
                % Value of single speed
                V_single = 250/3.6; % single value in m/s
                
                % Range of low and high payloads to be analized
                Wp_low = 10*1000; % low value in kg
                Wp_high = 30*1000; % high value in kg
                N_Wp_VAR_perf = 2;
                % Value of single payload
                W_single = 23*1000; % single value in kg
                Post_processing_PERFORMANCE = 1;

            case 6 % case_AC = 6 - CERVERA
                % Enter type of mission segments betwee brackets being
                type_missions = [5 13 13]; % mission for the weight fraction method
                num_missions = length(type_missions);
                %% User needs to define the inputs for each segment
                FlightMODE_var.climb_mode = [2 2 2];
                FlightMODE_var.turn_mode = [5 5 5];
                FlightMODE_var.descent_mode = [2 2 2];
                FlightMODE_var.cruise_mode = [2 2 2];
                
                prop_data(1) = 2; % Type of engine
                prop_data(2) = AC_CONFIGURATION.n_eng; % Number of engines
                prop_data(3) = 30; % Thrust (lbf) or Power (shp) per engine
                prop_data(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
                prop_data(5) = 2; % Normativa
                prop_data(6) = 0.82; % Prop efficiency
                prop_data(7) = 1; % By-pass
                Prop_sel = prop_data;
                
                % Range of low and high speed to be analized
                V_low = 250/3.6; % low value in m/s
                V_high = 450/3.6; % high value in m/s
                N_V_VAR_perf = 5; % number of iterations
                % Value of single speed
                V_single = 250/3.6; % singe value in m/s
                
                % Range of low and high payloads to be analized
                Wp_low = 10*1000; % low value in kg
                Wp_high = 30*1000; % high value in kg
                N_Wp_VAR_perf = 2;
                % Value of single payload
                W_single = 23*1000; % single value in kg
                
                % Conducts postprocessing Performences : CAPS, Fuel, Time,
                % etc...
                Post_processing_PERFORMANCE = 1;
        end
        
        % Creates the variable or single condictions
        if variable_speed_AP == 1
            V_VAR_perf = linspace(V_low,V_high,N_V_VAR_perf);
        else
            N_V_VAR_perf = 1; % number of iterations
            V_VAR_perf = V_single;
        end
        Plot_Options.V_VAR_perf = V_VAR_perf;
        Plot_Options.N_V_VAR_perf = N_V_VAR_perf;
        
        % Creates the variable or single condictions
        if variable_weight_AP == 1
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
                
                handles = caracteristicas_avanz_v2(Aero_TH,Prop_data,Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Aero,Prop_sel);
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
                
                Segments = Define_Segments(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_var,num_missions);
                Segments{1}.crucero(1) = h_inicial_cr;
                Segments{1}.crucero(2) = dist_final_cr;
                Segments{1}.crucero(4) = Mach_cr;
                CL_opt = sqrt(Aero_TH.CD0/Aero_TH.CD2);
                Segments{1}.crucero(5) = CL_opt;
                
                % Generates the different segments
                [seg] = Generation_Mission_Segments_v3(conv_UNITS,num_missions,type_missions,Performance,case_AC,Segments,FlightMODE_var);
                
                [Weights_AP,Total_Datos,datos] = procesar_mision_v3(seg,num_missions,handles,Segments,FlightMODE_var,Geo_tier);
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
conditions.x_XCG = XCG_data.x_XCG;

% Selects XCG depending if it's from desired stability conditions in Forward Flight (XCG_FF = 1) or
% AXIAL Flight XCG_FF = 0)
XCG_FF =1;

if STABILITY_STUDY == 1
    %% Dynamic Stability Formulation
    % Stability_Pamadi = 0 - Fran's Formulation
    % Stability_Pamadi = 1 - Pamadi
    % Stability_Pamadi = 2 - Roskam
    if STABILITY_STUDY_Trim == 1
        only_trim = 0;
        StabilityModel =1;
        [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts]  = Calculo_Stability_Derivatives_Abril2020...
            (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
            Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,XCG_data,AC_CONFIGURATION,only_trim);
        %% XCG range
        x_XCG_fwd = TRIM_RESULTS.x_XCG_fwd;
        x_XCG_rwd = TRIM_RESULTS.x_XCG_rwd;
        N_x_XCG_VAR = 100;
        x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);
        Plot_Options.x_XCG_VAR = x_XCG_VAR;
    end
        
    %% Speed range
    V_low = Performance.V_min;
    V_high = Performance.V_max;
    N_V_VAR = 20;
    V_VAR = linspace(V_low,V_high,N_V_VAR);
    Plot_Options.V_VAR = V_VAR;
    
    %% Weight range
    % m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_systems + Weight_tier.m_batteries/2 + Weight_tier.m_fuel/2);
    m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy);
    m_mid = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy)/2;
    m_high = Weight_tier.m_TOW;
    N_m_VAR = 20;
    m_VAR = linspace(m_low,m_high,N_m_VAR);
    Plot_Options.W_VAR = m_VAR;
    
    if STABILITY_STUDY_Trim_varV_m0 == 1
        % reduces the calculations since only for trim conditions
        only_trim = 1;
        %% Study of variation of Velocity and mass with XCG selected for trim conditions at desired Static Margin
        for i=1:N_V_VAR
            conditions.x_XCG = TRIM_RESULTS.x_XCG;
            conditions.V = V_VAR(i);
            conditions.x_XCG_var = TRIM_RESULTS.x_XCG;
            for j=1:N_m_VAR
                conditions.m_TOW = m_VAR(j);
                %             [TRIM_RESULTS_var{i,j},Trim_ITER_var{i,j},Stab_Der{i,j}] = Calculo_Stability_Derivatives_Feb2020...
                [TRIM_RESULTS_calc,Trim_ITER_calc,Stab_Der_calc,Stab_Der_parts_calc] = Calculo_Stability_Derivatives_Abril2020...
                    (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
                    Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,XCG_data,AC_CONFIGURATION,only_trim);
                TRIM_RESULTS_var_V_XCG{i,j} = TRIM_RESULTS_calc;
                Trim_ITER_var_V_XCG{i,j} = Trim_ITER_calc;
                Stab_Der_var_V_XCG{i,j} = Stab_Der_calc;
            end
        end
    end
    
    %% Study of variation of Velocity with XCG selected for trim conditions at desired Static Margin
    if STABILITY_STUDY_Trim_var_XCG
        conditions.study_var_xcg = 1;
        for i=1:N_x_XCG_VAR
            x_XCG_var = TRIM_RESULTS.x_XCG;
            %         conditions.V = 80; % - Velocity at which are analized with XCG variation
            conditions.V = Performance.V;
            conditions.x_XCG = x_XCG_VAR(i);
            conditions.x_XCG_var = TRIM_RESULTS.x_XCG;
            [TRIM_RESULTS_calc,Trim_ITER_calc,Stab_Der_calc,Stab_Der_parts_calc] = Calculo_Stability_Derivatives_Abril2020...
                (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
                Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,XCG_data,AC_CONFIGURATION,only_trim);
            TRIM_RESULTS_var_XCG{i} = TRIM_RESULTS_calc;
            Trim_ITER_var_XCG{i} = Trim_ITER_calc;
            Stab_Der_var_XCG{i} = Stab_Der_calc;
        end
    end
    
    % Estimation through Fran's Formulation  if Stability_Pamadi = 0;
    % Estimation through Pamadi if Stability_Pamadi = 1;
    % Estimation through ROSKAM if Stability_Pamadi = 2;
    StabilityModel = 1;
    if STABILITY_STUDY_Long_dyn == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Stab_Dyn_Long = longitudinal_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,Weight_tier,TRIM_RESULTS,Geo_tier);
    end
    
    if STABILITY_STUDY_LatDir_dyn == 1
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Stab_Dyn_LatDir = lateral_directional_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,Weight_tier,TRIM_RESULTS,Geo_tier);
    end
    
    if STABILITY_STUDY_Trim_lat == 1
        conditions.study_var_xcg = 0;
        % %%%%%%%%%%%%%%%%%%%%TRIM LATERAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        beta = 11.5*D2R; % vaue of beta at qhich wants to be analyze
        % variable beta
        beta_i = 0.1*D2R;
        beta_f = 25*D2R;
        Delta_beta = 30;
        beta_vec = linspace(beta_i,beta_f,Delta_beta);
        % Storing study conditions
        conditions_TRIM_lat.beta = beta;
        conditions_TRIM_lat.beta_vec = beta_vec;
        [Trim_ITER_LAT] = Calculo_Trim_ITER_LAT_v3(Performance,Geo_tier,conv_UNITS,Trim_ITER,...
            Stab_Der,Weight_tier,conditions_TRIM_lat);
    end
    
    if STABILITY_STUDY_Turning == 1
        % %%%%%%%%%%%%%%%%%%%% Turnng conditioms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % variable beta
        phi_i = 5*D2R;
        phi_f = 77*D2R;
        Delta_phi = 100;
        phi = 30*D2R;
        n_viraje = 4.4; % Load factor
        phi_vec = linspace(phi_i,phi_f,Delta_phi);
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
% Case 11 PERFORMANCE_STUDY

%% Defines options for the plots
[mark_Type COLOR_scheme, Plot_Options] = Generate_Plot_Options(prefix,mark_legend,VECTOR_XFLR5,Plot_Options);

for i=1:length(PLOTS)
    switch PLOTS(i)
        case 1 % compare = 1 elements: w1,w2, vtp... alone
            % Logic ensures that the studies have been conducted
            if AERODYNAMIC_STUDY == 1
                [Fig] = plot_XFLR5_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix);
            else
                Warning = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 2 % compare = 1 elements: w1,w2, vtp... alone
            % Logic ensures that the studies have been conducted
            if AERODYNAMIC_STUDY == 1
                [Fig] = plot_XFLR5_Polar_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix);
            else
                Warning = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 3
            % Logic ensures that the studies have been conducted
            if AERODYNAMIC_STUDY == 1
                [XAC,Fig] = plot_Generate_XAC(DATA_Ae,Aero,DATA_PL,CASOS,degrees_XAC,Fig,...
                    XFLR5_DATA,Performance,X_OC,Plot_Options,VECTOR_XFLR5);
                
            else
                Warning = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 4
            [Fig] = Generates_Plots_PropulsionModels(Prop_data,Data_Trim,V,alpha,D,Plot_Options,Fig);
            
        case 5
            [Fig] = plot_GEOMETRY_2018(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC,Performance);
            
        case 6
%             % Logic ensures that the studies have been conducted
%             if STABILITY_STUDY == 1
%                 if STABILITY_STUDY_Trim == 1
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
            if STABILITY_STUDY == 1
                if STABILITY_STUDY_Trim_var_XCG == 1
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
            if STABILITY_STUDY == 1
                if STABILITY_STUDY_Trim_varV_m0 == 1
                    [Fig] = Generates_Plots_Longitudinal_Trim_VAR(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var_V_XCG,Trim_ITER_var_V_XCG,Geo_tier,Plot_Options,Fig);
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
            if STABILITY_STUDY == 1
                if STABILITY_STUDY_Trim_lat == 1
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
            if STABILITY_STUDY == 1
                if STABILITY_STUDY_Turning == 1
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
            if PERFORMANCE_STUDY == 1
                if ANALYSIS_PERFORMANCE_AP_var == 1
                    [Fig] = Generates_Plots_Performance_v1(Geo_tier,Plot_Options,conv_UNITS,Fig,...
                        Segments, handles,Weights_AP_var,fuel_total_var,tiempo_total_var,distancia_total_var,W_var,datos_var,case_AC,Plots_performance);
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
    if Detailed_Profile
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