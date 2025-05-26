%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANAID - Stability ANAlysis Interactive Design Tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
% Defines Generic Path so that files can be organized in folders
% change for the different versions
addpath(genpath('../../MATLAB/src'))
addpath(genpath('../../MATLAB/XFLR5'))
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
case_AC = 5;

%% Studies
% Conducts Performance Study
AERODYNAMIC_STUDY = 1;
PERFORMANCE_STUDY = 1;
% Conducts study of the limits of the props
Study_Prop_Limits = 0;
% Conducts study of sensitivity analysis of the Performance
PERFORMANCE_seNsitivity = 0;
% Analizes poerformance integrated with Academic Performance Codes
ANALYSIS_PERFORMANCE_AP = 0;
% Analizes poerformance integrated with Academic Performance Codes varying
% conditions
ANALYSIS_PERFORMANCE_AP_var = 1;
%Stability Analysis
STABILITY_STUDY = 0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Geometry %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PERFORMANCE ANALYSIS %%%%%%%%%%%%%%%
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        % Scaling Factor 
        SF = 1.;
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
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_EMERGENTIA(SF,conv_UNITS,AC_type,case_AC);
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
        Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data,case_AC);        
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
        Geo_tier = Generation_Geometric_Data_v2(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data,case_AC);
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
end

%% identifies the aerodynamic surfaces being used
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;

% Defines S ref
S_ref = Geo_tier.S_ref;

%% Note
% Geo_tier = Generation_Geometric_Data_v3 doe sinclude component W3 (VTP)
% if no VTP then use Generation_Geometric_Data_v2

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
[mark_Type] = Generates_plot_Info; % Generates the style of the lines and the Legend

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
        h_init = 200; % (m)
        V_init = 30; % speed
        V_max_init = 1.25*V_init; 
        
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        h_init = 200; % (m)
        V_init = 30; % speed
        V_max_init = 1.25*V_init; 
        
    case 3 % case_AC = 3 - PEPIÑO XXL
        h_init = 3000; % (m)
        V_init = 80; % speed
        V_max_init = 1.25*V_init; 
        
    case 4 % case_AC = 4 - COMERCIAL
        h_init = 3000; % (m)
        V_init = 90; % speed
        V_max_init = 1.25*V_init; 

    case 5 % case_AC = 5 - WIG
        h_init = 9; % (m)
        V_init = 70; % speed
        V_max_init = 1.25*V_init; 
end

% Indentifies preliminary performance results
Performance_preliminar.h = h_init;
Performance_preliminar.V = V_init;
Performance_preliminar.V_max = V_max_init;

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

if AERODYNAMIC_STUDY == 1
    %% Generates the file for Aerodynamic Data
    switch case_AC
        case 1 % case_AC = 1 - EMERGENTIA 1:1
            [casos prefix mark_legend X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 2 % case_AC = 2 - EMERGENTIA 1:2
            [casos prefix mark_legend X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 3 % case_AC = 3 - PEPIÑO XXL
            [casos prefix mark_legend X_OC] = read_aero_files_June2019_v3_PEPINOXXL(Performance_preliminar);
        case 4 % case_AC = 4 - COMERCIAL
            [casos prefix mark_legend X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        case 5 % case_AC = 5 - WIGL
            [casos prefix mark_legend X_OC] = read_aero_files_Feb2020_v3_WIG(Performance_preliminar);
    end
    %% Read Aerodynamic Data
    [DATA_Ae] = Read_data_Aero(casos);
    
    %% Fusion Aerodynamic Propert1ies
    [Aero_TH] = Polar_Literatura(Performance_preliminar,Geo_tier,conv_UNITS,Body_Geo,AC_CONFIGURATION,Conf);
    
    % % Propulsion_CE =  Propulsion_Data_CE(AC_CONFIGURATION);
    % % Weight Fraction
    % % NOTE: Uses first stimation of Polar to determine the weights
    % Weights_FP = Weight_Fraction(segment,Propulsion_CE,conv_UNITS,Aero_TH);
    % % Updates Weight according to design - includding separate terms and
    % % reserve fuel
    % Weights = Weight_UPDATE(Weights);
    
    %% Initializes figures
    Fig = 0;
    
    % Premiliminary Aerodynamic Design
    [Aero,DATA_PL,Performance] = Generate_Aero_Data_2020_v1(DATA_Ae,Design_criteria,Performance_preliminar,Geo_tier,Weight_tier,conv_UNITS,AC_CONFIGURATION);
    % [Aero,DATA_PL,Fig] = Aerodynamic_Design_2020_v1(Geo_tier,...
    %     Weight_tier,conv_UNITS,Design_criteria,Performance,CASOS,degrees_XAC,Fig,XFLR5_DATA,...
    %     MAC_Estimation,DATA_Ae,X_OC,Plot_Options,VECTOR_XFLR5);
end


% actualizes the Prop geometrý
Prop_data = Generation_Propulsion_Data(SF,D_prop,AC_CONFIGURATION); % Defines Propulsion DATA

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
    
    %% Generates the file for Aerodynamic Data
    switch case_AC
        case 1 % case_AC = 1 - EMERGENTIA 1:1
        case 2 % case_AC = 2 - EMERGENTIA 1:2
        case 3 % case_AC = 3 - PEPIÑO XXL
        case 4 % case_AC = 4 - COMERCIAL
            % Enter type of mission segments betwee brackets being
            type_missions_WF = [1 2 3 5 8 12]; % mission for the weight fraction method
            num_missions_WF = length(type_missions_WF);
            
            %% User needs to define the inputs for each segment
            Segments = Define_Segments_AC(Aero,Aero_TH,case_AC,conv_UNITS);
        case 5 % case_AC = 5 - WIG
            % Enter type of mission segments betwee brackets being
            type_missions_WF = [5]; % mission for the weight fraction method
            num_missions_WF = length(type_missions_WF);
            
            %% User needs to define the inputs for each segment
            Segments = Define_Segments_WIG(Aero,Aero_TH,case_AC,conv_UNITS);
    end
    
    [seg_WF] = Generation_Mission_Segments_v1(conv_UNITS,num_missions_WF,type_missions_WF,Performance,case_AC,Segments);
    Weights = Weight_Fraction(seg_WF,Propulsion,conv_UNITS,Aero_TH,type_missions_WF,Weight_tier,Geo_tier,AC_CONFIGURATION);
    
    % Performance Analysis
    type_missions = [1 2 3 5 7 8 10 11 12];
    type_missions = [5 13 13];
    num_missions = length(type_missions);
    [seg] = Generation_Mission_Segments_v1(conv_UNITS,num_missions,type_missions,Performance,case_AC,Segments);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Determine Minim Prop Diameter %%%%%%%%%%%%%%%
    Prop_data = Generation_Propulsion_Data(SF,D_prop,AC_CONFIGURATION); % Defines Propulsion DATA
    
    if Study_Prop_Limits == 1
        % Determination
        switch case_AC
            case 1 % case_AC = 1 - EMERGENTIA 1:1
                pp = 0.80; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 0.01; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Selects to plot Get Vertical Flight Limits that determines the
                % minimum prop diameter for the different maneuvers
                PLOT_Get_Vertical_Flight_Limits = 0;
                [Fig] = Get_Vertical_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vv,Vv_vect,Fig,PLOT_Get_Vertical_Flight_Limits);
                
                % Range of horizontal speeds
                Vh_min = V_min; % Vertical climb speed
                Vh_max = V_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Selects to plot Get Vertical Flight Limits that determines the
                % minimum prop diameter for the different maneuvers
                PLOT_Get_Horizontal_Flight_Limits = 0;
                [Fig] = Get_Horizontal_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vh,Vh_vect,Fig,PLOT_Get_Horizontal_Flight_Limits);
                
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
                
                % Selects to plot Get Vertical Flight Limits that determines the
                % minimum prop diameter for the different maneuvers
                PLOT_Get_Vertical_Flight_Limits = 0;
                [Fig] = Get_Vertical_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vv,Vv_vect,Fig,PLOT_Get_Vertical_Flight_Limits);
                
                % Range of horizontal speeds
                Vh_min = V_min; % Vertical climb speed
                Vh_max = V_min*1.13; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Selects to plot Get Vertical Flight Limits that determines the
                % minimum prop diameter for the different maneuvers
                PLOT_Get_Horizontal_Flight_Limits = 0;
                [Fig] = Get_Horizontal_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vh,Vh_vect,Fig,PLOT_Get_Horizontal_Flight_Limits);
            case 3 % case_AC = 3 - PEPIÑO XXL
                pp = 0.95; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                % Range of horizontal speeds
                Vh_min = V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Selects to plot Get Vertical Flight Limits that determines the
                % minimum prop diameter for the different maneuvers
                PLOT_Get_Horizontal_Flight_Limits = 0;
                [Fig] = Get_Horizontal_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vh,Vh_vect,Fig,PLOT_Get_Horizontal_Flight_Limits);
        end
    end
    
    PLOTS_SENSITIVITY_VERTICAL = 0;
    PLOTS_SENSITIVITY_HORIZONTAL = 0;
    % Uses Fzero which takes longer time
    solve_w_fero = 0;
    
    plots_3d_sensitivity = 0; % represents 3D plots
    plots_contour_sensitivity = 0; % represents contour plots
    special_PLOTS = 0; % represents the plotys searchjing for minimum Energy
    
    if PERFORMANCE_seNsitivity == 1
        % Selects the Engine-Prop configuration according to the Limit Study
        N_prop = 100; % number of prop elements to analyze between the limits
        pp_D_prop_min = 0.5;
        pp_D_prop_max = 2;
        N_contour_lines = 15; % number of contour lines
        
        switch case_AC
            case 1 % case_AC = 1 - EMERGENTIA 1:1
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
            case 2 % case_AC = 1 - EMERGENTIA 1:2
                D_prop = 22*2.54/100;
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
            case 3 % case_AC = 3 - PEPIÑO XXL
                D_prop = 32*2.54/100;
                % Actualizes Prop Data
                SF_prop = 0.85;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
        end
        
        % actualizes the Prop geometrý
        Prop_data = Generation_Propulsion_Data(SF,D_prop,AC_CONFIGURATION); % Defines Propulsion DATA
        
        % Sensitivity analysis ofr the Vertical Performance
        [Fig] = Sensitivity_Study_Vertical_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,...
            Data_ATM,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_VERTICAL,Performance,solve_w_fero,...
            plots_3d_sensitivity,plots_contour_sensitivity,special_PLOTS,N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines);
        
        % Sensitivity analysis ofr the Horizontal Performance
        [Fig] = Sensitivity_Study_Horizontal_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,...
            Data_ATM,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_HORIZONTAL,Performance,solve_w_fero,...
            plots_3d_sensitivity,plots_contour_sensitivity,special_PLOTS,N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines);
    end
    
    %% Performance
    if ANALYSIS_PERFORMANCE_AP ==1
        output_caracteristics = caracteristicas_avanz_v1(Aero_TH,Prop_data,Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Aero);
        [Weights_AP,fuel_total,tiempo_total,distancia_total,W,datos] = procesar_mision_v1(seg,num_missions,output_caracteristics,Segments);
    end
    
    %% Performance
    if ANALYSIS_PERFORMANCE_AP_var ==1
        V_low = 250/3.6;
        V_high = 450/3.6;
        N_V_VAR_perf = 10;
        V_VAR_perf = linspace(V_low,V_high,N_V_VAR_perf);
        Plot_Options.V_VAR_perf = V_VAR_perf;
        Plot_Options.N_V_VAR_perf = N_V_VAR_perf;
        for i = 1:N_V_VAR_perf
            % Flight Conditions For Cruise
            h_inicial_cr = Segments{1}.crucero(1);% - [m] % Altura inicial
            dist_final_cr = Segments{1}.crucero(2);% - [m] % 1: DISTANCIA FINAL
            V_cr = V_VAR_perf(i);
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_cr);
            Mach_cr = V_cr/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            Segments = Define_Segments_WIG(Aero,Aero_TH,case_AC,conv_UNITS);
            Segments{1}.crucero(1) = h_inicial_cr;
            Segments{1}.crucero(2) = dist_final_cr;
            Segments{1}.crucero(3) = V_VAR_perf(i);
            Segments{1}.crucero(4) = Mach_cr;
            type_missions = [5 13 13];
            [seg] = Generation_Mission_Segments_v1(conv_UNITS,num_missions,type_missions,Performance,case_AC,Segments);
            output_caracteristics = caracteristicas_avanz_v1(Aero_TH,Prop_data,Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Aero);
            [Weights_AP,fuel_total,tiempo_total,distancia_total,W,datos] = procesar_mision_v1(seg,num_missions,output_caracteristics,Segments);
            Weights_AP_var{i} = Weights_AP;
            fuel_total_var{i} = fuel_total;
            tiempo_total_var{i} = tiempo_total;
            distancia_total_var{i} = distancia_total;
            W_var{i} = W;
            datos_var{i}.crucero = datos(1).crucero;
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
        Stability_Pamadi =1;
        [TRIM_RESULTS,Trim_ITER,Stab_Der]  = Calculo_Stability_Derivatives_Feb2020...
            (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
            Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,XCG_data,AC_CONFIGURATION,only_trim);
        
        %% XCG range
        x_XCG_fwd = TRIM_RESULTS.x_XCG_fwd;
        x_XCG_rwd = TRIM_RESULTS.x_XCG_rwd;
        N_x_XCG_VAR = 100;
        x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);
        Plot_Options.x_XCG_VAR = x_XCG_VAR;
    end
    
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
                [TRIM_RESULTS_calc,Trim_ITER_calc,Stab_Der_calc] = Calculo_Stability_Derivatives_Feb2020...
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
            [TRIM_RESULTS_calc,Trim_ITER_calc,Stab_Der_calc] = Calculo_Stability_Derivatives_Feb2020...
                (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
                Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,XCG_data,AC_CONFIGURATION,only_trim);
            TRIM_RESULTS_var_XCG{i} = TRIM_RESULTS_calc;
            Trim_ITER_var_XCG{i} = Trim_ITER_calc;
            Stab_Der_var_XCG{i} = Stab_Der_calc;
            conditions.study_var_xcg = 0;
        end
    end
    
    
    % Estimation through Fran's Formulation  if Stability_Pamadi = 0;
    % Estimation through Pamadi if Stability_Pamadi = 1;
    % Estimation through ROSKAM if Stability_Pamadi = 2;
    Stability_Pamadi = 1;
    if STABILITY_STUDY_Long_dyn == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Stab_Dyn_Long = longitudinal_analysis(Performance,Stab_Der,conv_UNITS,Stability_Pamadi,Weight_tier,TRIM_RESULTS,Geo_tier);
    end
    
    if STABILITY_STUDY_LatDir_dyn == 1
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Stab_Dyn_LatDir = lateral_directional_analysis(Performance,Stab_Der,conv_UNITS,Stability_Pamadi,Weight_tier,TRIM_RESULTS,Geo_tier);
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
% Case 11 PRINT_PLOTS_TURNING_LAT = 0; % prints plots Performance
PLOTS = [11];

    %% Defines options for the plots
    [mark_Type COLOR_scheme, Plot_Options] = Generate_Plot_Options(prefix,mark_legend,VECTOR_XFLR5,Plot_Options);

for i=1:length(PLOTS)
    switch PLOTS(i)
        case 1 % compare = 1 elements: w1,w2, vtp... alone
            [Fig] = plot_XFLR5_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix);
            
        case 2
            [Fig] = plot_XFLR5_Polar_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix);
            
        case 3
            [XAC,Fig] = plot_Generate_XAC(DATA_Ae,Aero,DATA_PL,CASOS,degrees_XAC,Fig,...
                XFLR5_DATA,Performance,X_OC,Plot_Options,VECTOR_XFLR5);
            
        case 4
            [Fig] = Generates_Plots_PropulsionModels(Prop_data,Data_Trim,V,alpha,D,Plot_Options,Fig);
            
        case 5
            [Fig] = plot_GEOMETRY_2018(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC);
            
        case 6
            [Fig] = Generates_Plots_Long_Stability(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var,Trim_ITER_var,Geo_tier,Plot_Options,Fig);
            
        case 7
            [Fig] = Generates_Plots_Longitudinal_Trim(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var_XCG,Trim_ITER_var_XCG,Geo_tier,Plot_Options,Fig);
            
        case 8
            [Fig] = Generates_Plots_Longitudinal_Trim_VAR(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var_V_XCG,Trim_ITER_var_V_XCG,Geo_tier,Plot_Options,Fig);
            
        case 9
            [Fig] = Generates_Plots_Lateral_Trim(TRIM_RESULTS,Trim_ITER,Trim_ITER_LAT,Geo_tier,...
                Plot_Options,conv_UNITS,conditions_TRIM_lat,Fig);
            
        case 10
            [Fig] = Generates_Plots_Lateral_Turning(TRIM_RESULTS,Trim_ITER,Trim_ITER_LAT_Viraje,Geo_tier,...
                Plot_Options,conv_UNITS,conditions_TRIM_turning,Fig);
        case 11
            [Fig] = Generates_Plots_Performance_v1(Geo_tier,Plot_Options,conv_UNITS,Fig,...        
            Segments, output_caracteristics,Weights_AP_var,fuel_total_var,tiempo_total_var,distancia_total_var,W_var,datos_var);
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
        
%% Gets Forces and Moments
% Wind Components
% Wx = 0;
% Wy = 0;
% Wz = 0;
% 
% Wind.Wx = Wx;
% Wind.Wy = Wy;
% Wind.Wz = Wz;

% [F, M,DynVar] = Get_ForcesMoments(x,delta,Data_Trim,Data_Der,h,V,Wind,Data_ATM,Propulsion,Geo_tier);