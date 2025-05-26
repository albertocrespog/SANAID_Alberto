%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANAID - Stability ANAlysis Interactive Design Tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
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
% case_AC = 1;

% prompt = {': 1- EMERGENTIA, 2 EMERGENTIA scaled, 3 PEPIÑOXXL'};
% uiwait(msgbox({'line1';'line2','line3','line4'},'CASE 1: EMERGENTIA','CASE 2: EMERGENTIA scaled','CASE 3: PEPIÑOXXL'));

prompt = sprintf('Enter the Type of Aircraft \n CASE 1: EMERGENTIA \n CASE 2: EMERGENTIA scaled \n CASE 3: PEPIÑOXXL');
dlgtitle = 'Input Aircraft Type';
dims = [1 35];
definput = {'1'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput);
case_AC = str2num(answer{1});

%% input data user defined
% Initial flight conditions (m)
h_in = 200; % altitude [m]
h_str = num2str(h_in);
% Range (km)
Range_in = 20;    %%% Distancia de la misión (km)
Range_str = num2str(Range_in)
% Cruise speed (m/s)
V_in = 26; % [m/s]
V_str = num2str(V_in);
% VTOL climb altitude (m)
h_climb_in = 200;
h_climb_str = num2str(h_climb_in);
% Scaling Factor
SF_in = 1;
SF_str = num2str(SF_in);
% Weight Estimation
Weight_in = 1;
Weight_str = num2str(Weight_in);
% Scaling Factor (in)
D_propin_in = 30;
D_propin_str = num2str(D_propin_in);

%% Asks user for the input data        
prompt = {'Enter flight altitude (m):','Enter Cruise Range (km)','Enter cruise speed (m/s):','Enter VTOL climb altitude (m)',...
    'Enter Prop Diameter Estimation (in):','Enter scale factor %:','Enter Weight Estimation Method: case 1 - composite Céfiro III, case 2 - wood Céfiro I'};
dlgtitle = 'Input';
dims = [1 35];
opts.Interpreter = 'tex';
definput = {h_str,Range_str,V_str,h_climb_str,D_propin_str,SF_str,Weight_str};
answer = inputdlg(prompt,dlgtitle,dims,definput)

%% Assigns values to variables
% Initial flight conditions
% h = 200; % altitude [m]
h = str2num(answer{1})
% Range
% Range = 20e3;    %%% Distancia de la misión
Range = str2num(answer{2})*1000
% Cruise speed
% V = 26; % [m/s]
V = str2num(answer{3})
% VTOL climb altitude
h_climb = str2num(answer{4})
% Scaling Factor
% SF = 1;
D_prop_in = str2num(answer{5})
D_prop = 100*D_prop_in/2.54;
SF = str2num(answer{6})
Weight_Estimation = str2num(answer{7})

V_max = 1.25*V; % Max Speed [m/s]
[Data_ATM V_Performance] = Flight_Conditions_2020(h,V,V_max);
% Range
V_Performance.Range = Range;    %%% Distancia de la misión
V_Performance.h_climb = h_climb; % Vertical climb altitude
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Geometry %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PERFORMANCE ANALYSIS %%%%%%%%%%%%%%%
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_EMERGENTIA(SF);
        Geo_input_tier = Generation_Input_Geometric_Data_EMERGENTIA(conv_UNITS,AC_CONFIGURATION,SF);
        [XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file] = Generation_data_fuselage_EMERGENTIA;
        % Defines the flag to determine method to Estimate weights
        % CASE1 = Weight_Estimation from composite - Cefiro III
        % CASE2 = Weight_Estimation from wood - Cefiro I
%         Weight_Estimation = 1;
        %% Estimation of prop diameter, just preliminary for geometric
%         % conditions
%         D_prop=30*2.54/100; 
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        % Scaling Factor 
        SF = 0.5;
        % Initial flight conditions
        h = 200; % altitude [m]
        % Cruise speed
        V = 25; % [m/s]
        V_max = 1.25*V; % Max Speed [m/s]
        [Data_ATM V_Performance] = Flight_Conditions_2020(h,V,V_max);
        % Range
        V_Performance.Range = 20e3;    %%% Distancia de la misión
        V_Performance.h_climb = 110; % Vertical climb altitude

        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_EMERGENTIA(SF);
        Geo_input_tier = Generation_Input_Geometric_Data_EMERGENTIA(conv_UNITS,AC_CONFIGURATION,SF);
        [XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file] = Generation_data_fuselage_EMERGENTIA;
        % Defines the flag to determine method to Estimate weights
        % CASE1 = Weight_Estimation from composite - Cefiro III
        % CASE2 = Weight_Estimation from wood - Cefiro I
        Weight_Estimation = 2;
        %% Estimation of prop diameter, just preliminary for geometric
        % conditions
        D_prop=22*2.54/100; 
    case 3 % case_AC = 3 - PEPIÑO XXL
        % Scaling Factor 
        SF = 1;
        % Initial flight conditions
        h = 1000; % altitude [m]
        % Cruise speed
        V = 80; % [m/s]
        V_max = 1.25*V; % Max Speed [m/s]
        [Data_ATM V_Performance] = Flight_Conditions_2020(h,V,V_max);
        % Range
        V_Performance.Range = 20e3;    %%% Distancia de la misión
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_PEPINOXXL(SF);
        Geo_input_tier = Generation_Input_Geometric_Data_PEPINOXXL(conv_UNITS,AC_CONFIGURATION,SF);
        [XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file] = Generation_data_fuselage_PEPINOXXL;
        % Defines the flag to determine method to Estimate weights
        % CASE1 = Weight_Estimation from composite - Cefiro III
        % CASE2 = Weight_Estimation from wood - Cefiro I
        Weight_Estimation = 1;
        %% Estimation of prop diameter, just preliminary for geometric
        % conditions
        D_prop=32*2.54/100; 
end


%

%% Generation of Geometry
% Geo_tier = Generation_Geometric_Data_v2(Geo_input_tier,Prop_data,conv_UNITS,AC_CONFIGURATION,XCG_data);
Geo_tier = Generation_Geometric_Data_v2(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data);
[Body_Geo,meshData] = Generation_Fuselage_Data(Geo_tier,XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file,SF); % Defines Propulsion DATA
%% Defines Estimation of Weights according to Cefiro III densities
Weight_tier = Generation_Weight_Data(Geo_tier,Body_Geo,AC_CONFIGURATION,conv_UNITS,Weight_Estimation,SF,case_AC);
m_TOW = Weight_tier.m_TOW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stores information for the plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mark_Type] = Generates_plot_Info; % Generates the style of the lines and the Legend

%% Generates the file for Aerodynamic Data
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        [casos prefix mark_legend X_OC] = read_aero_files_June2019_v3_EMERGENTIA(V_Performance);
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        [casos prefix mark_legend X_OC] = read_aero_files_June2019_v3_EMERGENTIA(V_Performance);
    case 3 % case_AC = 3 - PEPIÑO XXL
        [casos prefix mark_legend X_OC] = read_aero_files_June2019_v3_PEPINOXXL(V_Performance);
end

%% Read Aerodynamic Data
[DATA_Ae] = Read_data_Aero(casos); 

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
PLOTS = [0];

%% Defines options for the plots
[mark_Type COLOR_scheme, Plot_Options] = Generate_Plot_Options(prefix,mark_legend,VECTOR_XFLR5);

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
        Conf.w1 = 1; % wing
        Conf.h = 0; % HTP
        Conf.v = 0;% VTP
        Conf.v2 = 0;% double vertical
        Conf.vtail = 1;% V-tail
        Conf.fus = 1;% fuselage
        Conf.nac = 1;% nacelles
        Conf.landgear = 0;
        
         %% Theoretical Aerodynamics
        % Define how many cases have been analyzed to determine Xac
        Casos_XAC_w1 = 7;
        Casos_XAC_w2 = 9;
        
        % Defines Degres that are used to determine XAC
        degrees_XAC = [2.5 5 7.5 10 12.5 15 17.5];
        XAC_vec = [0 0.05 0.10 0.15 0.20 0.25];
        
        CASOS.Casos_XAC_w1 = Casos_XAC_w1;
        CASOS.Casos_XAC_w2 = Casos_XAC_w2;
        
        % mean aerodynamic chord according to XFLR5
        MAC_XFLR5 = 0.189;
        S_w1_XFLR5 = 0.429;
        XNP_XFLR5 = 0.030;
        XFLR5_DATA.MAC_XFLR5 = MAC_XFLR5;
        XFLR5_DATA.S_w1_XFLR5 = S_w1_XFLR5;
        XFLR5_DATA.XNP_XFLR5 = XNP_XFLR5;
        
        % Etimation of the MAC
        MAC_Estimation = 0;
        
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
        Conf.w1 = 1; % wing
        Conf.h = 0; % HTP
        Conf.v = 0;% VTP
        Conf.v2 = 0;% double vertical
        Conf.vtail = 1;% V-tail
        Conf.fus = 1;% fuselage
        Conf.nac = 1;% nacelles
        Conf.landgear = 0;
        
         %% Theoretical Aerodynamics
        % Define how many cases have been analyzed to determine Xac
        Casos_XAC_w1 = 7;
        Casos_XAC_w2 = 9;
        
        % Defines Degres that are used to determine XAC
        degrees_XAC = [2.5 5 7.5 10 12.5 15 17.5];
        XAC_vec = [0 0.05 0.10 0.15 0.20 0.25];
        
        CASOS.Casos_XAC_w1 = Casos_XAC_w1;
        CASOS.Casos_XAC_w2 = Casos_XAC_w2;
        
        % mean aerodynamic chord according to XFLR5
        MAC_XFLR5 = 0.189;
        S_w1_XFLR5 = 0.429;
        XNP_XFLR5 = 0.030;
        XFLR5_DATA.MAC_XFLR5 = MAC_XFLR5;
        XFLR5_DATA.S_w1_XFLR5 = S_w1_XFLR5;
        XFLR5_DATA.XNP_XFLR5 = XNP_XFLR5;
        
        % Etimation of the MAC
        MAC_Estimation = 0;

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
        Conf.w1 = 1; % wing
        Conf.h = 0; % HTP
        Conf.v = 0;% VTP
        Conf.v2 = 0;% double vertical
        Conf.vtail = 1;% V-tail
        Conf.fus = 1;% fuselage
        Conf.nac = 1;% nacelles
        Conf.landgear = 0;
        
        %% Theoretical Aerodynamics
        % Define how many cases have been analyzed to determine Xac
        Casos_XAC_w1 = 7;
        Casos_XAC_w2 = 9;
        
        % Defines Degres that are used to determine XAC
        degrees_XAC = [2.5 5 7.5 10 12.5 15 17.5];
        XAC_vec = [0 0.05 0.10 0.15 0.20 0.25];
        
        CASOS.Casos_XAC_w1 = Casos_XAC_w1;
        CASOS.Casos_XAC_w2 = Casos_XAC_w2;
        
        % mean aerodynamic chord according to XFLR5
        MAC_XFLR5 = 0.205;
        S_w1_XFLR5 = 0.452;
        XNP_XFLR5 = -0.031;
        XFLR5_DATA.MAC_XFLR5 = MAC_XFLR5;
        XFLR5_DATA.S_w1_XFLR5 = S_w1_XFLR5;
        XFLR5_DATA.XNP_XFLR5 = XNP_XFLR5;
        
        % Etimation of the MAC
        MAC_Estimation = 0;
end

%% Stores data
alpha_selected_w1 = i_w1; % (degs)
alpha_selected_w2 = i_w2; % (degs)
Design_criteria.index_w1 = index_w1;
Design_criteria.index_w2 = index_w2;
Design_criteria.alpha_selected_w1 = alpha_selected_w1;
Design_criteria.alpha_selected_w2 = alpha_selected_w2;
% Design choices for the incidence of the different surfaces
Design_criteria.i_w1 = i_w1*D2R; % incidence of Front Wing 
Design_criteria.i_w2 = i_w2*D2R; % incidence of Rear Wing
% Flight Safe Margin to calculate aerodynamic properties
Flight_SF = 1.2; % Stall Safe Margin Conditions
Design_criteria.Flight_SF = Flight_SF;
%% Propulsion Generation
% D = Prop_data.D_prop;
alpha_f = 0*D2R;
beta_f = 0*D2R;

%% Fusion Aerodynamic Propert1ies
[Aero_TH] = Polar_Literatura(V_Performance,Geo_tier,conv_UNITS,Body_Geo,AC_CONFIGURATION,Conf);

%% Initializes figures
Fig = 0;

% Premiliminary Aerodynamic Design
[Aero,DATA_PL,Fig] = Aerodynamic_Design_Aug2018(Geo_tier,...
    Weight_tier,conv_UNITS,Design_criteria,V_Performance,CASOS,degrees_XAC,Fig,XFLR5_DATA,...
    MAC_Estimation,DATA_Ae,X_OC,Plot_Options,VECTOR_XFLR5);
% Defines S ref
S_ref = Geo_tier.S_ref;

%% Determination of Stall conditions
alpha_max = Aero.alpha_max_w1_CR; 
CL_max_w1 = Aero.CL_w1_limit_max;
CL_max_w1_ope = CL_max_w1/(1.2^2);
alpha_max_w1_ope = interp1(DATA_Ae(index_w1).CL,DATA_Ae(index_w1).alpha,CL_max_w1_ope,'spline');
rho = V_Performance.rho;
V_stall = sqrt((2*m_TOW*g)/(rho*S_ref*CL_max_w1));
V_min = 1.2*V_stall;

V_Performance.alpha_max = alpha_max;
V_Performance.CL_max_w1 = CL_max_w1;
V_Performance.CL_max_w1_ope = CL_max_w1_ope;
V_Performance.alpha_max_w1_ope = alpha_max_w1_ope;
V_Performance.V_stall = V_stall;
V_Performance.V_min = V_min;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Determine Minim Prop Diameter %%%%%%%%%%%%%%%
Prop_data = Generation_Propulsion_Data(SF,D_prop,AC_CONFIGURATION); % Defines Propulsion DATA

%% Conducts study of the limits of the props
Study_Prop_Limits = 0;
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

PLOTS_SENSITIVITY_VERTICAL = 1;
PLOTS_SENSITIVITY_HORIZONTAL = 1;
% Uses Fzero which takes longer time
solve_w_fero = 0;
plots_3d_sensitivity = 0; % represents 3D plots
plots_contour_sensitivity = 1; % represents contour plots
special_PLOTS = 0; % represents the plotys searchjing for minimum Energy

%% Conducts study of sensitivity analysis of the Performance
PERFORMANCE_sensitivity = 1;
if PERFORMANCE_sensitivity == 1
    % This tool crates
    
    f = msgbox('Sensitivity Study To Determine the Propper Selection of Propeller Diameter');
            
    % Selects the Engine-Prop configuration according to the Limit Study
    N_prop = 100; % number of prop elements to analyze between the limits
    pp_D_prop_min = 0.5;
    pp_D_prop_max = 2;
    N_contour_lines = 15; % number of contour lines
    switch case_AC
        case 1 % case_AC = 1 - EMERGENTIA 1:1
            D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
            % Actualizes Prop Data
            SF_prop = 0.90;
            RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
            sensitivity_vertical = 1;
            sensitivity_horizontal = 1;
        case 2 % case_AC = 1 - EMERGENTIA 1:2
            D_prop = 22*2.54/100;
            % Actualizes Prop Data
            SF_prop = 0.90;
            RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
            % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
            sensitivity_vertical = 1;
            sensitivity_horizontal = 1;
        case 3 % case_AC = 3 - PEPIÑO XXL
            D_prop = 32*2.54/100;
            % Actualizes Prop Data
            SF_prop = 0.85;
            RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
            % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
            sensitivity_vertical = 0;
            sensitivity_horizontal = 1;
    end
    
    % actualizes the Prop geometry
    Prop_data = Generation_Propulsion_Data(SF,D_prop,AC_CONFIGURATION); % Defines Propulsion DATA
        
    if sensitivity_vertical == 1;
        % Sensitivity analysis ofr the Vertical Performance
        [X,Y,Z_Pe,Z_delta,Z_RPM,Z_Q,Z_E,Z_etha,Fig] = Sensitivity_Study_Vertical_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,...
            Data_ATM,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_VERTICAL,V_Performance,solve_w_fero,...
            plots_3d_sensitivity,plots_contour_sensitivity,special_PLOTS,N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines);
        
        f = msgbox('Go to the command window to select if you want to select a Prop Diameter and a velocity');
        prompt = 'Do you want more? Y/N [Y]: ';
        str = input(prompt,'s');
        if isempty(str)
            str = 'Y';
        end
        prompt = 'Select the prop diameter in inches?';
        D_propin_q = input(prompt); % query selection Prop
        prompt = 'Select the vertical speed ';
        Vv_q = input(prompt); % query selection Vertical Speed
        z_Pe = interp2(X,Y,Z_Pe,D_propin_q,Vv_q,'cubic');
        z_delta = interp2(X,Y,Z_delta,D_propin_q,Vv_q,'cubic');
        z_RPM = interp2(X,Y,Z_RPM,D_propin_q,Vv_q,'cubic');
        z_Q = interp2(X,Y,Z_Q,D_propin_q,Vv_q,'cubic');
        z_etha = interp2(X,Y,Z_etha,D_propin_q,Vv_q,'cubic');
        opts.Interpreter = 'tex';
        h = msgbox({sprintf('Electric Power: P_e = %g (KW)',z_Pe)},'Standard Values',CreateStruct);
        
    end
    
    if sensitivity_horizontal == 1
        % Sensitivity analysis ofr the Horizontal Performance
        [Fig] = Sensitivity_Study_Horizontal_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,...
            Data_ATM,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_HORIZONTAL,V_Performance,solve_w_fero,...
            plots_3d_sensitivity,plots_contour_sensitivity,special_PLOTS,N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines);
    end
    
    f = msgbox('Sensitivity Study Conducted - User required to selct Prop Diameter');
    pause
    
    prompt = 'Select the prop diameter in inches? ';
    D_propin = input(prompt);
    D_prop = 100*D_propin/2.54;
    prompt = 'Select the prop diameter in inches? ';
    D_propin = input(prompt);
    D_prop = 100*D_propin/2.54;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CONFIGURACIÓN PUSHER. MODELO DE HELICÓPTERO%%%%%%%%%%%%%%%
Posicion_Palanca = 1;
% Will be used in future version of propulsion
conditions.alpha_f = alpha_f;
conditions.beta_f = beta_f;
conditions.h = V_Performance.h;
conditions.V = V_Performance.V;
% Selects XCG depending if it's from desired stability conditions in Forward Flight (XCG_FF = 1) or
% AXIAL Flight XCG_FF = 0)
XCG_FF =1;

%% Dynamic Stability Formulation
% Stability_Pamadi = 0 - Fran's Formulation
% Stability_Pamadi = 1 - Pamadi
% Stability_Pamadi = 2 - Roskam
Stability_Pamadi =1;
[TRIM_RESULTS,Trim_ITER,Stab_Der]  = ...
     Calculo_Stability_Derivatives_Dec2019(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
    Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,V_Performance,Data_ATM,XCG_data,AC_CONFIGURATION);

%% XCG range
x_XCG_fwd = TRIM_RESULTS.x_XCG_fwd;
x_XCG_rwd = TRIM_RESULTS.x_XCG_rwd;
N_x_XCG_VAR = 100;
x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);
Plot_Options.x_XCG_VAR = x_XCG_VAR;

%% Speed range
V_low = V_min;
V_high = V_Performance.V_max;
N_V_VAR = 20;
V_VAR = linspace(V_low,V_high,N_V_VAR);
Plot_Options.V_VAR = V_VAR;

%% Weight range
m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_systems + Weight_tier.m_batteries/2 + Weight_tier.m_fuel/2);
m_high = Weight_tier.m_TOW;
N_m_VAR = 20;
m_VAR = linspace(m_low,m_high,N_m_VAR);
Plot_Options.W_VAR = m_VAR;

%% Study of variation of Velocity and mass with XCG selected for trim conditions at desired Static Margin
for i=1:N_V_VAR
    conditions.x_XCG = TRIM_RESULTS.x_XCG;
    conditions.V = V_VAR(i);
    conditions.x_XCG_var = TRIM_RESULTS.x_XCG;
    for j=1:N_m_VAR
        conditions.m_TOW = m_VAR(j);
        [TRIM_RESULTS_var{i,j},Trim_ITER_var{i,j}] = Calculo_Trimado_2019_var_v1...
            (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
    Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,V_Performance,Data_ATM,XCG_data,AC_CONFIGURATION);
    end
end

%% Study of variation of Velocity with XCG selected for trim conditions at desired Static Margin
for i=1:N_x_XCG_VAR
    x_XCG_var = TRIM_RESULTS.x_XCG;
    conditions.V = 80; % - Velocity at which are analized with XCG variation 
    conditions.x_XCG = x_XCG_VAR(i);
    conditions.x_XCG_var = TRIM_RESULTS.x_XCG;
        [TRIM_RESULTS_var{i,j},Trim_ITER_var{i,j}] = Calculo_Trimado_2019_var_v1...
            (conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
    Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,V_Performance,Data_ATM,XCG_data,AC_CONFIGURATION);
end

% Estimation through Fran's Formulation  if Stability_Pamadi = 0;
% Estimation through Pamadi if Stability_Pamadi = 1;
% Estimation through ROSKAM if Stability_Pamadi = 2;
Stability_Pamadi = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stab_Dyn_Long = longitudinal_analysis(V_Performance,Stab_Der,conv_UNITS,Stability_Pamadi,Weight_tier,TRIM_RESULTS,Geo_tier);
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stab_Dyn_LatDir = lateral_directional_analysis(V_Performance,Stab_Der,conv_UNITS,Stability_Pamadi,Weight_tier,TRIM_RESULTS,Geo_tier);

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
[Trim_ITER_LAT] = Calculo_Trim_ITER_LAT_v3(V_Performance,Geo_tier,conv_UNITS,Trim_ITER,...
    Stab_Der,Weight_tier,conditions_TRIM_lat);

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
[Trim_ITER_LAT_Viraje] = Calculo_Trim_ITER_LAT_Viraje_v3(V_Performance,Geo_tier,conv_UNITS,...
    Trim_ITER,Stab_Der,Weight_tier,conditions_TRIM_turning);

%% Plots graphics
for i=1:length(PLOTS)
    switch PLOTS(i)
        case 1 % compare = 1 elements: w1,w2, vtp... alone
            [Fig] = plot_XFLR5_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix);
            
        case 2
            [Fig] = plot_XFLR5_Polar_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix);
            
        case 3
            [XAC,Fig] = plot_Generate_XAC(DATA_Ae,Aero,DATA_PL,CASOS,degrees_XAC,Fig,...
                XFLR5_DATA,V_Performance,X_OC,Plot_Options,VECTOR_XFLR5);
            
        case 4
            [Fig] = Generates_Plots_PropulsionModels(Prop_data,Data_Trim,V,alpha,D,Plot_Options,Fig);
            
        case 5
            [Fig] = plot_GEOMETRY_2018(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC);
            
        case 6
            [Fig] = Generates_Plots_Long_Stability(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var,Trim_ITER_var,Geo_tier,Plot_Options,Fig);
            
        case 7
            [Fig] = Generates_Plots_Longitudinal_Trim(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var_XCG,Trim_ITER_var_XCG,Geo_tier,Plot_Options,Fig);
            
        case 8
            [Fig] = Generates_Plots_Longitudinal_Trim_VAR(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var,Trim_ITER_var,Geo_tier,Plot_Options,Fig);
            
        case 9
            [Fig] = Generates_Plots_Lateral_Trim(TRIM_RESULTS,Trim_ITER,Trim_ITER_LAT,Trim_ITER_var,Geo_tier,...
                Plot_Options,conv_UNITS,conditions_TRIM_lat,Fig);
            
        case 10
            [Fig] = Generates_Plots_Lateral_Turning(TRIM_RESULTS,Trim_ITER,Trim_ITER_LAT_Viraje,Trim_ITER_var,Geo_tier,...
                Plot_Options,conv_UNITS,conditions_TRIM_turning,Fig);
            
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