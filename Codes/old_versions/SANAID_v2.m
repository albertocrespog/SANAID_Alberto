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

%% Prompt USER 
% Defines de aircraft to be analyzed
% case_AC = 1 - EMERGENTIA 1:1
% case_AC = 2 - EMERGENTIA 1:2
% case_AC = 3 - PEPIÑOXXL
% case_AC = 1;

prompt_typeAC = sprintf('Enter the Type of Aircraft \n CASE 1: EMERGENTIA \n CASE 2: EMERGENTIA scaled \n CASE 3: PEPIÑOXXL');
dlgtitle_typeAC = 'Input Aircraft Type';
dims = [1 35];
definput_typeAC = {'1'};
opts.Interpreter = 'tex';
answer_typeAC = inputdlg(prompt_typeAC,dlgtitle_typeAC,dims,definput_typeAC);
case_AC = str2num(answer_typeAC{1});

%% Prompt USER 
% Asks user for the input data        
% Enter type of mission segments betwee brackets being
% - 1 TakeOff
% - 2 = Climb
% - 3 = VTOL Climb
% - 4 = Cruise (Range)
% - 5 = Cruise (Endurance)
% - 6 = Descent
% - 7 = Descent (VTOL)
prompt1 = sprintf('Enter the number of mission segments:');
prompt2 = sprintf('Enter type of mission segments betwee brackets being \n TakeOff = 1 \n Climb = 2 \n VTOL Climb = 3 \n Cruise (Range) = 4 \n Cruise (Endurance) = 5 \n Descent = 6 \n Descent (VTOL) = 7');
prompt_mission = {prompt1,prompt2};
dlgtitle = 'Input # Missions';
dims = [1 50];
opts.Interpreter = 'tex';
% Examples of initialization
num_missions = 2; 
num_missions_str = num2str(num_missions);
missions_str = '[3 4]';
definput_mission = {num_missions_str,missions_str};
answer_mission = inputdlg(prompt_mission,dlgtitle,dims,definput_mission);

%% Prompt USER 
% Asks user for the input data        
% Introduce data for the different segments
N_missions = str2num(answer_mission{1});
type_missions = str2num(answer_mission{2});
for i=1:N_missions
    type_mission = type_missions(i);
    switch type_mission
        case 1 % case 1 TakeOff
            mission_tex = 'TakeOff';
        case 2 % case 2 Climb
            mission_tex = 'Climb';
        case 3 % case 3 VTOL Climb
            mission_tex = 'VTOL Climb';
        case 4 % case 4 Cruise (Range)
            mission_tex = 'Cruise (Range)';
        case 5 % case 5 Cruise (Endurance)
            mission_tex = 'Cruise (Endurance)';
        case 6 % case 6 Descent
            mission_tex = 'Descent';
        case 7 % case 7 Descent (VTOL)
            mission_tex = 'Descent (VTOL)';
    end
    
    uiwait(msgbox({'Enter Properties for Mission ',num2str(type_mission),mission_tex},'input Mission Properties','modal'));
    
    switch type_mission
        case 1 % case 1 TakeOff
            prompt1 = sprintf('Enter the initial altitude (m):');
            prompt2 = sprintf('Enter the final altitude (m):');
            prompt3 = sprintf('Enter Take Off Speed (m/s) (It is Just an Estimae):');
            prompt_mission1 = {prompt1,prompt2,prompt3};
            dlgtitle_mission1 = 'Segment Type 1 - Take off';
            dims = [1 60];
            opts.Interpreter = 'tex';
            % Examples of initialization
            h_initial_str = '0';
            h_final_str = '0';
            V_TO_str = '15';
            definput_mission1 = {h_initial_str,h_final_str,V_TO_str};
            answer_mission1 = inputdlg(prompt_mission1,dlgtitle_mission1,dims,definput_mission1);
            segment{i}.data.mision = type_mission; 
            segment{i}.data.h_initial = str2num(answer_mission1{1});
            segment{i}.data.h_final = str2num(answer_mission1{2});
            segment{i}.data.V_TO = str2num(answer_mission1{3});
            [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_TO);
            segment{i}.data.Data_ATM = Data_ATM;
            segment{i}.data.Performance = Performance;
        case 2 % case 2 Climb
            prompt1 = sprintf('Enter the initial altitude (m):');
            prompt2 = sprintf('Enter the final altitude (m):');
            prompt3 = sprintf('Enter the climb speed (m/s):');
            prompt4 = sprintf('Enter the flight ascent path angle (deg):');
            prompt_mission2 = {prompt1,prompt2,prompt3,prompt4};
            dlgtitle_mission2 = 'Segment Type 2 - Climb';
            dims = [1 60];
            opts.Interpreter = 'tex';
            % Examples of initialization
            h_initial_str = '0';
            h_final_str = '0';
            V_cl_str = '22';
            gamma_cl_str = '3';
            definput_mission2 = {h_initial_str,h_final_str,V_cl_str,gamma_cl_str};
            answer_mission2 = inputdlg(prompt_mission2,dlgtitle_mission2,dims,definput_mission2);
            segment{i}.data.mision = type_mission; 
            segment{i}.data.h_initial = str2num(answer_mission2{1});
            segment{i}.data.h_final = str2num(answer_mission2{2});
            segment{i}.data.V_cl = str2num(answer_mission2{3});
            segment{i}.data.gamma_cl = str2num(answer_mission2{4});
            [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_cl);
            segment{i}.data.V_V = segment{i}.data.V_cl*sin(segment{i}.data.gamma_cl);
            segment{i}.data.Data_ATM = Data_ATM;
            segment{i}.data.Performance = Performance;
        case 3 % case 3 VTOL Climb
            prompt1 = sprintf('Enter the initial altitude (m) - VTOL Flight:');
            prompt2 = sprintf('Enter the final altitude (m) - VTOL Flight:');
            prompt3 = sprintf('Enter the climb speed (m/s) - VTOL Flight:');
            prompt_mission3 = {prompt1,prompt2,prompt3};
            dlgtitle_mission3 = 'Segment Type 3 - VTOL Climb';
            dims = [1 60];
            opts.Interpreter = 'tex';
            % Examples of initialization
            h_initial_str = '0';
            h_final_str = '200';
            V_VTOL_str = '5';
            definput_mission3 = {h_initial_str,h_final_str,V_VTOL_str};
            answer_mission3 = inputdlg(prompt_mission3,dlgtitle_mission3,dims,definput_mission3);
            segment{i}.data.mision = type_mission; 
            segment{i}.data.h_initial = str2num(answer_mission3{1});
            segment{i}.data.h_final = str2num(answer_mission3{2});
            segment{i}.data.V_VTOL = str2num(answer_mission3{3});
            [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_VTOL);
            segment{i}.data.Data_ATM = Data_ATM;
            segment{i}.data.Performance = Performance;
        case 4 % case 4 Cruise(Range)
            prompt1 = sprintf('Enter the initial altitude (m):');
            prompt2 = sprintf('Enter the final altitude (m):');
            prompt3 = sprintf('Enter the cruise speed (m/s):');
            prompt4 = sprintf('Enter the Range (Km):');
            prompt_mission4 = {prompt1,prompt2,prompt3,prompt4};
            dlgtitle_mission4 = 'Segment Type 4 - Cruise Range';
            dims = [1 60];
            opts.Interpreter = 'tex';
            % Examples of initialization
            h_initial_str = '200';
            h_final_str = '200';
            V_cr_str = '26';
            Range_str = '20';
            definput_mission4 = {h_initial_str,h_final_str,V_cr_str,Range_str};
            answer_mission4 = inputdlg(prompt_mission4,dlgtitle_mission4,dims,definput_mission4);
            segment{i}.data.mision = type_mission; 
            segment{i}.data.h_initial = str2num(answer_mission4{1});
            segment{i}.data.h_final = str2num(answer_mission4{2});
            segment{i}.data.V_cr = str2num(answer_mission4{3});
            segment{i}.data.Range = str2num(answer_mission4{4});
            [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_cr);
            segment{i}.data.Data_ATM = Data_ATM;
            segment{i}.data.Performance = Performance;
        case 5 % case 5 Cruise (Endurance)
            prompt1 = sprintf('Enter the initial altitude (m):');
            prompt2 = sprintf('Enter the final altitude (m):');
            prompt3 = sprintf('Enter the cruise speed (m/s):');
            prompt4 = sprintf('Enter the Endurance Time (min):');
            prompt_mission5 = {prompt1,prompt2,prompt3,prompt4};
            dlgtitle_mission5 = 'Segment Type 5 - Cruise Endurance';
            dims = [1 60];
            opts.Interpreter = 'tex';
            % Examples of initialization
            h_initial_str = '0';
            h_final_str = '0';
            V_cr_str = '5';
            Endurance_str = '10';
            definput_mission5 = {h_initial_str,h_final_str,V_cr_str,Endurance_str};
            answer_mission5 = inputdlg(prompt_mission5,dlgtitle_mission5,dims,definput_mission5);
            segment{i}.data.mision = type_mission; 
            segment{i}.data.h_initial = str2num(answer_mission5{1});
            segment{i}.data.h_final = str2num(answer_mission5{2});
            segment{i}.data.V_cr = str2num(answer_mission5{3});
            segment{i}.data.Endurance = str2num(answer_mission5{4});
            [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_cr);
            segment{i}.data.Data_ATM = Data_ATM;
            segment{i}.data.Performance = Performance;
        case 6 % case 6 Descent
            prompt1 = sprintf('Enter the initial altitude (m):');
            prompt2 = sprintf('Enter the final altitude (m):');
            prompt3 = sprintf('Enter the descent speed (m/s):');
            prompt4 = sprintf('Enter the flight path angle (min):');
            prompt_mission6 = {prompt1,prompt2,prompt3,prompt4};
            dlgtitle_mission6 = 'Segment Type 6 - Descent';
            dims = [1 60];
            opts.Interpreter = 'tex';
            % Examples of initialization
            h_initial_str = '0';
            h_final_str = '0';
            V_dsc_str = '5';
            gamma_dsc_str = '5';
            definput_mission6 = {h_initial_str,h_final_str,V_dsc_str,gamma_dsc_str};
            answer_mission6 = inputdlg(prompt_mission6,dlgtitle_mission6,dims,definput_mission6);
            segment{i}.data.mision = type_mission; 
            segment{i}.data.h_initial = str2num(answer_mission6{1});
            segment{i}.data.h_final = str2num(answer_mission6{2});
            segment{i}.data.V_dsc = str2num(answer_mission6{3});
            segment{i}.data.gamma_dsc = str2num(answer_mission6{4});
            [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_dsc);
            segment{i}.data.V_V = segment{i}.data.V_dsc*sin(segment{i}.data.gamma_dsc);
            segment{i}.data.Data_ATM = Data_ATM;
            segment{i}.data.Performance = Performance;
        case 7 % case 7 Descent (VTOL)
            prompt1 = sprintf('Enter the initial altitude (m):');
            prompt2 = sprintf('Enter the final altitude (m):');
            prompt3 = sprintf('Enter the descent speed (m/s) - VTOL:');
            prompt_mission7 = {prompt1,prompt2,prompt3};
            dlgtitle_mission7 = 'Segment Type 7 - Descent - VTOL';
            dims = [1 60];
            opts.Interpreter = 'tex';
            % Examples of initialization
            h_initial_str = '0';
            h_final_str = '0';
            V_VTOL_str = '-2';
            definput_mission7 = {h_initial_str,h_final_str,V_VTOL_str};
            answer_mission7 = inputdlg(prompt_mission7,dlgtitle_mission7,dims,definput_mission7);
            segment{i}.data.mision = type_mission; 
            segment{i}.data.h_initial = str2num(answer_mission5{1});
            segment{i}.data.h_final = str2num(answer_mission5{2});
            segment{i}.data.V_VTOL = str2num(answer_mission5{3});
            [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_VTOL);
            segment{i}.data.Data_ATM = Data_ATM;
            segment{i}.data.Performance = Performance;
    end
end

%% Prompt USER 
% Asks user for the input data        
prompt = {'Enter Prop Diameter Estimation (in):','Enter scale factor %:','Enter Weight Estimation Method: case 1 - composite Céfiro III, case 2 - wood Céfiro I'};
dlgtitle = 'Input';
dims = [1 35];
opts.Interpreter = 'tex';

% input data user defined that apears in the promt command window
% Initial flight conditions (m)
h_in = 200; % altitude [m]
h_str = num2str(h_in);
% Range (km)
Range_in = 20;    %%% Distancia de la misión (km)
Range_str = num2str(Range_in);
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

definput = {D_propin_str,SF_str,Weight_str};
answer = inputdlg(prompt,dlgtitle,dims,definput);

%% Assigns values to variables temporarily from adquisition of mission
% Initial flight conditions
% h = 200; % altitude [m]
h = str2num(answer_mission4{1});
% Range
Range = str2num(answer_mission4{4})*1000;
% Cruise speed
V = str2num(answer_mission4{3});
% VTOL climb altitude
h_climb = str2num(answer_mission3{1});

% Prop diameter in inches
D_prop_in = str2num(answer{1});
D_prop = D_prop_in*conv_UNITS.in2m;
% Scaling Factor
SF = str2num(answer{2});
% Weight Estimation Method
Weight_Estimation = str2num(answer{3});

%% Performance Analysis
V_max = 1.25*V; % Max Speed [m/s]
[Data_ATM V_Performance] = Flight_Conditions_2020(h,V,V_max);
% Range
V_Performance.Range = Range;    %%% Distancia de la misión
V_Performance.h_climb = h_climb; % Vertical climb altitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Geometry %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PERFORMANCE ANALYSIS %%%%%%%%%%%%%%%
    
% [file_fusealge,path_fuselage] = uigetfile('*.txt',...
%    'Select One or More Files for the Fuselage', ...
%    'MultiSelect', 'on');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Geometry %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PERFORMANCE ANALYSIS %%%%%%%%%%%%%%%
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_EMERGENTIA(SF);
        Geo_input_tier = Generation_Input_Geometric_Data_EMERGENTIA(conv_UNITS,AC_CONFIGURATION,SF);
        [XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file] = Generation_data_fuselage_EMERGENTIA;
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_EMERGENTIA(SF);
        Geo_input_tier = Generation_Input_Geometric_Data_EMERGENTIA(conv_UNITS,AC_CONFIGURATION,SF);
        [XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file] = Generation_data_fuselage_EMERGENTIA;
    case 3 % case_AC = 3 - PEPIÑO XXL
        [XCG_data AC_CONFIGURATION] = Generation_AC_configuration_PEPINOXXL(SF);
        Geo_input_tier = Generation_Input_Geometric_Data_PEPINOXXL(conv_UNITS,AC_CONFIGURATION,SF);
        [XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file] = Generation_data_fuselage_PEPINOXXL;
end

%% Generation of Geometry
% Geo_tier = Generation_Geometric_Data_v2(Geo_input_tier,Prop_data,conv_UNITS,AC_CONFIGURATION,XCG_data);
Geo_tier = Generation_Geometric_Data_v2(Geo_input_tier,D_prop,conv_UNITS,AC_CONFIGURATION,XCG_data);
[Body_Geo,meshData] = Generation_Fuselage_Data(Geo_tier,XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file,SF); % Defines Propulsion DATA
%% Defines Estimation of Weights according to Cefiro III densities
Weight_tier = Generation_Weight_Data(Geo_tier,Body_Geo,AC_CONFIGURATION,conv_UNITS,Weight_Estimation,SF,case_AC);
m_TOW = Weight_tier.m_TOW;

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

%% Asks user for the input data        
% prompt = {'Enter flight altitude (m):','Enter Cruise Range (km)','Enter cruise speed (m/s):','Enter VTOL climb altitude (m)',...
%     'Enter Prop Diameter Estimation (in):','Enter scale factor %:','Enter Weight Estimation Method: case 1 - composite Céfiro III, case 2 - wood Céfiro I'};
% dlgtitle = 'Input';
% dims = [1 35];
% opts.Interpreter = 'tex';
% definput = {h_str,Range_str,V_str,h_climb_str,D_propin_str,SF_str,Weight_str};
% answer = inputdlg(prompt,dlgtitle,dims,definput);

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

%% Lets user select if wants to conduct study for the propellers
answer = questdlg('Would you like to conduct study for the limits of the propellers?', ...
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

%% Conducts study of the limits of the props
% Study_Prop_Limits = 0;
% Uses Fzero which takes longer time
solve_w_fero = 0;
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
            [Fig] = Get_Vertical_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vv,Vv_vect,Fig,PLOT_Get_Vertical_Flight_Limits,solve_w_fero);
            
            % Range of horizontal speeds
            Vh_min = V_min; % Vertical climb speed
            Vh_max = V_min*1.2; % Vertical climb speed
            N_Vh = 100; % Number of points
            Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
            
            % Selects to plot Get Vertical Flight Limits that determines the
            % minimum prop diameter for the different maneuvers
            PLOT_Get_Horizontal_Flight_Limits = 0;
            [Fig] = Get_Horizontal_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vh,Vh_vect,Fig,PLOT_Get_Horizontal_Flight_Limits,solve_w_fero);
            
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
            [Fig] = Get_Vertical_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vv,Vv_vect,Fig,PLOT_Get_Vertical_Flight_Limits,solve_w_fero);
            
            % Range of horizontal speeds
            Vh_min = V_min; % Vertical climb speed
            Vh_max = V_min*1.13; % Vertical climb speed
            N_Vh = 100; % Number of points
            Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
            
            % Selects to plot Get Vertical Flight Limits that determines the
            % minimum prop diameter for the different maneuvers
            PLOT_Get_Horizontal_Flight_Limits = 0;
            [Fig] = Get_Horizontal_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vh,Vh_vect,Fig,PLOT_Get_Horizontal_Flight_Limits,solve_w_fero);
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
            [Fig] = Get_Horizontal_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Data_ATM,Geo_tier,AC_CONFIGURATION,N_Vh,Vh_vect,Fig,PLOT_Get_Horizontal_Flight_Limits,solve_w_fero);
    end
end

PLOTS_SENSITIVITY_VERTICAL = 1;
PLOTS_SENSITIVITY_HORIZONTAL = 1;
% plots_3d_sensitivity = 0; % represents 3D plots
% plots_contour_sensitivity = 1; % represents contour plots
special_PLOTS = 0; % represents the plotys searchjing for minimum Energy

%% Conducts study of sensitivity analysis of the Performance

%% Lets user select if wants to conduct study for the propellers
answer = questdlg('Would you like to conduct Sensitivity Study To Determine the Propper Selection of Propeller Diameter?', ...
	'Prop Limits', ...
	'Yes','No','No');
% Handle response
switch answer
    case 'Yes'
        % Selects if 3D Plots are shown
        answer1 = questdlg('Would you like to show 3D plots?', ...
            '3D plots', ...
            'Yes','No','No');
        % Handle response
        switch answer1
            case 'Yes'
                plots_3d_sensitivity = 1;
            case 'No'
                plots_3d_sensitivity = 0;
        end
        
        % Selects if Contour Plots are shown
        answer2 = questdlg('Would you like to show contour plots?', ...
            '3D plots', ...
            'Yes','No','No');
        % Handle response
        switch answer2
            case 'Yes'
                plots_contour_sensitivity = 1;
            case 'No'
                plots_contour_sensitivity = 0;
        end
        PERFORMANCE_sensitivity = 1;
    case 'No'
%         disp([answer ' coming right up.'])
        PERFORMANCE_sensitivity = 0;
end

if PERFORMANCE_sensitivity == 1
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
        [X,Y,Z_T,Z_Pe,Z_delta,Z_RPM,Z_Q,Z_E,Z_etha,Fig] = Sensitivity_Study_Vertical_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,...
            Weight_tier,conv_UNITS,...
            Data_ATM,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_VERTICAL,V_Performance,solve_w_fero,...
            plots_3d_sensitivity,plots_contour_sensitivity,special_PLOTS,N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines);
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
                    
                    z_T = interp2(X,Y,Z_T,D_propin_q,Vv_q,'cubic');
                    z_Pe = interp2(X,Y,Z_Pe,D_propin_q,Vv_q,'cubic');
                    z_Q = interp2(X,Y,Z_Q,D_propin_q,Vv_q,'cubic');
                    z_delta = interp2(X,Y,Z_delta,D_propin_q,Vv_q,'cubic');
                    z_RPM = interp2(X,Y,Z_RPM,D_propin_q,Vv_q,'cubic');
                    z_etha = interp2(X,Y,Z_etha,D_propin_q,Vv_q,'cubic');

                    uiwait(msgbox({['Thrust = ',num2str(z_T),' (N)']...
                         ['Electric Power = ',num2str(z_Pe),' (W)']...
                         ['Torque = ',num2str(z_Q),' (Nm)']...
                         ['Throttle = ',num2str(z_delta),' (%)']...
                         ['RPM = ',num2str(z_RPM),' (RPM)']...
                         ['\eta = ',num2str(z_etha),]}));
                     
                    query_properties = 1;
                case 'No'
                    query_properties = 0;
            end
        end
    end
    
   
    if sensitivity_horizontal == 1
 
        [X,Y,Z_T,Z_Pe,Z_delta,Z_RPM,Z_Q,Z_E,Z_etha,Fig] = Sensitivity_Study_Horizontal_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,...
            Data_ATM,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_HORIZONTAL,V_Performance,solve_w_fero,...
            plots_3d_sensitivity,plots_contour_sensitivity,special_PLOTS,N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines);

        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';        
        query_properties = 1;
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
                    definput = {'30','5'};
                    answer_sv = inputdlg(prompt,dlgtitle,dims,definput);
                    
                    D_propin_q = str2num(answer_sv{1});
                    Vv_q = str2num(answer_sv{2});
                    
                    z_T = interp2(X,Y,Z_T,D_propin_q,Vv_q,'cubic');
                    z_Pe = interp2(X,Y,Z_Pe,D_propin_q,Vv_q,'cubic');
                    z_Q = interp2(X,Y,Z_Q,D_propin_q,Vv_q,'cubic');
                    z_delta = interp2(X,Y,Z_delta,D_propin_q,Vv_q,'cubic');
                    z_RPM = interp2(X,Y,Z_RPM,D_propin_q,Vv_q,'cubic');
                    z_etha = interp2(X,Y,Z_etha,D_propin_q,Vv_q,'cubic');
                                        
                    uiwait( msgbox({['Thrust = ',num2str(z_T),' (N)']...
                         ['Electric Power = ',num2str(z_Pe),' (W)']...
                         ['Torque = ',num2str(z_Q),' (Nm)']...
                         ['Throttle = ',num2str(z_delta),' (%)']...
                         ['RPM = ',num2str(z_RPM),' (RPM)']...
                         ['\eta = ',num2str(z_etha),]}));
                     
                    query_properties = 1;
                case 'No'
                    query_properties = 0;
            end
        end
    end
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