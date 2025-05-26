%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANAID - Stability ANAlysis Interactive Design Tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
path(pathdef)
% clc
% Call function that Defines Generic Path so that files can be organized in folders
get_add_path

% Units conversion
conv_UNITS = conversion_UNITS;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
g = conv_UNITS.g;

%% Initializes figures
Fig = 0;

%% Prompts User for reading process
[case_AC OUTPUT_read_XLSX] = Initial_prompt_advanced;

%% MATLAB Compatibility
% Flag that determines MATLAB incompatibility (For example functions that do not work in MATLAB 2017 )
MATLAB_in = OUTPUT_read_XLSX.MATLAB_flags.MATLAB_in;
CHECK_Efficiency = OUTPUT_read_XLSX.MATLAB_flags.CHECK_Efficiency;
Detailed_Profile = OUTPUT_read_XLSX.MATLAB_flags.Detailed_Profile;

%% Defines options for the plots
[mark_Type COLOR_scheme, Plot_Options] = Generate_Plot_Options_init;
Plot_Options.MATLAB_in = MATLAB_in;
% Scaling Factor
% SF = OUTPUT_read_XLSX.AC_Data_flags.SF;

% Defines the flag to determine method to Estimate weights
Weight_Estimation = OUTPUT_read_XLSX.AC_Data_flags.Weight_Estimation;

%% Estimation of prop diameter, just preliminary for geometric conditions
% Stores the flags
% Propulsive_flags.type_battery = type_battery; %
% Propulsive_flags.Engine_loc = Engine_loc; %
% Propulsive_flags.Engine_conf = Engine_conf; %

prop_known = OUTPUT_read_XLSX.Propulsive_flags.prop_known; %
prop_properties = OUTPUT_read_XLSX.Propulsive_flags.prop_properties; %

if prop_known == 1
    number = prop_properties/100;
    integ = floor(number);
    fract = (number-integ)*10;
    D_prop = integ*2.54/100;
    OUTPUT_read_XLSX.Propulsive_flags.D_prop = D_prop;
else
    D_prop = OUTPUT_read_XLSX.Propulsive_flags.D_prop; %
end

%% Aircraft type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
% AC_type = OUTPUT_read_XLSX.AC_Data_flags.AC_type;
[AC_CONFIGURATION] = Generation_AC_configuration(conv_UNITS,OUTPUT_read_XLSX,case_AC);
AC_CONFIGURATION.type_battery = OUTPUT_read_XLSX.Propulsive_flags.type_battery;
Geo_input_tier = Generation_Input_Geometric_Data_v2(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX);
% Determines the fuselage than will be shown
% CASE_fuse = 1 - Just STL of fuselage
% CASE_fuse = 2 - STL of the entire aircraft
% CASE_fuse = 3 - Fuselage from XFLR5
% CASE_fuse = OUTPUT_read_XLSX.AC_Data_flags.CASE_fuse;
% Determine the type of scaling
% ESCALADO = 0; % no scaling Ratio =1
% ESCALADO = 1; % scaling acording to desired for each 3 axix
% ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
ESCALADO = OUTPUT_read_XLSX.AC_Data_flags.ESCALADO;
% [XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_EMERGENTIA;
[XFLR5_DATA,XFLR5_file,STL_files] = Extraction_fuselage_data(OUTPUT_read_XLSX);

Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX,case_AC);

if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    save('data/Geo_tier.mat', 'Geo_tier')
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

% %% Generation of Geometry
% [Body_Geo,meshData] = Read_Fuselage_Data(Geo_tier,XFLR5_DATA,ESCALADO,XFLR5_file,STL_files,OUTPUT_read_XLSX); % Defines Fuselage DATA
% [Body_Geo,meshData] = Generation_Fuselage_Data(Geo_tier,OUTPUT_read_XLSX,XFLR5_DATA,XFLR5_file,STL_file,conversion_UNITS); % Defines Fuselage DATA
[Fig,Body_Geo,meshData] = Generation_Fuselage_Data_old(Geo_tier,XFLR5_DATA,ESCALADO,XFLR5_file,STL_files,OUTPUT_read_XLSX,Fig);% Defines Propulsion DATA 


%% Defines Estimation of Weights according to Cefiro III densities
Weight_tier = Generation_Weight_Data(Geo_tier,Body_Geo,AC_CONFIGURATION,conv_UNITS,Weight_Estimation,OUTPUT_read_XLSX);
m_TOW = Weight_tier.m_TOW;

% Save MAT Structure
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    save('data/Weight_tier.mat', 'Weight_tier')
end

% Propulsion data for combustion Engine
Propulsion = Propulsion_Data_CE(AC_CONFIGURATION);

% Desired Static Margin
SM_des = OUTPUT_read_XLSX.Stability_flags.SM_des;

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
Performance_preliminar.Flight_cruise = OUTPUT_read_XLSX.Performance_pre_flags.Flight_cruise; %
Performance_preliminar.Flight_takeoff = OUTPUT_read_XLSX.Performance_pre_flags.Flight_takeoff; %
Performance_preliminar.Flight_climb = OUTPUT_read_XLSX.Performance_pre_flags.Flight_climb; %

%% Defines the Reading files    
if OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY == 1 
    
    % Data for analysis of XFLR5 DATA
    % Generates the file for Aerodynamic Data
    if OUTPUT_read_XLSX.PLOT_flags.plot(1) == 1 || OUTPUT_read_XLSX.PLOT_flags.plot(2) == 1 % Decided to plot Aero and Aero Polar
        VECTOR_XFLR5 = Get_Comparing_vectors_Aerodynamic(case_AC,OUTPUT_read_XLSX);
    else
       DUMMY = 1;
       VECTOR_XFLR5 = DUMMY;
    end
        
    % Selection of areodynamic properties (incidence angles)
    Design_criteria = Generate_aerodynamic_data(OUTPUT_read_XLSX,AC_CONFIGURATION,conv_UNITS);
    
    % Read Aerodynamic Data
    [DATA_Ae,casos,prefix,mark_legend,X_OC] = Read_data_Aero(Performance_preliminar,case_AC,OUTPUT_read_XLSX);
    
    % Fusion Aerodynamic Propert1ies
    [Aero_TH] = Polar_Literatura_v1(Performance_preliminar,Geo_tier,conv_UNITS,Body_Geo,AC_CONFIGURATION,OUTPUT_read_XLSX);
    
    % Premiliminary Aerodynamic Design
    [Aero,DATA_PL,Performance] = Generate_Aero_Data_2021_v1(DATA_Ae,Design_criteria,Performance_preliminar,Geo_tier...
        ,Weight_tier,conv_UNITS,AC_CONFIGURATION, Body_Geo,OUTPUT_read_XLSX);
    % [Aero,DATA_PL,Fig] = Aerodynamic_Design_2020_v1(Geo_tier,...
    %     Weight_tier,conv_UNITS,Design_criteria,Performance,CASOS,degrees_XAC,Fig,XFLR5_DATA,...
    %     MAC_Estimation,DATA_Ae,X_OC,Plot_Options,VECTOR_XFLR5);
    
    %   Xac_cre_W_B = get_xac_cre_wing(Geo_tier.lambda_w1, Geo_tier.AR_w1,Geo_tier.Lambda_LE_w1,Performance.Mach);
    %   Geo_tier.xbar_w1 = Xac_cre_W_B*cR_w1;  %Consultar!!!
    %   Geo_tier.xbar_w1_e = Xac_cre_W_B*cR_w1;  %Consultar!!!    

    % Selects the Polar estimates that will be used including configuration
    % of Cruise, TakeOff or Climb
    Polar = get_Aero_Polar(Aero,Aero_TH,OUTPUT_read_XLSX,Performance);
    Aero.Polar = Polar;

    %% Saves the data to a mat file
    if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
        Storing_AERO_DATA_1 = Saving_data_Aero(VECTOR_XFLR5,Design_criteria,DATA_Ae,casos,prefix,mark_legend,X_OC,Aero_TH,Aero,DATA_PL,Performance);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_AERO_DATA_1 = dummy;
    end
    
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_AERO_DATA_1 = dummy;
    % Warning
    warning_areo = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';;
    disp(warning_areo)
    pause
end

%% Data for analysis of Prop DATA
VECTOR_Prop = Get_Comparing_vectors_Propulsion(case_AC,OUTPUT_read_XLSX);

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

[Prop_data] = Generation_Propulsion_Data(AC_CONFIGURATION,Prop_selection,OUTPUT_read_XLSX);

% Defines Propulsion DATA
% [Fig] = plot_prop_APC(Data_P,Plot_Options,Fig,prefix,VECTOR_Prop)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Conducts Propulsion optimization %
if OUTPUT_read_XLSX.STUDY_flags.PROP_STUDY == 1
    Conduct_Prop_Optimzation
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Conducts Performance optimization %
if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
    Storing_PERFORMANCE_DATA_1 = Conduct_Performance_Optimization;
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_PERFORMANCE_DATA_1 = dummy;
% 
%         % Dummy variables
%     Storing_AERO_DATA_1
% Segments = 0;
%     handles = 0;
%     Weights_AP_var = 0;
%     fuel_total_var = 0;
%     tiempo_total_var = 0;
%     distancia_total_var = 0;
%     W_var = 0;
%     datos_var = 0;
%     Plots_performance = 0;
%     Post_processing_PERFORMANCE = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CONFIGURACIÓN PUSHER. MODELO DE HELICÓPTERO%%%%%%%%%%%%%%%
Posicion_Palanca = 1;
% Will be used in future version of propulsion
conditions.alpha_f = alpha_f;
conditions.beta_f = beta_f;
conditions.h = Performance.h;
conditions.V = Performance.V;
conditions.study_var_xcg = 0;
conditions.x_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
% conditions.x_XCG = 1.0013;
% conditions.V = 30;

% Selects XCG depending if it's from desired stability conditions in Forward Flight (XCG_FF = 1) or
% AXIAL Flight XCG_FF = 0)
XCG_FF =1;
XCG_FF = OUTPUT_read_XLSX.Stability_flags.XCG_FF;

if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dynamic Stability Formulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Trim Studies
%     if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim == 1
%         conditions.m_TOW = m_TOW;
%         % only_trim = OUTPUT_read_XLSX.Stability_flags.only_trim;
%         % Forces only trim so that obtains preliminary resutls
%         only_trim = 1;
%         StabilityModel = OUTPUT_read_XLSX.Stability_flags;
%         [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts,Propulsion_values] = ...
%             Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
%             Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,...
%             Performance_preliminar);
% 
% %         %% XCG range
% %         x_XCG_fwd = TRIM_RESULTS.x_XCG_fwd;
% %         x_XCG_rwd = TRIM_RESULTS.x_XCG_rwd;
% %         N_x_XCG_VAR = 100;
% %         x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);
% %         Plot_Options.x_XCG_VAR = x_XCG_VAR;
%     end

    %% Regular Stability Study
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Regular == 1
        conditions.m_TOW = m_TOW;
        only_trim = OUTPUT_read_XLSX.Stability_flags.only_trim;
        StabilityModel = OUTPUT_read_XLSX.Stability_flags;
        [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts,Propulsion_values] = ...
                Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
            Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar);

        StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
        Variable_Study = 0;
        %ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
            show_messages_screen = OUTPUT_read_XLSX.MATLAB_flags.show_messages_screen;
            Stab_Dyn_Long = longitudinal_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
                Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
        else
            Dummy = 0;
            Stab_Dyn_Long = Dummy;
        end
        %ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
           show_messages_screen = OUTPUT_read_XLSX.MATLAB_flags.show_messages_screen;
           Stab_Dyn_LatDir = lateral_directional_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
                Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
        else
            Dummy = 0;
            Stab_Dyn_LatDir = Dummy;
        end
        
        % Saving variables
        if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
            Storing_STABILITY_DATA_1 = Saving_data_Stability_TRIM(TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts,Stab_Dyn_Long,Stab_Dyn_LatDir,Propulsion_values);
        else
            % Stores dummy variable for continuity
            dummy = 1;
            Storing_STABILITY_DATA_1 = dummy;
        end       
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA_1 = dummy;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stability Sensitivity Study for Velocity and Mass
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
%     m_low = Weight_tier.m_TOW - (Weight_tier.m_energy);
%     m_mid = Weight_tier.m_TOW - (Weight_tier.m_energy)/2;
    m_high = Weight_tier.m_TOW;
    N_m_VAR = OUTPUT_read_XLSX.Stability_flags.N_m_VAR;
    m_VAR = linspace(m_low,m_high,N_m_VAR);
    Plot_Options.W_VAR = m_VAR;
    
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
        % reduces the calculations since only for trim conditions
        only_trim = OUTPUT_read_XLSX.Stability_flags.only_trim;
        %% Study of variation of Velocity and mass with XCG selected for trim conditions at desired Static Margin
        for i=1:N_V_VAR
%             conditions.x_XCG = TRIM_RESULTS.x_XCG_des;
            conditions.V = V_VAR(i);
%             conditions.x_XCG_var = TRIM_RESULTS.x_XCG_des;
            
            % actul values 
            conditions.x_XCG = conditions.x_XCG;
            conditions.x_XCG_var = conditions.x_XCG;
            
            for j=1:N_m_VAR
                conditions.m_TOW = m_VAR(j);
                V_stall = sqrt(2*m_VAR(j)*g/(Geo_tier.S_w1*Performance_preliminar.rho*Aero.CL_max_w1_CR));
                V_TO = 1.2*V_stall;
                V_min_ope = 1.3*V_stall;
                Restrictions.V_stall = V_stall;
                Restrictions.V_TO = V_TO;
                Restrictions.V_min_ope = V_min_ope;
                [TRIM_RESULTS_calc_var_V_XCGv,Trim_ITER_calc_var_V_XCGv,Stab_Der_calc_var_V_XCGv,Stab_Der_parts_calc_var_V_XCGv,Propulsion_calc_var_V_XCGv] = ...
                    Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
                    Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar);
                Restrictions_var_V_XCG{i,j} = Restrictions;
                TRIM_RESULTS_var_V_XCG{i,j} = TRIM_RESULTS_calc_var_V_XCGv;
                Trim_ITER_var_V_XCG{i,j} = Trim_ITER_calc_var_V_XCGv;
                Stab_Der_var_V_XCG{i,j} = Stab_Der_calc_var_V_XCGv;
                Stab_Der_parts_V_XCG{i,j} = Stab_Der_parts_calc_var_V_XCGv;
                Propulsion_var_V_XCG{i,j} =  Propulsion_calc_var_V_XCGv;
                
                StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
                Variable_Study = 1;
                
                %ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
                    show_messages_screen = 0;
                    Stab_Dyn_Long_var_V_XCG{i,j} = longitudinal_analysis(Performance,Stab_Der_calc_var_V_XCGv,conv_UNITS,StabilityModel,...
                        Weight_tier,TRIM_RESULTS_calc_var_V_XCGv,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
                else
                    Dummy = 0;
                    Stab_Dyn_Long_var_V_XCG = Dummy;
                end
                %ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
                    show_messages_screen = 0;
                    Stab_Dyn_LatDir_var_V_XCG{i,j} = lateral_directional_analysis(Performance,Stab_Der_calc_var_V_XCGv,conv_UNITS,StabilityModel,...
                        Weight_tier,TRIM_RESULTS_calc_var_V_XCGv,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
                else
                    Dummy = 0;
                    Stab_Dyn_LatDir_var_V_XCG = Dummy;
                end
            end
        end
        
        for i=1:N_V_VAR
            for j=1:N_m_VAR                
                TRIM_RESULTS_var_V_XCG{i,j}
            end
        end
        %% Saves the data to a mat file
        if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
            Storing_STABILITY_DATA_2 = Saving_data_Stability_varV(TRIM_RESULTS_var_V_XCG,Trim_ITER_var_V_XCG,Stab_Der_var_V_XCG,Stab_Der_parts_V_XCG,Stab_Dyn_Long_var_V_XCG,Stab_Dyn_LatDir_var_V_XCG,case_AC,Propulsion_var_V_XCG,Restrictions_var_V_XCG);
        else
            % Stores dummy variable for continuity
            dummy = 1;
            Storing_STABILITY_DATA_2 = dummy;
        end
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA_2 = dummy;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stability Sensitivity for varying Xcg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Study of variation of Velocity with XCG selected for trim conditions at desired Static Margin
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_var_XCG
        
        %% XCG range
        x_XCG_fwd = TRIM_RESULTS.x_XCG_fwd;
        x_XCG_rwd = TRIM_RESULTS.x_XCG_rwd;
        N_x_XCG_VAR = 100;
        x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);
        Plot_Options.x_XCG_VAR = x_XCG_VAR;
        
        conditions.study_var_xcg = 1;
        
        for i=1:N_x_XCG_VAR
            conditions.V = Performance.V;
            conditions.x_XCG = x_XCG_VAR(i);
            
            % Modification of the propulsion arms in X direction
            x_eng_xbar = Geo_tier.x_eng_xbar;
            % Determine the propulsive arms - Positive for engine behind Xcg
            x_d_T = x_eng_xbar - conditions.x_XCG; % Positive for engine behind Xcg
            % Storing DATA
            Geo_tier.x_d_T = x_d_T;            
            
            conditions.m_TOW = m_TOW;
            [TRIM_RESULTS_calc_var_xcg,Trim_ITER_calc_var_xcg,Stab_Der_calc_var_xcg,Stab_Der_parts_calc_var_xcg,Propulsion_calc_var_xcg] = ...
                Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
                Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar);
            TRIM_RESULTS_var_XCG{i} = TRIM_RESULTS_calc_var_xcg;
            Trim_ITER_var_XCG{i} = Trim_ITER_calc_var_xcg;
            Stab_Der_var_XCG{i} = Stab_Der_calc_var_xcg;
            Stab_Der_parts_var_XCG{i} = Stab_Der_parts_calc_var_xcg;
            Propulsion_var_xcg{i} = Propulsion_calc_var_xcg;
            %% Dynamic Stability Analysis
            StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
            Variable_Study = 1;
            %ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
                show_messages_screen = 0;
                Stab_Dyn_Long_var_XCG{i} = longitudinal_analysis(Performance,Stab_Der_calc_var_xcg,conv_UNITS,StabilityModel,...
                    Weight_tier,TRIM_RESULTS_calc_var_xcg,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
            else
                Dummy = 0;
                Stab_Dyn_Long_var_XCG = Dummy;
            end
            %ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
                show_messages_screen = 0;
                Stab_Dyn_LatDir_var_XCG{i} = lateral_directional_analysis(Performance,Stab_Der_calc_var_xcg,conv_UNITS,StabilityModel,...
                    Weight_tier,TRIM_RESULTS_calc_var_xcg,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
            else
                Dummy = 0;
                Stab_Dyn_LatDir_var_XCG = Dummy;
            end
        end
        %% Saves the data to a mat file
        if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
            Storing_STABILITY_DATA_3 = Saving_data_Stability_varXCG(TRIM_RESULTS_var_XCG,Trim_ITER_var_XCG,Stab_Der_var_XCG,Stab_Der_parts_var_XCG,...
                Stab_Dyn_Long_var_XCG,Stab_Dyn_LatDir_var_XCG,case_AC,Propulsion_var_xcg);
        else
            % Stores dummy variable for continuity
            dummy = 1;
            Storing_STABILITY_DATA_3 = dummy;
        end
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA_3 = dummy;
    end
    
    StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
%     Variable_Study = 0;
%     if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         show_messages_screen = OUTPUT_read_XLSX.MATLAB_flags.show_messages_screen;
%         Stab_Dyn_Long = longitudinal_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
%             Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
%     end
%     
%     if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
%         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % %%%%%%%%%%%%%ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
%         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         show_messages_screen = OUTPUT_read_XLSX.MATLAB_flags.show_messages_screen;
%         Stab_Dyn_LatDir = lateral_directional_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
%             Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
%     end
    
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
        %% Saves the data to a mat file
        if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
%             save Study_Trim_ITER_LAT.mat Trim_ITER_LAT conditions_TRIM_lat
            save('data/Study_Trim_ITER_LAT.mat', 'Trim_ITER_LAT','conditions_TRIM_lat')

            Storing_STABILITY_DATA_4.Trim_ITER_LAT = Trim_ITER_LAT;
            Storing_STABILITY_DATA_4.conditions_TRIM_lat = conditions_TRIM_lat;
        else
            % Stores dummy variable for continuity
            dummy = 1;
            Storing_STABILITY_DATA_4 = dummy;
        end
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA_4 = dummy;
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
        %% Saves the data to a mat file
        if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
%             save Study_Trim_ITER_LAT_Turn.mat Trim_ITER_LAT_Viraje conditions_TRIM_turning
            save('data/Study_Trim_ITER_LAT_Turn.mat', 'Trim_ITER_LAT_Viraje','conditions_TRIM_turning')
            Storing_STABILITY_DATA_5.Trim_ITER_LAT_Viraje = Trim_ITER_LAT_Viraje;
            Storing_STABILITY_DATA_5.conditions_TRIM_turning = conditions_TRIM_turning;
        else
            % Stores dummy variable for continuity
            dummy = 1;
            Storing_STABILITY_DATA_5 = dummy;
        end
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA_5 = dummy;
    end
    
else
    Plot_Options.dummy = 1;
    dummy = 1;
    Storing_STABILITY_DATA_1 = dummy;
    Storing_STABILITY_DATA_2 = dummy;
    Storing_STABILITY_DATA_3 = dummy;
    Storing_STABILITY_DATA_4 = dummy;
    Storing_STABILITY_DATA_5 = dummy;
end

%% Defines options for the plots
%% Generate Plots
GENERATE_PLOTS(OUTPUT_read_XLSX,Geo_tier,Plot_Options,Fig,Storing_STABILITY_DATA_1,Storing_STABILITY_DATA_2,Storing_STABILITY_DATA_3,...
    Storing_STABILITY_DATA_4,Storing_STABILITY_DATA_5,Storing_AERO_DATA_1,Performance,Prop_data,Storing_PERFORMANCE_DATA_1,conv_UNITS,Body_Geo,...
    meshData,AC_CONFIGURATION,case_AC);

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