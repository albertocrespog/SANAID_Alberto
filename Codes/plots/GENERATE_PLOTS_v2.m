function [M_alpha_cero,V_alpha_cero] = GENERATE_PLOTS_v2(OUTPUT_read_XLSX,Storing_DATA,Plot_Options,Fig,conv_UNITS,AC_CONFIGURATION,case_AC,filenameS)

Storing_GEO_DATA_1 = Storing_DATA.Storing_GEO_DATA_1;
Storing_WEIGHT_DATA_1 = Storing_DATA.Storing_WEIGHT_DATA_1;
Storing_AERO_DATA_1 = Storing_DATA.Storing_AERO_DATA_1;
Storing_PROPULSION_DATA_1 = Storing_DATA.Storing_PROPULSION_DATA_1;
% Performance
Storing_PERFORMANCE_DATA_1 = Storing_DATA.Storing_PERFORMANCE_DATA_1;
Storing_PERFORMANCE_DATA_21 = Storing_DATA.Storing_PERFORMANCE_DATA_21;
Storing_PERFORMANCE_DATA_22 = Storing_DATA.Storing_PERFORMANCE_DATA_22;
Storing_PERFORMANCE_DATA_23 = Storing_DATA.Storing_PERFORMANCE_DATA_23;
Storing_PERFORMANCE_DATA_24 = Storing_DATA.Storing_PERFORMANCE_DATA_24;
Storing_PERFORMANCE_DATA_25 = Storing_DATA.Storing_PERFORMANCE_DATA_25;
Storing_PERFORMANCE_DATA_26 = Storing_DATA.Storing_PERFORMANCE_DATA_26;
Storing_PERFORMANCE_DATA_27 = Storing_DATA.Storing_PERFORMANCE_DATA_27;
% Stability
Storing_STABILITY_DATA_1 = Storing_DATA.Storing_STABILITY_DATA_1;
Storing_STABILITY_DATA_2 = Storing_DATA.Storing_STABILITY_DATA_2;
Storing_STABILITY_DATA_2B = Storing_DATA.Storing_STABILITY_DATA_2B;
Storing_STABILITY_DATA_2C = Storing_DATA.Storing_STABILITY_DATA_2C;
Storing_STABILITY_DATA_3 = Storing_DATA.Storing_STABILITY_DATA_3;
Storing_STABILITY_DATA_4A = Storing_DATA.Storing_STABILITY_DATA_4A;
Storing_STABILITY_DATA_4B = Storing_DATA.Storing_STABILITY_DATA_4B;
Storing_STABILITY_DATA_4C = Storing_DATA.Storing_STABILITY_DATA_4C;
Storing_STABILITY_DATA_4D = Storing_DATA.Storing_STABILITY_DATA_4D;
Storing_STABILITY_DATA_5 = Storing_DATA.Storing_STABILITY_DATA_5;

Geo_tier = Storing_GEO_DATA_1.Geo_tier;
Body_Geo = Storing_GEO_DATA_1.Body_Geo;
meshData = Storing_GEO_DATA_1.meshData;
Performance = Storing_AERO_DATA_1.Performance;
Weight_tier = Storing_WEIGHT_DATA_1.Weight_tier;
Prop_data = Storing_PROPULSION_DATA_1.Prop_data;

if OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY == 1
    VECTOR_XFLR5 = Storing_AERO_DATA_1.VECTOR_XFLR5;
    Design_criteria = Storing_AERO_DATA_1.Design_criteria;
    DATA_Ae = Storing_AERO_DATA_1.DATA_Ae;
    casos = Storing_AERO_DATA_1.casos;
    prefix = Storing_AERO_DATA_1.prefix;
    mark_legend = Storing_AERO_DATA_1.mark_legend;
    X_OC = Storing_AERO_DATA_1.X_OC;
    Aero_TH = Storing_AERO_DATA_1.Aero_TH;
    Aero = Storing_AERO_DATA_1.Aero;
    DATA_PL = Storing_AERO_DATA_1.DATA_PL;
    Performance = Storing_AERO_DATA_1.Performance;
end

if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim == 1
        TRIM_RESULTS = Storing_STABILITY_DATA_1.TRIM_RESULTS;
        Trim_ITER = Storing_STABILITY_DATA_1.Trim_ITER;
        Stab_Der = Storing_STABILITY_DATA_1.Stab_Der;
        Stab_Der_parts = Storing_STABILITY_DATA_1.Stab_Der_parts;
        Stab_Dyn_Long = Storing_STABILITY_DATA_1.Stab_Dyn_Long;
        Stab_Dyn_LatDir = Storing_STABILITY_DATA_1.Stab_Dyn_LatDir;
    end
    
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
        TRIM_RESULTS_var_V_XCG = Storing_STABILITY_DATA_2.TRIM_RESULTS_var_V_XCG;
        Trim_ITER_var_V_XCG = Storing_STABILITY_DATA_2.Trim_ITER_var_V_XCG;
        Restrictions_var_V_XCG = Storing_STABILITY_DATA_2.Restrictions_var_V_XCG;
        Stab_Der_var_V_XCG = Storing_STABILITY_DATA_2.Stab_Der_var_V_XCG;
        Stab_Der_parts_V_XCG = Storing_STABILITY_DATA_2.Stab_Der_parts_V_XCG;
        Stab_Dyn_Long_var_V_XCG = Storing_STABILITY_DATA_2.Stab_Dyn_Long_var_V_XCG;
        Stab_Dyn_LatDir_var_V_XCG = Storing_STABILITY_DATA_2.Stab_Dyn_LatDir_var_V_XCG;
        x_XCG_VAR = Storing_STABILITY_DATA_2.x_XCG_VAR;
        % VN diagram
        TRIM_RESULTS_var_V_XCG2 = Storing_STABILITY_DATA_2B.TRIM_RESULTS_var_V_XCG;
        Trim_ITER_var_V_XCG2 = Storing_STABILITY_DATA_2B.Trim_ITER_var_V_XCG;
        Restrictions_var_V_XCG2 = Storing_STABILITY_DATA_2B.Restrictions_var_V_XCG;
        Stab_Der_var_V_XCG2 = Storing_STABILITY_DATA_2B.Stab_Der_var_V_XCG;
        Stab_Der_parts_V_XCG2 = Storing_STABILITY_DATA_2B.Stab_Der_parts_V_XCG;
        Stab_Dyn_Long_var_V_XCG2 = Storing_STABILITY_DATA_2B.Stab_Dyn_Long_var_V_XCG;
        Stab_Dyn_LatDir_var_V_XCG2 = Storing_STABILITY_DATA_2B.Stab_Dyn_LatDir_var_V_XCG;
        Plot_options2 = Storing_STABILITY_DATA_2B.Plot_Options;
        Aero2 = Storing_AERO_DATA_1.Aero;        
        
    end

      if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_gamma == 1
        TRIM_RESULTS_var_V_gamma = Storing_STABILITY_DATA_2C.TRIM_RESULTS_var_V_gamma;
        Trim_ITER_var_V_gamma = Storing_STABILITY_DATA_2C.Trim_ITER_var_V_gamma;
        Restrictions_var_V_gamma = Storing_STABILITY_DATA_2C.Restrictions_var_V_gamma;
        Stab_Der_var_V_gamma = Storing_STABILITY_DATA_2C.Stab_Der_var_V_gamma;
        Stab_Der_parts_V_gamma = Storing_STABILITY_DATA_2C.Stab_Der_parts_V_gamma;
        Stab_Dyn_Long_var_V_gamma = Storing_STABILITY_DATA_2C.Stab_Dyn_Long_var_V_gamma;
        Stab_Dyn_LatDir_var_V_gamma = Storing_STABILITY_DATA_2C.Stab_Dyn_LatDir_var_V_gamma;
        x_XCG_VAR = Storing_STABILITY_DATA_2C.x_XCG_VAR;
        
      end


    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_var_XCG
        TRIM_RESULTS_var_XCG = Storing_STABILITY_DATA_3.TRIM_RESULTS_var_XCG;
        Trim_ITER_var_XCG = Storing_STABILITY_DATA_3.Trim_ITER_var_XCG;
        Stab_Der_var_XCG = Storing_STABILITY_DATA_3.Stab_Der_var_XCG;
        Stab_Der_parts_var_XCG = Storing_STABILITY_DATA_3.Stab_Der_parts_var_XCG;
        Stab_Dyn_Long_var_XCG = Storing_STABILITY_DATA_3.Stab_Dyn_Long_var_XCG;
        Stab_Dyn_LatDir_var_XCG = Storing_STABILITY_DATA_3.Stab_Dyn_LatDir_var_XCG;
    end

    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat == 1
        Trim_ITER_LAT = Storing_STABILITY_DATA_4A.Trim_ITER_LAT;
        conditions_TRIM_lat = Storing_STABILITY_DATA_4A.conditions_TRIM_lat;
    end
    
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat_asymmetries == 1
        Trim_ITER_LAT2 = Storing_STABILITY_DATA_4B.Trim_ITER_LAT2;
        conditions_TRIM_lat2 = Storing_STABILITY_DATA_4B.conditions_TRIM_lat;
    end

    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat_accelerations == 1
        Trim_ITER_LAT3 = Storing_STABILITY_DATA_4C.Trim_ITER_LAT3;
        conditions_TRIM_lat3 = Storing_STABILITY_DATA_4C.conditions_TRIM_lat;
    end

    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat_Trim_TAB == 1
        Trim_ITER_LAT4 = Storing_STABILITY_DATA_4D.Trim_ITER_LAT4;
        Trim_ITER_LAT4B = Storing_STABILITY_DATA_4D.Trim_ITER_LAT4B;
        Trim_ITER_LAT4C = Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C;
        Trim_ITER_LAT4D = Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D;
        conditions_TRIM_lat4 = Storing_STABILITY_DATA_4D.conditions_TRIM_lat;
    end
    
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Turning == 1
        Trim_ITER_LAT_Viraje = Storing_STABILITY_DATA_5.Trim_ITER_LAT_Viraje;
        conditions_TRIM_turning = Storing_STABILITY_DATA_5.conditions_TRIM_turning;
    end
   
end

if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
    if OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP == 1
        % Segments = Storing_PERFORMANCE_DATA_1{1};
        % handles = Storing_PERFORMANCE_DATA_1{2};
        % seg = Storing_PERFORMANCE_DATA_1{3};
        % Weights_AP = Storing_PERFORMANCE_DATA_1{4};
        % Total_Datos = Storing_PERFORMANCE_DATA_1{5};
        % datos = Storing_PERFORMANCE_DATA_1{6};
        Segments = Storing_PERFORMANCE_DATA_1.Segments;
        handles = Storing_PERFORMANCE_DATA_1.handles;
        seg = Storing_PERFORMANCE_DATA_1.seg;
        Weights_AP = Storing_PERFORMANCE_DATA_1.Weights_AP;
        Total_Datos = Storing_PERFORMANCE_DATA_1.Total_Datos;
        datos = Storing_PERFORMANCE_DATA_1.datos;
    end
    
    % Performance Analysis integrated with AP Codes varying Conditions
    if OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var == 1
        if OUTPUT_read_XLSX.STUDY_flags.Perfo1 == 1
            Storing_PERFORMANCE_21 = Storing_PERFORMANCE_DATA_21.Storing_PERFORMANCE_21;
            Plot_Options_21 = Storing_PERFORMANCE_DATA_21.Plot_Options;
            Restrictions_var_V_m = Storing_PERFORMANCE_DATA_21.Restrictions_var_V_m;
        end
        if OUTPUT_read_XLSX.STUDY_flags.Perfo7 == 1
            Storing_PERFORMANCE_27 = Storing_PERFORMANCE_DATA_27.Storing_PERFORMANCE_27;
            Plot_Options5 = Storing_PERFORMANCE_DATA_21.Plot_Options;
            Restrictions_var_V_m = Storing_PERFORMANCE_DATA_27.Restrictions_var_V_m;
        end
        % if  OUTPUT_read_XLSX.STUDY_flags.variable_speed_altitude_AP_Glide == 1
        %     Storing_PERFORMANCE_4 = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_4;
        %     Plot_Options6 = Storing_PERFORMANCE_4.Plot_Options;
        %     Restrictions_var_h_V = Storing_PERFORMANCE_4.Restrictions_var_h_V;
        % end
    end
end

%% Plots graphics
%% Defines the plots that will be plotted, and in which order
% Prints plots - Aero	plot_aero1
% Prints plots - Aero Polar	plot_aero2
% Print Plots for Estimation ofXAC	plot_aero3
% Print Plots for Propulsive Models	plot_prop1
% Prints plots for 3D	plot_3D
% Prints plots of SM analysis	plot_stability1
% Prints plots of longitudinal Trim	plot_stability2
% Prints plots of longitudinal Trim with Variable mass & Variable Speed	plot_stability3
% Prints plots of longitudinal Trim with Variable gamma & Variable Speed	plot_stability4
% Prints plots of lateral Trim	plot_stability5
% Prints plots of lateral Trim - asymmetries	plot_stability6
% Prints plots of lateral Trim - accelerations	plot_stability7
% Prints plots of lateral Trim - Trim Tabs	plot_stability8
% Prints plots of lateral Turning	plot_stability9
% Prints PLOTS STABILITY DERIVATIVES FOR VAR MASS AND VELOCITY	plot_stability10
% Prints PLOTS STABILITY ANALYSIS FOR VAR MASS AND VELOCITY	plot_stability11
% Prints PLOTS DYNAMIC STABILITY ANALYSIS AFTER IMPULSE	plot_stability12
% Prints plots of longitudinal Trim with V-n diagram n & Variable Speed	plot_stability13
% Prints PLOTS PERFORMANCE STUDY plot_perfo1
% Prints plots of Performance for Variable V and mass - Electric	plot_perfo2
% Prints plots of Performance for Variable h and V	plot_perfo3
% Prints plots of Performance Glide for Variable h and V	plot_perfo4
% Prints plots of Performance Glide for Variable h Max	plot_perfo5

% plot = OUTPUT_read_XLSX.PLOT_flags.plot;

%% Generates the PLOTS 
a1 = 1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_aero1; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_aero2; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_aero3; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_prop1; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_3D; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability1; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability2; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability3; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability4; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability5; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability6; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability7; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability8; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability9; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability10; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability11; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability12; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_stability13; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_perfo1; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_perfo2; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_perfo3; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_perfo4; a1=a1+1;
PLOTS(a1) = OUTPUT_read_XLSX.PLOT_flags.plot_perfo5; a1=a1+1;

PLOT_Case = zeros(1,length(PLOTS)); % dummy
for i=1:length(PLOTS)
    if PLOTS(i) == 1
        PLOT_Case(i) = i;
    end
end

% Initialize Variables
dummy = 0;
M_alpha_cero = 0;
V_alpha_cero = 0;

[mark_Type COLOR_scheme, Plot_Options] = Generate_Plot_Options(prefix,mark_legend,VECTOR_XFLR5,Plot_Options,OUTPUT_read_XLSX);

% Creates Folder For Plots
st1 = filenameS.filename_Plots;
st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st3 = strcat(st1,st2,'\');
% Verificar si la carpeta para `name` existe y crearla si no
folder = fileparts(st3); % Obtiene solo la ruta de la carpeta
if ~exist(folder, 'dir') % Verifica si la carpeta no existe
    mkdir(folder); % Crea la carpeta
end

filenameS.plots = st3;

for i=1:length(PLOTS)
    switch PLOT_Case(i) % PLOT 1
        case 1 % compare = 1 elements: w1,w2, vtp... alone
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY == 1
                [Fig] = plot_XFLR5_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix,OUTPUT_read_XLSX,AC_CONFIGURATION,filenameS);
            else
                Warning = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 2  % PLOT 2 % compare = 1 elements: w1,w2, vtp... alone
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY == 1
                [Fig] = plot_XFLR5_Polar_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix,OUTPUT_read_XLSX,AC_CONFIGURATION,filenameS);
            else
                Warning = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 3  % PLOT 3
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY == 1
%                 [XAC,Fig] = plot_Generate_XAC(DATA_Ae,Aero,DATA_PL,CASOS,degrees_XAC,Fig,...
%                     XFLR5_DATA,Performance,X_OC,Plot_Options,VECTOR_XFLR5);
                
            else
                Warning = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
         case 4  % PLOT 4
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
            
        case 5  % PLOT 5
             % [Fig] = plot_GEOMETRY_2022(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC,Performance,OUTPUT_read_XLSX,filenameS);
             [Fig] = plot_GEOMETRY_2022_v2(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC,Performance,OUTPUT_read_XLSX,filenameS);
            

        case 6 % PLOT 6
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
        case 7  % PLOT 7
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_var_XCG == 1
                    [Fig] = Generates_Plots_Longitudinal_Trim(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var_XCG,Trim_ITER_var_XCG,Geo_tier,Plot_Options,OUTPUT_read_XLSX,Fig,filenameS);
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
            
        case 8  % PLOT 8
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
                    [Fig,M_alpha_cero,V_alpha_cero] = Generates_Plots_Longitudinal_Trim_VAR(TRIM_RESULTS_var_V_XCG,Trim_ITER_var_V_XCG,...
                        Restrictions_var_V_XCG,Stab_Der_var_V_XCG,Geo_tier,Plot_Options,Fig,OUTPUT_read_XLSX,x_XCG_VAR,filenameS);
                    % [Fig,M_alpha_cero,V_alpha_cero] = Generates_Plots_Longitudinal_Trim_VAR2(TRIM_RESULTS_var_V_XCG2,Trim_ITER_var_V_XCG2,...
                        % Restrictions_var_V_XCG2,Stab_Der_var_V_XCG2,Geo_tier,Plot_Options,Plot_options2,Fig,OUTPUT_read_XLSX,Aero2);
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

        case 9  % PLOT 9
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_gamma == 1
                    [Fig,M_alpha_cero,V_alpha_cero] = Generates_Plots_Longitudinal_Trim_VAR_gamma(TRIM_RESULTS_var_V_gamma,Trim_ITER_var_V_gamma,...
                        Restrictions_var_V_gamma,Stab_Der_var_V_gamma,Geo_tier,Plot_Options,Fig,OUTPUT_read_XLSX,x_XCG_VAR,filenameS,conv_UNITS);
                    % [Fig,M_alpha_cero,V_alpha_cero] = Generates_Plots_Longitudinal_Trim_VAR2(TRIM_RESULTS_var_V_XCG2,Trim_ITER_var_V_XCG2,...
                    % Restrictions_var_V_XCG2,Stab_Der_var_V_XCG2,Geo_tier,Plot_Options,Plot_options2,Fig,OUTPUT_read_XLSX,Aero2);
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
            
        case 10  % PLOT 10
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat == 1
                    [Fig] = Generates_Plots_Lateral_Trim(Trim_ITER_LAT,Geo_tier,...
                        Plot_Options,conv_UNITS,conditions_TRIM_lat,OUTPUT_read_XLSX,Fig,filenameS);
%                 elseif OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat_asymmetries_accelerations == 1
%                     % If Special cases, asymmetry and accelerations
%                     [Fig] = Generates_Plots_Lateral_Trim2(Trim_ITER_LAT2,Trim_ITER_LAT3,Geo_tier,...
%                         Plot_Options,conv_UNITS,conditions_TRIM_lat,OUTPUT_read_XLSX,Fig);
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

        case 11  % PLOT 11
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat_asymmetries == 1
                    % If Special cases, asymmetry and accelerations
                    [Fig] = Generates_Plots_Lateral_Trim2(Trim_ITER_LAT2,Geo_tier,...
                        Plot_Options,conv_UNITS,conditions_TRIM_lat2,OUTPUT_read_XLSX,Fig,filenameS);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Trim_lat_asymmetries == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end

        case 12  % PLOT 12
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat_accelerations == 1
                    % If Special cases, asymmetry and accelerations
                    [Fig] = Generates_Plots_Lateral_Trim3(Trim_ITER_LAT3,Geo_tier,...
                        Plot_Options,conv_UNITS,conditions_TRIM_lat3,OUTPUT_read_XLSX,Fig,filenameS);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Trim_lat_asymmetries == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end

        case 13  % PLOT 13
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat_Trim_TAB == 1
                    % If Special cases, asymmetry and accelerations
                    % [Fig] = Generates_Plots_Lateral_Trim4(Trim_ITER_LAT4,Geo_tier,...
                    %     Plot_Options,conv_UNITS,conditions_TRIM_lat4,OUTPUT_read_XLSX,Fig);
                    % [Fig] = Generates_Plots_Lateral_Trim4B(Trim_ITER_LAT4B,Geo_tier,...
                    %     Plot_Options,conv_UNITS,conditions_TRIM_lat4,OUTPUT_read_XLSX,Fig);
                    % [Fig] = Generates_Plots_Lateral_Trim4C(Trim_ITER_LAT4C,Geo_tier,...
                    %     Plot_Options,conv_UNITS,conditions_TRIM_lat4,OUTPUT_read_XLSX,Fig);
                    [Fig] = Generates_Plots_Lateral_Trim4D(Trim_ITER_LAT4D,Geo_tier,...
                        Plot_Options,conv_UNITS,conditions_TRIM_lat4,OUTPUT_read_XLSX,Fig,filenameS);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Trim_lat_asymmetries == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end
            
        case 14  % PLOT 14
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Turning == 1
                    [Fig] = Generates_Plots_Lateral_Turning(Trim_ITER_LAT_Viraje,Geo_tier,...
                        Plot_Options,conv_UNITS,conditions_TRIM_turning,OUTPUT_read_XLSX,Fig,filenameS);
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
            
        case 15  % PLOT 15
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Stability_Derivatives_varV_m0 == 1
                    [Fig] = Generates_Plots_Derivatives_VAR(Stab_Der_var_V_XCG,Geo_tier,Plot_Options,OUTPUT_read_XLSX,Fig,filenameS);
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

        case 16 % PLOT 16
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
                    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Stability_Analysis_varV_m0_long == 1
                        [Fig] = Generates_Plots_StabilityAnalysis_long_VAR(Geo_tier,Plot_Options,Stab_Dyn_Long_var_V_XCG,OUTPUT_read_XLSX,Fig,filenameS);
                    end
                    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Stability_Analysis_varV_m0_lat == 1
                        [Fig] = Generates_Plots_StabilityAnalysis_lat_VAR(Geo_tier,Plot_Options,Stab_Dyn_LatDir_var_V_XCG,OUTPUT_read_XLSX,Fig,filenameS);
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

        case 17  % PLOT 17
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                %                 if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Stability_Analysis_dyna_impulse_long == 1
                    [Fig] = Generate_Plots_longitudinal_dynamic(Plot_Options,OUTPUT_read_XLSX,conv_UNITS,Fig,Stab_Dyn_Long,filenameS);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Stability_Analysis_dyna_impulse_long == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Stability_Analysis_dyna_impulse_lat == 1
                    [Fig] = Generate_Plots_lateral_dynamic(Plot_Options,OUTPUT_read_XLSX,conv_UNITS,Fig,Stab_Dyn_LatDir,filenameS);
                else
                    Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY_Stability_Analysis_dyna_impulse_lat == 1  - Code in PAUSE';
                    disp(Warning)
                    pause
                end
                %             end
                
%             end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that STABILITY_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end

          case 18 % PLOT 18
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
                if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
                    [Fig,M_alpha_cero,V_alpha_cero] = Generates_Plots_Longitudinal_Trim_VAR2(TRIM_RESULTS_var_V_XCG2,Trim_ITER_var_V_XCG2,...
                        Restrictions_var_V_XCG2,Stab_Der_var_V_XCG2,Geo_tier,Plot_Options,Plot_options2,Fig,OUTPUT_read_XLSX,Aero2,filenameS,conv_UNITS);

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

        case 19  % PLOT 19
            % Logic ensures that the studies have been conducted
            if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
                SEGMENTS = OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF;
                SEGMENTS_STYPE = OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF;

                if SEGMENTS == 5 && SEGMENTS_STYPE == 7
                
                else
                    % Generates_Plots_Performance_v3(datos,OUTPUT_read_XLSX,Fig,Plot_Options,filenameS);
                    Generates_Plots_Performance_v3B(datos,OUTPUT_read_XLSX,Fig,Plot_Options,filenameS);
                end

            else
                Warning = 'WARNING!!! No results generated, please Make sure that PERFORMANCE_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end

        case 20 % PLOT 21
            if OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var == 1
                if OUTPUT_read_XLSX.STUDY_flags.Perfo1 == 1
                    [Fig] = Generates_Plots_Performance_v4(Storing_PERFORMANCE_21,Plot_Options,Plot_Options_21,OUTPUT_read_XLSX,Restrictions_var_V_m,Fig,filenameS);
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that PERFORMANCE_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end

        case 21 % PLOT 21
            if OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var == 1
                if OUTPUT_read_XLSX.STUDY_flags.Perfo3 == 1
                    Storing_PERFORMANCE_23 = Storing_DATA.Storing_PERFORMANCE_DATA_23;
                    Restrictions_var_h_V = Storing_DATA.Storing_PERFORMANCE_DATA_23.Restrictions_var_h_V;
                    [Fig] = Generates_Plots_Performance_v5B(Storing_PERFORMANCE_23,Plot_Options,OUTPUT_read_XLSX,Restrictions_var_h_V,Fig,filenameS);
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that PERFORMANCE_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end

        case 22  % PLOT 22
            if OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var == 1
                if OUTPUT_read_XLSX.STUDY_flags.Perfo4 == 1
                    Storing_PERFORMANCE_24 = Storing_DATA.Storing_PERFORMANCE_DATA_24;
                    Restrictions_var_h_V = Storing_DATA.Storing_PERFORMANCE_DATA_24.Restrictions_var_h_V;
                    [Fig] = Generates_Plots_Performance_v5(Storing_PERFORMANCE_24,Plot_Options,OUTPUT_read_XLSX,Restrictions_var_h_V,Fig,filenameS);
                end
            else
                Warning = 'WARNING!!! No results generated, please Make sure that PERFORMANCE_STUDY == 1  - Code in PAUSE';
                disp(Warning)
                pause
            end

        case 23  % PLOT 23
               
    end
end