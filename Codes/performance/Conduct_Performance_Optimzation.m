function Storing_PERFORMANCE_DATA_1 = Conduct_Performance_Optimization

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
    
% Storing Fata
Storing_PERFORMANCE_DATA_1{1} = Segments;
Storing_PERFORMANCE_DATA_1{2} = handles;
Storing_PERFORMANCE_DATA_1{3} = Weights_AP_var;
Storing_PERFORMANCE_DATA_1{4} = fuel_total_var;
Storing_PERFORMANCE_DATA_1{5} = tiempo_total_var;
Storing_PERFORMANCE_DATA_1{6} = distancia_total_var;
Storing_PERFORMANCE_DATA_1{7} = W_var;
Storing_PERFORMANCE_DATA_1{8} = datos_var;
Storing_PERFORMANCE_DATA_1{9} = Plots_performance;
Storing_PERFORMANCE_DATA_1{1} = Post_processing_PERFORMANCE;



