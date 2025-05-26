function Storing_PERFORMANCE_DATA_1 = Conduct_Performance_Optimization_v2(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
    OUTPUT_read_XLSX,Storing_AERO_DATA_1,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_PROPULSION_DATA,filenameS)

Aero = Storing_AERO_DATA_1.Aero;
Aero_TH = Storing_AERO_DATA_1.Aero_TH;
Performance = Storing_AERO_DATA_1.Performance;
Weight_tier = Storing_WEIGHT_DATA_1.Weight_tier;
Geo_tier = Storing_GEO_DATA_1.Geo_tier;

type_missions_WF        = OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF;
sub_type_missions_WF    = OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF;
num_missions_WF         = OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF;
FlightMODE_WF = sub_type_missions_WF;

Segments = Define_Segments_v2(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF,OUTPUT_read_XLSX);

% Weight_tier
if OUTPUT_read_XLSX.STUDY_flags.STUDY_Weight_Fraction == 1
    [seg_WF] = Generation_Mission_Segments_v3(conv_UNITS,num_missions_WF,type_missions_WF,Performance,case_AC,Segments,FlightMODE_WF);
%     Weights = Weight_Fraction(seg_WF,Propulsion,conv_UNITS,Aero_TH,type_missions_WF,Weight_tier,Geo_tier,AC_CONFIGURATION);
end

type_missions        = OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF;
sub_type_missions    = OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF;
num_missions         = OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF;

FlightMODE_var = sub_type_missions;

%% Propulsive selected data

propul = OUTPUT_read_XLSX.Propulsive_flags.propul;
prop_data(1) = propul(1); % Type of engine
prop_data(2) = propul(2); % Number of engines
prop_data(3) = propul(3); % Thrust (lbf) or Power (shp) per engine
prop_data(4) = propul(4); % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
prop_data(5) = propul(7); % Normativa
prop_data(6) = propul(6); % Prop efficiency
prop_data(7) = propul(5); % By-pass
Prop_sel = prop_data;

% Conducts Performance studies
handles = caracteristicas_avanz_v3(Aero_TH,Prop_data,Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Aero,OUTPUT_read_XLSX);

% Generates the different segments
[seg] = Generation_Mission_Segments_v3(conv_UNITS,num_missions,type_missions,case_AC,Segments,FlightMODE_var);
% Analyse the mission
% [Weights_AP,Total_Datos,datos] = procesar_mision_v3(seg,num_missions,handles,Segments,FlightMODE_var,Geo_tier,Prop_data,AC_CONFIGURATION,conv_UNITS,OUTPUT_read_XLSX);
[Weights_AP,Total_Datos,datos] = procesar_mision_v4(seg,num_missions,handles,FlightMODE_WF,Geo_tier,Prop_data,AC_CONFIGURATION);

OUTPUT_read_XLSX.Propulsive_flags.seg = seg;

CI = OUTPUT_read_XLSX.Propulsive_flags.CI; %0.3; % kg/s
density_fuel = OUTPUT_read_XLSX.Propulsive_flags.CI; %  0.809; % kg/l
cost_fuel = OUTPUT_read_XLSX.Propulsive_flags.cost_fuel; %161.69;%  cts/gal - 7 Feb 2020
l2gal = conv_UNITS.l2gal;
m2km = conv_UNITS.m2km;
kg2Tm = conv_UNITS.kg2Tm;
nm2m = conv_UNITS.nm2m; %1852; % nautical mile/meter
C_fuel = cost_fuel*l2gal/density_fuel; % coste combustible centavos/kg
%         cost_fuel_cent/density_fuel/lb2kg; % coste combustible centavos/kg
%         Cost_fuel = (Total_Datos.fuel_total/density_fuel)*l2gal*cost_fuel;
DOC_PLOT = (Total_Datos.tiempo_total*CI + Total_Datos.fuel_total)*C_fuel;
%         ASM = Weights_AP.m_p*kg2Tm * Total_Datos.distancia_total * (1/nm2m);
ASM_PLOT = (Total_Datos.distancia_total*(1/nm2m))*(Weights_AP.m_p*kg2Tm);
CAPM_PLOT = DOC_PLOT/ASM_PLOT;

Storing_PERFORMANCE_DATA_1.datos.TOTAL.C_fuel = C_fuel;
Storing_PERFORMANCE_DATA_1.datos.TOTAL.DOC_PLOT = DOC_PLOT;
Storing_PERFORMANCE_DATA_1.datos.TOTAL.ASM_PLOT = ASM_PLOT;
Storing_PERFORMANCE_DATA_1.datos.TOTAL.CAPM_PLOT = CAPM_PLOT;

% Does not rpresnt results in screen if multimission
MISSIONS_STUDY = OUTPUT_read_XLSX.STUDY_flags.MISSIONS_STUDY; % Conducts Mission Studies
if MISSIONS_STUDY == 1
    % Represents on screen

    Time_msg = strcat('Total Time = ',num2str(Total_Datos.tiempo_total/60),' (min)');
    Distance_msg = strcat('Total Distance = ',num2str(Total_Datos.distancia_total/1000),' (km)');
    Wempty_msg = strcat('Initial Weight = ',num2str(Weights_AP.m_TOW),' (kg)');
    Battery_msg = strcat('Bateries Mass = ',num2str(Weights_AP.m_bat),' (kg)');
    W_final_msg = strcat('Final Weight = ',num2str(Weights_AP.m_F),' (kg)');
    Speed_msg = strcat('Cruise speed = ',num2str(OUTPUT_read_XLSX.IPP_flags.V_cr),' (m/s)');
    separate_msg = strcat('-----------------------------------------------------------------');
   
    disp(Time_msg)
    disp(Distance_msg)
    disp(Wempty_msg)
    disp(Battery_msg)
    disp(W_final_msg)
    disp(Speed_msg)
    disp(separate_msg)

else
    % Resultados globales:
    if handles.propul(1) == 4
        uiwait(msgbox({['Total Time = ',num2str(Total_Datos.tiempo_total/60),' (min)']...
            ['Total Distance = ',num2str(Total_Datos.distancia_total/1000),' (km)']...
            ['Initial Weight = ',num2str(Weights_AP.m_TOW),' (kg)']...
            ['Total Fuel = ','N/A']...
            ['Bateries Mass = ',num2str(Weights_AP.m_bat),' (kg)']...
            ['Final Weight = ',num2str(Weights_AP.m_F),' (kg)']...
            ['Fuel Cost = ','N/A']...
            ['DOC = ',num2str(DOC_PLOT),' (cts$)']...
            ['CAPM = ',num2str(CAPM_PLOT),' (cts$/km ton)']}));
    else
        uiwait(msgbox({['Total Time = ',num2str(Total_Datos.tiempo_total/60),' (min)']...
            ['Total Distance = ',num2str(Total_Datos.distancia_total/1000),' (km)']...
            ['Total Distance = ',num2str(Total_Datos.distancia_total/1000),' (km)']...
            ['Initial Weight = ',num2str(Weights_AP.m_TOW),' (kg)']...
            ['Total Fuel = ',num2str(Weights_AP.m_f),' (kg)']...
            ['Bateries Mass = ','N/A']...
            ['Final Weight = ',num2str(Weights_AP.m_F),' (kg)']...
            ['Fuel Cost = ',num2str(cost_fuel),' (cts$)']...
            ['DOC = ',num2str(DOC_PLOT),' (cts$)']...
            ['CAPM = ',num2str(CAPM_PLOT),' (cts$/km ton)']}));
    end
end


%         if Post_processing_PERFORMANCE == 1
%             l2gal = conv_UNITS.l2gal;
%             m2km = conv_UNITS.m2km;
%             kg2Tm = conv_UNITS.kg2Tm;
%             % Generation of variable to plot
%             for i=1:N_V_VAR_perf
%                 for j=1:N_Wp_VAR_perf
%
%                     m_f_PLOT(i,j) = Weights_AP_var{i,j}.m_f;
%                     tiempo_total_PLOT(i,j) = tiempo_total_var{i,j};
%                     distancia_total_PLOT(i,j) = distancia_total_var{i,j};
%                     Cost_fuel(i,j) = (m_f_PLOT(i,j)/density_fuel)*l2gal*cost_fuel;
%                     DOC_PLOT(i,j) = (tiempo_total_PLOT(i,j)*CI + Cost_fuel(i,j));
%                     CAPM_PLOT(i,j) = DOC_PLOT(i,j)/((distancia_total_PLOT(i,j)*m2km)*(Wp_VAR_perf(j)*kg2Tm));
%                 end
%             end
%
%             for j=1:length(datos_var{i}.tiempo)
%                 distancia_PLOT(j) = datos_var{i}.distancia(j);
%             end
%
%             % Generation of variable to plot for Cruise
%             for i=1:N_V_VAR_perf
%                 for j=1:length(datos_var{i,N_Wp_plot}.L_D)
%                     L_D_PLOT(i,j) = datos_var{i,N_Wp_plot}.L_D(j);
%                     tiempo_PLOT(i,j) = datos_var{i,N_Wp_plot}.tiempo(j);
% %                     CL_PLOT(i,j)  = datos_var{i,N_Wp_plot}.subida.CL(j);
%                     palanca_PLOT(i,j) = datos_var{i,N_Wp_plot}.palanca(j);
%                 end
%             end
%
%             Plots_performance.m_f_PLOT = m_f_PLOT;
%             Plots_performance.tiempo_total_PLOT = tiempo_total_PLOT;
%             Plots_performance.distancia_total_PLOT = distancia_total_PLOT;
%             Plots_performance.Cost_fuel = Cost_fuel;
%             Plots_performance.DOC_PLOT = DOC_PLOT;
%             Plots_performance.CAPM_PLOT = CAPM_PLOT;
%             Plots_performance.distancia_PLOT = distancia_PLOT;
%             Plots_performance.L_D_PLOT = L_D_PLOT;
%             Plots_performance.tiempo_PLOT = tiempo_PLOT;
% %             Plots_performance.CL_PLOT = CL_PLOT;
%             Plots_performance.palanca_PLOT = palanca_PLOT;
%             Plots_performance.CI = CI;
%         else
%             % Dummy Variable
%             Plots_performance = 0;
%         end
%     end

% Storing Data
% Storing_PERFORMANCE_DATA_1.Segments = Segments;
% Storing_PERFORMANCE_DATA_1.handles = handles;
% Storing_PERFORMANCE_DATA_1.Weights_AP = Weights_AP;
% Storing_PERFORMANCE_DATA_1.Total_Datos = Total_Datos;
% Storing_PERFORMANCE_DATA_1.datos = datos;
% Storing_PERFORMANCE_DATA_1.Storing_PROPULSION_DATA = Storing_PROPULSION_DATA;
% Storing_PERFORMANCE_DATA_1.seg = seg;


if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_PERFORMANCE_DATA_1 = Saving_data_Performance_Optimization_v2(Segments,handles,Weights_AP,Total_Datos,datos,Storing_PROPULSION_DATA,seg,filenameS,OUTPUT_read_XLSX);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_PERFORMANCE_DATA_1 = dummy;
end

% Storing_PERFORMANCE_DATA_1{1} = Segments;
% Storing_PERFORMANCE_DATA_1{2} = handles;
% Storing_PERFORMANCE_DATA_1{3} = seg;
% Storing_PERFORMANCE_DATA_1{4} = Weights_AP;
% Storing_PERFORMANCE_DATA_1{5} = Total_Datos;
% Storing_PERFORMANCE_DATA_1{6} = datos;

% Storing_PERFORMANCE_DATA_1{3} = Weights_AP_var;
% Storing_PERFORMANCE_DATA_1{4} = fuel_total_var;
% Storing_PERFORMANCE_DATA_1{5} = tiempo_total_var;
% Storing_PERFORMANCE_DATA_1{6} = distancia_total_var;
% Storing_PERFORMANCE_DATA_1{7} = W_var;
% Storing_PERFORMANCE_DATA_1{8} = datos_var;
% Storing_PERFORMANCE_DATA_1{9} = Plots_performance;
% Storing_PERFORMANCE_DATA_1{1} = Post_processing_PERFORMANCE;



