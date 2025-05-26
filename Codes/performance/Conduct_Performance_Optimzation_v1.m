function Storing_PERFORMANCE_DATA_1 = Conduct_Performance_Optimzation_v1(Aero,Aero_TH,case_AC,conv_UNITS,OUTPUT_read_XLSX,...
    Performance,Propulsion,Weight_tier,Geo_tier,AC_CONFIGURATION,Prop_data)

%% Weight Fraction
% Enter the number of mission segments
type_missions_WF        = OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF;
sub_type_missions_WF    = OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF;
num_missions_WF         = OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF;
%     climb_mode = OUTPUT_read_XLSX.PerforMisionSelection_flags.climb_mode;
%     cruise_mode = OUTPUT_read_XLSX.PerforMisionSelection_flags.cruise_mode;
%     turn_mode = OUTPUT_read_XLSX.PerforMisionSelection_flags.turn_mode;
%     descent_mode = OUTPUT_read_XLSX.PerforMisionSelection_flags.descent_mode;

%     FlightMODE_WF.climb_mode = climb_mode;
%     FlightMODE_WF.cruise_mode = cruise_mode;
%     FlightMODE_WF.turn_mode = turn_mode;
%     FlightMODE_WF.descent_mode = descent_mode;
FlightMODE_WF = sub_type_missions_WF;

Segments = Define_Segments_v2(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE_WF,OUTPUT_read_XLSX);

% Weight_tier
if OUTPUT_read_XLSX.STUDY_flags.STUDY_Weight_Fraction == 1
    [seg_WF] = Generation_Mission_Segments_v3(conv_UNITS,num_missions_WF,type_missions_WF,Performance,case_AC,Segments,FlightMODE_WF);
    Weights = Weight_Fraction(seg_WF,Propulsion,conv_UNITS,Aero_TH,type_missions_WF,Weight_tier,Geo_tier,AC_CONFIGURATION);
end

%% Performance Variable
type_missions        = OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF;
sub_type_missions    = OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF;
num_missions         = OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF;

FlightMODE_var = sub_type_missions_WF;

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

handles = caracteristicas_avanz_v3(Aero_TH,Prop_data,Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Aero,OUTPUT_read_XLSX);

% Generates the different segments
[seg] = Generation_Mission_Segments_v3(conv_UNITS,num_missions,type_missions,case_AC,Segments,FlightMODE_var);
% Analyse the mission
[Weights_AP,Total_Datos,datos] = procesar_mision_v3(seg,num_missions,handles,Segments,FlightMODE_var,Geo_tier,Prop_data,AC_CONFIGURATION);

% Cálculo de Datos de Performance
CI = OUTPUT_read_XLSX.Propulsive_flags.CI; % kg/s
density_fuel = OUTPUT_read_XLSX.Propulsive_flags.density_fuel; % kg/l
cost_fuel = OUTPUT_read_XLSX.Propulsive_flags.cost_fuel;%  cts/gal - 7 Feb 2020
l2gal = conv_UNITS.l2gal;
m2km = conv_UNITS.m2km;
kg2Tm = conv_UNITS.kg2Tm;
Cost_fuel = (Total_Datos.fuel_total/density_fuel)*l2gal*cost_fuel;
DOC_PLOT = (Total_Datos.tiempo_total*CI + Cost_fuel);
CAPM_PLOT = DOC_PLOT/((Total_Datos.distancia_total*m2km)*(Weights_AP.m_p*kg2Tm));

% Resultados globales:
uiwait(msgbox({['Total Time = ',num2str(Total_Datos.tiempo_total/60),' (min)']...
    ['Total Distance = ',num2str(Total_Datos.distancia_total/1000),' (km)']...
    ['Initial Weight = ',num2str(Weights_AP.m_TOW),' (kg)']...
    ['Total Fuel = ',num2str(Weights_AP.m_f),' (kg)']...
    ['Final Weight = ',num2str(Weights_AP.m_F),' (kg)']...
    ['Fuel Cost = ',num2str(Cost_fuel),' (cts$)']...
    ['DOC = ',num2str(DOC_PLOT),' (cts$)']...
    ['CAPM = ',num2str(CAPM_PLOT),' (cts$/km ton)']}));

% Storing Fata
Storing_PERFORMANCE_DATA_1{1} = Segments;
Storing_PERFORMANCE_DATA_1{2} = handles;
Storing_PERFORMANCE_DATA_1{3} = Weights_AP;
Storing_PERFORMANCE_DATA_1{4} = Total_Datos;
Storing_PERFORMANCE_DATA_1{5} = datos;
Storing_PERFORMANCE_DATA_1{6} = Cost_fuel;
Storing_PERFORMANCE_DATA_1{7} = DOC_PLOT;
Storing_PERFORMANCE_DATA_1{8} = CAPM_PLOT;
