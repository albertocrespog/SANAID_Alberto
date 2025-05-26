function [conditions OUTPUT_read_XLSX] = mission_changes_PEPINO_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;
% Folder where to stores de Figures
m2ft = conv_UNITS.m2ft;

dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

%% Selects engine on or off
% Posicion_Palanca = 0; % Engine OFF
% Posicion_Palanca = 1; % Engine ON
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

%% Initialize Stability Studies
conditions = initialization_STABILITY_Analysis(conditions,conv_UNITS);

OUTPUT_read_XLSX.Fuselage_flags.Use_Storing_DATA == 1;

switch mission_actual
    case 1 % ID_01
        % Geometry
        % Aerodynamics
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_01_PEPINO_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 2 % ID_01
        % Geometry
        % Aerodynamics
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr*1.25;
        % Plots prefix
        prefix = 'CASO_ID_02_PEPINO_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
end