function [conditions OUTPUT_read_XLSX] = mission_changes_Pepino(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);   
           
%Loadsa de geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

filename = '../Results/PEPINO/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

%% Initialize Stability Studies
conditions = initialization_STABILITY_Analysis(conditions,conv_UNITS)

switch mission_actual
    case 1 % ID_01
        % Geometry
        % Aerodynamics

        
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_01_PEPINO_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 2 % ID_01
        % Geometry
        % Aerodynamics
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr*1.25;
        % Plots prefix
        prefix = 'CASO_ID_02_PEPINO_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
end