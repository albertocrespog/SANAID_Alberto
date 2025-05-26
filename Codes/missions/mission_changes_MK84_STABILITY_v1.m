function [conditions OUTPUT_read_XLSX] = mission_changes_MK84_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
             
%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;
% Folder where to stores de Figures
m2ft = conv_UNITS.m2ft;

%% Selects engine on or off
% Posicion_Palanca = 0; % Engine OFF
% Posicion_Palanca = 1; % Engine ON
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

%% Initialize Stability Studies
conditions = initialization_STABILITY_Analysis(conditions,conv_UNITS);

switch mission_actual
    case 1 % ID_01

        %% Aerodynamic

        %% Plots prefix
        prefix = 'CASEID_01_AC_24_MK84_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Geometry

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45; 
end

