function [conditions OUTPUT_read_XLSX] = mission_changes_FR_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
             
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
        
        %% Plots prefix
        prefix = 'CASEID_01_AC_27_FR_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Geometry
        i_w1 = -2;
        i_vee = -5;

        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = i_w1; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_vee = i_vee; %

        % Performance
        % Porpulsion
        Thrust = 176*0.5;
        OUTPUT_read_XLSX.Propulsive_flags.propul(3) = Thrust; % EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine

        %% Stability
        % OUTPUT_read_XLSX.Stability_flags.V_high = 45;

    case 2 % ID_01
        
        %% Plots prefix
        prefix = 'CASEID_02_AC_27_FR_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Geometry
        i_w1 = 2;
        i_vee = -2;

        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = i_w1; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_vee = i_vee; %

        % Performance
        % Porpulsion
        Thrust = 176*0.5;
        OUTPUT_read_XLSX.Propulsive_flags.propul(3) = Thrust; % EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine

        %% Stability
        % OUTPUT_read_XLSX.Stability_flags.V_high = 45;
      
end


