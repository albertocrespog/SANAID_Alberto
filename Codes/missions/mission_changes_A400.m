function [conditions OUTPUT_read_XLSX] = mission_changes_A400(OUTPUT_read_XLSX,mission_actual,conv_UNITS);   
           
%Loadsa de geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

switch mission_actual
    case 1 % ID_01
        % Geometry
        % Aerodynamics
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_01_A400_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        OUTPUT_read_XLSX.PLOT_flags.fname = 'results\A400';

    case 2 % ID_01
        % Geometry
        % Aerodynamics
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr*1.25;
        % Plots prefix
        prefix = 'CASO_ID_01_A4000_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        OUTPUT_read_XLSX.PLOT_flags.fname = 'results\A400';
end