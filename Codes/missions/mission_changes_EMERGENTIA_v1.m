function [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS)
             
%Loadsa de geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

%% Selects engione on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

filename = '../Results/EMERGENTIA/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

%% Initialize Stability Studies
conditions = initialization_STABILITY_Analysis(conditions,conv_UNITS)

switch mission_actual
    case 0 % 
        % Geometry
        % Aerodynamics
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_00_EMERGENTIA_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
                
    case 1 % ID_01
        % Geometry
%          OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.06; %1.0716
%          OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453; %
        % Aerodynamics
%         OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 3; %
%         OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 3; %
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27; 
        % Plots prefix
        
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = 0*pi/180;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45*pi/180;        
        prefix = 'CASO_ID_01_EMERGENTIA_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
          
    case 2 % ID_02
        % Geometry
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.06; %1.0716
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = -0.1; %1.0716
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 + 0.005; %
        % Aerodynamics
%         OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 3; %
%         OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 3; %
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_02_EMERGENTIA_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

    case 3 % ID_03
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.06; %1.0716
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 + 0.01; %
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 3; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 3; %
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;        
        % Plots prefix
        prefix = 'CASO_ID_03_EMERGENTIA_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 4 % ID_04
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.06; %1.0716
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 + 0.02; %
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 3; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 3; %
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;      
        % Plots prefix
        prefix = 'CASO_ID_04_EMERGENTIA_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
end



