function [conditions OUTPUT_read_XLSX] = mission_changes_MILVUS2_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
m2ft = conv_UNITS.m2ft;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;

% Identifies folder where to store Figures
filename = '../Results/29_MILVUS2/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

%% Selects engine on or off
% Posicion_Palanca = 0; % Engine OFF
% Posicion_Palanca = 1; % Engine ON
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

%% Initialize Stability Studies
conditions = initialization_STABILITY_Analysis(conditions,conv_UNITS)

%% Studies
OUTPUT_read_XLSX.STUDY_flags.AERODYNAMIC_STUDY = 1;
OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY = 0; 
OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY = 1;
OUTPUT_read_XLSX.STUDY_flags.MISSIONS_STUDY = 1;
OUTPUT_read_XLSX.Fuselage_flags.Use_Storing_DATA == 1;

%% Switching CASES
switch mission_actual
    case 1 % ID_01
        
        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 9;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 8;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 0;

        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_vee = -2;
        
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 1;
        % Uses estimation of POLAR using CD0 from CFD and K1 and K2 from FLOW5
        fuse_aero_FLOW_and_CBM = 2;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %

        %% Plots prefix
        prefix = 'CASO_ID_01_MILVUS2_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        % OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        % Fuselage contribution 
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 1;
        
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.481211305490165;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.125;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e = (45*D2R);
end