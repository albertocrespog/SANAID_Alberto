function [conditions OUTPUT_read_XLSX] = mission_changes_BAT_July2023(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

% Weight_mission = T120_weightconfiguration;
% conditions.m_TOW = Weight_mission.W_vec_f(1);

% Percentage where bayoneta is located
pp_bayoneta = (1100.55-966.57)/469;
cR_w1_T75 = 0.421;
x_loc_1R_y1_w1_CAD_T75 = 1.45766;
x_loc_1R_y1_w1_CAD_T75_bayoneta = x_loc_1R_y1_w1_CAD_T75 + pp_bayoneta*cR_w1_T75;
x_loc_1R_y1_w2_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;


%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;


% Conversion factor depending on the units of the CG_TARSIS120 (1 for m 1000 for mm)
% conv_factor = 1 old version
% conv_factor = 1000 new version
conv_factor = 1000;
OUTPUT_read_XLSX.Weights_flags.conv_factor = conv_factor;

% Extra manual modification of Xcg
Dxcg = Modifications.Dxcg;

%% Correction of geometry according to CATIA - AERTEC - Junio 2022
% Wing Corrections
% OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 1527.565/1000;
% 
% OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 1527.565/1000;
% OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 1527.565/1000;
% 
% OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_w1_CAD = -60.865/1000;
% 
% % HTP Corrections
% OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD = 3479.874/1000;
% OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_w2_CAD = 96.47/1000;
% % VTP1 Corrections
% OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = 3392.884/1000;
% OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = -92.192/1000;
% % VTP2 Corrections
% OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP2_CAD = 3391.59/1000;
% OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP2_CAD = -95.315/1000;

% Geometry
% OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 50/100;
OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.469;
OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.469;
% DeltaCR = OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 - cR_w1_T75;
% location_XAC_wing = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD + OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1/4;
% OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD_T75_bayoneta - pp_bayoneta*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;


%% Correction of geometry according to CATIA - AERTEC - 28 Junio 2022
% Wing Corrections
OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 1447.697/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD = 180.678/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 3109.99/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_w1_CAD = -62.394/1000;
% HTP Corrections
OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD = 3559.712/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = 0/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = 662.382/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_w2_CAD = 90.467/1000;
% VTP1 Corrections
OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = 674.343/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = -96.953/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = 381.213/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = 3472.882/1000;
% VTP2 Corrections
OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = -674.343/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = -96.953/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = 381.213/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = 3472.882/1000;

% Engine configuration SP210
OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1 = 2672.835/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1 = 0.014/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1 = 34.41/1000;

OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar1 = 2672.835/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar1 = 0.014/1000;
OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar1 = 34.41/1000;

% % Engine configuration DA215
% OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1 =  -2648.012/1000;
% OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1 = 0.014/1000;
% OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1 = 34.41/1000;

%  of control surface aileron (chordwise)	cf_ail
OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_ail = 0.255;
%  of control surface flap (chordwise)	
OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_flap = 0.255;

y_loc_1R_y2_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD;
y_loc_1R_y1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD;
b_w1_e2 =  y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CAD;
y_offset_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD;

y_1R_y2_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD;
y_1R_y1_ail =  y_1R_y2_ail - 671/1000;
K_y1_ail_w1 = (y_1R_y1_ail - y_offset_w1)/(b_w1_e2);
K_y2_ail_w1 = (y_1R_y2_ail - y_offset_w1)/(b_w1_e2);

OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_ail_w1 = K_y1_ail_w1;
OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_ail_w1 = K_y2_ail_w1;

% fracción de cuerda que ocupan flaps/alerones: 0.255
% posición en evergadura del aleron (midiendo desde punta de ala): desde 0 hasta 671 [mm]
% posición en evergadura del flaps (midiendo desde punta de ala): desde 909 hasta 2304 [mm]


% Identifies folder where to store Figures
filename = '../Results/BAT/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;


switch mission_actual
    case 1 % ID_05. c_a, b+50 - Case 1 Final
        % CASE 1
        %% Modification that include interference with RACK
        % Used to elimnate a portion of the wing assuming that does not generate
        % lift to take into account the loss of effective area due to the racks
        % Racks_wing_span = 0.00; no wing aerodynamic interference considered
        % Racks_wing_span = 0.30; wing aerodynamic interference considered
        Racks_wing_span = 0.00;
        
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD + Racks_wing_span;
        y_offset_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD;
        b_w1_e2 = (OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD - OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD);
        K_y1_ail_w1 = (y_1R_y1_ail - y_offset_w1)/(b_w1_e2);
        K_y2_ail_w1 = (y_1R_y2_ail - y_offset_w1)/(b_w1_e2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_ail_w1 = K_y1_ail_w1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_ail_w1 = K_y2_ail_w1;
        
        % Mass and Inertia: Assumed to be those of T120 with missiles
        % Geometry
        % Weight Configuration from AERTEC
        conditions.MOTOR = 1; % 1 para SP210, 2 para DA215
        conditions.RACK = 1; % 1 == rack, 0 == sin rack
        conditions.DEPOSITO = 1; % 1 == deposito original, 2 == depósito húmedo
        conditions.FUEL = 1;  % en tanto por 1
        conditions.W_PL = 3.75;
%         conditions.n_MSL = 0; %0,1,2,3,4
        % Condiciones que definen en que situación se ubica el KSA
%         conditions.n_MSL = 4; % se han sustituido los datos de la configuración de 4 micros por los de la configuración KSA BAT+BOLA, sin racks a analizar, de forma que no es necesaria la modificación de los datos de entrada.
        conditions.n_MSL = 3; % se han sustituido los datos de la configuración de 3 micros por los de la configuración KSA SÓLO BOLA, sin racks.

        OUTPUT = CG_TARSIS120_KSA(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
               
        conditions.m_TOW = OUTPUT(1);
        conditions.x_XCG = OUTPUT(2)/conv_factor  + Dxcg;
        
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = OUTPUT(2)/conv_factor + Dxcg;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = OUTPUT(3)/conv_factor;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = OUTPUT(4)/conv_factor;

        OUTPUT_read_XLSX.Weights_flags.I_xx = OUTPUT(5); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = OUTPUT(6); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = OUTPUT(7); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = OUTPUT(8); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = OUTPUT(9); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = OUTPUT(10);

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        % Flag that does use thesimplificationof wing to span
        OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam = 1; 

        % Drag of missiles and pod simulaton each rack
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 0;% 0 no missile, 1 missile, 2 KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.n_missile = 2; % asume que hay resistencia igual a 2 KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.pod = 0;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.KSA = 1;

        % Plots prefix
        prefix = 'CASO_ID_05_1_TARSIS_Final_KSA_case1_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
        %% Rolling Moments definition
        % CAse 1: NO rolling movement
        % CAse 2: Rolling movement asymetric weight 1 different in each wing
        % CAse 3: Rolling movement asymetric weight 2 - only in 1 wing 1
        % CAse 4: Rolling movement asymetric weight 2 - only in 1 wing 2
        conditions.conditions_TRIM_lat.Rolling_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir = -1; 
        
        % BAT + Bola
        % BAT
        conditions.conditions_TRIM_lat.sign1 = 1;
        conditions.conditions_TRIM_lat.L_arm1 = 1359.676/1000;
        conditions.conditions_TRIM_lat.L_weight1 = 12*9.81; % 12 Kg
        % Bola
        conditions.conditions_TRIM_lat.sign2 = -1;
        conditions.conditions_TRIM_lat.L_arm2 = 1359.676/1000;
        conditions.conditions_TRIM_lat.L_weight2 = 4.5*9.81;% 4.5 Kg
        
        %% Yawing Moments definition
        % CAse 1: NO Yawing movement
        % CAse 2: Yawing movement drag asymetry
        % CAse 3: Yawing movement engine asymetry
        % change yawing moment direction
        conditions.conditions_TRIM_lat.Yawing_maneuver = 2;
        N_moment = b_w1/2; % Worse case scenarios
        N_moment = (1119.227/1000); % Worse case scenarios
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        % To correcto older Version  - Oct 2024
        conditions.conditions_TRIM_lat.Z_arm1 = 0;
        conditions.conditions_TRIM_lat.Z_arm2 = 0;
        conditions.conditions_TRIM_lat.direction_beta = 1; % defines the directión of sideslip
        
end
