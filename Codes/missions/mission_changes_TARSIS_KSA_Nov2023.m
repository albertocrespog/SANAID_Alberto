function [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS_KSA_Nov2023(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

% Weight_mission = T120_weightconfiguration;
%
% conditions.m_TOW = Weight_mission.W_vec_f(1);


%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

% Percentage where bayoneta is located
pp_bayoneta = (1100.55-966.57)/469;
cR_w1_T75 = 0.421;
x_loc_1R_y1_w1_CAD_T75 = 1.45766;
x_loc_1R_y1_w1_CAD_T75_bayoneta = x_loc_1R_y1_w1_CAD_T75 + pp_bayoneta*cR_w1_T75;
x_loc_1R_y1_w2_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

% conditions.T120 = 1; % Tarsis 120 with KSA
% if conditions.T120 == 1
%     CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
% elseif conditions.T120 == 2
%     CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
% end

% CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
% CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
%OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;

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
filename = '../Results/TARSIS_120/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;


switch mission_actual
    case 1 % CASE 5
        %% Modification that include interference with RACK
        % Used to elimnate a portion of the wing assuming that does not generate
        % lift to take into account the loss of effective area due to the racks
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
        % For KSA 3 micros (solo bola) y 4 micros (KSA completo: BAT+bola).
        conditions.MOTOR = 1; % 1 para SP210, 2 para DA215
        conditions.RACK = 0; % 1 == rack, 0 == sin rack
        conditions.DEPOSITO = 1; % 1 == deposito original, 2 == depósito húmedo
        conditions.FUEL = 1;  % en tanto por 1
        conditions.W_PL = 3.75;
        conditions.n_MSL = 0; %0,1,2,3,4
        
        % conditions.T120 = 1; % Tarsis 120 with FOX
        % conditions.T120 = 2; % Tarsis 120 with KSA
        % conditions.T120 = 3; % Tarsis 120 with Experimental values (Wrong)
        % conditions.T120 = 3; % Tarsis 120 with Experimental values (Nov 30th 2023)
        conditions.T120 = 1; % Tarsis 120 with KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1; % Missile Configuration
        if conditions.T120 == 1
            CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            %Lectura datos CG depósito
            CG_DEP = xlsread('Tabla cdg e inercias deposito T120.xlsx');
            OUTPUT_read_XLSX.Weights_flags.CG_DEP = CG_DEP;
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA_V5(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120,CG_DEP);
        elseif conditions.T120 == 4
           % Lectura datos T120
           conditions.Tabla_cg_inercias_T120 = xlsread('Tabla cdg e inercias T120 KSA 2.xlsx', 'D4:M13');
           % Lectura datos CG depósito
           FUEL = xlsread('Tabla cdg  e inercias depósito v2.xlsx', 'B3:N17');
           [percentage_FUEL OUTPUT] = CG_TARSIS120_KSA_V6(conditions.W_PL, conditions.Tabla_cg_inercias_T120, FUEL, conditions.n_MSL);
           conditions.FUEL = percentage_FUEL;
        end
        
        OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
             
        conditions.FUEL = percentage_FUEL;
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
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 0;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.n_missile = 0;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.pod = 0;

        % Plots prefix
        prefix = 'CASO_ID_05_1_TARSIS_Final_case5_KSA_subcase1_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
             
        %%Stability
        OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;
        
        % FOX/ BAT in right wing
        conditions.conditions_TRIM_lat.sign1 = 1;
        conditions.conditions_TRIM_lat.L_arm1 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm1 = -210.997/1000;
        conditions.conditions_TRIM_lat.L_weight1 = 0*9.81; % 12 Kg
        % Bola in left wing
        conditions.conditions_TRIM_lat.sign2 = -1;
        conditions.conditions_TRIM_lat.L_arm2 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm2 = -166.826/1000;
        conditions.conditions_TRIM_lat.L_weight2 = 0*9.81;% 4.5 Kg
        
        %% Rolling Moments definition
        % CAse 1: NO rolling movement
        % CAse 2: Rolling moment asymetric weight 1 different in each wing
        % CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
        % CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
        conditions.conditions_TRIM_lat.Rolling_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 

        %% Yawing Moments definition
        % CAse 1: NO Yawing moment
        % CAse 2: Yawing movement drag asymetry different in each wing
        % CAse 3: Yawing movement drag asymetry  - only in wing 1
        % CAse 4: Yawing movement drag asymetry  - only in wing 2
        % CAse 5: Yawing movement engine asymetry
        % CAse 6: Yawing movement engine asymetry and drag asymmetry
        % change yawing moment direction
        conditions.conditions_TRIM_lat.Yawing_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_yaw = 1; 
        conditions.conditions_TRIM_lat.direction_beta = 1; % defines the directión of sideslip
        
%         b_w1_e2 = (OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD - OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD);
%         N_moment = b_w1_e2/2; % Worse case scenarios
%         N_moment = (1119.227/1000); % Worse case scenarios
        N_moment = (1359.676/1000); % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

    case 2 % CASE 5
        %% Modification that include interference with RACK
        % Used to elimnate a portion of the wing assuming that does not generate
        % lift to take into account the loss of effective area due to the racks
        Racks_wing_span = 0.30;
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
%         conditions.FUEL = 1;  % en tanto por 1
        conditions.FUEL = 0.94229;  % en tanto por 1
        conditions.W_PL = 3.75;
        % Condiciones que definen en que situación se ubica el KSA
        % conditions.n_MSL = 3; % se han sustituido los datos de la configuración de 3 micros por los de la configuración KSA SÓLO BOLA, sin racks.
        % conditions.n_MSL = 4; % se han sustituido los datos de la configuración de 4 micros por los de la configuración KSA BAT+BOLA,
        % sin racks a analizar, de forma que no es necesaria la modificación de los datos de entrada.
        conditions.n_MSL = 3;
        
        % conditions.T120 = 1; % Tarsis 120 with FOX
        % conditions.T120 = 2; % Tarsis 120 with KSA
        % conditions.T120 = 3; % Tarsis 120 with Experimental values
        conditions.T120 = 3; % Tarsis 120 with KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1; % Missile Configuration
        if conditions.T120 == 1
            CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            %Lectura datos CG depósito
            CG_DEP = xlsread('Tabla cdg e inercias deposito T120.xlsx')
            OUTPUT_read_XLSX.Weights_flags.CG_DEP = CG_DEP;
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA_V5(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120,CG_DEP)
            pause
        end
        
        % Temporal ya que códigos from AERTC incorrect
        solve_ERROR_AERTEC = 2;
        if solve_ERROR_AERTEC == 1
            CATIAT120_temp = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT2 = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120_temp);
            
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
            conditions.m_TOW = OUTPUT(1);
            conditions.x_XCG = OUTPUT(2)/conv_factor  + Dxcg;
            
            OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = OUTPUT(2)/conv_factor + Dxcg;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = OUTPUT(3)/conv_factor;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = OUTPUT(4)/conv_factor;
            
            OUTPUT_read_XLSX.Weights_flags.I_xx = OUTPUT2(5); %
            OUTPUT_read_XLSX.Weights_flags.I_yy = OUTPUT2(6); %
            OUTPUT_read_XLSX.Weights_flags.I_zz = OUTPUT2(7); %
            OUTPUT_read_XLSX.Weights_flags.I_xy = OUTPUT2(8); %
            OUTPUT_read_XLSX.Weights_flags.I_xz = OUTPUT2(9); %
            OUTPUT_read_XLSX.Weights_flags.I_yz = OUTPUT2(10);
        else
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
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
        end

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        % Flag that does use thesimplificationof wing to span
        OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam = 1; 

        % Drag of missiles and pod simulaton each rack
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.n_missile = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.pod = 1;

        % Plots prefix
        prefix = 'CASO_ID_05_1_TARSIS_Final_case5_KSA_subcase2_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
             
        %%Stability
        OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;
        
        % FOX/ BAT in right wing
        conditions.conditions_TRIM_lat.sign1 = 1;
        conditions.conditions_TRIM_lat.L_arm1 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm1 = -210.997/1000;
        conditions.conditions_TRIM_lat.L_weight1 = 0*9.81; % 12 Kg
        % Bola in left wing
        conditions.conditions_TRIM_lat.sign2 = -1;
        conditions.conditions_TRIM_lat.L_arm2 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm2 = -166.826/1000;
        conditions.conditions_TRIM_lat.L_weight2 = 4.5*9.81;% 4.5 Kg
        
        %% Rolling Moments definition
        % CAse 1: NO rolling movement
        % CAse 2: Rolling moment asymetric weight 1 different in each wing
        % CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
        % CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
        conditions.conditions_TRIM_lat.Rolling_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 

        %% Yawing Moments definition
        % CAse 1: NO Yawing moment
        % CAse 2: Yawing movement drag asymetry different in each wing
        % CAse 3: Yawing movement drag asymetry  - only in wing 1
        % CAse 4: Yawing movement drag asymetry  - only in wing 2
        % CAse 5: Yawing movement engine asymetry
        % CAse 6: Yawing movement engine asymetry and drag asymmetry
        % change yawing moment direction
        conditions.conditions_TRIM_lat.Yawing_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_yaw = 1; 
        conditions.conditions_TRIM_lat.direction_beta = 1; % defines the directión of sideslip
        
%         b_w1_e2 = (OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD - OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD);
%         N_moment = b_w1_e2/2; % Worse case scenarios
%         N_moment = (1119.227/1000); % Worse case scenarios
        N_moment = (1359.676/1000); % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;
        
    case 3 % CASE 5
        %% Modification that include interference with RACK
        % Used to elimnate a portion of the wing assuming that does not generate
        % lift to take into account the loss of effective area due to the racks
        % Racks_wing_span = 0.00; no wing aerodynamic interference considered
        % Racks_wing_span = 0.30; wing aerodynamic interference considered
        Racks_wing_span = 0.30;
        
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
%         conditions.FUEL = 1;  % en tanto por 1
        conditions.FUEL = 0.3519;  % en tanto por 1
        conditions.W_PL = 3.75;
        % Condiciones que definen en que situación se ubica el KSA
        % conditions.n_MSL = 3; % se han sustituido los datos de la configuración de 3 micros por los de la configuración KSA SÓLO BOLA, sin racks.
        % conditions.n_MSL = 4; % se han sustituido los datos de la configuración de 4 micros por los de la configuración KSA BAT+BOLA,
        % sin racks a analizar, de forma que no es necesaria la modificación de los datos de entrada.
        conditions.n_MSL = 4;
        
        % conditions.T120 = 1; % Tarsis 120 with FOX
        % conditions.T120 = 2; % Tarsis 120 with KSA
        % conditions.T120 = 3; % Tarsis 120 with Experimental values
        conditions.T120 = 3; % Tarsis 120 with KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1; % Missile Configuration
        if conditions.T120 == 1
            CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            %Lectura datos CG depósito
            CG_DEP = xlsread('Tabla cdg e inercias deposito T120.xlsx');
            OUTPUT_read_XLSX.Weights_flags.CG_DEP = CG_DEP;
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA_V5(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120,CG_DEP);
        end
        
        % Temporal ya que códigos from AERTC incorrect
        solve_ERROR_AERTEC = 2;
        if solve_ERROR_AERTEC == 1
            CATIAT120_temp = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT2 = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120_temp);
            
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
            conditions.m_TOW = OUTPUT(1);
            conditions.x_XCG = OUTPUT(2)/conv_factor  + Dxcg;
            
            OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = OUTPUT(2)/conv_factor + Dxcg;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = OUTPUT(3)/conv_factor;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = OUTPUT(4)/conv_factor;
            
            OUTPUT_read_XLSX.Weights_flags.I_xx = OUTPUT2(5); %
            OUTPUT_read_XLSX.Weights_flags.I_yy = OUTPUT2(6); %
            OUTPUT_read_XLSX.Weights_flags.I_zz = OUTPUT2(7); %
            OUTPUT_read_XLSX.Weights_flags.I_xy = OUTPUT2(8); %
            OUTPUT_read_XLSX.Weights_flags.I_xz = OUTPUT2(9); %
            OUTPUT_read_XLSX.Weights_flags.I_yz = OUTPUT2(10);
        else
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
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
        end

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        % Flag that does use thesimplificationof wing to span
        OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam = 1; 

        % Drag of missiles and pod simulaton each rack
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 0;% 0 no missile, 1 missile, 2 KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.n_missile = 2; % asume que hay resistencia igual a 2 KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.pod = 0;

        % Plots prefix
        prefix = 'CASO_ID_05_1_TARSIS_Final_case5_KSA_subcase3_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
             
        %%Stability
        OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;
        
        % FOX/ BAT in right wing
        conditions.conditions_TRIM_lat.sign1 = 1;
        conditions.conditions_TRIM_lat.L_arm1 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm1 = -210.997/1000;
        conditions.conditions_TRIM_lat.L_weight1 = 12*9.81; % 12 Kg
        % Bola in left wing
        conditions.conditions_TRIM_lat.sign2 = -1;
        conditions.conditions_TRIM_lat.L_arm2 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm2 = -166.826/1000;
        conditions.conditions_TRIM_lat.L_weight2 = 4.5*9.81;% 4.5 Kg
        
        %% Rolling Moments definition
        % CAse 1: NO rolling movement
        % CAse 2: Rolling moment asymetric weight 1 different in each wing
        % CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
        % CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
        conditions.conditions_TRIM_lat.Rolling_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 

        %% Yawing Moments definition
        % CAse 1: NO Yawing moment
        % CAse 2: Yawing movement drag asymetry different in each wing
        % CAse 3: Yawing movement drag asymetry  - only in wing 1
        % CAse 4: Yawing movement drag asymetry  - only in wing 2
        % CAse 5: Yawing movement engine asymetry
        % CAse 6: Yawing movement engine asymetry and drag asymmetry
        % change yawing moment direction
        conditions.conditions_TRIM_lat.Yawing_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_yaw = 1; 
        conditions.conditions_TRIM_lat.direction_beta = 1; % defines the directión of sideslip
        
%         b_w1_e2 = (OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD - OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD);
%         N_moment = b_w1_e2/2; % Worse case scenarios
%         N_moment = (1119.227/1000); % Worse case scenarios
        N_moment = (1359.676/1000); % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;    
        
    case 4 % CASE 5
        %% Modification that include interference with RACK
        % Used to elimnate a portion of the wing assuming that does not generate
        % lift to take into account the loss of effective area due to the racks
        Racks_wing_span = 0.30;
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
%         conditions.FUEL = 1;  % en tanto por 1
        conditions.FUEL = 0.94229;  % en tanto por 1
        conditions.W_PL = 3.75;
        % Condiciones que definen en que situación se ubica el KSA
        % conditions.n_MSL = 3; % se han sustituido los datos de la configuración de 3 micros por los de la configuración KSA SÓLO BOLA, sin racks.
        % conditions.n_MSL = 4; % se han sustituido los datos de la configuración de 4 micros por los de la configuración KSA BAT+BOLA,
        % sin racks a analizar, de forma que no es necesaria la modificación de los datos de entrada.
        conditions.n_MSL = 3;
        
        % conditions.T120 = 1; % Tarsis 120 with FOX
        % conditions.T120 = 2; % Tarsis 120 with KSA
        % conditions.T120 = 3; % Tarsis 120 with Experimental values
        conditions.T120 = 3; % Tarsis 120 with KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1; % Missile Configuration
        if conditions.T120 == 1
            CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            %Lectura datos CG depósito
            CG_DEP = xlsread('Tabla cdg e inercias deposito T120.xlsx');
            OUTPUT_read_XLSX.Weights_flags.CG_DEP = CG_DEP;
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA_V5(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120,CG_DEP);
        end
                
        % Temporal ya que códigos from AERTC incorrect
        solve_ERROR_AERTEC = 2;
        if solve_ERROR_AERTEC == 1
            CATIAT120_temp = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT2 = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120_temp);
            
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
            conditions.m_TOW = OUTPUT(1);
            conditions.x_XCG = OUTPUT(2)/conv_factor  + Dxcg;
            
            OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = OUTPUT(2)/conv_factor + Dxcg;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = OUTPUT(3)/conv_factor;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = OUTPUT(4)/conv_factor;
            
            OUTPUT_read_XLSX.Weights_flags.I_xx = OUTPUT2(5); %
            OUTPUT_read_XLSX.Weights_flags.I_yy = OUTPUT2(6); %
            OUTPUT_read_XLSX.Weights_flags.I_zz = OUTPUT2(7); %
            OUTPUT_read_XLSX.Weights_flags.I_xy = OUTPUT2(8); %
            OUTPUT_read_XLSX.Weights_flags.I_xz = OUTPUT2(9); %
            OUTPUT_read_XLSX.Weights_flags.I_yz = OUTPUT2(10);
        else
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
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
        end

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        % Flag that does use thesimplificationof wing to span
        OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam = 1; 

        % Drag of missiles and pod simulaton each rack
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.n_missile = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.pod = 1;

        % Plots prefix
        prefix = 'CASO_ID_05_1_TARSIS_Final_case5_KSA_subcase4_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
             
        %%Stability
        OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;
        
        % FOX/ BAT in right wing
        conditions.conditions_TRIM_lat.sign1 = 1;
        conditions.conditions_TRIM_lat.L_arm1 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm1 = -210.997/1000;
        conditions.conditions_TRIM_lat.L_weight1 = 0*9.81; % 12 Kg
        % Bola in left wing
        conditions.conditions_TRIM_lat.sign2 = -1;
        conditions.conditions_TRIM_lat.L_arm2 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm2 = -166.826/1000;
        conditions.conditions_TRIM_lat.L_weight2 = 4.5*9.81;% 4.5 Kg
        
        %% Rolling Moments definition
        % CAse 1: NO rolling movement
        % CAse 2: Rolling moment asymetric weight 1 different in each wing
        % CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
        % CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
        conditions.conditions_TRIM_lat.Rolling_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 

        %% Yawing Moments definition
        % CAse 1: NO Yawing moment
        % CAse 2: Yawing movement drag asymetry different in each wing
        % CAse 3: Yawing movement drag asymetry  - only in wing 1
        % CAse 4: Yawing movement drag asymetry  - only in wing 2
        % CAse 5: Yawing movement engine asymetry
        % CAse 6: Yawing movement engine asymetry and drag asymmetry
        % change yawing moment direction
        conditions.conditions_TRIM_lat.Yawing_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_yaw = 1; 
        conditions.conditions_TRIM_lat.direction_beta = -1; % defines the directión of sideslip
        
%         b_w1_e2 = (OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD - OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD);
%         N_moment = b_w1_e2/2; % Worse case scenarios
%         N_moment = (1119.227/1000); % Worse case scenarios
        N_moment = (1359.676/1000); % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;
      
    case 5 % CASE 5
        %% Modification that include interference with RACK
        % Used to elimnate a portion of the wing assuming that does not generate
        % lift to take into account the loss of effective area due to the racks
        % Racks_wing_span = 0.00; no wing aerodynamic interference considered
        % Racks_wing_span = 0.30; wing aerodynamic interference considered
        Racks_wing_span = 0.30;
        
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
%         conditions.FUEL = 1;  % en tanto por 1
        conditions.FUEL = 0.3519;  % en tanto por 1
        conditions.W_PL = 3.75;
        % Condiciones que definen en que situación se ubica el KSA
        % conditions.n_MSL = 3; % se han sustituido los datos de la configuración de 3 micros por los de la configuración KSA SÓLO BOLA, sin racks.
        % conditions.n_MSL = 4; % se han sustituido los datos de la configuración de 4 micros por los de la configuración KSA BAT+BOLA,
        % sin racks a analizar, de forma que no es necesaria la modificación de los datos de entrada.
        conditions.n_MSL = 4;
        
        % conditions.T120 = 1; % Tarsis 120 with FOX
        % conditions.T120 = 2; % Tarsis 120 with KSA
        % conditions.T120 = 3; % Tarsis 120 with Experimental values
        conditions.T120 = 3; % Tarsis 120 with KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1; % Missile Configuration
        if conditions.T120 == 1
            CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            %Lectura datos CG depósito
            CG_DEP = xlsread('Tabla cdg e inercias deposito T120.xlsx');
            OUTPUT_read_XLSX.Weights_flags.CG_DEP = CG_DEP;
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA_V5(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120,CG_DEP);
        end
        
        % Temporal ya que códigos from AERTC incorrect
        solve_ERROR_AERTEC = 2;
        if solve_ERROR_AERTEC == 1
            CATIAT120_temp = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT2 = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120_temp);
            
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
            conditions.m_TOW = OUTPUT(1);
            conditions.x_XCG = OUTPUT(2)/conv_factor  + Dxcg;
            
            OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = OUTPUT(2)/conv_factor + Dxcg;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = OUTPUT(3)/conv_factor;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = OUTPUT(4)/conv_factor;
            
            OUTPUT_read_XLSX.Weights_flags.I_xx = OUTPUT2(5); %
            OUTPUT_read_XLSX.Weights_flags.I_yy = OUTPUT2(6); %
            OUTPUT_read_XLSX.Weights_flags.I_zz = OUTPUT2(7); %
            OUTPUT_read_XLSX.Weights_flags.I_xy = OUTPUT2(8); %
            OUTPUT_read_XLSX.Weights_flags.I_xz = OUTPUT2(9); %
            OUTPUT_read_XLSX.Weights_flags.I_yz = OUTPUT2(10);
        else
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
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
        end

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        % Flag that does use thesimplificationof wing to span
        OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam = 1; 

        % Drag of missiles and pod simulaton each rack
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 0;% 0 no missile, 1 missile, 2 KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.n_missile = 2; % asume que hay resistencia igual a 2 KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.pod = 0;

        % Plots prefix
        prefix = 'CASO_ID_05_1_TARSIS_Final_case5_KSA_subcase5_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
             
        %%Stability
        OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;
        
        % FOX/ BAT in right wing
        conditions.conditions_TRIM_lat.sign1 = 1;
        conditions.conditions_TRIM_lat.L_arm1 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm1 = -210.997/1000;
        conditions.conditions_TRIM_lat.L_weight1 = 12*9.81; % 12 Kg
        % Bola in left wing
        conditions.conditions_TRIM_lat.sign2 = -1;
        conditions.conditions_TRIM_lat.L_arm2 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm2 = -166.826/1000;
        conditions.conditions_TRIM_lat.L_weight2 = 4.5*9.81;% 4.5 Kg
        
        %% Rolling Moments definition
        % CAse 1: NO rolling movement
        % CAse 2: Rolling moment asymetric weight 1 different in each wing
        % CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
        % CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
        conditions.conditions_TRIM_lat.Rolling_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 

        %% Yawing Moments definition
        % CAse 1: NO Yawing moment
        % CAse 2: Yawing movement drag asymetry different in each wing
        % CAse 3: Yawing movement drag asymetry  - only in wing 1
        % CAse 4: Yawing movement drag asymetry  - only in wing 2
        % CAse 5: Yawing movement engine asymetry
        % CAse 6: Yawing movement engine asymetry and drag asymmetry
        % change yawing moment direction
        conditions.conditions_TRIM_lat.Yawing_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_yaw = 1; 
        conditions.conditions_TRIM_lat.direction_beta = -1; % defines the directión of sideslip
        
%         b_w1_e2 = (OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD - OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD);
%         N_moment = b_w1_e2/2; % Worse case scenarios
%         N_moment = (1119.227/1000); % Worse case scenarios
        N_moment = (1359.676/1000); % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;
        
    case 6 % CASE 5
        %% Modification that include interference with RACK
        % Used to elimnate a portion of the wing assuming that does not generate
        % lift to take into account the loss of effective area due to the racks
        Racks_wing_span = 0.30;
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
%         conditions.FUEL = 1;  % en tanto por 1
        conditions.FUEL = 0.1;  % en tanto por 1
        conditions.W_PL = 3.75;
        % Condiciones que definen en que situación se ubica el KSA
        % conditions.n_MSL = 3; % se han sustituido los datos de la configuración de 3 micros por los de la configuración KSA SÓLO BOLA, sin racks.
        % conditions.n_MSL = 4; % se han sustituido los datos de la configuración de 4 micros por los de la configuración KSA BAT+BOLA,
        % sin racks a analizar, de forma que no es necesaria la modificación de los datos de entrada.
        conditions.n_MSL = 3;
        
        % conditions.T120 = 1; % Tarsis 120 with FOX
        % conditions.T120 = 2; % Tarsis 120 with KSA
        % conditions.T120 = 3; % Tarsis 120 with Experimental values
        conditions.T120 = 3; % Tarsis 120 with KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1; % Missile Configuration
        if conditions.T120 == 1
            CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            %Lectura datos CG depósito
            CG_DEP = xlsread('Tabla cdg e inercias deposito T120.xlsx');
            OUTPUT_read_XLSX.Weights_flags.CG_DEP = CG_DEP;
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA_V5(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120,CG_DEP);
        end
                
        % Temporal ya que códigos from AERTC incorrect
        solve_ERROR_AERTEC = 2;
        if solve_ERROR_AERTEC == 1
            CATIAT120_temp = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT2 = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120_temp);
            
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
            conditions.m_TOW = OUTPUT(1);
            conditions.x_XCG = OUTPUT(2)/conv_factor  + Dxcg;
            
            OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = OUTPUT(2)/conv_factor + Dxcg;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = OUTPUT(3)/conv_factor;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = OUTPUT(4)/conv_factor;
            
            OUTPUT_read_XLSX.Weights_flags.I_xx = OUTPUT2(5); %
            OUTPUT_read_XLSX.Weights_flags.I_yy = OUTPUT2(6); %
            OUTPUT_read_XLSX.Weights_flags.I_zz = OUTPUT2(7); %
            OUTPUT_read_XLSX.Weights_flags.I_xy = OUTPUT2(8); %
            OUTPUT_read_XLSX.Weights_flags.I_xz = OUTPUT2(9); %
            OUTPUT_read_XLSX.Weights_flags.I_yz = OUTPUT2(10);
        else
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
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
        end

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        % Flag that does use thesimplificationof wing to span
        OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam = 1; 

        % Drag of missiles and pod simulaton each rack
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.n_missile = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.pod = 1;

        % Plots prefix
        prefix = 'CASO_ID_05_1_TARSIS_Final_case5_KSA_subcase6_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
             
        %%Stability
        OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;
        
        % FOX/ BAT in right wing
        conditions.conditions_TRIM_lat.sign1 = 1;
        conditions.conditions_TRIM_lat.L_arm1 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm1 = -210.997/1000;
        conditions.conditions_TRIM_lat.L_weight1 = 0*9.81; % 12 Kg
        % Bola in left wing
        conditions.conditions_TRIM_lat.sign2 = -1;
        conditions.conditions_TRIM_lat.L_arm2 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm2 = -166.826/1000;
        conditions.conditions_TRIM_lat.L_weight2 = 4.5*9.81;% 4.5 Kg
        
        %% Rolling Moments definition
        % CAse 1: NO rolling movement
        % CAse 2: Rolling moment asymetric weight 1 different in each wing
        % CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
        % CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
        conditions.conditions_TRIM_lat.Rolling_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 

        %% Yawing Moments definition
        % CAse 1: NO Yawing moment
        % CAse 2: Yawing movement drag asymetry different in each wing
        % CAse 3: Yawing movement drag asymetry  - only in wing 1
        % CAse 4: Yawing movement drag asymetry  - only in wing 2
        % CAse 5: Yawing movement engine asymetry
        % CAse 6: Yawing movement engine asymetry and drag asymmetry
        % change yawing moment direction
        conditions.conditions_TRIM_lat.Yawing_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_yaw = 1; 
        conditions.conditions_TRIM_lat.direction_beta = -1; % defines the directión of sideslip
        
%         b_w1_e2 = (OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD - OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD);
%         N_moment = b_w1_e2/2; % Worse case scenarios
%         N_moment = (1119.227/1000); % Worse case scenarios
        N_moment = (1359.676/1000); % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;
    case 7 % CASE 6
        %% Modification that include interference with RACK
        % Used to elimnate a portion of the wing assuming that does not generate
        % lift to take into account the loss of effective area due to the racks
        % Racks_wing_span = 0.00; no wing aerodynamic interference considered
        % Racks_wing_span = 0.30; wing aerodynamic interference considered
        Racks_wing_span = 0.30;
        
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
%         conditions.FUEL = 1;  % en tanto por 1
        conditions.FUEL = 0.1;  % en tanto por 1
        conditions.W_PL = 3.75;
        % Condiciones que definen en que situación se ubica el KSA
        % conditions.n_MSL = 3; % se han sustituido los datos de la configuración de 3 micros por los de la configuración KSA SÓLO BOLA, sin racks.
        % conditions.n_MSL = 4; % se han sustituido los datos de la configuración de 4 micros por los de la configuración KSA BAT+BOLA,
        % sin racks a analizar, de forma que no es necesaria la modificación de los datos de entrada.
        conditions.n_MSL = 4;
        
        % conditions.T120 = 1; % Tarsis 120 with FOX
        % conditions.T120 = 2; % Tarsis 120 with KSA
        % conditions.T120 = 3; % Tarsis 120 with Experimental values
        conditions.T120 = 3; % Tarsis 120 with KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 1; % Missile Configuration
        if conditions.T120 == 1
            CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            %Lectura datos CG depósito
            CG_DEP = xlsread('Tabla cdg e inercias deposito T120.xlsx');
            OUTPUT_read_XLSX.Weights_flags.CG_DEP = CG_DEP;
            CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');
            OUTPUT = CG_TARSIS120_KSA_V5(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120,CG_DEP);
        end
        
        % Temporal ya que códigos from AERTC incorrect
        solve_ERROR_AERTEC = 2;
        if solve_ERROR_AERTEC == 1
            CATIAT120_temp = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
            OUTPUT2 = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120_temp);
            
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
            conditions.m_TOW = OUTPUT(1);
            conditions.x_XCG = OUTPUT(2)/conv_factor  + Dxcg;
            
            OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = OUTPUT(2)/conv_factor + Dxcg;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = OUTPUT(3)/conv_factor;
            OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = OUTPUT(4)/conv_factor;
            
            OUTPUT_read_XLSX.Weights_flags.I_xx = OUTPUT2(5); %
            OUTPUT_read_XLSX.Weights_flags.I_yy = OUTPUT2(6); %
            OUTPUT_read_XLSX.Weights_flags.I_zz = OUTPUT2(7); %
            OUTPUT_read_XLSX.Weights_flags.I_xy = OUTPUT2(8); %
            OUTPUT_read_XLSX.Weights_flags.I_xz = OUTPUT2(9); %
            OUTPUT_read_XLSX.Weights_flags.I_yz = OUTPUT2(10);
        else
            OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;
            
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
        end

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        % Flag that does use thesimplificationof wing to span
        OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam = 1; 

        % Drag of missiles and pod simulaton each rack
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.missile = 0;% 0 no missile, 1 missile, 2 KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.n_missile = 2; % asume que hay resistencia igual a 2 KSA
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf.pod = 0;

        % Plots prefix
        prefix = 'CASO_ID_05_1_TARSIS_Final_case5_KSA_subcase7_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
             
        %%Stability
        OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;
        
        % FOX/ BAT in right wing
        conditions.conditions_TRIM_lat.sign1 = 1;
        conditions.conditions_TRIM_lat.L_arm1 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm1 = -210.997/1000;
        conditions.conditions_TRIM_lat.L_weight1 = 12*9.81; % 12 Kg
        % Bola in left wing
        conditions.conditions_TRIM_lat.sign2 = -1;
        conditions.conditions_TRIM_lat.L_arm2 = 1359.676/1000;
        conditions.conditions_TRIM_lat.Z_arm2 = -166.826/1000;
        conditions.conditions_TRIM_lat.L_weight2 = 4.5*9.81;% 4.5 Kg
        
        %% Rolling Moments definition
        % CAse 1: NO rolling movement
        % CAse 2: Rolling moment asymetric weight 1 different in each wing
        % CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
        % CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
        conditions.conditions_TRIM_lat.Rolling_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 

        %% Yawing Moments definition
        % CAse 1: NO Yawing moment
        % CAse 2: Yawing movement drag asymetry different in each wing
        % CAse 3: Yawing movement drag asymetry  - only in wing 1
        % CAse 4: Yawing movement drag asymetry  - only in wing 2
        % CAse 5: Yawing movement engine asymetry
        % CAse 6: Yawing movement engine asymetry and drag asymmetry
        % change yawing moment direction
        conditions.conditions_TRIM_lat.Yawing_maneuver = 2;
        conditions.conditions_TRIM_lat.change_dir_yaw = 1; 
        conditions.conditions_TRIM_lat.direction_beta = -1; % defines the directión of sideslip
        
%         b_w1_e2 = (OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD - OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD);
%         N_moment = b_w1_e2/2; % Worse case scenarios
%         N_moment = (1119.227/1000); % Worse case scenarios
        N_moment = (1359.676/1000); % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;
end
