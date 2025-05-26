function Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX,case_AC)

% Assigns the value
Geo_tier = Geo_input_tier;

% Replaces the input data into the geometric file
% % Generation of XCG data
% Geometric_Data
% Legend of elements
% Wing 1 Front Wing (FW) (Canard)
% _w1 - wing 1
% Wing 2 Rear Wing (RW)
% _w2 - wing 2
% _v1 - VTP #1
% _v1 - VTP #2
% _fus - fuselage
% Legend of properties
% _LE - Leading Edge
% _TE - Trailing Edge
% b_ - wing span
% AR_ - Aspect Ratio
% S_ - Area
% c_ - chord
% _e - effective (eliminating covered area for example fuselage)
% _s - surface (disatance meassured along the surface taking into account dihedral
% lambda_ - taper ratio
% Lambda_LE - sweep of leading edge
% Lambda_TE - sweep of trailing edge
% Lambda_c4 - sweep of 1/4 chord
% Lambda_c2 - sweep of 1/2 chord
% cmac_ - Mean Aerodynamic Chord
% ybar_ - y-location of the cmac
% xbar_ - x-location of the cmac
% cR_ - root chord
% cT_ - tip chord
% c_ele_R  - chord elevator in root
% _ele - elevator
% _rud - rudder
% _ail - aileron
% Distances
% x_ - distance from reference point in the x direction. Positive rearwars for example:
% x_w1_LE - distance to LE of w1
% x_w1_xbar - distance to cmac of w1
% z_ - distance from reference point (origin set al origin in fuselage.
% Positive up
% y_ - distance from reference point (origin set al origin in fuselage)
% positive right

% Determines the ammount of control surface in the w2 aileron & elevator
% For Elevons
% Control_surface 0 = ELEVON;
% Control_surface 1 = W1 AND W2: FLAPS, AILERONS + RUDDERVATOR;
% Control_surface 2 = 2-SURFACE;
% Control_surface 3 = 3-SURFACE;
% Aircraft type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
% Engine location
% Engine_loc = 1 - under wings
% Engine_loc = 2 - fuselage
% Engine_loc = 3 - wingtips
% Engine Configuration
% Engine_conf = 1 - pusher prop
% Engine_conf = 2 - puller prop
% Engine_conf = 3 - turbine

% Stores Aircraft Configuration
CASE = AC_CONFIGURATION.CASE;
Control_surface = AC_CONFIGURATION.Control_surface;
% case_AC = AC_CONFIGURATION.case_AC;
Engine_loc = AC_CONFIGURATION.Engine_loc;
Engine_conf = AC_CONFIGURATION.Engine_conf;

% Propeller
D_prop = OUTPUT_read_XLSX.Propulsive_flags.D_prop; %

%% identifies the aerodynamic surfaces being used
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
HTP2 = AC_CONFIGURATION.HTP2;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
twin_VTP = AC_CONFIGURATION.twin_VTP;
d_ail = AC_CONFIGURATION.d_ail;
d_ele = AC_CONFIGURATION.d_ele;
d_elevon = AC_CONFIGURATION.d_elevon;
d_flap = AC_CONFIGURATION.d_flap;
d_rudder = AC_CONFIGURATION.d_rudder;
d_rudvtr = AC_CONFIGURATION.d_rudvtr;
d_rudvtr2 = AC_CONFIGURATION.d_rudvtr2;
d_can = AC_CONFIGURATION.d_can;

% Conversion
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;
% Propulsion Data
% D_prop = Prop_data.D_prop;

% Distance from origin in CATIa to Fuselage (offset value)
x_offset_CAD = Geo_tier.x_offset_CAD;
z_offset_CAD = Geo_tier.z_offset_CAD;

% Geometry

% For multiple Fuselages
if case_AC == 5 
        % case 5 % WIG - multiple fuselage
        % Fuselage geomtry (Units in m)
        w_fus = Geo_tier.w_fus1;
        h_fus = Geo_tier.h_fus1;
        l_fus = Geo_tier.l_fus1;
        d_fus = Geo_tier.d_fus1;
else
        w_fus = Geo_tier.w_fus;
        h_fus = Geo_tier.h_fus;
        l_fus = Geo_tier.l_fus;
        d_fus = Geo_tier.d_fus;   
end

% Center of gravity
% x_XCG = XCG_data.x_XCG;
% y_XCG = XCG_data.y_XCG;
% z_XCG = XCG_data.z_XCG;
x_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
y_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG;
z_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;

% if OUTPUT_read_XLSX.Fuselage_flags.CAD_Kink_Wing == 1
    % 
    % 
    % % Wingspan
    %     y_loc_1R_y1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD;
    %     y_loc_1R_yB1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_w1_CAD;
    %     y_loc_1R_yB2_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_w1_CAD;
    %     % Sweep
    %     Lambda_LE_w1_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e;
    %     Lambda_LE_w1_k1_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k1_e;
    %     Lambda_LE_w1_k2_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k2_e;
    %     % Dihedral
    %     dihedral_w1_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e;
    %     dihedral_w1_k1_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_k1_e;
    %     dihedral_w1_k2_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_k2_e;
    %     % Chrod
    %     cR_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;
    %     cB_k1_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_w1;
    %     cB_k2_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_w1;
    % 
    % 
    %         %% Aerodynamic surfaces location
    % % Location of LE w1
    % % 1R identifies LE
    % % 2R identifies TE
    % % y1 identifies inner position
    % % y2 identifies outer position
    % x_loc_1R_y1_w1_CAD = Geo_tier.x_loc_1R_y1_w1_CAD;
    % y_loc_1R_y1_w1_CAD = Geo_tier.y_loc_1R_y1_w1_CAD;
    % z_loc_1R_y1_w1_CAD = Geo_tier.z_loc_1R_y1_w1_CAD;
    % 
    % x_loc_1R_y2_w1_CAD = Geo_tier.x_loc_1R_y2_w1_CAD;
    % y_loc_1R_y2_w1_CAD = Geo_tier.y_loc_1R_y2_w1_CAD;
    % z_loc_1R_y2_w1_CAD = Geo_tier.z_loc_1R_y2_w1_CAD;
    % 
    % y_loc_1R_yB1_w1_CAD = Geo_tier.y_loc_1R_yB1_w1_CAD;
    % y_loc_1R_yB2_w1_CAD = Geo_tier.y_loc_1R_yB2_w1_CAD;
    % 
    % x_y1_y2_w1_CAD = Geo_tier.x_y1_y2_w1_CAD;
    % y_y1_y2_w1_CAD = Geo_tier.y_y1_y2_w1_CAD;
    % z_y1_y2_w1_CAD = Geo_tier.z_y1_y2_w1_CAD; 
    % 
    % % Wingspan - w1
    % b_w1 = Geo_tier.b_w1;
    % b_w1_e = Geo_tier.b_w1_e;
    % b_w1_fus = Geo_tier.b_w1_fus;
    % % Location of LE w1
    % x_loc_LE_w1_CAD = x_loc_1R_y1_w1_CAD; % distance from CAD refference point
    % y_loc_LE_w1_CAD = y_loc_1R_y1_w1_CAD; % distance from CAD refference point
    % z_loc_LE_w1_CAD = z_loc_1R_y1_w1_CAD; % distance from CAD refference point
    % % Correction from CAD identifying the nose as new Origin
    % x_loc_LE_w1 = x_loc_LE_w1_CAD + x_offset_CAD; % corrected to the nose of the aircraft
    % y_loc_LE_w1 = b_w1_fus/2;
    % z_loc_LE_w1 = z_loc_LE_w1_CAD + z_offset_CAD; % corrected to the nose of the aircraft
    % % Storing DATA
    % Geo_tier.x_loc_LE_w1 = x_loc_LE_w1;
    % Geo_tier.y_loc_LE_w1 = y_loc_LE_w1;
    % Geo_tier.z_loc_LE_w1 = z_loc_LE_w1;
    % Geo_tier.x_loc_LE_w1_CAD = x_loc_LE_w1_CAD;
    % Geo_tier.y_loc_LE_w1_CAD = y_loc_LE_w1_CAD;
    % Geo_tier.z_loc_LE_w1_CAD = z_loc_LE_w1_CAD;
    % 
    % % Location of LE w1
    % x_w1_LE = x_loc_LE_w1;
    % y_w1_LE = y_loc_LE_w1;
    % z_w1_LE = z_loc_LE_w1;
    % % Storing DATA
    % Geo_tier.x_w1_LE = x_w1_LE;
    % Geo_tier.y_w1_LE = y_w1_LE;
    % Geo_tier.z_w1_LE = z_w1_LE;
    % 
    % % Angles
    % Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
    % Lambda_LE_w1_k1_e = Geo_tier.Lambda_LE_w1_k1_e;
    % Lambda_LE_w1_k2_e = Geo_tier.Lambda_LE_w1_k2_e;
    % dihedral_w1 = Geo_tier.dihedral_w1;
    % dihedral_w1_e = Geo_tier.dihedral_w1_e;
    % dihedral_w1_k1_e = Geo_tier.dihedral_w1_k1_e;
    % dihedral_w1_k2_e = Geo_tier.dihedral_w1_k2_e;
    % % CAD meassurements
    % % Chord
    % cR_w1 = Geo_tier.cR_w1;
    % cB_k1_w1 = Geo_tier.cB_k1_w1;
    % cB_k2_w1 = Geo_tier.cB_k2_w1;
    % cT_w1 = Geo_tier.cT_w1;
    % 
    % 
    % 
    %         %% Geometric location of the w1 in order to define the corners of the surface
    % % 1R identifies LE
    % % 2R identifies TE
    % % y1 identifies inner position
    % % y2 identifies outer position
    % x_1R_y1_w1 = 0; % defines inner position of wing chord (aileron section) LE
    % x_2R_y1_w1 = x_1R_y1_w1 + cR_w1; % defines inner position of wing chord TE
    % % chord at each location
    % y_offset_w1 = b_w1_fus/2;
    % x_1R_y2_w1 =  x_1R_y1_w1 + (b_w1/2 - y_offset_w1)*tan(Lambda_LE_w1); % defines outter position of wing chord (aileron section) LE
    % x_2R_y2_w1 = x_1R_y2_w1 + cT_w1; % defines outter position of wing chord TE
    % % spanwise location
    % y_1R_y1_w1 = b_w1/2 - b_w1_e/2; % inner position from the center line
    % y_2R_y1_w1 = y_1R_y1_w1; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    % y_1R_y2_w1 = b_w1/2; % outter position from the center line
    % y_2R_y2_w1 = y_1R_y2_w1; % outter position from the center line (same as LE assumes chord paralel to x-axis)
    % % z-position
    % z_1R_y1_w1 = 0; % inner position from the center line
    % z_2R_y1_w1 = 0; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    % z_1R_y2_w1 = z_1R_y1_w1 + (b_w1/2 - y_offset_w1)*tan(dihedral_w1);
    % z_2R_y2_w1 = z_2R_y1_w1 + (b_w1/2 - y_offset_w1)*tan(dihedral_w1);
    % 
    % end


%% Wing Aerodynamic Surface
if W1 == 1
    %% Aerodynamic surfaces location
    % Location of LE w1
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_loc_1R_y1_w1_CAD = Geo_tier.x_loc_1R_y1_w1_CAD;
    y_loc_1R_y1_w1_CAD = Geo_tier.y_loc_1R_y1_w1_CAD;
    z_loc_1R_y1_w1_CAD = Geo_tier.z_loc_1R_y1_w1_CAD;
    x_loc_1R_y2_w1_CAD = Geo_tier.x_loc_1R_y2_w1_CAD;
    y_loc_1R_y2_w1_CAD = Geo_tier.y_loc_1R_y2_w1_CAD;
    z_loc_1R_y2_w1_CAD = Geo_tier.z_loc_1R_y2_w1_CAD;
    y_loc_1R_yB1_w1_CAD = Geo_tier.y_loc_1R_yB1_w1_CAD;
    y_loc_1R_yB2_w1_CAD = Geo_tier.y_loc_1R_yB2_w1_CAD;
    x_y1_y2_w1_CAD = Geo_tier.x_y1_y2_w1_CAD;
    y_y1_y2_w1_CAD = Geo_tier.y_y1_y2_w1_CAD;
    z_y1_y2_w1_CAD = Geo_tier.z_y1_y2_w1_CAD; 
    
    % Wingspan - w1
    b_w1 = Geo_tier.b_w1;
    b_w1_e = Geo_tier.b_w1_e;
    b_w1_fus = Geo_tier.b_w1_fus;
    % Location of LE w1
    x_loc_LE_w1_CAD = x_loc_1R_y1_w1_CAD; % distance from CAD refference point
    y_loc_LE_w1_CAD = y_loc_1R_y1_w1_CAD; % distance from CAD refference point
    z_loc_LE_w1_CAD = z_loc_1R_y1_w1_CAD; % distance from CAD refference point
    % Correction from CAD identifying the nose as new Origin
    x_loc_LE_w1 = x_loc_LE_w1_CAD + x_offset_CAD; % corrected to the nose of the aircraft
    y_loc_LE_w1 = b_w1_fus/2;
    z_loc_LE_w1 = z_loc_LE_w1_CAD + z_offset_CAD; % corrected to the nose of the aircraft
    % Storing DATA
    Geo_tier.x_loc_LE_w1 = x_loc_LE_w1;
    Geo_tier.y_loc_LE_w1 = y_loc_LE_w1;
    Geo_tier.z_loc_LE_w1 = z_loc_LE_w1;
    Geo_tier.x_loc_LE_w1_CAD = x_loc_LE_w1_CAD;
    Geo_tier.y_loc_LE_w1_CAD = y_loc_LE_w1_CAD;
    Geo_tier.z_loc_LE_w1_CAD = z_loc_LE_w1_CAD;
    
    % Location of LE w1
    x_w1_LE = x_loc_LE_w1;
    y_w1_LE = y_loc_LE_w1;
    z_w1_LE = z_loc_LE_w1;
    % Storing DATA
    Geo_tier.x_w1_LE = x_w1_LE;
    Geo_tier.y_w1_LE = y_w1_LE;
    Geo_tier.z_w1_LE = z_w1_LE;
    
    % Angles
    Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
    Lambda_LE_w1_k1_e = Geo_tier.Lambda_LE_w1_k1_e;
    Lambda_LE_w1_k2_e = Geo_tier.Lambda_LE_w1_k2_e;
    dihedral_w1 = Geo_tier.dihedral_w1;
    dihedral_w1_e = Geo_tier.dihedral_w1_e;
    dihedral_w1_k1_e = Geo_tier.dihedral_w1_k1_e;
    dihedral_w1_k2_e = Geo_tier.dihedral_w1_k2_e;
    % CAD meassurements
    % Chord
    cR_w1 = Geo_tier.cR_w1;
    cB_k1_w1 = Geo_tier.cB_k1_w1;
    cB_k2_w1 = Geo_tier.cB_k2_w1;
    cT_w1 = Geo_tier.cT_w1;
    
    % Wingspan - w1
    b_w1_s = (b_w1_e/cos(dihedral_w1_e)); % effective wingspan in y direction
    b_w1_pv = b_w1_s*sin(dihedral_w1_e); % effective wingspan in vertical direction direction
    
    % Storing Data
    Geo_tier.b_w1_s = b_w1_s;
    Geo_tier.b_w1_pv = b_w1_pv;
    
    % Wing area 1
    S_w1_e = b_w1_e*(cR_w1+cT_w1)/2; % effective wing area: proyected in the y-plane
    S_w1_s = b_w1_s*(cR_w1+cT_w1)/2; % real wing area: proyected along the surface of wing
    S_w1_fus = cR_w1*b_w1_fus; % w1 area within the fuselage
    S_w1 = S_w1_e + S_w1_fus; % w1 Reference area
    S_w1_ph = S_w1_s*cos(dihedral_w1_e); % Proyected area of the w1 ino the horizontal plane
    S_w1_pv = S_w1_s*sin(dihedral_w1_e); % Proyected area of the w1 ino the vertical plane
    
    % Storing DATA
    Geo_tier.S_w1 = S_w1;
    Geo_tier.S_w1_e  = S_w1_e;
    Geo_tier.S_w1_s  = S_w1_s;
    Geo_tier.S_w1_fus  = S_w1_fus;
    Geo_tier.S_w1_ph  = S_w1_ph;
    Geo_tier.S_w1_pv  = S_w1_pv;
    
    %%%%% IMPORTANT %%%%
    % The user select the reference wing area that will be used
    % Reference area
    S_ref = S_w1;
    Geo_tier.S_ref = S_ref;
    
    % Aspect Ratio
    AR_w1 = b_w1^2/S_w1;
    AR_w1_e = b_w1_e^2/S_w1_e;
    AR_w1_s = b_w1_s^2/S_w1_s;
    
    % Storing DATA
    Geo_tier.AR_w1 = AR_w1;
    Geo_tier.AR_w1_e = AR_w1_e;
    Geo_tier.AR_w1_s = AR_w1_s;
    
    % Taper Ratio
    lambda_w1 = cT_w1/cR_w1;
    lambda_w1_e = cT_w1/cR_w1;
    % Storing DATA
    Geo_tier.lambda_w1 = lambda_w1;
    Geo_tier.lambda_w1_e = lambda_w1_e;
    
    %% Calculation of Sweep angles
    % Quarter Chord
    Lambda_c4_w1 = Get_Nth_Lambda(1/4,AR_w1_e,Lambda_LE_w1,lambda_w1_e);
    % Half chord
    Lambda_c2_w1 = Get_Nth_Lambda(1/2,AR_w1_e,Lambda_LE_w1,lambda_w1_e);
    % Trailing Edge
    Lambda_TE_w1 = Get_Nth_Lambda(1,AR_w1_e,Lambda_LE_w1,lambda_w1_e);
    % 3/4 of chord
    %% 
    Lambda_c34_w1 = Get_Nth_Lambda(3/4,AR_w1_e,Lambda_LE_w1,lambda_w1_e);

    % Storing DATA
    Geo_tier.Lambda_TE_w1  = Lambda_TE_w1;
    Geo_tier.Lambda_TE_w1  = Lambda_TE_w1;
    Geo_tier.Lambda_c4_w1  = Lambda_c4_w1;
    Geo_tier.Lambda_c2_w1  = Lambda_c2_w1;
    Geo_tier.Lambda_c34_w1  = Lambda_c34_w1;

    %% Geometric location of the w1 in order to define the cornes of the surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_1R_y1_w1 = 0; % defines inner position of wing chord (aileron section) LE
    x_2R_y1_w1 = x_1R_y1_w1 + cR_w1; % defines inner position of wing chord TE
    % chord at each location
    y_offset_w1 = b_w1_fus/2;
    x_1R_y2_w1 =  x_1R_y1_w1 + (b_w1/2 - y_offset_w1)*tan(Lambda_LE_w1); % defines outter position of wing chord (aileron section) LE
    x_2R_y2_w1 = x_1R_y2_w1 + cT_w1; % defines outter position of wing chord TE
    % spanwise location
    y_1R_y1_w1 = b_w1/2 - b_w1_e/2; % inner position from the center line
    y_2R_y1_w1 = y_1R_y1_w1; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    y_1R_y2_w1 = b_w1/2; % outter position from the center line
    y_2R_y2_w1 = y_1R_y2_w1; % outter position from the center line (same as LE assumes chord paralel to x-axis)
    % z-position
    z_1R_y1_w1 = 0; % inner position from the center line
    z_2R_y1_w1 = 0; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    z_1R_y2_w1 = z_1R_y1_w1 + (b_w1/2 - y_offset_w1)*tan(dihedral_w1);
    z_2R_y2_w1 = z_2R_y1_w1 + (b_w1/2 - y_offset_w1)*tan(dihedral_w1);
    % Storing DATA
    Geo_tier.x_1R_y1_w1 = x_1R_y1_w1;
    Geo_tier.x_2R_y1_w1 = x_2R_y1_w1;
    Geo_tier.x_1R_y2_w1 = x_1R_y2_w1;
    Geo_tier.x_2R_y2_w1 = x_2R_y2_w1;
    Geo_tier.y_1R_y1_w1 = y_1R_y1_w1;
    Geo_tier.y_2R_y1_w1 = y_2R_y1_w1;
    Geo_tier.y_1R_y2_w1 = y_1R_y2_w1;
    Geo_tier.y_2R_y2_w1 = y_2R_y2_w1;
    Geo_tier.z_1R_y1_w1 = z_1R_y1_w1;
    Geo_tier.z_2R_y1_w1 = z_2R_y1_w1;
    Geo_tier.z_1R_y2_w1 = z_1R_y2_w1;
    Geo_tier.z_2R_y2_w1 = z_2R_y2_w1;
    Geo_tier.y_offset_w1 = y_offset_w1;
    
    % Position relative to the Origin
    % X-position
    x_cR_w1_LE = x_w1_LE + x_1R_y1_w1;
    x_cR_w1_TE = x_w1_LE + x_2R_y1_w1;
    x_cT_w1_LE = x_w1_LE + x_1R_y2_w1;
    x_cT_w1_TE = x_w1_LE + x_2R_y2_w1;
    % Y-position
    y_cR_w1_LE = y_offset_w1;
    y_cR_w1_TE = y_cR_w1_LE;
    y_cT_w1_LE = y_cR_w1_LE + b_w1_e/2;
    y_cT_w1_TE = y_cT_w1_LE;
    % Z-position
    z_cR_w1_LE = z_loc_LE_w1 + z_1R_y1_w1;
    z_cT_w1_LE = z_loc_LE_w1 + z_1R_y2_w1;
    z_cR_w1_TE = z_loc_LE_w1 + z_2R_y1_w1;
    z_cT_w1_TE = z_loc_LE_w1 + z_2R_y2_w1;
    % Storing DATA
    Geo_tier.x_cR_w1_LE = x_cR_w1_LE;
    Geo_tier.x_cR_w1_TE = x_cR_w1_TE;
    Geo_tier.x_cT_w1_LE = x_cT_w1_LE;
    Geo_tier.x_cT_w1_TE = x_cT_w1_TE;
    Geo_tier.y_cR_w1_LE = y_cR_w1_LE;
    Geo_tier.y_cR_w1_TE = y_cR_w1_TE;
    Geo_tier.y_cT_w1_LE = y_cT_w1_LE;
    Geo_tier.y_cT_w1_TE = y_cT_w1_TE;
    Geo_tier.z_cR_w1_LE = z_cR_w1_LE;
    Geo_tier.z_cT_w1_LE = z_cT_w1_LE;
    Geo_tier.z_cR_w1_TE = z_cR_w1_TE;
    Geo_tier.z_cT_w1_TE = z_cT_w1_TE;
    
    % Geometric position of w1 MAC with respect to wing 1R
    XYZ_MAC = Get_MAC_Coordinates(b_w1,lambda_w1,cR_w1,dihedral_w1,Lambda_c4_w1);
    Geo_tier.cmac_w1 = XYZ_MAC.cbar;
    Geo_tier.xbar_w1 = XYZ_MAC.xbar_w;
    Geo_tier.ybar_w1 = XYZ_MAC.ybar_w;
    Geo_tier.zbar_w1 = XYZ_MAC.zbar_w;
    
    % Geometric position of w1 MAC with respect to wing 1R
    XYZ_MAC = Get_MAC_Coordinates(b_w1_e,lambda_w1_e,cR_w1,dihedral_w1,Lambda_c4_w1);
    Geo_tier.cmac_w1_e = XYZ_MAC.cbar;
    Geo_tier.xbar_w1_e = XYZ_MAC.xbar_w;
    Geo_tier.ybar_w1_e = XYZ_MAC.ybar_w;
    Geo_tier.zbar_w1_e = XYZ_MAC.zbar_w;
    
    [S_e, S, S_s, S_fus, S_ph, S_pv, AR, AR_e, AR_s, XYZ_MAC, XYZ_MAC_e] = get_crankedwing_geom_airbus_method(cR_w1, cB_k1_w1, cB_k2_w1, cT_w1,...
    y_loc_LE_w1_CAD, y_loc_1R_yB1_w1_CAD, y_loc_1R_yB2_w1_CAD, y_loc_1R_y2_w1_CAD, ...
    Lambda_LE_w1, Lambda_LE_w1_k1_e , Lambda_LE_w1_k2_e, dihedral_w1_e, dihedral_w1_k1_e, dihedral_w1_k2_e);

    Geo_tier.S_w1 = S;
    Geo_tier.S_w1_e = S_e;
    Geo_tier.S_w1_s = S_s;
    Geo_tier.S_w1_fus = S_fus;
    Geo_tier.S_w1_ph = S_ph;
    Geo_tier.S_w1_pv = S_pv;
    S_ref = S; 
    Geo_tier.S_ref = S_ref;

    Geo_tier.AR_w1 = AR;
    Geo_tier.AR_w1_e = AR_e;
    Geo_tier.AR_w1_s = AR_s;
    
    Geo_tier.cmac_w1 = XYZ_MAC.cbar_w;
    Geo_tier.xbar_w1 = XYZ_MAC.xbar_w;
    Geo_tier.ybar_w1 = XYZ_MAC.ybar_w;
    Geo_tier.zbar_w1 = XYZ_MAC.zbar_w;
    
    Geo_tier.cmac_w1_e = XYZ_MAC_e.cbar_w;
    Geo_tier.xbar_w1_e = XYZ_MAC_e.xbar_w;
    Geo_tier.ybar_w1_e = XYZ_MAC_e.ybar_w;
    Geo_tier.zbar_w1_e = XYZ_MAC_e.zbar_w;
    
    % Distances relative to the origin
    x_xbar_w1 = x_loc_LE_w1 + Geo_tier.xbar_w1;
    y_ybar_w1 = y_loc_LE_w1 + Geo_tier.ybar_w1;
    z_zbar_w1 = z_loc_LE_w1 + Geo_tier.zbar_w1;
    
    x_xbar_w1_e = x_loc_LE_w1 + Geo_tier.xbar_w1_e;
    y_ybar_w1_e = y_loc_LE_w1 + Geo_tier.ybar_w1_e;
    z_zbar_w1_e = z_loc_LE_w1 + Geo_tier.zbar_w1_e;
    
    % Storing DATA
    Geo_tier.x_xbar_w1 = x_xbar_w1;
    Geo_tier.y_ybar_w1 = y_ybar_w1;
    Geo_tier.z_zbar_w1 = z_zbar_w1;
    Geo_tier.x_xbar_w1_e = x_xbar_w1_e;
    Geo_tier.y_ybar_w1_e = y_ybar_w1_e;
    Geo_tier.z_zbar_w1_e = z_zbar_w1_e;

    % fuselage length
    lfus_b_w1 = l_fus/b_w1;
    Geo_tier.lfus_b_w1 = lfus_b_w1;
end

%% Canard Aerodynamic Surface
if Can == 1
    %% Aerodynamic surfaces location
    % Location of LE can
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_loc_1R_y1_can_CAD = Geo_tier.x_loc_1R_y1_can_CAD;
    y_loc_1R_y1_can_CAD = Geo_tier.y_loc_1R_y1_can_CAD;
    z_loc_1R_y1_can_CAD = Geo_tier.z_loc_1R_y1_can_CAD;
    x_loc_1R_y2_can_CAD = Geo_tier.x_loc_1R_y2_can_CAD;
    y_loc_1R_y2_can_CAD = Geo_tier.y_loc_1R_y2_can_CAD;
    z_loc_1R_y2_can_CAD = Geo_tier.z_loc_1R_y2_can_CAD;
    y_loc_1R_yB1_can_CAD = Geo_tier.y_loc_1R_yB1_can_CAD;
    y_loc_1R_yB2_can_CAD = Geo_tier.y_loc_1R_yB2_can_CAD;
    x_y1_y2_can_CAD = Geo_tier.x_y1_y2_can_CAD;
    y_y1_y2_can_CAD = Geo_tier.y_y1_y2_can_CAD;
    z_y1_y2_can_CAD = Geo_tier.z_y1_y2_can_CAD;
    
    % Wingspan - can
    b_can = Geo_tier.b_can;
    b_can_e = Geo_tier.b_can_e;
    b_can_fus = Geo_tier.b_can_fus;
    % Location of LE can
    x_loc_LE_can_CAD = x_loc_1R_y1_can_CAD; % distance from CAD refference point
    y_loc_LE_can_CAD = y_loc_1R_y1_can_CAD; % distance from CAD refference point
    z_loc_LE_can_CAD = z_loc_1R_y1_can_CAD; % distance from CAD refference point
    % Correction from CAD identifying the nose as new Origin
    x_loc_LE_can = x_loc_LE_can_CAD + x_offset_CAD; % corrected to the nose of the aircraft
    y_loc_LE_can = b_can_fus/2;
    z_loc_LE_can = z_loc_LE_can_CAD + z_offset_CAD; % corrected to the nose of the aircraft
    % Storing DATA
    Geo_tier.x_loc_LE_can = x_loc_LE_can;
    Geo_tier.y_loc_LE_can = y_loc_LE_can;
    Geo_tier.z_loc_LE_can = z_loc_LE_can;
    Geo_tier.x_loc_LE_can_CAD = x_loc_LE_can_CAD;
    Geo_tier.y_loc_LE_can_CAD = y_loc_LE_can_CAD;
    Geo_tier.z_loc_LE_can_CAD = z_loc_LE_can_CAD;
    
    % Location of LE can
    x_can_LE = x_loc_LE_can;
    y_can_LE = y_loc_LE_can;
    z_can_LE = z_loc_LE_can;
    
    % Storing DATA
    Geo_tier.x_can_LE = x_can_LE;
    Geo_tier.y_can_LE = y_can_LE;
    Geo_tier.z_can_LE = z_can_LE;
    
    % Angles
    Lambda_LE_can = Geo_tier.Lambda_LE_can;
    Lambda_LE_can_k1_e = Geo_tier.Lambda_LE_can_k1_e;
    Lambda_LE_can_k2_e = Geo_tier.Lambda_LE_can_k2_e;
    dihedral_can = Geo_tier.dihedral_can;
    dihedral_can_e = Geo_tier.dihedral_can_e;
    dihedral_can_k1_e = Geo_tier.dihedral_can_k1_e;
    dihedral_can_k2_e = Geo_tier.dihedral_can_k2_e;
    
    % CAD meassurements
    % Chord
    cR_can = Geo_tier.cR_can;
    cB_k1_can = Geo_tier.cB_k1_can;
    cB_k2_can = Geo_tier.cB_k2_can;
    cT_can = Geo_tier.cT_can;
    
    % Wingspan - can
    b_can_s = (b_can_e/cos(dihedral_can_e)); % effective wingspan in y direction
    b_can_pv = b_can_s*sin(dihedral_can_e); % effective wingspan in vertical direction direction
    
    % Storing Data
    Geo_tier.b_can_s = b_can_s;
    Geo_tier.b_can_pv = b_can_pv;
    
    % Wing area 1
    S_can_e = b_can_e*(cR_can+cT_can)/2; % effective wing area: proyected in teh y-plane
    S_can_s = b_can_s*(cR_can+cT_can)/2; % real wing area: proyected along the surface of wing
    S_can_fus = cR_can*b_can_fus; % can area within the fuselage
    S_can = S_can_e + S_can_fus; % can Reference area
    S_can_ph = S_can_s*cos(dihedral_can_e); % Proyected area of the can ino the horizontal plane
    S_can_pv = S_can_s*sin(dihedral_can_e); % Proyected area of the can ino the vertical plane
    
    % Storing DATA
    Geo_tier.S_can = S_can;
    Geo_tier.S_can_e  = S_can_e;
    Geo_tier.S_can_s  = S_can_s;
    Geo_tier.S_can_fus  = S_can_fus;
    Geo_tier.S_can_ph  = S_can_ph;
    Geo_tier.S_can_pv  = S_can_pv;
        
    % Aspect Ratio
    AR_can = b_can^2/S_can;
    AR_can_e = b_can_e^2/S_can_e;
    AR_can_s = b_can_s^2/S_can_s;
    
    % Storing DATA
    Geo_tier.AR_can = AR_can;
    Geo_tier.AR_can_e = AR_can_e;
    Geo_tier.AR_can_s = AR_can_s;
    
    % Taper Ratio
    lambda_can = cT_can/cR_can;
    lambda_can_e = cT_can/cR_can;
    % Storing DATA
    Geo_tier.lambda_can = lambda_can;
    Geo_tier.lambda_can_e = lambda_can_e;
    
    %% Calculation of Sweep angles
    % Quarter Chord
    Lambda_c4_can = Get_Nth_Lambda(1/4,AR_can_e,Lambda_LE_can,lambda_can_e);
    % Half chord
    Lambda_c2_can = Get_Nth_Lambda(1/2,AR_can_e,Lambda_LE_can,lambda_can_e);
    % Trailing Edge
    Lambda_TE_can = Get_Nth_Lambda(1,AR_can_e,Lambda_LE_can,lambda_can_e);
    % 3/4 of chord
    Lambda_c34_can = Get_Nth_Lambda(3/4,AR_can_e,Lambda_LE_can,lambda_can_e);
    
    % Storing DATA
    Geo_tier.Lambda_TE_can  = Lambda_TE_can;
    Geo_tier.Lambda_TE_can  = Lambda_TE_can;
    Geo_tier.Lambda_c4_can  = Lambda_c4_can;
    Geo_tier.Lambda_c2_can  = Lambda_c2_can;
    Geo_tier.Lambda_c34_can  = Lambda_c34_can;
    
    %% Geometric location of the can in order to define the cornes of the surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_1R_y1_can = 0; % defines inner position of wing chord (aileron section) LE
    x_2R_y1_can = x_1R_y1_can + cR_can; % defines inner position of wing chord TE
    % chord at each location
    y_offset_can = b_can_fus/2;
    x_1R_y2_can =  x_1R_y1_can + (b_can/2 - y_offset_can)*tan(Lambda_LE_can); % defines outter position of wing chord (aileron section) LE
    x_2R_y2_can = x_1R_y2_can + cT_can; % defines outter position of wing chord TE
    % spanwise location
    y_1R_y1_can = b_can/2 - b_can_e/2; % inner position from the center line
    y_2R_y1_can = y_1R_y1_can; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    y_1R_y2_can = b_can/2; % outter position from the center line
    y_2R_y2_can = y_1R_y2_can; % outter position from the center line (same as LE assumes chord paralel to x-axis)
    % z-position
    z_1R_y1_can = 0; % inner position from the center line
    z_2R_y1_can = 0; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    z_1R_y2_can = z_1R_y1_can + (b_can/2 - y_offset_can)*tan(dihedral_can);
    z_2R_y2_can = z_2R_y1_can + (b_can/2 - y_offset_can)*tan(dihedral_can);
    % Storing DATA
    Geo_tier.x_1R_y1_can = x_1R_y1_can;
    Geo_tier.x_2R_y1_can = x_2R_y1_can;
    Geo_tier.x_1R_y2_can = x_1R_y2_can;
    Geo_tier.x_2R_y2_can = x_2R_y2_can;
    Geo_tier.y_1R_y1_can = y_1R_y1_can;
    Geo_tier.y_2R_y1_can = y_2R_y1_can;
    Geo_tier.y_1R_y2_can = y_1R_y2_can;
    Geo_tier.y_2R_y2_can = y_2R_y2_can;
    Geo_tier.z_1R_y1_can = z_1R_y1_can;
    Geo_tier.z_2R_y1_can = z_2R_y1_can;
    Geo_tier.z_1R_y2_can = z_1R_y2_can;
    Geo_tier.z_2R_y2_can = z_2R_y2_can;
    Geo_tier.y_offset_can = y_offset_can;
    
    % Position relative to the Origin
    % X-position
    x_cR_can_LE = x_can_LE + x_1R_y1_can;
    x_cR_can_TE = x_can_LE + x_2R_y1_can;
    x_cT_can_LE = x_can_LE + x_1R_y2_can;
    x_cT_can_TE = x_can_LE + x_2R_y2_can;
    % Y-position
    y_cR_can_LE = y_offset_can;
    y_cR_can_TE = y_cR_can_LE;
    y_cT_can_LE = y_cR_can_LE + b_can_e/2;
    y_cT_can_TE = y_cT_can_LE;
    % Z-position
    z_cR_can_LE = z_loc_LE_can + z_1R_y1_can;
    z_cT_can_LE = z_loc_LE_can + z_1R_y2_can;
    z_cR_can_TE = z_loc_LE_can + z_2R_y1_can;
    z_cT_can_TE = z_loc_LE_can + z_2R_y2_can;
    % Storing DATA
    Geo_tier.x_cR_can_LE = x_cR_can_LE;
    Geo_tier.x_cR_can_TE = x_cR_can_TE;
    Geo_tier.x_cT_can_LE = x_cT_can_LE;
    Geo_tier.x_cT_can_TE = x_cT_can_TE;
    Geo_tier.y_cR_can_LE = y_cR_can_LE;
    Geo_tier.y_cR_can_TE = y_cR_can_TE;
    Geo_tier.y_cT_can_LE = y_cT_can_LE;
    Geo_tier.y_cT_can_TE = y_cT_can_TE;
    Geo_tier.z_cR_can_LE = z_cR_can_LE;
    Geo_tier.z_cT_can_LE = z_cT_can_LE;
    Geo_tier.z_cR_can_TE = z_cR_can_TE;
    Geo_tier.z_cT_can_TE = z_cT_can_TE;
    
    % Geometric position of can MAC with respect to wing 1R
    XYZ_MAC = Get_MAC_Coordinates(b_can,lambda_can,cR_can,dihedral_can,Lambda_c4_can);
    Geo_tier.cmac_can = XYZ_MAC.cbar;
    Geo_tier.xbar_can = XYZ_MAC.xbar_w;
    Geo_tier.ybar_can = XYZ_MAC.ybar_w;
    Geo_tier.zbar_can = XYZ_MAC.zbar_w;
    
    % Geometric position of can MAC with respect to wing 1R
    XYZ_MAC = Get_MAC_Coordinates(b_can_e,lambda_can_e,cR_can,dihedral_can,Lambda_c4_can);
    Geo_tier.cmac_can_e = XYZ_MAC.cbar;
    Geo_tier.xbar_can_e = XYZ_MAC.xbar_w;
    Geo_tier.ybar_can_e = XYZ_MAC.ybar_w;
    Geo_tier.zbar_can_e = XYZ_MAC.zbar_w;
    
    [S_e_can, S_can, S_s_can, S_fus_can, S_ph_can, S_pv_can, AR_can, AR_e_can, AR_s_can, XYZ_MAC_can, ...
    XYZ_MAC_e_can] = get_crankedwing_geom_airbus_method(cR_can, cB_k1_can, cB_k2_can, cT_can,...
    y_loc_LE_can_CAD, y_loc_1R_yB1_can_CAD, y_loc_1R_yB2_can_CAD, y_loc_1R_y2_can_CAD, ...
    Lambda_LE_can, Lambda_LE_can_k1_e , Lambda_LE_can_k2_e, dihedral_can, dihedral_can_k1_e, dihedral_can_k2_e);

    Geo_tier.S_can = S_can;
    Geo_tier.S_can_e = S_e_can;
    Geo_tier.S_can_s = S_s_can;
    Geo_tier.S_can_fus = S_fus_can;
    Geo_tier.S_can_ph = S_ph_can;
    Geo_tier.S_can_pv = S_pv_can;

    Geo_tier.AR_can = AR_can;
    Geo_tier.AR_can_e = AR_e_can;
    Geo_tier.AR_can_s = AR_s_can;
    
    Geo_tier.cmac_can = XYZ_MAC_can.cbar_w;
    Geo_tier.xbar_can = XYZ_MAC_can.xbar_w;
    Geo_tier.ybar_can = XYZ_MAC_can.ybar_w;
    Geo_tier.zbar_can = XYZ_MAC_can.zbar_w;
    
    Geo_tier.cmac_can_e = XYZ_MAC_e_can.cbar_w;
    Geo_tier.xbar_can_e = XYZ_MAC_e_can.xbar_w;
    Geo_tier.ybar_can_e = XYZ_MAC_e_can.ybar_w;
    Geo_tier.zbar_can_e = XYZ_MAC_e_can.zbar_w;
    
    % Distances relative to the origin
    x_xbar_can = x_loc_LE_can + Geo_tier.xbar_can;
    y_ybar_can = y_loc_LE_can + Geo_tier.ybar_can;
    z_zbar_can = z_loc_LE_can + Geo_tier.zbar_can;
    x_xbar_can_e = x_loc_LE_can + Geo_tier.xbar_can_e;
    y_ybar_can_e = y_loc_LE_can + Geo_tier.ybar_can_e;
    z_zbar_can_e = z_loc_LE_can + Geo_tier.zbar_can_e;
    
    % Storing DATA
    Geo_tier.x_xbar_can = x_xbar_can;
    Geo_tier.y_ybar_can = y_ybar_can;
    Geo_tier.z_zbar_can = z_zbar_can;
    Geo_tier.x_xbar_can_e = x_xbar_can_e;
    Geo_tier.y_ybar_can_e = y_ybar_can_e;
    Geo_tier.z_zbar_can_e = z_zbar_can_e;

    % fuselage length
    lfus_b_can = l_fus/b_can;
    Geo_tier.lfus_b_can = lfus_b_can;
end


%% HTP
if HTP == 1 | Vee == 1
    %% Aerodynamic surfaces location
    % Location of LE w2
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_loc_1R_y1_w2_CAD = Geo_tier.x_loc_1R_y1_w2_CAD;
    y_loc_1R_y1_w2_CAD = Geo_tier.y_loc_1R_y1_w2_CAD;
    z_loc_1R_y1_w2_CAD = Geo_tier.z_loc_1R_y1_w2_CAD;
    x_loc_1R_y2_w2_CAD = Geo_tier.x_loc_1R_y2_w2_CAD;
    y_loc_1R_y2_w2_CAD = Geo_tier.y_loc_1R_y2_w2_CAD;
    z_loc_1R_y2_w2_CAD = Geo_tier.z_loc_1R_y2_w2_CAD;
    if HTP == 1
        y_loc_1R_yB1_w2_CAD = Geo_tier.y_loc_1R_yB1_w2_CAD;
        y_loc_1R_yB2_w2_CAD = Geo_tier.y_loc_1R_yB2_w2_CAD;
    end
    x_y1_y2_w2_CAD = Geo_tier.x_y1_y2_w2_CAD;
    y_y1_y2_w2_CAD = Geo_tier.y_y1_y2_w2_CAD;
    z_y1_y2_w2_CAD = Geo_tier.z_y1_y2_w2_CAD;
    % Wingspan - w2
    b_w2 = Geo_tier.b_w2;
    b_w2_e = Geo_tier.b_w2_e;
    b_w2_fus = Geo_tier.b_w2_fus;
    
    % Location of LE w2
    x_loc_LE_w2_CAD = x_loc_1R_y1_w2_CAD; % distance from CAD refference point
    y_loc_LE_w2_CAD = y_loc_1R_y1_w2_CAD; % distance from CAD refference point
    z_loc_LE_w2_CAD = z_loc_1R_y1_w2_CAD; % distance from CAD refference point
    x_loc_LE_w2 = x_loc_LE_w2_CAD + x_offset_CAD;
    y_loc_LE_w2 = b_w2_fus/2;
    z_loc_LE_w2 = z_loc_LE_w2_CAD + z_offset_CAD;
    % Storing DATA
    Geo_tier.x_loc_LE_w2 = x_loc_LE_w2;
    Geo_tier.y_loc_LE_w2 = y_loc_LE_w2;
    Geo_tier.z_loc_LE_w2 = z_loc_LE_w2;
    Geo_tier.x_loc_LE_w2_CAD = x_loc_LE_w2_CAD;
    Geo_tier.y_loc_LE_w2_CAD = y_loc_LE_w2_CAD;
    Geo_tier.z_loc_LE_w2_CAD = z_loc_LE_w2_CAD;
    
    % Location of LE w1 and w2
    x_w2_LE = x_loc_LE_w2;
    y_w2_LE = y_loc_LE_w2;
    z_w2_LE = z_loc_LE_w2;
    % Storing DATA
    Geo_tier.x_w2_LE = x_w2_LE;
    Geo_tier.y_w2_LE = y_w2_LE;
    Geo_tier.z_w2_LE = z_w2_LE;
    
    % Angles
    Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
    dihedral_w2 = Geo_tier.dihedral_w2;
    dihedral_w2_e = Geo_tier.dihedral_w2_e;
    if HTP == 1
        Lambda_LE_w2_k1_e = Geo_tier.Lambda_LE_w2_k1_e;
        Lambda_LE_w2_k2_e = Geo_tier.Lambda_LE_w2_k2_e;
        dihedral_w2_k1_e = Geo_tier.dihedral_w2_k1_e;
        dihedral_w2_k2_e = Geo_tier.dihedral_w2_k2_e;
    end
    
    % CAD meassurements
    % Chord
    cR_w2 = Geo_tier.cR_w2;
    if HTP == 1
        cB_k1_w2 = Geo_tier.cB_k1_w2;
        cB_k2_w2 = Geo_tier.cB_k2_w2;
    end
    cT_w2 = Geo_tier.cT_w2;
    
    % Wingspan - w2
    b_w2_s = (b_w2_e/cos(dihedral_w2_e));
    % effective wingspan in y direction
    b_w2_pv = b_w2_s*sin(dihedral_w2_e); % effective wingspan in vertical direction direction
    
    % Storing Data
    Geo_tier.b_w2_s = b_w2_s;
    Geo_tier.b_w2_pv = b_w2_pv;
        
    % Wing area 2
    S_w2_e = b_w2_e*(cR_w2+cT_w2)/2; % effective wing area: proyected in teh y-plane
    S_w2_s = b_w2_s*(cR_w2+cT_w2)/2; % real wing area: proyected along the surface of wing
    S_w2_fus = cR_w2*b_w2_fus; % w2 area within the fuselage
%     S_w2 = S_w2_e; % w2 Reference area % Cola en V with no fuselage interference
    S_w2 = S_w2_e + S_w2_fus; % w2 Reference area
    S_w2_ph = S_w2_s*cos(dihedral_w2_e); % Proyected area of the w1 ino the horizontal plane
    S_w2_pv = S_w2_s*sin(dihedral_w2_e); % Proyected area of the w1 ino the vertical plane
        
    % Storing DATA
    Geo_tier.S_w2 = S_w2;
    Geo_tier.S_w2_e  = S_w2_e;
    Geo_tier.S_w2_s  = S_w2_s;
    Geo_tier.S_w2_fus  = S_w2_fus;
    Geo_tier.S_w2_ph  = S_w2_ph;
    Geo_tier.S_w2_pv  = S_w2_pv;
    
    % Aspect Ratio
    AR_w2 = b_w2^2/S_w2;
    AR_w2_e = b_w2_e^2/S_w2_e;
    AR_w2_s = b_w2_s^2/S_w2_s;
    
    % Storing DATA
    Geo_tier.AR_w2 = AR_w2;
    Geo_tier.AR_w2_e = AR_w2_e;
    Geo_tier.AR_w2_s = AR_w2_s;
    
    % Taper Ratio
    lambda_w2 = cT_w2/cR_w2;
    lambda_w2_e = cT_w2/cR_w2;
    % Storing DATA
    Geo_tier.lambda_w2 = lambda_w2;
    Geo_tier.lambda_w2_e = lambda_w2_e;
    
    %% Calculation of Sweep angles
    % Quarter Chord
    Lambda_c4_w2 = Get_Nth_Lambda(1/4,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
    % Half chord
    Lambda_c2_w2 = Get_Nth_Lambda(1/2,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
    % Trailing Edge
    Lambda_TE_w2 = Get_Nth_Lambda(1,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
    % 3/4 of chord
    Lambda_c34_w2 = Get_Nth_Lambda(3/4,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
    
    % Storing DATA
    Geo_tier.Lambda_TE_w2  = Lambda_TE_w2;
    Geo_tier.Lambda_TE_w2  = Lambda_TE_w2;
    Geo_tier.Lambda_c4_w2  = Lambda_c4_w2;
    Geo_tier.Lambda_c2_w2  = Lambda_c2_w2;
    Geo_tier.Lambda_c34_w2  = Lambda_c34_w2;
        
    %% Geometric location of the w2 in order to define the corners of the surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_1R_y1_w2 = 0; % defines inner position of wing chord (aileron section) LE
    x_2R_y1_w2 = x_1R_y1_w2 + cR_w2; % defines inner position of wing chord TE
    % chord at each location
    y_offset_w2 = b_w2_fus/2;
    x_1R_y2_w2 = x_1R_y1_w2 + (b_w2/2 - y_offset_w2)*tan(Lambda_LE_w2); % defines inner position of wing chord (ruddervator section) LE
    x_2R_y2_w2 = x_1R_y2_w2 + cT_w2; % defines inner position of wing chord TE
    % spanwise location
    y_1R_y1_w2 = b_w2/2 - b_w2_e/2; % inner position from the center line
    y_2R_y1_w2 = y_1R_y1_w2; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    y_1R_y2_w2 = b_w2/2; % outter position from the center line
    y_2R_y2_w2 = y_1R_y2_w2; % outter position from the center line (same as LE assumes chord paralel to x-axis)
    % z-position
    z_1R_y1_w2 = 0; % inner position from the center line
    z_2R_y1_w2 = 0; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    z_1R_y2_w2 = z_1R_y1_w2 + (b_w2/2 - y_offset_w2)*tan(dihedral_w2);
    z_2R_y2_w2 = z_2R_y1_w2 + (b_w2/2 - y_offset_w2)*tan(dihedral_w2);
    
    % Storing DATA
    Geo_tier.x_1R_y1_w2 = x_1R_y1_w2;
    Geo_tier.x_2R_y1_w2 = x_2R_y1_w2;
    Geo_tier.x_1R_y2_w2 = x_1R_y2_w2;
    Geo_tier.x_2R_y2_w2 = x_2R_y2_w2;
    Geo_tier.y_1R_y1_w2 = y_1R_y1_w2;
    Geo_tier.y_2R_y1_w2 = y_2R_y1_w2;
    Geo_tier.y_1R_y2_w2 = y_1R_y2_w2;
    Geo_tier.y_2R_y2_w2 = y_2R_y2_w2;
    Geo_tier.z_1R_y1_w2 = z_1R_y1_w2;
    Geo_tier.z_2R_y1_w2 = z_2R_y1_w2;
    Geo_tier.z_1R_y2_w2 = z_1R_y2_w2;
    Geo_tier.z_2R_y2_w2 = z_2R_y2_w2;
    Geo_tier.y_offset_w2 = y_offset_w2;
    
    % Position relative to the Origin
    % X-position
    x_cR_w2_LE = x_w2_LE + x_1R_y1_w2;
    x_cR_w2_TE = x_w2_LE + x_2R_y1_w2;
    x_cT_w2_LE = x_w2_LE + x_1R_y2_w2;
    x_cT_w2_TE = x_w2_LE + x_2R_y2_w2;
    % Y-position
    y_cR_w2_LE = y_offset_w2;
    y_cR_w2_TE = y_cR_w2_LE;
    y_cT_w2_LE = y_cR_w2_LE + b_w2_e/2;
    y_cT_w2_TE = y_cT_w2_LE;
    % Z-position
    
    z_cR_w2_LE = z_loc_LE_w2 + z_1R_y1_w2; % inner position from the center line
    z_cT_w2_LE = z_loc_LE_w2 + z_1R_y2_w2; % inner position from the center line
    z_cR_w2_TE = z_loc_LE_w2 + z_2R_y1_w2;
    z_cT_w2_TE = z_loc_LE_w2 + z_2R_y2_w2;
    
    % Storing DATA
    Geo_tier.x_cR_w2_LE = x_cR_w2_LE;
    Geo_tier.x_cR_w2_TE = x_cR_w2_TE;
    Geo_tier.x_cT_w2_LE = x_cT_w2_LE;
    Geo_tier.x_cT_w2_TE = x_cT_w2_TE;
    Geo_tier.y_cR_w2_LE = y_cR_w2_LE;
    Geo_tier.y_cT_w2_LE = y_cT_w2_LE;
    Geo_tier.y_cR_w2_TE = y_cR_w2_TE;
    Geo_tier.y_cT_w2_TE = y_cT_w2_TE;
    Geo_tier.z_cR_w2_LE = z_cR_w2_LE;
    Geo_tier.z_cT_w2_LE = z_cT_w2_LE;
    Geo_tier.z_cR_w2_TE = z_cR_w2_TE;
    Geo_tier.z_cT_w2_TE = z_cT_w2_TE;

    % Geometric position of w2 MAC
    XYZ_MAC = Get_MAC_Coordinates(b_w2,lambda_w2,cR_w2,dihedral_w2,Lambda_c4_w2);
    Geo_tier.cmac_w2 = XYZ_MAC.cbar;
    Geo_tier.xbar_w2 = XYZ_MAC.xbar_w;
    Geo_tier.ybar_w2 = XYZ_MAC.ybar_w;
    Geo_tier.zbar_w2 = XYZ_MAC.zbar_w;
    
    % Geometric position of w2 MAC
    XYZ_MAC = Get_MAC_Coordinates(b_w2_e,lambda_w2_e,cR_w2,dihedral_w2,Lambda_c4_w2);
    Geo_tier.cmac_w2_e = XYZ_MAC.cbar;
    Geo_tier.xbar_w2_e = XYZ_MAC.xbar_w;
    Geo_tier.ybar_w2_e = XYZ_MAC.ybar_w;
    Geo_tier.zbar_w2_e = XYZ_MAC.zbar_w;
    
    if HTP == 1
        [S_e_w2, S_w2, S_s_w2, S_fus_w2, S_ph_w2, S_pv_w2, AR_w2, AR_e_w2, AR_s_w2, XYZ_MAC_w2, ...
            XYZ_MAC_e_w2] = get_crankedwing_geom_airbus_method(cR_w2, cB_k1_w2, cB_k2_w2, cT_w2,...
            y_loc_LE_w2_CAD, y_loc_1R_yB1_w2_CAD, y_loc_1R_yB2_w2_CAD, y_loc_1R_y2_w2_CAD, ...
            Lambda_LE_w2, Lambda_LE_w2_k1_e , Lambda_LE_w2_k2_e, dihedral_w2, dihedral_w2_k1_e, dihedral_w2_k2_e);
        
        Geo_tier.S_w2 = S_w2;
        Geo_tier.S_w2_e = S_e_w2;
        Geo_tier.S_w2_s = S_s_w2;
        Geo_tier.S_w2_fus = S_fus_w2;
        Geo_tier.S_w2_ph = S_ph_w2;
        Geo_tier.S_w2_pv = S_pv_w2;
        
        Geo_tier.AR_w2 = AR_w2;
        Geo_tier.AR_w2_e = AR_e_w2;
        Geo_tier.AR_w2_s = AR_s_w2;
        
        Geo_tier.cmac_w2 = XYZ_MAC_w2.cbar_w;
        Geo_tier.xbar_w2 = XYZ_MAC_w2.xbar_w;
        Geo_tier.ybar_w2 = XYZ_MAC_w2.ybar_w;
        Geo_tier.zbar_w2 = XYZ_MAC_w2.zbar_w;
        
        Geo_tier.cmac_w2_e = XYZ_MAC_e_w2.cbar_w;
        Geo_tier.xbar_w2_e = XYZ_MAC_e_w2.xbar_w;
        Geo_tier.ybar_w2_e = XYZ_MAC_e_w2.ybar_w;
        Geo_tier.zbar_w2_e = XYZ_MAC_e_w2.zbar_w;
        
        
    end
    % Distances relative to the origin
    x_xbar_w2 = x_loc_LE_w2 + Geo_tier.xbar_w2;
    y_ybar_w2 = y_loc_LE_w2 + Geo_tier.ybar_w2;
    z_zbar_w2 = z_loc_LE_w2 + Geo_tier.zbar_w2;
    x_xbar_w2_e = x_loc_LE_w2 + Geo_tier.xbar_w2_e;
    y_ybar_w2_e = y_loc_LE_w2 + Geo_tier.ybar_w2_e;
    z_zbar_w2_e = z_loc_LE_w2 + Geo_tier.zbar_w2_e;
    
    % Storing DATA
    Geo_tier.x_xbar_w2 = x_xbar_w2;
    Geo_tier.y_ybar_w2 = y_ybar_w2;
    Geo_tier.z_zbar_w2 = z_zbar_w2;
    Geo_tier.x_xbar_w2_e = x_xbar_w2_e;
    Geo_tier.y_ybar_w2_e = y_ybar_w2_e;
    Geo_tier.z_zbar_w2_e = z_zbar_w2_e;
    

    %% Tail Volume Coefficients
    l_xac_w1w2 = x_xbar_w2 - x_xbar_w1; % from xac_wing1 to xac_wing2
    % tail volume coefficient
    Cw2 = S_w2*l_xac_w1w2/(Geo_tier.cmac_w2*S_ref);
    
    lfus_b_HTP = l_fus/b_w2;
    Geo_tier.lfus_b_w2 = lfus_b_HTP;
    
    % Storing DATA
    Geo_tier.l_xac_w1w2 = l_xac_w1w2;
    Geo_tier.Cw2 = Cw2;
end

%% HTP
if HTP2 == 1 | Vee2 == 1
    %% Aerodynamic surfaces location
    % Location of LE w2
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_loc_1R_y1_w2_CAD = Geo_tier.x_loc_1R_y1_w2_CAD;
    y_loc_1R_y1_w2_CAD = Geo_tier.y_loc_1R_y1_w2_CAD;
    z_loc_1R_y1_w2_CAD = Geo_tier.z_loc_1R_y1_w2_CAD;
    x_loc_1R_y2_w2_CAD = Geo_tier.x_loc_1R_y2_w2_CAD;
    y_loc_1R_y2_w2_CAD = Geo_tier.y_loc_1R_y2_w2_CAD;
    z_loc_1R_y2_w2_CAD = Geo_tier.z_loc_1R_y2_w2_CAD;
    if HTP == 1
        y_loc_1R_yB1_w2_CAD = Geo_tier.y_loc_1R_yB1_w2_CAD;
        y_loc_1R_yB2_w2_CAD = Geo_tier.y_loc_1R_yB2_w2_CAD;
    end
    x_y1_y2_w2_CAD = Geo_tier.x_y1_y2_w2_CAD;
    y_y1_y2_w2_CAD = Geo_tier.y_y1_y2_w2_CAD;
    z_y1_y2_w2_CAD = Geo_tier.z_y1_y2_w2_CAD;
    % Wingspan - w2
    b_w2 = Geo_tier.b_w2;
    b_w2_e = Geo_tier.b_w2_e;
    b_w2_fus = Geo_tier.b_w2_fus;
    
    % Location of LE w2
    x_loc_LE_w2_CAD = x_loc_1R_y1_w2_CAD; % distance from CAD refference point
    y_loc_LE_w2_CAD = y_loc_1R_y1_w2_CAD; % distance from CAD refference point
    z_loc_LE_w2_CAD = z_loc_1R_y1_w2_CAD; % distance from CAD refference point
    x_loc_LE_w2 = x_loc_LE_w2_CAD + x_offset_CAD;
    y_loc_LE_w2 = b_w2_fus/2;
    z_loc_LE_w2 = z_loc_LE_w2_CAD + z_offset_CAD;
    % Storing DATA
    Geo_tier.x_loc_LE_w2 = x_loc_LE_w2;
    Geo_tier.y_loc_LE_w2 = y_loc_LE_w2;
    Geo_tier.z_loc_LE_w2 = z_loc_LE_w2;
    Geo_tier.x_loc_LE_w2_CAD = x_loc_LE_w2_CAD;
    Geo_tier.y_loc_LE_w2_CAD = y_loc_LE_w2_CAD;
    Geo_tier.z_loc_LE_w2_CAD = z_loc_LE_w2_CAD;
    
    % Location of LE w1 and w2
    x_w2_LE = x_loc_LE_w2;
    y_w2_LE = y_loc_LE_w2;
    z_w2_LE = z_loc_LE_w2;
    % Storing DATA
    Geo_tier.x_w2_LE = x_w2_LE;
    Geo_tier.y_w2_LE = y_w2_LE;
    Geo_tier.z_w2_LE = z_w2_LE;
    
    % Angles
    Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
    dihedral_w2 = Geo_tier.dihedral_w2;
    dihedral_w2_e = Geo_tier.dihedral_w2_e;
    if HTP == 1
        Lambda_LE_w2_k1_e = Geo_tier.Lambda_LE_w2_k1_e;
        Lambda_LE_w2_k2_e = Geo_tier.Lambda_LE_w2_k2_e;
        dihedral_w2_k1_e = Geo_tier.dihedral_w2_k1_e;
        dihedral_w2_k2_e = Geo_tier.dihedral_w2_k2_e;
    end
    
    % CAD meassurements
    % Chord
    cR_w2 = Geo_tier.cR_w2;
    if HTP == 1
        cB_k1_w2 = Geo_tier.cB_k1_w2;
        cB_k2_w2 = Geo_tier.cB_k2_w2;
    end
    cT_w2 = Geo_tier.cT_w2;
    
    % Wingspan - w2
    b_w2_s = (b_w2_e/cos(dihedral_w2_e));
    % effective wingspan in y direction
    b_w2_pv = b_w2_s*sin(dihedral_w2_e); % effective wingspan in vertical direction direction
    
    % Storing Data
    Geo_tier.b_w2_s = b_w2_s;
    Geo_tier.b_w2_pv = b_w2_pv;
        
    % Wing area 2
    S_w2_e = b_w2_e*(cR_w2+cT_w2)/2; % effective wing area: proyected in teh y-plane
    S_w2_s = b_w2_s*(cR_w2+cT_w2)/2; % real wing area: proyected along the surface of wing
    S_w2_fus = cR_w2*b_w2_fus; % w2 area within the fuselage
%     S_w2 = S_w2_e; % w2 Reference area % Cola en V with no fuselage interference
    S_w2 = S_w2_e + S_w2_fus; % w2 Reference area
    S_w2_ph = S_w2_s*cos(dihedral_w2_e); % Proyected area of the w1 ino the horizontal plane
    S_w2_pv = S_w2_s*sin(dihedral_w2_e); % Proyected area of the w1 ino the vertical plane
        
    % Storing DATA
    Geo_tier.S_w2 = S_w2;
    Geo_tier.S_w2_e  = S_w2_e;
    Geo_tier.S_w2_s  = S_w2_s;
    Geo_tier.S_w2_fus  = S_w2_fus;
    Geo_tier.S_w2_ph  = S_w2_ph;
    Geo_tier.S_w2_pv  = S_w2_pv;
    
    % Aspect Ratio
    AR_w2 = b_w2^2/S_w2;
    AR_w2_e = b_w2_e^2/S_w2_e;
    AR_w2_s = b_w2_s^2/S_w2_s;
    
    % Storing DATA
    Geo_tier.AR_w2 = AR_w2;
    Geo_tier.AR_w2_e = AR_w2_e;
    Geo_tier.AR_w2_s = AR_w2_s;
    
    % Taper Ratio
    lambda_w2 = cT_w2/cR_w2;
    lambda_w2_e = cT_w2/cR_w2;
    % Storing DATA
    Geo_tier.lambda_w2 = lambda_w2;
    Geo_tier.lambda_w2_e = lambda_w2_e;
    
    %% Calculation of Sweep angles
    % Quarter Chord
    Lambda_c4_w2 = Get_Nth_Lambda(1/4,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
    % Half chord
    Lambda_c2_w2 = Get_Nth_Lambda(1/2,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
    % Trailing Edge
    Lambda_TE_w2 = Get_Nth_Lambda(1,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
    % 3/4 of chord
    Lambda_c34_w2 = Get_Nth_Lambda(3/4,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
    
    % Storing DATA
    Geo_tier.Lambda_TE_w2  = Lambda_TE_w2;
    Geo_tier.Lambda_TE_w2  = Lambda_TE_w2;
    Geo_tier.Lambda_c4_w2  = Lambda_c4_w2;
    Geo_tier.Lambda_c2_w2  = Lambda_c2_w2;
    Geo_tier.Lambda_c34_w2  = Lambda_c34_w2;
        
    %% Geometric location of the w2 in order to define the corners of the surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_1R_y1_w2 = 0; % defines inner position of wing chord (aileron section) LE
    x_2R_y1_w2 = x_1R_y1_w2 + cR_w2; % defines inner position of wing chord TE
    % chord at each location
    y_offset_w2 = b_w2_fus/2;
    x_1R_y2_w2 = x_1R_y1_w2 + (b_w2/2 - y_offset_w2)*tan(Lambda_LE_w2); % defines inner position of wing chord (ruddervator section) LE
    x_2R_y2_w2 = x_1R_y2_w2 + cT_w2; % defines inner position of wing chord TE
    % spanwise location
    y_1R_y1_w2 = b_w2/2 - b_w2_e/2; % inner position from the center line
    y_2R_y1_w2 = y_1R_y1_w2; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    y_1R_y2_w2 = b_w2/2; % outter position from the center line
    y_2R_y2_w2 = y_1R_y2_w2; % outter position from the center line (same as LE assumes chord paralel to x-axis)
    % z-position
    z_1R_y1_w2 = 0; % inner position from the center line
    z_2R_y1_w2 = 0; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    z_1R_y2_w2 = z_1R_y1_w2 + (b_w2/2 - y_offset_w2)*tan(dihedral_w2);
    z_2R_y2_w2 = z_2R_y1_w2 + (b_w2/2 - y_offset_w2)*tan(dihedral_w2);
    
    % Storing DATA
    Geo_tier.x_1R_y1_w2 = x_1R_y1_w2;
    Geo_tier.x_2R_y1_w2 = x_2R_y1_w2;
    Geo_tier.x_1R_y2_w2 = x_1R_y2_w2;
    Geo_tier.x_2R_y2_w2 = x_2R_y2_w2;
    Geo_tier.y_1R_y1_w2 = y_1R_y1_w2;
    Geo_tier.y_2R_y1_w2 = y_2R_y1_w2;
    Geo_tier.y_1R_y2_w2 = y_1R_y2_w2;
    Geo_tier.y_2R_y2_w2 = y_2R_y2_w2;
    Geo_tier.z_1R_y1_w2 = z_1R_y1_w2;
    Geo_tier.z_2R_y1_w2 = z_2R_y1_w2;
    Geo_tier.z_1R_y2_w2 = z_1R_y2_w2;
    Geo_tier.z_2R_y2_w2 = z_2R_y2_w2;
    Geo_tier.y_offset_w2 = y_offset_w2;
    
    % Position relative to the Origin
    % X-position
    x_cR_w2_LE = x_w2_LE + x_1R_y1_w2;
    x_cR_w2_TE = x_w2_LE + x_2R_y1_w2;
    x_cT_w2_LE = x_w2_LE + x_1R_y2_w2;
    x_cT_w2_TE = x_w2_LE + x_2R_y2_w2;
    % Y-position
    y_cR_w2_LE = y_offset_w2;
    y_cR_w2_TE = y_cR_w2_LE;
    y_cT_w2_LE = y_cR_w2_LE + b_w2_e/2;
    y_cT_w2_TE = y_cT_w2_LE;
    % Z-position
    
    z_cR_w2_LE = z_loc_LE_w2 + z_1R_y1_w2; % inner position from the center line
    z_cT_w2_LE = z_loc_LE_w2 + z_1R_y2_w2; % inner position from the center line
    z_cR_w2_TE = z_loc_LE_w2 + z_2R_y1_w2;
    z_cT_w2_TE = z_loc_LE_w2 + z_2R_y2_w2;
    
    % Storing DATA
    Geo_tier.x_cR_w2_LE = x_cR_w2_LE;
    Geo_tier.x_cR_w2_TE = x_cR_w2_TE;
    Geo_tier.x_cT_w2_LE = x_cT_w2_LE;
    Geo_tier.x_cT_w2_TE = x_cT_w2_TE;
    Geo_tier.y_cR_w2_LE = y_cR_w2_LE;
    Geo_tier.y_cT_w2_LE = y_cT_w2_LE;
    Geo_tier.y_cR_w2_TE = y_cR_w2_TE;
    Geo_tier.y_cT_w2_TE = y_cT_w2_TE;
    Geo_tier.z_cR_w2_LE = z_cR_w2_LE;
    Geo_tier.z_cT_w2_LE = z_cT_w2_LE;
    Geo_tier.z_cR_w2_TE = z_cR_w2_TE;
    Geo_tier.z_cT_w2_TE = z_cT_w2_TE;

    % Geometric position of w2 MAC
    XYZ_MAC = Get_MAC_Coordinates(b_w2,lambda_w2,cR_w2,dihedral_w2,Lambda_c4_w2);
    Geo_tier.cmac_w2 = XYZ_MAC.cbar;
    Geo_tier.xbar_w2 = XYZ_MAC.xbar_w;
    Geo_tier.ybar_w2 = XYZ_MAC.ybar_w;
    Geo_tier.zbar_w2 = XYZ_MAC.zbar_w;
    
    % Geometric position of w2 MAC
    XYZ_MAC = Get_MAC_Coordinates(b_w2_e,lambda_w2_e,cR_w2,dihedral_w2,Lambda_c4_w2);
    Geo_tier.cmac_w2_e = XYZ_MAC.cbar;
    Geo_tier.xbar_w2_e = XYZ_MAC.xbar_w;
    Geo_tier.ybar_w2_e = XYZ_MAC.ybar_w;
    Geo_tier.zbar_w2_e = XYZ_MAC.zbar_w;
    
    if HTP == 1
        [S_e_w2, S_w2, S_s_w2, S_fus_w2, S_ph_w2, S_pv_w2, AR_w2, AR_e_w2, AR_s_w2, XYZ_MAC_w2, ...
            XYZ_MAC_e_w2] = get_crankedwing_geom_airbus_method(cR_w2, cB_k1_w2, cB_k2_w2, cT_w2,...
            y_loc_LE_w2_CAD, y_loc_1R_yB1_w2_CAD, y_loc_1R_yB2_w2_CAD, y_loc_1R_y2_w2_CAD, ...
            Lambda_LE_w2, Lambda_LE_w2_k1_e , Lambda_LE_w2_k2_e, dihedral_w2, dihedral_w2_k1_e, dihedral_w2_k2_e);
        
        Geo_tier.S_w2 = S_w2;
        Geo_tier.S_w2_e = S_e_w2;
        Geo_tier.S_w2_s = S_s_w2;
        Geo_tier.S_w2_fus = S_fus_w2;
        Geo_tier.S_w2_ph = S_ph_w2;
        Geo_tier.S_w2_pv = S_pv_w2;
        
        Geo_tier.AR_w2 = AR_w2;
        Geo_tier.AR_w2_e = AR_e_w2;
        Geo_tier.AR_w2_s = AR_s_w2;
        
        Geo_tier.cmac_w2 = XYZ_MAC_w2.cbar_w;
        Geo_tier.xbar_w2 = XYZ_MAC_w2.xbar_w;
        Geo_tier.ybar_w2 = XYZ_MAC_w2.ybar_w;
        Geo_tier.zbar_w2 = XYZ_MAC_w2.zbar_w;
        
        Geo_tier.cmac_w2_e = XYZ_MAC_e_w2.cbar_w;
        Geo_tier.xbar_w2_e = XYZ_MAC_e_w2.xbar_w;
        Geo_tier.ybar_w2_e = XYZ_MAC_e_w2.ybar_w;
        Geo_tier.zbar_w2_e = XYZ_MAC_e_w2.zbar_w;
        
        
    end
    % Distances relative to the origin
    x_xbar_w2 = x_loc_LE_w2 + Geo_tier.xbar_w2;
    y_ybar_w2 = y_loc_LE_w2 + Geo_tier.ybar_w2;
    z_zbar_w2 = z_loc_LE_w2 + Geo_tier.zbar_w2;
    x_xbar_w2_e = x_loc_LE_w2 + Geo_tier.xbar_w2_e;
    y_ybar_w2_e = y_loc_LE_w2 + Geo_tier.ybar_w2_e;
    z_zbar_w2_e = z_loc_LE_w2 + Geo_tier.zbar_w2_e;
    
    % Storing DATA
    Geo_tier.x_xbar_w2 = x_xbar_w2;
    Geo_tier.y_ybar_w2 = y_ybar_w2;
    Geo_tier.z_zbar_w2 = z_zbar_w2;
    Geo_tier.x_xbar_w2_e = x_xbar_w2_e;
    Geo_tier.y_ybar_w2_e = y_ybar_w2_e;
    Geo_tier.z_zbar_w2_e = z_zbar_w2_e;
    

    %% Tail Volume Coefficients
    l_xac_w1w2 = x_xbar_w2 - x_xbar_w1; % from xac_wing1 to xac_wing2
    % tail volume coefficient
    Cw2 = S_w2*l_xac_w1w2/(Geo_tier.cmac_w2*S_ref);
    
    lfus_b_HTP = l_fus/b_w2;
    Geo_tier.lfus_b_w2 = lfus_b_HTP;
    
    % Storing DATA
    Geo_tier.l_xac_w1w2 = l_xac_w1w2;
    Geo_tier.Cw2 = Cw2;
end
%% VTP
if VTP == 1
    %% Aerodynamic surfaces location
    % Location of LE VTP
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_loc_1R_y1_VTP_CAD = Geo_tier.x_loc_1R_y1_VTP_CAD;
    y_loc_1R_y1_VTP_CAD = Geo_tier.y_loc_1R_y1_VTP_CAD;
    z_loc_1R_y1_VTP_CAD = Geo_tier.z_loc_1R_y1_VTP_CAD;
    x_loc_1R_y2_VTP_CAD = Geo_tier.x_loc_1R_y2_VTP_CAD;
    y_loc_1R_y2_VTP_CAD = Geo_tier.y_loc_1R_y2_VTP_CAD;
    z_loc_1R_y2_VTP_CAD = Geo_tier.z_loc_1R_y2_VTP_CAD;
    x_y1_y2_VTP_CAD = Geo_tier.x_y1_y2_VTP_CAD;
    y_y1_y2_VTP_CAD = Geo_tier.y_y1_y2_VTP_CAD;
    z_y1_y2_VTP_CAD = Geo_tier.z_y1_y2_VTP_CAD;
    % Wingspan - VTP
    b_VTP = Geo_tier.b_VTP;
    b_VTP_e = Geo_tier.b_VTP_e;
    b_VTP_fus = Geo_tier.b_VTP_fus;
    
    % Location of LE VTP
    x_loc_LE_VTP_CAD = x_loc_1R_y1_VTP_CAD; % distance from CAD refference point
    y_loc_LE_VTP_CAD = y_loc_1R_y1_VTP_CAD; % distance from CAD refference point
    z_loc_LE_VTP_CAD = z_loc_1R_y1_VTP_CAD; % distance from CAD refference point
    % Correction from CAD identifying the nose as new Origin
    x_loc_LE_VTP = x_loc_LE_VTP_CAD + x_offset_CAD;
    y_loc_LE_VTP = 0;
    z_loc_LE_VTP = z_loc_LE_VTP_CAD + z_offset_CAD;
    % Storing DATA
    Geo_tier.x_loc_LE_VTP = x_loc_LE_VTP;
    Geo_tier.y_loc_LE_VTP = y_loc_LE_VTP;
    Geo_tier.z_loc_LE_VTP = z_loc_LE_VTP;
    Geo_tier.x_loc_LE_VTP_CAD = x_loc_LE_VTP_CAD;
    Geo_tier.y_loc_LE_VTP_CAD = y_loc_LE_VTP_CAD;
    Geo_tier.z_loc_LE_VTP_CAD = z_loc_LE_VTP_CAD;
    
    % Location of LE w1 and w2
    x_VTP_LE = x_loc_LE_VTP;
    y_VTP_LE = y_loc_LE_VTP;
    z_VTP_LE = z_loc_LE_VTP;
    % Storing DATA
    Geo_tier.x_VTP_LE = x_VTP_LE;
    Geo_tier.y_VTP_LE = y_VTP_LE;
    Geo_tier.z_VTP_LE = z_VTP_LE;
    
    % Angles
    Lambda_LE_VTP = Geo_tier.Lambda_LE_VTP;
    dihedral_VTP = Geo_tier.dihedral_VTP;
    dihedral_VTP_e = Geo_tier.dihedral_VTP_e;
    
    % CAD meassurements
    % Chord
    cR_VTP = Geo_tier.cR_VTP;
    cT_VTP = Geo_tier.cT_VTP;
    
    % Wingspan - VTP
    b_VTP_s = (b_VTP*cos(dihedral_VTP_e)); % effective wingspan along the surface direction
    b_VTP_pv = b_VTP*cos(dihedral_VTP_e); % effective wingspan in vertical direction direction
    b_VTP_ph = b_VTP*sin(dihedral_VTP_e); % effective wingspan in vertical direction direction

    % Storing Data
    Geo_tier.b_VTP_s = b_VTP_s;
    Geo_tier.b_VTP_pv = b_VTP_pv;
    
    if twin_VTP == 1
        S_VTP_e = b_VTP_e*(cR_VTP + cT_VTP)/2; % effective wing area: proyected in teh y-plane
        S_VTP_s = b_VTP_s*(cR_VTP + cT_VTP)/2; % real wing area: proyected along the surface of wing
        S_VTP_fus = cR_VTP*b_VTP_fus; % w2 area within the fuselage
        
        S_VTP = S_VTP_e; % w2 Reference area % Cola en V with no fuselage interference
        S_VTP_ph = S_VTP_s*sin(dihedral_VTP_e); % Proyected area of the w1 ino the horizontal plane
        S_VTP_pv = S_VTP_s*cos(dihedral_VTP_e); % Proyected area of the w1 ino the vertical plane
        
        % Storing the value for just one of the two VTP
        S_VTP_e1 = b_VTP_e*(cR_VTP+cT_VTP)/2; % effective wing area: proyected in teh y-plane
        S_VTP_s1 = b_VTP_s*(cR_VTP+cT_VTP)/2; % real wing area: proyected along the surface of wing
        S_VTP_fus1 = cR_VTP*b_VTP_fus; % w2 area within the fuselage
        
        S_VTP1 = S_VTP_e1; % w2 Reference area % Cola en V with no fuselage interference
        S_VTP_ph1 = S_VTP_s*sin(dihedral_VTP_e); % Proyected area of the w1 ino the horizontal plane
        S_VTP_pv1 = S_VTP_s*cos(dihedral_VTP_e); % Proyected area of the w1 ino the vertical plane

        % Storing DATA
        Geo_tier.S_VTP = S_VTP;
        Geo_tier.S_VTP_e  = S_VTP_e;
        Geo_tier.S_VTP_s  = S_VTP_s;
        Geo_tier.S_VTP_fus  = S_VTP_fus;
        Geo_tier.S_VTP_ph  = S_VTP_ph;
        Geo_tier.S_VTP_pv  = S_VTP_pv;
        
        % Storing DATA
        Geo_tier.S_VTP1 = S_VTP1;
        Geo_tier.S_VTP_e1  = S_VTP_e1;
        Geo_tier.S_VTP_s1  = S_VTP_s1;
        Geo_tier.S_VTP_fus1  = S_VTP_fus1;
        Geo_tier.S_VTP_ph1  = S_VTP_ph1;
        Geo_tier.S_VTP_pv1  = S_VTP_pv1;

    else
        S_VTP_e = b_VTP_e*(cR_VTP + cT_VTP)/2; % effective wing area: proyected in teh y-plane
        S_VTP_s = b_VTP_s*(cR_VTP + cT_VTP)/2; % real wing area: proyected along the surface of wing
        S_VTP_fus = 2*cR_VTP*b_VTP_fus; % w2 area within the fuselage
        
        S_VTP = S_VTP_e; % w2 Reference area % Cola en V with no fuselage interference
        S_VTP_ph = S_VTP_s*sin(dihedral_VTP_e); % Proyected area of the w1 ino the horizontal plane
        S_VTP_pv = S_VTP_s*cos(dihedral_VTP_e); % Proyected area of the w1 ino the vertical plane

        % Storing DATA
        Geo_tier.S_VTP = S_VTP;
        Geo_tier.S_VTP_e  = S_VTP_e;
        Geo_tier.S_VTP_s  = S_VTP_s;
        Geo_tier.S_VTP_fus  = S_VTP_fus;
        Geo_tier.S_VTP_ph  = S_VTP_ph;
        Geo_tier.S_VTP_pv  = S_VTP_pv;
    end
    
    if twin_VTP == 1
        % Aspect Ratio
        AR_VTP = b_VTP^2/S_VTP;
        AR_VTP_e = b_VTP_e^2/S_VTP_e;
        AR_VTP_s = b_VTP_s^2/S_VTP_s;
    else
        % Aspect Ratio
        AR_VTP = b_VTP^2/S_VTP;
        AR_VTP_e = b_VTP_e^2/S_VTP_e;
        AR_VTP_s = b_VTP_s^2/S_VTP_s;
    end
    
    % Storing DATA
    Geo_tier.AR_VTP = AR_VTP;
    Geo_tier.AR_VTP_e = AR_VTP_e;
    Geo_tier.AR_VTP_s = AR_VTP_s;
    
    % Taper Ratio
    lambda_VTP = cT_VTP/cR_VTP;
    lambda_VTP_e = cT_VTP/cR_VTP;
    % Storing DATA
    Geo_tier.lambda_VTP = lambda_VTP;
    Geo_tier.lambda_VTP_e = lambda_VTP_e;
    
    %% Calculation of Sweep angles
    % Quarter Chord
    Lambda_c4_VTP = Get_Nth_Lambda(1/4,AR_VTP_e,Lambda_LE_VTP,lambda_VTP_e);
    % Half chord
    Lambda_c2_VTP = Get_Nth_Lambda(1/2,AR_VTP_e,Lambda_LE_VTP,lambda_VTP_e);
    % Trailing Edge
    Lambda_TE_VTP = Get_Nth_Lambda(1,AR_VTP_e,Lambda_LE_VTP,lambda_VTP_e);
    % 3/4 of chord
    Lambda_c34_VTP = Get_Nth_Lambda(3/4,AR_VTP_e,Lambda_LE_VTP,lambda_VTP_e);
    
    % Storing DATA
    Geo_tier.Lambda_TE_VTP  = Lambda_TE_VTP;
    Geo_tier.Lambda_TE_VTP  = Lambda_TE_VTP;
    Geo_tier.Lambda_c4_VTP  = Lambda_c4_VTP;
    Geo_tier.Lambda_c2_VTP  = Lambda_c2_VTP;
    Geo_tier.Lambda_c34_VTP  = Lambda_c34_VTP;
        
    %% Geometric location of the VTP in order to define the corners of the surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    x_1R_y1_VTP = 0; % defines inner position of wing chord (aileron section) LE
    x_2R_y1_VTP = x_1R_y1_VTP + cR_VTP; % defines inner position of wing chord TE
    % chord at each location
%     y_offset_VTP = b_VTP_fus/2;
    y_offset_VTP = Geo_tier.y_loc_1R_y1_2VTP_CAD;
    z_offset_VTP = 0;
    x_1R_y2_VTP = x_1R_y1_VTP + (b_VTP)*tan(Lambda_LE_VTP);

    % defines inner position of wing chord (rudder section) LE
    x_2R_y2_VTP = x_1R_y2_VTP + cT_VTP; % defines inner position of wing chord TE
    % spanwise location 
    y_1R_y1_VTP = b_VTP - b_VTP_e; % inner position from the center line
    y_1R_y1_VTP = 0; % inner position from the center line
    y_2R_y1_VTP = y_1R_y1_VTP; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    y_1R_y2_VTP = b_VTP; % outter position from the center line
    y_1R_y2_VTP = b_VTP_s*sin(dihedral_VTP_e); % inner position from the center line
    y_2R_y2_VTP = y_1R_y2_VTP; % outter position from the center line (same as LE assumes chord paralel to x-axis)
    % z-position
    z_1R_y1_VTP = 0; % inner position from the center line
    z_2R_y1_VTP = 0; % inner position from the center line (same as LE assumes chord paralel to x-axis)
    z_1R_y2_VTP = z_1R_y1_VTP + (b_VTP - z_offset_VTP)*cos(dihedral_VTP);
    z_2R_y2_VTP = z_2R_y1_VTP + (b_VTP - z_offset_VTP)*cos(dihedral_VTP);
    

    % Storing DATA
    Geo_tier.x_1R_y1_VTP = x_1R_y1_VTP;
    Geo_tier.x_2R_y1_VTP = x_2R_y1_VTP;
    Geo_tier.x_1R_y2_VTP = x_1R_y2_VTP;
    Geo_tier.x_2R_y2_VTP = x_2R_y2_VTP;
    Geo_tier.y_1R_y1_VTP = y_1R_y1_VTP;
    Geo_tier.y_2R_y1_VTP = y_2R_y1_VTP;
    Geo_tier.y_1R_y2_VTP = y_1R_y2_VTP;
    Geo_tier.y_2R_y2_VTP = y_2R_y2_VTP;
    Geo_tier.z_1R_y1_VTP = z_1R_y1_VTP;
    Geo_tier.z_2R_y1_VTP = z_2R_y1_VTP;
    Geo_tier.z_1R_y2_VTP = z_1R_y2_VTP;
    Geo_tier.z_2R_y2_VTP = z_2R_y2_VTP;
    Geo_tier.y_offset_VTP = y_offset_VTP;
    Geo_tier.z_offset_VTP = z_offset_VTP;
    
    % Position relative to the Origin
    % X-position
    x_cR_VTP_LE = x_VTP_LE + x_1R_y1_VTP;
    x_cR_VTP_TE = x_VTP_LE + x_2R_y1_VTP;
    x_cT_VTP_LE = x_VTP_LE + x_1R_y2_VTP;
    x_cT_VTP_TE = x_VTP_LE + x_2R_y2_VTP;
    % Y-position
    y_cR_VTP_LE = y_offset_VTP;
    y_cR_VTP_TE = y_cR_VTP_LE;
    y_cT_VTP_LE = y_cR_VTP_LE + b_VTP_e*sin(dihedral_VTP);
    y_cT_VTP_TE = y_cT_VTP_LE;
    % Z-position
    z_cR_VTP_LE = z_loc_LE_VTP + z_1R_y1_VTP; % inner position from the center line
    z_cT_VTP_LE = z_loc_LE_VTP + z_1R_y2_VTP; % inner position from the center line
    z_cR_VTP_TE = z_loc_LE_VTP + z_2R_y1_VTP;
    z_cT_VTP_TE = z_loc_LE_VTP + z_2R_y2_VTP;
    % Storing DATA
    Geo_tier.x_cR_VTP_LE = x_cR_VTP_LE;
    Geo_tier.x_cR_VTP_TE = x_cR_VTP_TE;
    Geo_tier.x_cT_VTP_LE = x_cT_VTP_LE;
    Geo_tier.x_cT_VTP_TE = x_cT_VTP_TE;
    Geo_tier.y_cR_VTP_LE = y_cR_VTP_LE;
    Geo_tier.y_cT_VTP_LE = y_cT_VTP_LE;
    Geo_tier.y_cR_VTP_TE = y_cR_VTP_TE;
    Geo_tier.y_cT_VTP_TE = y_cT_VTP_TE;
    Geo_tier.z_cR_VTP_LE = z_cR_VTP_LE;
    Geo_tier.z_cT_VTP_LE = z_cT_VTP_LE;
    Geo_tier.z_cR_VTP_TE = z_cR_VTP_TE;
    Geo_tier.z_cT_VTP_TE = z_cT_VTP_TE;
    
    % Geometric position of VTP MAC
    if VTP == 1
        XYZ_MAC = Get_MAC_Coordinates_VTP(b_VTP,lambda_VTP,cR_VTP,dihedral_VTP,Lambda_c4_VTP);
        Geo_tier.cmac_VTP = XYZ_MAC.cbar;
        Geo_tier.xbar_VTP = XYZ_MAC.xbar_w;
        Geo_tier.zbar_VTP = XYZ_MAC.zbar_w;
        Geo_tier.ybar_VTP = XYZ_MAC.ybar_w;
        % Geometric position of VTP MAC
        XYZ_MAC = Get_MAC_Coordinates_VTP(b_VTP_e,lambda_VTP_e,cR_VTP,dihedral_VTP,Lambda_c4_VTP);
        Geo_tier.cmac_VTP_e = XYZ_MAC.cbar;
        Geo_tier.xbar_VTP_e = XYZ_MAC.xbar_w;
        Geo_tier.ybar_VTP_e = XYZ_MAC.ybar_w;
        Geo_tier.zbar_VTP_e = XYZ_MAC.zbar_w;
    end
    
    % Distances relative to the origin
    x_xbar_VTP = x_loc_LE_VTP + Geo_tier.xbar_VTP;
    y_ybar_VTP = y_loc_LE_VTP + Geo_tier.ybar_VTP;
    z_zbar_VTP = z_loc_LE_VTP + Geo_tier.zbar_VTP;
    x_xbar_VTP_e = x_loc_LE_VTP + Geo_tier.xbar_VTP_e;
    y_ybar_VTP_e = y_loc_LE_VTP + Geo_tier.ybar_VTP_e;
    z_zbar_VTP_e = z_loc_LE_VTP + Geo_tier.zbar_VTP_e;
    
    % Storing DATA
    Geo_tier.x_xbar_VTP = x_xbar_VTP;
    Geo_tier.y_ybar_VTP = y_ybar_VTP;
    Geo_tier.z_zbar_VTP = z_zbar_VTP;
    Geo_tier.x_xbar_VTP_e = x_xbar_VTP_e;
    Geo_tier.y_ybar_VTP_e = y_ybar_VTP_e;
    Geo_tier.z_zbar_VTP_e = z_zbar_VTP_e;
    
    %% Tail Volume Coefficients
    l_xac_w1VTP = x_xbar_VTP - x_xbar_w1; % from xac_wing1 to xac_VTP
    
    if twin_VTP == 1
        % tail volume coefficient
        CVTP = 2*S_VTP*l_xac_w1VTP/(Geo_tier.cmac_VTP*S_ref);
    else
        % tail volume coefficient
        CVTP = S_VTP*l_xac_w1VTP/(Geo_tier.cmac_VTP*S_ref);
    end
    lfus_b_VTP = l_fus/b_VTP;
    Geo_tier.lfus_b_w2 = lfus_b_VTP;
    
    % Storing DATA
    Geo_tier.l_xac_w1VTP = l_xac_w1VTP;
    Geo_tier.CVTP = CVTP;
end

%% Control surfaces
% Determines the ammount of control surface in the w2 aileron & elevator
% For Elevons
% Control_surface 0 = ELEVON;
% Control_surface 1 = RUDDERVATOR;
% Control_surface 2 = 2-SURFACE;
% Control_surface 3 = 3-SURFACE;
% Control_surface = 1;

if d_ail ==1
    K_y1_ail_w1 = Geo_tier.K_y1_ail_w1;
    K_y2_ail_w1 = Geo_tier.K_y2_ail_w1;
    K_ail = Geo_tier.K_ail;
    % Flag that determines if it is a VTP control surface
    VTP_cs = 0;    

    %% AILERON
    cf_ail = Geo_tier.cf_ail;
    t_c_ail = Geo_tier.t_c_ail; %
    b_ail = b_w1_e*K_ail; % length of aileron's (both surfaces)
    % inner and outter location of the control surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    y_1R_y1_ail = y_offset_w1 + (b_w1_e/2)*K_y1_ail_w1; % inner position from the center line
    y_2R_y1_ail = y_1R_y1_ail; % inner position from the center line
    y_1R_y2_ail = y_offset_w1 + (b_w1_e/2)*K_y2_ail_w1; % outter position from the center line
    y_2R_y2_ail = y_1R_y2_ail; % outter position from the center line
    
    % Calculates geometric information regarding the control surfaces
    CS_geo = Get_ControlSurface_Coordinates(y_1R_y1_ail,y_2R_y1_ail,y_1R_y2_ail,y_2R_y2_ail,...
        x_1R_y1_w1,x_2R_y1_w1,y_offset_w1,Lambda_LE_w1,Lambda_TE_w1,Lambda_c4_w1,dihedral_w1,cf_ail,AC_CONFIGURATION,VTP_cs);
    
    Geo_tier.cmac_ail = CS_geo.cmac_cs;
    Geo_tier.xbar_ail = CS_geo.xbar_cs;
    Geo_tier.ybar_ail = CS_geo.ybar_cs;
    Geo_tier.zbar_ail = CS_geo.zbar_cs;
    
    % Distances relative to the origin
    x_xbar_ail = x_loc_LE_w1 + CS_geo.xbar_cs;
    y_ybar_ail = y_loc_LE_w1 + CS_geo.ybar_cs;
    z_zbar_ail = z_loc_LE_w1 + CS_geo.ybar_cs;
    
    % Store DATA
    Geo_tier.x_xbar_ail = x_xbar_ail;
    Geo_tier.y_ybar_ail = y_ybar_ail;
    Geo_tier.z_zbar_ail = z_zbar_ail;
    
    Geo_tier.x_1R_y1_ail = CS_geo.x_1R_y1_cs;
    Geo_tier.x_2R_y1_ail = CS_geo.x_2R_y1_cs;
    Geo_tier.x_1R_y2_ail = CS_geo.x_1R_y2_cs;
    Geo_tier.x_2R_y2_ail = CS_geo.x_2R_y2_cs;
    Geo_tier.y_1R_y1_ail = CS_geo.y_1R_y1_cs;
    Geo_tier.y_2R_y1_ail = CS_geo.y_2R_y1_cs;
    Geo_tier.y_1R_y2_ail = CS_geo.y_1R_y2_cs;
    Geo_tier.y_2R_y2_ail = CS_geo.y_2R_y2_cs;
    Geo_tier.z_1R_y1_ail = CS_geo.z_1R_y1_cs;
    Geo_tier.z_2R_y1_ail = CS_geo.z_2R_y1_cs;
    Geo_tier.z_1R_y2_ail = CS_geo.z_1R_y2_cs;
    Geo_tier.z_2R_y2_ail = CS_geo.z_2R_y2_cs;
    Geo_tier.cR_ail = CS_geo.c_y1_cs;
    Geo_tier.cT_ail = CS_geo.c_y2_cs;
    Geo_tier.b_ail = CS_geo.b_cs;
    Geo_tier.S_ail = CS_geo.S_cs;
end

if d_ele == 1
%     if AC_type == 1 % Flying wing
% 
%     else       
    K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
    K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;
    K_ele = Geo_tier.K_ele;
    % Flag that determines if it is a VTP control surface
    VTP_cs = 0;

    %% ELEVATOR
    cf_ele = Geo_tier.cf_ele;
    t_c_ele = Geo_tier.t_c_ele; %
    b_ele = b_w1_e*K_ele; % length of aileron's (both surfaces)
    % inner and outter location of the control surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    y_1R_y1_ele = y_offset_w2 + (b_w2_e/2)*K_y1_ele_w2; % inner position from the center line
    y_2R_y1_ele = y_1R_y1_ele; % inner position from the center line
    y_1R_y2_ele = y_offset_w2 + (b_w2_e/2)*K_y2_ele_w2; % outter position from the center line
    y_2R_y2_ele = y_1R_y2_ele; % outter position from the center line
    
    % Calculates geometric information regarding the control surfaces
    CS_geo = Get_ControlSurface_Coordinates(y_1R_y1_ele,y_2R_y1_ele,y_1R_y2_ele,y_2R_y2_ele,...
        x_1R_y1_w2,x_2R_y1_w2,y_offset_w2,Lambda_LE_w2,Lambda_TE_w2,Lambda_c4_w2,dihedral_w2,cf_ele,AC_CONFIGURATION,VTP_cs);
    
    Geo_tier.cmac_ele = CS_geo.cmac_cs;
    Geo_tier.xbar_ele = CS_geo.xbar_cs;
    Geo_tier.ybar_ele = CS_geo.ybar_cs;
    Geo_tier.zbar_ele = CS_geo.zbar_cs;
    
    % Distances relative to the origin
    x_xbar_ele = x_loc_LE_w2 + CS_geo.xbar_cs;
    y_ybar_ele = y_loc_LE_w2 + CS_geo.ybar_cs;
    z_zbar_ele = z_loc_LE_w2 + CS_geo.ybar_cs;
    % Store DATA
    Geo_tier.x_xbar_ele = x_xbar_ele;
    Geo_tier.y_ybar_ele = y_ybar_ele;
    Geo_tier.z_zbar_ele = z_zbar_ele;
    Geo_tier.x_1R_y1_ele = CS_geo.x_1R_y1_cs;
    Geo_tier.x_2R_y1_ele = CS_geo.x_2R_y1_cs;
    Geo_tier.x_1R_y2_ele = CS_geo.x_1R_y2_cs;
    Geo_tier.x_2R_y2_ele = CS_geo.x_2R_y2_cs;
    Geo_tier.y_1R_y1_ele = CS_geo.y_1R_y1_cs;
    Geo_tier.y_2R_y1_ele = CS_geo.y_2R_y1_cs;
    Geo_tier.y_1R_y2_ele = CS_geo.y_1R_y2_cs;
    Geo_tier.y_2R_y2_ele = CS_geo.y_2R_y2_cs;
    Geo_tier.z_1R_y1_ele = CS_geo.z_1R_y1_cs;
    Geo_tier.z_2R_y1_ele = CS_geo.z_2R_y1_cs;
    Geo_tier.z_1R_y2_ele = CS_geo.z_1R_y2_cs;
    Geo_tier.z_2R_y2_ele = CS_geo.z_2R_y2_cs;
    Geo_tier.cR_ele = CS_geo.c_y1_cs;
    Geo_tier.cT_ele = CS_geo.c_y2_cs;
    Geo_tier.b_ele = CS_geo.b_cs;
    Geo_tier.S_ele = CS_geo.S_cs;
    
end
if d_elevon == 1
    K_y1_elevon_w1 = Geo_tier.K_y1_elevon_w1;
    K_y2_elevon_w1 = Geo_tier.K_y2_elevon_w1;
    K_elevon = Geo_tier.K_elevon;
    % Flag that determines if it is a VTP control surface
    VTP_cs = 0;
    
    %% ELEVON
    cf_elevon = Geo_tier.cf_elevon;
    t_c_elevon = Geo_tier.t_c_elevon; %
    b_elevon = b_w1_e*K_elevon; % length of aileron's (both surfaces)
    % inner and outter location of the control surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    y_1R_y1_elevon = y_offset_w1 + (b_w1_e/2)*K_y1_elevon_w1; % inner position from the center line
    y_2R_y1_elevon = y_1R_y1_elevon; % inner position from the center line
    y_1R_y2_elevon = y_offset_w1 + (b_w1_e/2)*K_y2_elevon_w1; % outter position from the center line
    y_2R_y2_elevon = y_1R_y2_elevon; % outter position from the center line
    
    % Calculates geometric information regarding the control surfaces
    CS_geo = Get_ControlSurface_Coordinates(y_1R_y1_elevon,y_2R_y1_elevon,y_1R_y2_elevon,y_2R_y2_elevon,...
        x_1R_y1_w1,x_2R_y1_w1,y_offset_w1,Lambda_LE_w1,Lambda_TE_w1,Lambda_c4_w1,dihedral_w1,cf_elevon,AC_CONFIGURATION,VTP_cs);
    
    Geo_tier.cmac_elevon = CS_geo.cmac_cs;
    Geo_tier.xbar_elevon = CS_geo.xbar_cs;
    Geo_tier.ybar_elevon = CS_geo.ybar_cs;
    Geo_tier.zbar_elevon = CS_geo.zbar_cs;
    
    % Distances relative to the origin
    x_xbar_elevon = x_loc_LE_w1 + CS_geo.xbar_cs;
    y_ybar_elevon = y_loc_LE_w1 + CS_geo.ybar_cs;
    z_zbar_elevon = z_loc_LE_w1 + CS_geo.ybar_cs;
    % Store DATA
    Geo_tier.x_xbar_elevon = x_xbar_elevon;
    Geo_tier.y_ybar_elevon = y_ybar_elevon;
    Geo_tier.z_zbar_elevon = z_zbar_elevon;
    Geo_tier.x_1R_y1_elevon = CS_geo.x_1R_y1_cs;
    Geo_tier.x_2R_y1_elevon = CS_geo.x_2R_y1_cs;
    Geo_tier.x_1R_y2_elevon = CS_geo.x_1R_y2_cs;
    Geo_tier.x_2R_y2_elevon = CS_geo.x_2R_y2_cs;
    Geo_tier.y_1R_y1_elevon = CS_geo.y_1R_y1_cs;
    Geo_tier.y_2R_y1_elevon = CS_geo.y_2R_y1_cs;
    Geo_tier.y_1R_y2_elevon = CS_geo.y_1R_y2_cs;
    Geo_tier.y_2R_y2_elevon = CS_geo.y_2R_y2_cs;
    Geo_tier.z_1R_y1_elevon = CS_geo.z_1R_y1_cs;
    Geo_tier.z_2R_y1_elevon = CS_geo.z_2R_y1_cs;
    Geo_tier.z_1R_y2_elevon = CS_geo.z_1R_y2_cs;
    Geo_tier.z_2R_y2_elevon = CS_geo.z_2R_y2_cs;
    Geo_tier.cR_elevon = CS_geo.c_y1_cs;
    Geo_tier.cT_elevon = CS_geo.c_y2_cs;
    Geo_tier.b_elevon = CS_geo.b_cs;
    Geo_tier.S_elevon = CS_geo.S_cs;
    
end
if d_flap == 1
    K_y1_flap_w1 = Geo_tier.K_y1_flap_w1;
    K_y2_flap_w1 = Geo_tier.K_y2_flap_w1;
    K_flap = Geo_tier.K_flap;
    % Flag that determines if it is a VTP control surface
    VTP_cs = 0;
    
    %% FLAP
    cf_flap = Geo_tier.cf_flap;
    t_c_flap = Geo_tier.t_c_flap;%
    b_flap = b_w1_e*K_flap; % length of flap (both surfaces)
    % inner and outter location of the control surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    y_1R_y1_flap = y_offset_w1 + (b_w1_e/2)*K_y1_flap_w1; % inner position from the center line
    y_2R_y1_flap = y_1R_y1_flap; % inner position from the center line
    y_1R_y2_flap = y_offset_w1 + (b_w1_e/2)*K_y2_flap_w1; % outter position from the center line
    y_2R_y2_flap = y_1R_y2_flap; % outter position from the center line
    
    % Calculates geometric information regarding the control surfaces
    CS_geo = Get_ControlSurface_Coordinates(y_1R_y1_flap,y_2R_y1_flap,y_1R_y2_flap,y_2R_y2_flap,...
        x_1R_y1_w1,x_2R_y1_w1,y_offset_w1,Lambda_LE_w1,Lambda_TE_w1,Lambda_c4_w1,dihedral_w1,cf_flap,AC_CONFIGURATION,VTP_cs);
    
    Geo_tier.cmac_flap = CS_geo.cmac_cs;
    Geo_tier.xbar_flap = CS_geo.xbar_cs;
    Geo_tier.ybar_flap = CS_geo.ybar_cs;
    Geo_tier.zbar_flap = CS_geo.zbar_cs;
    
    % Distances relative to the origin
    x_xbar_flap = x_loc_LE_w1 + CS_geo.xbar_cs;
    y_ybar_flap = y_loc_LE_w1 + CS_geo.ybar_cs;
    z_zbar_flap = z_loc_LE_w1 + CS_geo.ybar_cs;
    
    % Store DATA
    Geo_tier.x_xbar_flap = x_xbar_flap;
    Geo_tier.y_ybar_flap = y_ybar_flap;
    Geo_tier.z_zbar_flap = z_zbar_flap;
    Geo_tier.x_1R_y1_flap = CS_geo.x_1R_y1_cs;
    Geo_tier.x_2R_y1_flap = CS_geo.x_2R_y1_cs;
    Geo_tier.x_1R_y2_flap = CS_geo.x_1R_y2_cs;
    Geo_tier.x_2R_y2_flap = CS_geo.x_2R_y2_cs;
    Geo_tier.y_1R_y1_flap = CS_geo.y_1R_y1_cs;
    Geo_tier.y_2R_y1_flap = CS_geo.y_2R_y1_cs;
    Geo_tier.y_1R_y2_flap = CS_geo.y_1R_y2_cs;
    Geo_tier.y_2R_y2_flap = CS_geo.y_2R_y2_cs;
    Geo_tier.z_1R_y1_flap = CS_geo.z_1R_y1_cs;
    Geo_tier.z_2R_y1_flap = CS_geo.z_2R_y1_cs;
    Geo_tier.z_1R_y2_flap = CS_geo.z_1R_y2_cs;
    Geo_tier.z_2R_y2_flap = CS_geo.z_2R_y2_cs;
    Geo_tier.cR_flap = CS_geo.c_y1_cs;
    Geo_tier.cT_flap = CS_geo.c_y2_cs;
    Geo_tier.b_flap = CS_geo.b_cs;
    Geo_tier.S_flap = CS_geo.S_cs;
end
if d_rudder == 1
    K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
    K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
    K_rudder = Geo_tier.K_rudder;
    % Flag that determines if it is a VTP control surface
    VTP_cs = 1;
    %% RUDDER
    cf_rudder = Geo_tier.cf_rudder;
    t_c_rudder = Geo_tier.t_c_rudder; %
    b_rudder = b_VTP_e*K_rudder; % length of flap (both surfaces)
    % inner and outter location of the control surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    z_1R_y1_rudder = z_offset_VTP + (b_VTP_e)*K_y1_rudder_VTP; % inner position from the center line
    z_2R_y1_rudder = z_1R_y1_rudder; % inner position from the center line
    z_1R_y2_rudder = z_offset_VTP + (b_VTP_e)*K_y2_rudder_VTP; % outter position from the center line
    z_2R_y2_rudder = z_1R_y2_rudder; % outter position from the center line
    
    % Calculates geometric information regarding the control surfaces
    CS_geo = Get_ControlSurface_Coordinates(z_1R_y1_rudder,z_2R_y1_rudder,z_1R_y2_rudder,z_2R_y2_rudder,...
        x_1R_y1_VTP,x_2R_y1_VTP,z_offset_VTP,Lambda_LE_VTP,Lambda_TE_VTP,Lambda_c4_VTP,dihedral_VTP,cf_rudder,AC_CONFIGURATION,VTP_cs);
    
    Geo_tier.cmac_rudder = CS_geo.cmac_cs;
    Geo_tier.xbar_rudder = CS_geo.xbar_cs;
    Geo_tier.ybar_rudder = CS_geo.ybar_cs;
    Geo_tier.zbar_rudder = CS_geo.zbar_cs;
    
    % Distances relative to the origin
    x_xbar_rudder = x_loc_LE_VTP + CS_geo.xbar_cs;
    y_ybar_rudder = y_loc_LE_VTP + CS_geo.ybar_cs;
    z_zbar_rudder = z_loc_LE_VTP + CS_geo.ybar_cs;
    
    % Store DATA
    Geo_tier.x_xbar_rudder = x_xbar_rudder;
    Geo_tier.y_ybar_rudder = y_ybar_rudder;
    Geo_tier.z_zbar_rudder = z_zbar_rudder;
    Geo_tier.x_1R_y1_rudder = CS_geo.x_1R_y1_cs;
    Geo_tier.x_2R_y1_rudder = CS_geo.x_2R_y1_cs;
    Geo_tier.x_1R_y2_rudder = CS_geo.x_1R_y2_cs;
    Geo_tier.x_2R_y2_rudder = CS_geo.x_2R_y2_cs;
    Geo_tier.y_1R_y1_rudder = CS_geo.y_1R_y1_cs;
    Geo_tier.y_2R_y1_rudder = CS_geo.y_2R_y1_cs;
    Geo_tier.y_1R_y2_rudder = CS_geo.y_1R_y2_cs;
    Geo_tier.y_2R_y2_rudder = CS_geo.y_2R_y2_cs;
    Geo_tier.z_1R_y1_rudder = CS_geo.z_1R_y1_cs;
    Geo_tier.z_2R_y1_rudder = CS_geo.z_2R_y1_cs;
    Geo_tier.z_1R_y2_rudder = CS_geo.z_1R_y2_cs;
    Geo_tier.z_2R_y2_rudder = CS_geo.z_2R_y2_cs;
    Geo_tier.cR_rudder = CS_geo.c_y1_cs;
    Geo_tier.cT_rudder = CS_geo.c_y2_cs;
    Geo_tier.b_rudder = CS_geo.b_cs;
    Geo_tier.S_rudder = CS_geo.S_cs;
    Geo_tier.S_rudder = CS_geo.S_cs;
end
if d_rudvtr == 1
    K_y1_rudvtr_w2 = Geo_tier.K_y1_rudvtr_w2;
    K_y2_rudvtr_w2 = Geo_tier.K_y2_rudvtr_w2;
    K_rudvtr = Geo_tier.K_rudvtr;
    % Flag that determines if it is a VTP control surface
    VTP_cs = 0;
    
    %% RUDDERVATOR
    cf_rudvtr = Geo_tier.cf_rudvtr;
    t_c_rudvtr = Geo_tier.t_c_rudvtr; %
    b_rudvtr = b_w2_e*K_rudvtr; % length of flap (both surfaces)
    % inner and outter location of the control surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    y_1R_y1_rudvtr = y_offset_w2 + (b_w2_e/2)*K_y1_rudvtr_w2; % inner position from the center line
    y_2R_y1_rudvtr = y_1R_y1_rudvtr; % inner position from the center line
    y_1R_y2_rudvtr = y_offset_w2 + (b_w2_e/2)*K_y2_rudvtr_w2; % outter position from the center line
    y_2R_y2_rudvtr = y_1R_y2_rudvtr; % outter position from the center line
    
    % Calculates geometric information regarding the control surfaces
    CS_geo = Get_ControlSurface_Coordinates(y_1R_y1_rudvtr,y_2R_y1_rudvtr,y_1R_y2_rudvtr,y_2R_y2_rudvtr,...
        x_1R_y1_w2,x_2R_y1_w2,y_offset_w2,Lambda_LE_w2,Lambda_TE_w2,Lambda_c4_w2,dihedral_w2,cf_rudvtr,AC_CONFIGURATION,VTP_cs);
    
    Geo_tier.cmac_rudvtr = CS_geo.cmac_cs;
    Geo_tier.xbar_rudvtr = CS_geo.xbar_cs;
    Geo_tier.ybar_rudvtr = CS_geo.ybar_cs;
    Geo_tier.zbar_rudvtr = CS_geo.zbar_cs;
    
    % Distances relative to the origin
    x_xbar_rudvtr = x_loc_LE_w2 + CS_geo.xbar_cs;
    y_ybar_rudvtr = y_loc_LE_w2 + CS_geo.ybar_cs;
    z_zbar_rudvtr = z_loc_LE_w2 + CS_geo.ybar_cs;
    
    % Store DATA
    Geo_tier.x_xbar_rudvtr = x_xbar_rudvtr;
    Geo_tier.y_ybar_rudvtr = y_ybar_rudvtr;
    Geo_tier.z_zbar_rudvtr = z_zbar_rudvtr;
    Geo_tier.x_1R_y1_rudvtr = CS_geo.x_1R_y1_cs;
    Geo_tier.x_2R_y1_rudvtr = CS_geo.x_2R_y1_cs;
    Geo_tier.x_1R_y2_rudvtr = CS_geo.x_1R_y2_cs;
    Geo_tier.x_2R_y2_rudvtr = CS_geo.x_2R_y2_cs;
    Geo_tier.y_1R_y1_rudvtr = CS_geo.y_1R_y1_cs;
    Geo_tier.y_2R_y1_rudvtr = CS_geo.y_2R_y1_cs;
    Geo_tier.y_1R_y2_rudvtr = CS_geo.y_1R_y2_cs;
    Geo_tier.y_2R_y2_rudvtr = CS_geo.y_2R_y2_cs;
    Geo_tier.z_1R_y1_rudvtr = CS_geo.z_1R_y1_cs;
    Geo_tier.z_2R_y1_rudvtr = CS_geo.z_2R_y1_cs;
    Geo_tier.z_1R_y2_rudvtr = CS_geo.z_1R_y2_cs;
    Geo_tier.z_2R_y2_rudvtr = CS_geo.z_2R_y2_cs;
    Geo_tier.cR_rudvtr = CS_geo.c_y1_cs;
    Geo_tier.cT_rudvtr = CS_geo.c_y2_cs;
    Geo_tier.b_rudvtr = CS_geo.b_cs;
    Geo_tier.S_rudvtr = CS_geo.S_cs;
    Geo_tier.S_rudvtr = CS_geo.S_cs;
end
if d_can == 1
    K_y1_canard_can = Geo_tier.K_y1_canard_can;
    K_y2_canard_can = Geo_tier.K_y2_canard_can;
    K_can = Geo_tier.K_can;
    % Flag that determines if it is a VTP control surface
    VTP_cs = 0;
    
    %% CANARD
    cf_canard = Geo_tier.cf_canard;
    t_c_canard = Geo_tier.t_c_canard; %
    b_canard = b_can_e*K_can; % length of flap (both surfaces)
    % inner and outter location of the control surface
    % 1R identifies LE
    % 2R identifies TE
    % y1 identifies inner position
    % y2 identifies outer position
    y_1R_y1_canard = y_offset_can + (b_can_e/2)*K_y1_canard_can; % inner position from the center line
    y_2R_y1_canard = y_1R_y1_canard; % inner position from the center line
    y_1R_y2_canard = y_offset_can + (b_can_e/2)*K_y2_canard_can; % outter position from the center line
    y_2R_y2_canard = y_1R_y2_canard; % outter position from the center line
    
    % Calculates geometric information regarding the control surfaces
    CS_geo = Get_ControlSurface_Coordinates(y_1R_y1_canard,y_2R_y1_canard,y_1R_y2_canard,y_2R_y2_canard,...
        x_1R_y1_can,x_2R_y1_can,y_offset_can,Lambda_LE_can,Lambda_TE_can,Lambda_c4_can,dihedral_can,cf_canard,AC_CONFIGURATION,VTP_cs);
    
    Geo_tier.cmac_canard = CS_geo.cmac_cs;
    Geo_tier.xbar_canard = CS_geo.xbar_cs;
    Geo_tier.ybar_canard = CS_geo.ybar_cs;
    Geo_tier.zbar_canard = CS_geo.zbar_cs;
    
    % Distances relative to the origin
    x_xbar_canard = x_loc_LE_can + CS_geo.xbar_cs;
    y_ybar_canard = y_loc_LE_can + CS_geo.ybar_cs;
    z_zbar_canard = z_loc_LE_can + CS_geo.ybar_cs;
    
    % Store DATA
    Geo_tier.x_xbar_canard = x_xbar_canard;
    Geo_tier.y_ybar_canard = y_ybar_canard;
    Geo_tier.z_zbar_canard = z_zbar_canard;
    Geo_tier.x_1R_y1_canard = CS_geo.x_1R_y1_cs;
    Geo_tier.x_2R_y1_canard = CS_geo.x_2R_y1_cs;
    Geo_tier.x_1R_y2_canard = CS_geo.x_1R_y2_cs;
    Geo_tier.x_2R_y2_canard = CS_geo.x_2R_y2_cs;
    Geo_tier.y_1R_y1_canard = CS_geo.y_1R_y1_cs;
    Geo_tier.y_2R_y1_canard = CS_geo.y_2R_y1_cs;
    Geo_tier.y_1R_y2_canard = CS_geo.y_1R_y2_cs;
    Geo_tier.y_2R_y2_canard = CS_geo.y_2R_y2_cs;
    Geo_tier.z_1R_y1_canard = CS_geo.z_1R_y1_cs;
    Geo_tier.z_2R_y1_canard = CS_geo.z_2R_y1_cs;
    Geo_tier.z_1R_y2_canard = CS_geo.z_1R_y2_cs;
    Geo_tier.z_2R_y2_canard = CS_geo.z_2R_y2_cs;
    Geo_tier.cR_canard = CS_geo.c_y1_cs;
    Geo_tier.cT_canard = CS_geo.c_y2_cs;
    Geo_tier.b_canard = CS_geo.b_cs;
    Geo_tier.S_canard = CS_geo.S_cs;
    Geo_tier.S_canard = CS_geo.S_cs;
end

%% Engine and propeller DATA
% Engine location
% Engine_loc = 1 - under wings
% Engine_loc = 2 - fuselage
% Engine_loc = 3 - wingtips
% Engine Configuration
% Engine_conf = 1 - pusher prop
% Engine_conf = 2 - puller prop
% Engine_conf = 3 - turbine
% Geometry of nacelle
l_nc = AC_CONFIGURATION.l_nc; % Length nacelle
d_nc = AC_CONFIGURATION.d_nc; % depth nacelle

switch Engine_loc
    case 1 % Engine_loc = 1 - under wings - n-engines symetrical
        % Assumes engines located at 1/4 of cT_w1 and rotates along the centroid of
        % the nacelle
        x_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
        y_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
        z_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;
        
%         x_eng_xbar = x_loc_LE_w1 + x_1R_y2_w1/2;
%         y_eng_ybar = (b_w1/2 + d_nc/2)/2; % at half wing
%         z_eng_zbar = z_loc_LE_w1 + (y_eng_ybar - y_offset_w1)*tan(dihedral_w1);
        x_eng_xbar = x_eng_ybar1;
        y_eng_ybar = y_eng_ybar1;
        z_eng_zbar = z_eng_ybar1;
        % Storing DATA
        Geo_tier.x_eng_xbar = x_eng_xbar;
        Geo_tier.y_eng_xbar = y_eng_ybar;
        Geo_tier.z_eng_xbar = z_eng_zbar;
        
        % Prop Downwash area
        % Location of engines at the cT_w1/4
        x_prop_cF = x_eng_xbar - l_nc/2; % Location of the prop disk
        y_prop_cF = y_eng_ybar; % Location of the center of the Force prop disk
        z_prop_cF = z_eng_zbar; % Location of the center of the Force prop disk
        y_y1_prop_dw = y_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        y_y2_prop_dw = y_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y1_prop_dw = z_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y2_prop_dw = z_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
        % Storing DATA
        Geo_tier.x_prop_cF = x_prop_cF;
        Geo_tier.y_prop_cF = y_prop_cF;
        Geo_tier.y_y1_prop_dw = y_y1_prop_dw;
        Geo_tier.y_y2_prop_dw = y_y2_prop_dw;  
        Geo_tier.z_y1_prop_dw = z_y1_prop_dw;
        Geo_tier.z_y2_prop_dw = z_y2_prop_dw;  
        
        % Determine the propulsive arms - Positive for engine behind Xcg
        x_d_T = x_eng_xbar - x_XCG; % Positive for engine behind Xcg
        y_d_T = y_eng_ybar - y_XCG; % Positive for engine Ycg
        z_d_T = z_eng_zbar - z_XCG; % Positive for engine Zcg
        
        % Storing DATA
        Geo_tier.x_d_T = x_d_T;
        Geo_tier.y_d_T = y_d_T;
        Geo_tier.z_d_T = z_d_T;
        
    case 2 % Engine_loc = 2 - fuselage - in front
        % the nacelle
        % Engine locations
        x_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
        y_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
        z_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;
        x_eng_xbar = x_eng_ybar1;
        y_eng_ybar = y_eng_ybar1;
        z_eng_zbar = z_eng_ybar1;
        % Storing DATA
        Geo_tier.x_eng_xbar = x_eng_xbar;
        Geo_tier.y_eng_xbar = y_eng_ybar;
        Geo_tier.z_eng_xbar = z_eng_zbar;
        
        % Prop Downwash area
        % Location of engines at the cT_w1/4
        x_prop_cF = x_eng_xbar - l_nc/2; % Location of the prop disk
        y_prop_cF = y_eng_ybar; % Location of the center of the Force prop disk
        z_prop_cF = z_eng_zbar; % Location of the center of the Force prop disk
        y_y1_prop_dw = y_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        y_y2_prop_dw = y_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y1_prop_dw = z_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y2_prop_dw = z_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
        % Storing DATA
        Geo_tier.x_prop_cF = x_prop_cF;
        Geo_tier.y_prop_cF = y_prop_cF;
        Geo_tier.y_y1_prop_dw = y_y1_prop_dw;
        Geo_tier.y_y2_prop_dw = y_y2_prop_dw;  
        Geo_tier.z_y1_prop_dw = z_y1_prop_dw;
        Geo_tier.z_y2_prop_dw = z_y2_prop_dw;  
        
        % Determine the propulsive arms - Positive for engine behind Xcg
        x_d_T = x_eng_xbar - x_XCG; % Positive for engine behind Xcg
        y_d_T = y_eng_ybar - y_XCG; % Positive for engine Ycg
        z_d_T = z_eng_zbar - z_XCG; % Positive for engine Zcg
        % Storing DATA
        Geo_tier.x_d_T = x_d_T;
        Geo_tier.y_d_T = y_d_T;
        Geo_tier.z_d_T = z_d_T;
    case 3 % Engine_loc = 3 - fuselage - rear
        % the nacelle
        % Engine locations
        x_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
        y_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
        z_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;
%         x_eng_xbar = Geo_input_tier.l_fus; % At the rear of the fuselage
%         y_eng_ybar = 0;
%         z_eng_zbar = 0;
        x_eng_xbar = x_eng_ybar1;
        y_eng_ybar = y_eng_ybar1;
        z_eng_zbar = z_eng_ybar1;
        % Storing DATA
        Geo_tier.x_eng_xbar = x_eng_xbar;
        Geo_tier.y_eng_xbar = y_eng_ybar;
        Geo_tier.z_eng_xbar = z_eng_zbar;
        
        % Prop Downwash area
        % Location of engines at the cT_w1/4
        x_prop_cF = x_eng_xbar - l_nc/2; % Location of the prop disk
        y_prop_cF = y_eng_ybar; % Location of the center of the Force prop disk
        z_prop_cF = z_eng_zbar; % Location of the center of the Force prop disk
        y_y1_prop_dw = y_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        y_y2_prop_dw = y_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)

        % Prop in rear fuselage
        y_y1_prop_dw = 0; % inner portion of the prop downwash (no tube correction)
        y_y2_prop_dw = D_prop/2; % inner portion of the prop downwash (no tube correction)

%         if y_eng_ybar == 0
%             y_y1_prop_dw = - D_prop/2; % inner portion of the prop downwash (no tube correction)
%             y_y2_prop_dw = + D_prop/2; % inner portion of the prop downwash (no tube correction)
%         else
%             y_y1_prop_dw = y_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
%             y_y2_prop_dw = y_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
%         end

        z_y1_prop_dw = z_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y2_prop_dw = z_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
        % Storing DATA
        Geo_tier.x_prop_cF = x_prop_cF;
        Geo_tier.y_prop_cF = y_prop_cF;
        Geo_tier.y_y1_prop_dw = y_y1_prop_dw;
        Geo_tier.y_y2_prop_dw = y_y2_prop_dw;  
        Geo_tier.z_y1_prop_dw = z_y1_prop_dw;
        Geo_tier.z_y2_prop_dw = z_y2_prop_dw;  
        
        % Determine the propulsive arms - Positive for engine behind Xcg
        x_d_T = x_eng_xbar - x_XCG; % Positive for engine behind Xcg
        y_d_T = y_eng_ybar - y_XCG; % Positive for engine Ycg
        z_d_T = z_eng_zbar - z_XCG; % Positive for engine Zcg
        % Storing DATA
        Geo_tier.x_d_T = x_d_T;
        Geo_tier.y_d_T = y_d_T;
        Geo_tier.z_d_T = z_d_T;

    case 4 % Engine_loc = 4 - wingtips - n_eng at each side
        % Assumes engines located at 1/4 of cT_w1 and rotates along the centroid of
        % the nacelle
%         x_eng_xbar = x_loc_LE_w1 + x_1R_y2_w1 + cT_w1/4;
%         y_eng_ybar = b_w1/2 + d_nc/2;
%         z_eng_zbar = z_loc_LE_w1 + (y_eng_ybar - y_offset_w1)*tan(dihedral_w1);

        % Engine locations
        x_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
        y_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
        z_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;
        x_eng_xbar = x_eng_ybar1;
        y_eng_ybar = y_eng_ybar1;
        z_eng_zbar = z_eng_ybar1;
        % Storing DATA
        Geo_tier.x_eng_xbar = x_eng_xbar;
        Geo_tier.y_eng_xbar = y_eng_ybar;
        Geo_tier.z_eng_xbar = z_eng_zbar;
        
        % Prop Downwash area
        % Location of engines at the cT_w1/4
        x_prop_cF = x_eng_xbar - l_nc/2; % Location of the prop disk
        y_prop_cF = y_eng_ybar; % Location of the center of the Force prop disk
        z_prop_cF = z_eng_zbar; % Location of the center of the Force prop disk
        y_y1_prop_dw = y_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        y_y2_prop_dw = y_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y1_prop_dw = z_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y2_prop_dw = z_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
        % Storing DATA
        Geo_tier.x_prop_cF = x_prop_cF;
        Geo_tier.y_prop_cF = y_prop_cF;
        Geo_tier.y_y1_prop_dw = y_y1_prop_dw;
        Geo_tier.y_y2_prop_dw = y_y2_prop_dw;  
        Geo_tier.z_y1_prop_dw = z_y1_prop_dw;
        Geo_tier.z_y2_prop_dw = z_y2_prop_dw;  
        
        % Determine the propulsive arms - Positive for engine behind Xcg
        x_d_T = x_eng_xbar - x_XCG; % Positive for engine behind Xcg
        y_d_T = y_eng_ybar - y_XCG; % Positive for engine Ycg
        z_d_T = z_eng_zbar - z_XCG; % Positive for engine Zcg
        % Storing DATA
        Geo_tier.x_d_T = x_d_T;
        Geo_tier.y_d_T = y_d_T;
        Geo_tier.z_d_T = z_d_T;
    case 5 % Engine_loc = 5 - wingtips for wing and canard configuration n_eng at each side
        
        % -------------------- Wing Engines ------------------------------
        % Assumes engines located at 1/4 of cT_w1 and rotates along the centroid of
        % the nacelle
        % Engine locations
        x_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
        y_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
        z_eng_ybar1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;
        x_eng_xbar = x_eng_ybar1;
        y_eng_ybar = y_eng_ybar1;
        z_eng_zbar = z_eng_ybar1;
        x_eng_xbar = x_loc_LE_w1 + x_1R_y2_w1 + cT_w1/4;
        y_eng_ybar = b_w1/2 + d_nc/2;
        z_eng_zbar = z_loc_LE_w1 + (y_eng_ybar - y_offset_w1)*tan(dihedral_w1);
        % Storing DATA
        Geo_tier.x_eng_xbar = x_eng_xbar;
        Geo_tier.y_eng_xbar = y_eng_ybar;
        Geo_tier.z_eng_xbar = z_eng_zbar;

        % Prop Downwash area
        % Location of engines at the cT_w1/4
        x_prop_cF = x_eng_xbar - l_nc/2; % Location of the prop disk
        y_prop_cF = y_eng_ybar; % Location of the center of the Force prop disk
        z_prop_cF = z_eng_zbar; % Location of the center of the Force prop disk
        y_y1_prop_dw = y_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        y_y2_prop_dw = y_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y1_prop_dw = z_prop_cF - D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y2_prop_dw = z_prop_cF + D_prop/2; % inner portion of the prop downwash (no tube correction)
        % Storing DATA
        Geo_tier.x_prop_cF = x_prop_cF;
        Geo_tier.y_prop_cF = y_prop_cF;
        Geo_tier.y_y1_prop_dw = y_y1_prop_dw;
        Geo_tier.y_y2_prop_dw = y_y2_prop_dw;  
        Geo_tier.z_y1_prop_dw = z_y1_prop_dw;
        Geo_tier.z_y2_prop_dw = z_y2_prop_dw;  
        
        % Determine the propulsive arms - Positive for engine behind Xcg
        x_d_T = x_eng_xbar - x_XCG; % Positive for engine behind Xcg
        y_d_T = y_eng_ybar - y_XCG; % Positive for engine Ycg
        z_d_T = z_eng_zbar - z_XCG; % Positive for engine Zcg
        % Storing DATA
        Geo_tier.x_d_T = x_d_T;
        Geo_tier.y_d_T = y_d_T;
        Geo_tier.z_d_T = z_d_T;
        
        % -------------------- Canard Engines ------------------------------
        % Assumes engines located at 1/4 of cT_w1 and rotates along the centroid of
        % the nacelle
        x_eng_ybar2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar2;
        y_eng_ybar2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar2;
        z_eng_ybar2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar2;
        x_eng_xbar_2 = x_eng_ybar2;
        y_eng_ybar_2 = y_eng_ybar2;
        z_eng_zbar_2 = z_eng_ybar2;
        x_eng_xbar_2 = x_loc_LE_w1 + x_1R_y2_w1 + cT_w1/4;
        y_eng_ybar_2 = b_w1/2 + d_nc/2;
        z_eng_zbar_2 = z_loc_LE_w1 + (y_eng_ybar_2 - y_offset_w1)*tan(dihedral_w1);
        % Storing DATA
        Geo_tier.x_eng_xbar2 = x_eng_xbar_2;
        Geo_tier.y_eng_xbar2 = y_eng_ybar_2;
        Geo_tier.z_eng_xbar2 = z_eng_zbar_2;
        
        % Prop Downwash area
        % Location of engines at the cT_w1/4
        x_prop_cF_2 = x_eng_xbar_2 - l_nc/2; % Location of the prop disk
        y_prop_cF_2 = y_eng_ybar_2; % Location of the center of the Force prop disk
        z_prop_cF_2 = z_eng_zbar_2; % Location of the center of the Force prop disk
        y_y1_prop_dw_2 = y_prop_cF_2 - D_prop/2; % inner portion of the prop downwash (no tube correction)
        y_y2_prop_dw_2 = y_prop_cF_2 + D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y1_prop_dw_2 = z_prop_cF_2 - D_prop/2; % inner portion of the prop downwash (no tube correction)
        z_y2_prop_dw_2 = z_prop_cF_2 + D_prop/2; % inner portion of the prop downwash (no tube correction)
        
        % Storing DATA
        Geo_tier.x_prop_cF_2 = x_prop_cF_2;
        Geo_tier.y_prop_cF_2 = y_prop_cF_2;
        Geo_tier.y_y1_prop_dw_2 = y_y1_prop_dw_2;
        Geo_tier.y_y2_prop_dw_2 = y_y2_prop_dw_2;  
        Geo_tier.z_y1_prop_dw_2 = z_y1_prop_dw_2;
        Geo_tier.z_y2_prop_dw_2 = z_y2_prop_dw_2;  
        
        % Determine the propulsive arms - Positive for engine behind Xcg
        x_d_T_2 = x_eng_xbar_2 - x_XCG; % Positive for engine behind Xcg
        y_d_T_2 = y_eng_ybar_2 - y_XCG; % Positive for engine Ycg
        z_d_T_2 = z_eng_zbar_2 - z_XCG; % Positive for engine Zcg
        % Storing DATA
        Geo_tier.x_d_T_2 = x_d_T_2;
        Geo_tier.y_d_T_2 = y_d_T_2;
        Geo_tier.z_d_T_2 = z_d_T_2;
end

%% Dynamic Pressure Corrections
% Checks if there is prop wash affecting w1
        
%% Geometry
% Initial Geometry
% Propulsion
% Defines the Number of engines according to the type of aircraft defined
switch Engine_conf
    case 1 % Engine_conf = 1 - pusher prop
        % Zero propwash affecting aerodynamic surfaces
        S_w1_pw = 0;
        Geo_tier.S_w1_pw = S_w1_pw;
        if HTP == 1 | Vee == 1
            S_w2_pw = 0;
            Geo_tier.S_w2_pw = S_w2_pw;
            
            % % Checks if there is prop wash affecting w2
            if y_y1_prop_dw < y_cT_w2_LE
                
                % Estimation of the corners of the affected w2 surface affected by prop
                % wash
                x_1R_y1_w2_pw = x_1R_y1_w2 + (y_y1_prop_dw - y_offset_w2)*tan(Lambda_LE_w2); % defines inner position of wing affeted by proopwash
                x_2R_y1_w2_pw = x_2R_y1_w2 + (y_y1_prop_dw - y_offset_w2)*tan(Lambda_TE_w2); % defines inner position of wing affeted by proopwash
                x_1R_y2_w2_pw = x_cT_w2_LE;
                x_2R_y2_w2_pw = x_cT_w2_TE;
                % Chord
                c_y1_w2_pw = x_2R_y1_w2_pw - x_1R_y1_w2_pw; % defines chord of wing at beginning of prop wash surface
                c_y2_w2_pw = x_2R_y2_w2_pw - x_1R_y2_w2_pw; % defines chord of wing at beginning of prop wash surface
                % Storing DATA
                Geo_tier.x_1R_y1_w2_pw = x_1R_y1_w2_pw;
                Geo_tier.x_2R_y1_w2_pw = x_2R_y1_w2_pw;
                Geo_tier.x_1R_y2_w2_pw = x_1R_y2_w2_pw;
                Geo_tier.x_2R_y2_w2_pw = x_2R_y2_w2_pw;
                Geo_tier.c_y1_w2_pw = c_y1_w2_pw;
                Geo_tier.c_y2_w2_pw = c_y2_w2_pw;
                % y-location
                y_1R_y1_w2_pw = y_y1_prop_dw;
                y_2R_y1_w2_pw = y_y1_prop_dw;
                y_1R_y2_w2_pw = y_cT_w2_LE;
                y_2R_y2_w2_pw = y_cT_w2_TE;
                % Storing DATA
                Geo_tier.y_1R_y1_w2_pw = y_1R_y1_w2_pw;
                Geo_tier.y_2R_y1_w2_pw = y_2R_y1_w2_pw;
                Geo_tier.y_1R_y2_w2_pw = y_1R_y2_w2_pw;
                Geo_tier.y_2R_y2_w2_pw = y_2R_y2_w2_pw;
                
                % Geometry
                lambda_w2_pw = c_y2_w2_pw/c_y1_w2_pw; % w2 prop wash taper ratio
                b_w2_pw = y_1R_y2_w2_pw - y_1R_y1_w2_pw; % same span
                S_w2_pw = b_w2_pw*((c_y1_w2_pw + c_y2_w2_pw)/2); % total area effective surface
                
                b_w2_pw = D_prop; % same span
                S_w2_pw = b_w2_pw*((c_y1_w2_pw + c_y2_w2_pw)/2); % total area effective surface
                
                % Saturates surface associated to propwash
                if S_w2_pw > S_w2
                    S_w2_pw = S_w2;
                end
                
                AR_w2_pw = b_w2_pw^2/S_w2_pw; % Aspect Ratio
                % Storing DATA
                Geo_tier.lambda_w2_pw = lambda_w2_pw;
                Geo_tier.b_w2_pw = b_w2_pw;
                Geo_tier.S_w2_pw = S_w2_pw;
                Geo_tier.AR_w2_pw = AR_w2_pw;
                
            else
                S_w2_pw = 0;
                Geo_tier.S_w2_pw = S_w2_pw;
            end
            
        end
        if VTP == 1
            S_VTP_pw = 0;
            Geo_tier.S_VTP_pw = S_VTP_pw;
        end
        
    case 2 % Engine_conf = 2 - puller prop
        
        % For wing
        if x_prop_cF < x_cT_w1_LE && y_y1_prop_dw < y_cT_w1_LE
            % Estimation of the corners of the affected w1 surface affected by prop
            % wash
            x_1R_y1_w1_pw = x_1R_y1_w1 + (y_y1_prop_dw - y_offset_w1)*tan(Lambda_LE_w1); % defines inner position of wing affeted by proopwash
            x_2R_y1_w1_pw = x_2R_y1_w1 + (y_y1_prop_dw - y_offset_w1)*tan(Lambda_TE_w1); % defines inner position of wing affeted by proopwash
            x_1R_y2_w1_pw = x_cT_w1_LE;
            x_2R_y2_w1_pw = x_cT_w1_TE;
            % Chord
            c_y1_w1_pw = x_2R_y1_w1_pw - x_1R_y1_w1_pw; % defines chord of wing at beginning of prop wash surface
            c_y2_w1_pw = x_2R_y2_w1_pw - x_1R_y2_w1_pw; % defines chord of wing at beginning of prop wash surface
            % Storing DATA
            Geo_tier.x_1R_y1_w1_pw = x_1R_y1_w1_pw;
            Geo_tier.x_2R_y1_w1_pw = x_2R_y1_w1_pw;
            Geo_tier.x_1R_y2_w1_pw = x_1R_y2_w1_pw;
            Geo_tier.x_2R_y2_w1_pw = x_2R_y2_w1_pw;
            Geo_tier.c_y1_w1_pw = c_y1_w1_pw;
            Geo_tier.c_y2_w1_pw = c_y2_w1_pw;
            % y-location
            y_1R_y1_w1_pw = y_y1_prop_dw;
            y_2R_y1_w1_pw = y_y1_prop_dw;
            y_1R_y2_w1_pw = y_cT_w1_LE;
            y_2R_y2_w1_pw = y_cT_w1_TE;
            % Storing DATA
            Geo_tier.y_1R_y1_w1_pw = y_1R_y1_w1_pw;
            Geo_tier.y_2R_y1_w1_pw = y_2R_y1_w1_pw;
            Geo_tier.y_1R_y2_w1_pw = y_1R_y2_w1_pw;
            Geo_tier.y_2R_y2_w1_pw = y_2R_y2_w1_pw;
            
            % Geometry
            lambda_w1_pw = c_y2_w1_pw/c_y1_w1_pw; % w1 prop wash taper ratio
            b_w1_pw = y_1R_y2_w1_pw - y_1R_y1_w1_pw; % same span
            S_w1_pw = b_w1_pw*((c_y1_w1_pw + c_y2_w1_pw)/2); % total area effective surface
            AR_w1_pw = b_w1_pw^2/S_w1_pw; % Aspect Ratio
            % Storing DATA
            Geo_tier.lambda_w1_pw = lambda_w1_pw;
            Geo_tier.b_w1_pw = b_w1_pw;
            Geo_tier.S_w1_pw = S_w1_pw;
            Geo_tier.AR_w1_pw = AR_w1_pw;
        else
            S_w1_pw = 0;
            Geo_tier.S_w1_pw = S_w1_pw;
        end
        
        if HTP ==1 | Vee == 1
            % % Checks if there is prop wash affecting w2
            if y_y1_prop_dw < y_cT_w2_LE
                
                % Estimation of the corners of the affected w2 surface affected by prop
                % wash
                x_1R_y1_w2_pw = x_1R_y1_w2 + (y_y1_prop_dw - y_offset_w2)*tan(Lambda_LE_w2); % defines inner position of wing affeted by proopwash
                x_2R_y1_w2_pw = x_2R_y1_w2 + (y_y1_prop_dw - y_offset_w2)*tan(Lambda_TE_w2); % defines inner position of wing affeted by proopwash
                x_1R_y2_w2_pw = x_cT_w2_LE;
                x_2R_y2_w2_pw = x_cT_w2_TE;
                % Chord
                c_y1_w2_pw = x_2R_y1_w2_pw - x_1R_y1_w2_pw; % defines chord of wing at beginning of prop wash surface
                c_y2_w2_pw = x_2R_y2_w2_pw - x_1R_y2_w2_pw; % defines chord of wing at beginning of prop wash surface
                % Storing DATA
                Geo_tier.x_1R_y1_w2_pw = x_1R_y1_w2_pw;
                Geo_tier.x_2R_y1_w2_pw = x_2R_y1_w2_pw;
                Geo_tier.x_1R_y2_w2_pw = x_1R_y2_w2_pw;
                Geo_tier.x_2R_y2_w2_pw = x_2R_y2_w2_pw;
                Geo_tier.c_y1_w2_pw = c_y1_w2_pw;
                Geo_tier.c_y2_w2_pw = c_y2_w2_pw;
                % y-location
                y_1R_y1_w2_pw = y_y1_prop_dw;
                y_2R_y1_w2_pw = y_y1_prop_dw;
                y_1R_y2_w2_pw = y_cT_w2_LE;
                y_2R_y2_w2_pw = y_cT_w2_TE;
                % Storing DATA
                Geo_tier.y_1R_y1_w2_pw = y_1R_y1_w2_pw;
                Geo_tier.y_2R_y1_w2_pw = y_2R_y1_w2_pw;
                Geo_tier.y_1R_y2_w2_pw = y_1R_y2_w2_pw;
                Geo_tier.y_2R_y2_w2_pw = y_2R_y2_w2_pw;
                
                % Geometry
                lambda_w2_pw = c_y2_w2_pw/c_y1_w2_pw; % w2 prop wash taper ratio
                b_w2_pw = y_1R_y2_w2_pw - y_1R_y1_w2_pw; % same span
                S_w2_pw = b_w2_pw*((c_y1_w2_pw + c_y2_w2_pw)/2); % total area effective surface
                
                b_w2_pw = D_prop; % same span
                S_w2_pw = b_w2_pw*((c_y1_w2_pw + c_y2_w2_pw)/2); % total area effective surface
                
                % Saturates surface associated to propwash
                if S_w2_pw > S_w2
                    S_w2_pw = S_w2;
                end
                
                AR_w2_pw = b_w2_pw^2/S_w2_pw; % Aspect Ratio
                % Storing DATA
                Geo_tier.lambda_w2_pw = lambda_w2_pw;
                Geo_tier.b_w2_pw = b_w2_pw;
                Geo_tier.S_w2_pw = S_w2_pw;
                Geo_tier.AR_w2_pw = AR_w2_pw;
                
            else
                S_w2_pw = 0;
                Geo_tier.S_w2_pw = S_w2_pw;
            end
        end
        
        if VTP ==1
            % % Checks if there is prop wash affecting w1
            if y_y1_prop_dw < y_cT_VTP_LE
                
                % Estimation of the corners of the affected VTP surface affected by prop
                % wash
                x_1R_y1_VTP_pw = x_1R_y1_VTP + (y_y1_prop_dw - y_offset_VTP)*tan(Lambda_LE_VTP); % defines inner position of wing affeted by proopwash
                x_2R_y1_VTP_pw = x_2R_y1_VTP + (y_y1_prop_dw - y_offset_VTP)*tan(Lambda_TE_VTP); % defines inner position of wing affeted by proopwash
                x_1R_y2_VTP_pw = x_cT_VTP_LE;
                x_2R_y2_VTP_pw = x_cT_VTP_TE;
                % Chord
                c_y1_VTP_pw = x_2R_y1_VTP_pw - x_1R_y1_VTP_pw; % defines chord of wing at beginning of prop wash surface
                c_y2_VTP_pw = x_2R_y2_VTP_pw - x_1R_y2_VTP_pw; % defines chord of wing at beginning of prop wash surface
                % Storing DATA
                Geo_tier.x_1R_y1_VTP_pw = x_1R_y1_VTP_pw;
                Geo_tier.x_2R_y1_VTP_pw = x_2R_y1_VTP_pw;
                Geo_tier.x_1R_y2_VTP_pw = x_1R_y2_VTP_pw;
                Geo_tier.x_2R_y2_VTP_pw = x_2R_y2_VTP_pw;
                Geo_tier.c_y1_VTP_pw = c_y1_VTP_pw;
                Geo_tier.c_y2_VTP_pw = c_y2_VTP_pw;
                % y-location
                y_1R_y1_VTP_pw = y_y1_prop_dw;
                y_2R_y1_VTP_pw = y_y1_prop_dw;
                y_1R_y2_VTP_pw = y_cT_VTP_LE;
                y_2R_y2_VTP_pw = y_cT_VTP_TE;
                % Storing DATA
                Geo_tier.y_1R_y1_VTP_pw = y_1R_y1_VTP_pw;
                Geo_tier.y_2R_y1_VTP_pw = y_2R_y1_VTP_pw;
                Geo_tier.y_1R_y2_VTP_pw = y_1R_y2_VTP_pw;
                Geo_tier.y_2R_y2_VTP_pw = y_2R_y2_VTP_pw;
                
                % Geometry
                lambda_VTP_pw = c_y2_VTP_pw/c_y1_VTP_pw; % VTP prop wash taper ratio
                b_VTP_pw = y_1R_y2_VTP_pw - y_1R_y1_VTP_pw; % same span
                S_VTP_pw = b_VTP_pw*((c_y1_VTP_pw + c_y2_VTP_pw)/2); % total area effective surface

                if twin_VTP == 1
                    % Saturates surface associated to propwash
                    if S_VTP_pw > S_VTP1
                        S_VTP_pw = S_VTP1;
                    end
                else
                    % Saturates surface associated to propwash
                    if S_VTP_pw > S_VTP
                        S_VTP_pw = S_VTP;
                    end
                end

                AR_VTP_pw = b_VTP_pw^2/S_VTP_pw; % Aspect Ratio
                % Storing DATA
                Geo_tier.lambda_VTP_pw = lambda_VTP_pw;
                Geo_tier.b_VTP_pw = b_VTP_pw;
                Geo_tier.S_VTP_pw = S_VTP_pw;
                Geo_tier.AR_VTP_pw = AR_VTP_pw;
                
            else
                S_VTP_pw = 0;
                Geo_tier.S_VTP_pw = S_VTP_pw;
            end
        end
    case 3 % Engine_conf = 3 - turbine
        % Zero propwash affecting aerodynamic surfaces
        %% TEMPORARY
        % Zero propwash affecting aerodynamic surfaces
        if W1 == 1
            S_w1_pw = 0;
            Geo_tier.S_w1_pw = S_w1_pw;
        end
        if HTP == 1 | Vee == 1
            S_w2_pw = 0;
            Geo_tier.S_w2_pw = S_w2_pw;
        end
        if VTP == 1
            S_VTP_pw = 0;
            Geo_tier.S_VTP_pw = S_VTP_pw;
        end
end