function Geo_input_tier = Generation_Input_Geometric_Data_v3(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX)

SF = OUTPUT_read_XLSX.AC_Data_flags.SF;

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
AC_type = AC_CONFIGURATION.AC_type;
Engine_loc = AC_CONFIGURATION.Engine_loc;
Engine_conf = AC_CONFIGURATION.Engine_conf;

% Conversion
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;

% Propulsion Data
% D_prop = Prop_data.D_prop;

% Distance from origin in CATIa to Fuselage (offset value)
% x_offset_CAD = SF*-18.829/1000;
% z_offset_CAD = SF*+8.425/1000;

x_offset_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.x_offset_CAD; %
z_offset_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_offset_CAD; %

% Soting DATA
Geo_input_tier.x_offset_CAD = x_offset_CAD;
Geo_input_tier.z_offset_CAD = z_offset_CAD;

% Geometry
% Fuselage geomtry (Units in m)
w_fus = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.w_fus;
h_fus = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.h_fus;
l_fus = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.l_fus;
d_fus = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.d_fus;
% w_fus = SF*0.26; % width
% h_fus = SF*0.25; % height
% l_fus = SF*1.679595; % length
% d_fus = (w_fus + h_fus)/2; % diameter
% Storing Data
Geo_input_tier.w_fus = w_fus;
Geo_input_tier.h_fus = h_fus;
Geo_input_tier.l_fus = l_fus;
Geo_input_tier.d_fus = d_fus;

%% Aerodynamic surfaces location
% Location of LE w1 
% 1R identifies LE
% 2R identifies TE
% y1 identifies inner position
% y2 identifies outer position

%% Wing 1
% y_loc_1R_y1_w1_CAD = SF*105.453/1000; % distance from CAD refference point
% y_loc_1R_y2_w1_CAD = SF*1256.6/1000; % distance from CAD refference point
% z_loc_1R_y1_w1_CAD = SF*152.999/1000;% distance from CAD refference point
y_loc_1R_y1_w1_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD; %
y_loc_1R_yB1_w1_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_w1_CAD; %
y_loc_1R_yB2_w1_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_w1_CAD; %
y_loc_1R_y2_w1_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD; %
z_loc_1R_y1_w1_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_w1_CAD; %
x_loc_1R_y1_w1_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD; %
Lambda_LE_w1_e      = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e; %
Lambda_LE_w1_k1_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k1_e; %
Lambda_LE_w1_k2_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k2_e; %
dihedral_w1_e       = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e; %
dihedral_w1_k1_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_k1_e; %
dihedral_w1_k2_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_k2_e; %

b_w1 = y_loc_1R_y2_w1_CAD*2;
b_w1_e = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CAD)*2;
x_loc_1R_y2_w1_CAD = x_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(Lambda_LE_w1_e); % distance from CAD refference point
z_loc_1R_y2_w1_CAD = z_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(dihedral_w1_e); % distance from CAD refference point

%% Wing 2
% referenced with rtespect ot the inner point (y1) and outter point (y2)
% y_loc_1R_y1_w2_CAD = SF*105.453/1000; % distance from CAD refference point
% y_loc_1R_y2_w2_CAD =SF*650.663/1000; % distance from CAD refference point
% z_loc_1R_y1_w2_CAD = SF*162.796/1000;% distance from CAD refference point
y_loc_1R_y1_w2_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD; %
y_loc_1R_yB1_w2_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_w2_CAD; %
y_loc_1R_yB2_w2_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_w2_CAD; %
y_loc_1R_y2_w2_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD; %
z_loc_1R_y1_w2_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_w2_CAD; %
x_loc_1R_y1_w2_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD; %
Lambda_LE_w2_e      = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e; %
Lambda_LE_w2_k1_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_k1_e; %
Lambda_LE_w2_k2_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_k2_e; %
dihedral_w2_e       = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e; %
dihedral_w2_k1_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_k1_e; %
dihedral_w2_k2_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_k2_e; %

b_w2 = y_loc_1R_y2_w2_CAD*2;
b_w2_e = (y_loc_1R_y2_w2_CAD - y_loc_1R_y1_w2_CAD)*2;
x_loc_1R_y2_w2_CAD = x_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(Lambda_LE_w2_e); % distance from CAD refference point
z_loc_1R_y2_w2_CAD = z_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(dihedral_w2_e); % distance from CAD refference point
        
%% Canard
y_loc_1R_y1_can_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD; %
y_loc_1R_yB1_can_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_can_CAD; %
y_loc_1R_yB2_can_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_can_CAD; %
y_loc_1R_y2_can_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_can_CAD; %
z_loc_1R_y1_can_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_can_CAD; %
x_loc_1R_y1_can_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_can_CAD; %
Lambda_LE_can_e      = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_can_e; %
Lambda_LE_can_k1_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_can_k1_e; %
Lambda_LE_can_k2_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_can_k2_e; %
dihedral_can_e       = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_can_e; %
dihedral_can_k1_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_can_k1_e; %
dihedral_can_k2_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_can_k2_e; %

b_can = y_loc_1R_y2_can_CAD*2;
b_can_e = (y_loc_1R_y2_can_CAD - y_loc_1R_y1_can_CAD)*2;
x_loc_1R_y2_can_CAD = x_loc_1R_y1_can_CAD + (b_can_e/2)*tan(Lambda_LE_can_e); % distance from CAD refference point
z_loc_1R_y2_can_CAD = z_loc_1R_y1_can_CAD + (b_can_e/2)*tan(dihedral_can_e); % distance from CAD refference point

%% HTP
y_loc_1R_y1_HTP_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_HTP_CAD; %
y_loc_1R_yB1_HTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_HTP_CAD; %
y_loc_1R_yB2_HTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_HTP_CAD; %
y_loc_1R_y2_HTP_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_HTP_CAD; %
z_loc_1R_y1_HTP_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_HTP_CAD; %
x_loc_1R_y1_HTP_CAD  = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_HTP_CAD; %
Lambda_LE_HTP_e      = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_HTP_e; %
Lambda_LE_HTP_k1_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_HTP_k1_e; %
Lambda_LE_HTP_k2_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_HTP_k2_e; %
dihedral_HTP_e       = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_HTP_e; %
dihedral_HTP_k1_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_HTP_k1_e; %
dihedral_HTP_k2_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_HTP_k2_e; %

b_HTP = y_loc_1R_y2_HTP_CAD*2;
b_HTP_e = (y_loc_1R_y2_HTP_CAD - y_loc_1R_y1_HTP_CAD)*2;
x_loc_1R_y2_HTP_CAD = x_loc_1R_y1_HTP_CAD + (b_HTP_e/2)*tan(Lambda_LE_HTP_e); % distance from CAD refference point
z_loc_1R_y2_HTP_CAD = z_loc_1R_y1_HTP_CAD + (b_HTP_e/2)*tan(dihedral_HTP_e); % distance from CAD refference point

%% VTP
y_loc_1R_y1_VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD; %
z_loc_1R_y1_VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD; %
z_loc_1R_yB1_VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_yB1_VTP_CAD; %
z_loc_1R_yB2_VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_yB2_VTP_CAD; %
z_loc_1R_y2_VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD; %
x_loc_1R_y1_VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD; %
Lambda_LE_VTP_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_VTP_e; %
Lambda_LE_VTP_k1_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_VTP_k1_e; %
Lambda_LE_VTP_k2_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_VTP_k2_e; %
dihedral_VTP_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_VTP_e; %
dihedral_VTP_k1_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_VTP_k1_e; %
dihedral_VTP_k2_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_VTP_k2_e; %

b_VTP = z_loc_1R_y2_VTP_CAD - z_loc_1R_y1_VTP_CAD;
b_VTP_e = z_loc_1R_y2_VTP_CAD - z_loc_1R_y1_VTP_CAD;
x_loc_1R_y2_VTP_CAD = x_loc_1R_y1_VTP_CAD + (b_VTP_e)*tan(Lambda_LE_VTP_e); % distance from CAD refference point
% z_loc_1R_y2_VTP_CAD = z_loc_1R_y1_VTP_CAD + (b_VTP_e)*tan(dihedral_VTP_e); % distance from CAD refference point
y_loc_1R_y2_VTP_CAD = y_loc_1R_y1_VTP_CAD + b_VTP_e*sin(dihedral_VTP_e);

%% Vtail
y_loc_1R_y1_vee_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_vee_CAD; %
y_loc_1R_yB1_vee_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_vee_CAD; %
y_loc_1R_yB2_vee_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_vee_CAD; %
y_loc_1R_y2_vee_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_vee_CAD; %
z_loc_1R_y1_vee_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_vee_CAD; %
x_loc_1R_y1_vee_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_vee_CAD; %
Lambda_LE_vee_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_vee_e; %
Lambda_LE_vee_k1_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_vee_k1_e; %
Lambda_LE_vee_k2_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_vee_k2_e; %
dihedral_vee_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e; %
dihedral_vee_k1_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_k1_e; %
dihedral_vee_k2_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_k2_e; %

b_vee = y_loc_1R_y2_vee_CAD*2;
b_vee_e = (y_loc_1R_y2_vee_CAD - y_loc_1R_y1_vee_CAD)*2;
x_loc_1R_y2_vee_CAD = x_loc_1R_y1_vee_CAD + (b_vee_e/2)*tan(Lambda_LE_vee_e); % distance from CAD refference point
z_loc_1R_y2_vee_CAD = z_loc_1R_y1_vee_CAD + (b_vee_e/2)*tan(dihedral_vee_e); % distance from CAD refference point

%% Vtail2
y_loc_1R_y1_vee2_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_vee2_CAD; %
y_loc_1R_yB1_vee2_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_vee2_CAD; %
y_loc_1R_yB2_vee2_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_vee2_CAD; %
y_loc_1R_y2_vee2_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_vee2_CAD; %
z_loc_1R_y1_vee2_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_vee2_CAD; %
x_loc_1R_y1_vee2_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_vee2_CAD; %
Lambda_LE_vee2_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_vee2_e; %
Lambda_LE_vee2_k1_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_vee2_k1_e; %
Lambda_LE_vee2_k2_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_vee2_k2_e; %
dihedral_vee2_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee2_e; %
dihedral_vee2_k1_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee2_k1_e; %
dihedral_vee2_k2_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee2_k2_e; %

b_vee2 = y_loc_1R_y2_vee2_CAD*2;
b_vee2_e = (y_loc_1R_y2_vee2_CAD - y_loc_1R_y1_vee2_CAD)*2;
x_loc_1R_y2_vee2_CAD = x_loc_1R_y1_vee2_CAD + (b_vee2_e/2)*tan(Lambda_LE_vee2_e); % distance from CAD refference point
z_loc_1R_y2_vee2_CAD = z_loc_1R_y1_vee2_CAD + (b_vee2_e/2)*tan(dihedral_vee2_e); % distance from CAD refference point

%% Twin VTP
y_loc_1R_y1_2VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD; %
z_loc_1R_y1_2VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD; %
z_loc_1R_yB1_2VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_yB1_2VTP_CAD; %
z_loc_1R_yB2_2VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_yB2_2VTP_CAD; %
z_loc_1R_y2_2VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD; %
x_loc_1R_y1_2VTP_CAD = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD; %
Lambda_LE_2VTP_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_2VTP_e; %
Lambda_LE_2VTP_k1_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_2VTP_k1_e; %
Lambda_LE_2VTP_k2_e   = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_2VTP_k2_e; %
dihedral_2VTP_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_2VTP_e; %
dihedral_2VTP_k1_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_2VTP_k1_e; %
dihedral_2VTP_k2_e    = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_2VTP_k2_e; %

b_2VTP = z_loc_1R_y2_2VTP_CAD - z_loc_1R_y1_2VTP_CAD;
b_2VTP_e = z_loc_1R_y2_2VTP_CAD - z_loc_1R_y1_2VTP_CAD;
x_loc_1R_y2_2VTP_CAD = x_loc_1R_y1_2VTP_CAD + (b_2VTP_e)*tan(Lambda_LE_2VTP_e); % distance from CAD refference point
% z_loc_1R_y2_2VTP_CAD = z_loc_1R_y1_2VTP_CAD + (b_2VTP_e)*tan(dihedral_2VTP_e); % distance from CAD refference point
y_loc_1R_y2_2VTP_CAD = y_loc_1R_y1_2VTP_CAD + b_2VTP_e*cos(dihedral_2VTP_e);

%% w1
% Entry Data, Surface LE sweep (in RAD)
% Leading Edge Sweept
Lambda_LE_w1 = Lambda_LE_w1_e;
dihedral_w1 = dihedral_w1_e;
% Relative distances between the inner and outter portion of w1
x_y1_y2_w1_CAD = x_loc_1R_y2_w1_CAD - x_loc_1R_y1_w1_CAD;
y_y1_y2_w1_CAD = y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CAD;
z_y1_y2_w1_CAD = z_loc_1R_y2_w1_CAD - z_loc_1R_y1_w1_CAD;
% Storing DATA
Geo_input_tier.Lambda_LE_w1  = Lambda_LE_w1;
Geo_input_tier.Lambda_LE_w1_e  = Lambda_LE_w1_e;
Geo_input_tier.Lambda_LE_w1_k1_e  = Lambda_LE_w1_k1_e;
Geo_input_tier.Lambda_LE_w1_k2_e  = Lambda_LE_w1_k2_e;
Geo_input_tier.dihedral_w1 = dihedral_w1;
Geo_input_tier.dihedral_w1_e = dihedral_w1_e;
Geo_input_tier.dihedral_w1_k1_e = dihedral_w1_k1_e;
Geo_input_tier.dihedral_w1_k2_e = dihedral_w1_k2_e;


% Storing DATA
Geo_input_tier.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD;
Geo_input_tier.y_loc_1R_y1_w1_CAD = y_loc_1R_y1_w1_CAD;
Geo_input_tier.z_loc_1R_y1_w1_CAD = z_loc_1R_y1_w1_CAD;
Geo_input_tier.x_loc_1R_y2_w1_CAD = x_loc_1R_y2_w1_CAD;
Geo_input_tier.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD;
Geo_input_tier.z_loc_1R_y2_w1_CAD = z_loc_1R_y2_w1_CAD;
Geo_input_tier.x_y1_y2_w1_CAD = x_y1_y2_w1_CAD;
Geo_input_tier.y_y1_y2_w1_CAD = y_y1_y2_w1_CAD;
Geo_input_tier.z_y1_y2_w1_CAD = z_y1_y2_w1_CAD;
Geo_input_tier.y_loc_1R_yB1_w1_CAD = y_loc_1R_yB1_w1_CAD;
Geo_input_tier.y_loc_1R_yB2_w1_CAD = y_loc_1R_yB2_w1_CAD;

% Wingspan - w1
% Central Section Section of fuselage with no wing
b_w1_fus = b_w1 - b_w1_e;
% Storing Data
Geo_input_tier.b_w1 = b_w1;
Geo_input_tier.b_w1_e = b_w1_e;
Geo_input_tier.b_w1_fus = b_w1_fus; 

%% w2
% Entry Data, Surface LE sweep (in RAD)
% Leading Edge Sweept
Lambda_LE_w2 = Lambda_LE_w2_e;
dihedral_w2 = dihedral_w2_e;
% Relative distances between the inner and outter portion of w2
x_y1_y2_w2_CAD = x_loc_1R_y2_w2_CAD - x_loc_1R_y1_w2_CAD;
y_y1_y2_w2_CAD = y_loc_1R_y2_w2_CAD - y_loc_1R_y1_w2_CAD;
z_y1_y2_w2_CAD = z_loc_1R_y2_w2_CAD - z_loc_1R_y1_w2_CAD;
% Storing DATA
Geo_input_tier.Lambda_LE_w2  = Lambda_LE_w2;
Geo_input_tier.Lambda_LE_w2_e  = Lambda_LE_w2_e;
Geo_input_tier.Lambda_LE_w2_k1_e  = Lambda_LE_w2_k1_e;
Geo_input_tier.Lambda_LE_w2_k2_e  = Lambda_LE_w2_k2_e;
Geo_input_tier.dihedral_w2 = dihedral_w2;
Geo_input_tier.dihedral_w2_e = dihedral_w2_e;
Geo_input_tier.dihedral_w2_k1_e = dihedral_w2_k1_e;
Geo_input_tier.dihedral_w2_k2_e = dihedral_w2_k2_e;
% Storing DATA
Geo_input_tier.x_loc_1R_y1_w2_CAD = x_loc_1R_y1_w2_CAD;
Geo_input_tier.y_loc_1R_y1_w2_CAD = y_loc_1R_y1_w2_CAD;
Geo_input_tier.z_loc_1R_y1_w2_CAD = z_loc_1R_y1_w2_CAD;
Geo_input_tier.x_loc_1R_y2_w2_CAD = x_loc_1R_y2_w2_CAD;
Geo_input_tier.y_loc_1R_y2_w2_CAD = y_loc_1R_y2_w2_CAD;
Geo_input_tier.z_loc_1R_y2_w2_CAD = z_loc_1R_y2_w2_CAD;
Geo_input_tier.x_y1_y2_w2_CAD = x_y1_y2_w2_CAD;
Geo_input_tier.y_y1_y2_w2_CAD = y_y1_y2_w2_CAD;
Geo_input_tier.z_y1_y2_w2_CAD = z_y1_y2_w2_CAD;
Geo_input_tier.y_loc_1R_yB1_w2_CAD = y_loc_1R_yB1_w2_CAD;
Geo_input_tier.y_loc_1R_yB2_w2_CAD = y_loc_1R_yB2_w2_CAD;

% Wingspan - w2
% Central Section Section of fuselage with no wing
b_w2_fus = b_w2 - b_w2_e;
% Storing Data
Geo_input_tier.b_w2 = b_w2;
Geo_input_tier.b_w2_e = b_w2_e;
Geo_input_tier.b_w2_fus = b_w2_fus; 

% Angles v
% dihedral_v = 0*D2R;

% % Storing DATA
% Geo_input_tier.Lambda_LE_w2  = Lambda_LE_w2;
% Geo_input_tier.Lambda_LE_w2_e  = Lambda_LE_w2_e;
% Geo_input_tier.dihedral_w2 = dihedral_w2;
% Geo_input_tier.dihedral_w2_e = dihedral_w2_e;
% Geo_input_tier.dihedral_v = dihedral_v;
% 
% % Relative distances between the inner and outter portion of w1
% x_y1_y2_w2_CAD = x_loc_1R_y2_w2_CAD - x_loc_1R_y1_w2_CAD;
% y_y1_y2_w2_CAD = y_loc_1R_y2_w2_CAD - y_loc_1R_y1_w2_CAD;
% z_y1_y2_w2_CAD = z_loc_1R_y2_w2_CAD - z_loc_1R_y1_w2_CAD;
% 
% % Storing DATA
% Geo_input_tier.x_loc_1R_y1_w2_CAD = x_loc_1R_y1_w2_CAD;
% Geo_input_tier.y_loc_1R_y1_w2_CAD = y_loc_1R_y1_w2_CAD;
% Geo_input_tier.z_loc_1R_y1_w2_CAD = z_loc_1R_y1_w2_CAD;
% Geo_input_tier.x_loc_1R_y2_w2_CAD = x_loc_1R_y2_w2_CAD;
% Geo_input_tier.y_loc_1R_y2_w2_CAD = y_loc_1R_y2_w2_CAD;
% Geo_input_tier.z_loc_1R_y2_w2_CAD = z_loc_1R_y2_w2_CAD;
% Geo_input_tier.x_y1_y2_w2_CAD = x_y1_y2_w2_CAD;
% Geo_input_tier.y_y1_y2_w2_CAD = y_y1_y2_w2_CAD;
% Geo_input_tier.z_y1_y2_w2_CAD = z_y1_y2_w2_CAD;
% 
% % Wingspan - w2
% % Central Section Section of fuselage with no wing
% b_w2_fus = b_w2 - b_w2_e; % distance from CAD refference point;
% % Storing Data
% Geo_input_tier.b_w2 = b_w2;
% Geo_input_tier.b_w2_e = b_w2_e;
% Geo_input_tier.b_w2_fus = b_w2_fus;

%% can
% Entry Data, Surface LE sweep (in RAD)
% Leading Edge Sweept
Lambda_LE_can = Lambda_LE_can_e;
dihedral_can = dihedral_can_e;
% Relative distances between the inner and outter portion of can
x_y1_y2_can_CAD = x_loc_1R_y2_can_CAD - x_loc_1R_y1_can_CAD;
y_y1_y2_can_CAD = y_loc_1R_y2_can_CAD - y_loc_1R_y1_can_CAD;
z_y1_y2_can_CAD = z_loc_1R_y2_can_CAD - z_loc_1R_y1_can_CAD;
% Storing DATA
Geo_input_tier.Lambda_LE_can  = Lambda_LE_can;
Geo_input_tier.Lambda_LE_can_e  = Lambda_LE_can_e;
Geo_input_tier.Lambda_LE_can_k1_e  = Lambda_LE_can_k1_e;
Geo_input_tier.Lambda_LE_can_k2_e  = Lambda_LE_can_k2_e;
Geo_input_tier.dihedral_can = dihedral_can;
Geo_input_tier.dihedral_can_e = dihedral_can_e;
Geo_input_tier.dihedral_can_k1_e = dihedral_can_k1_e;
Geo_input_tier.dihedral_can_k2_e = dihedral_can_k2_e;

% Storing DATA
Geo_input_tier.x_loc_1R_y1_can_CAD = x_loc_1R_y1_can_CAD;
Geo_input_tier.y_loc_1R_y1_can_CAD = y_loc_1R_y1_can_CAD;
Geo_input_tier.z_loc_1R_y1_can_CAD = z_loc_1R_y1_can_CAD;
Geo_input_tier.x_loc_1R_y2_can_CAD = x_loc_1R_y2_can_CAD;
Geo_input_tier.y_loc_1R_y2_can_CAD = y_loc_1R_y2_can_CAD;
Geo_input_tier.z_loc_1R_y2_can_CAD = z_loc_1R_y2_can_CAD;
Geo_input_tier.x_y1_y2_can_CAD = x_y1_y2_can_CAD;
Geo_input_tier.y_y1_y2_can_CAD = y_y1_y2_can_CAD;
Geo_input_tier.z_y1_y2_can_CAD = z_y1_y2_can_CAD;
Geo_input_tier.y_loc_1R_yB1_can_CAD = y_loc_1R_yB1_can_CAD;
Geo_input_tier.y_loc_1R_yB2_can_CAD = y_loc_1R_yB2_can_CAD;


% Wingspan - can
% Central Section Section of fuselage with no wing
b_can_fus = b_can - b_can_e;
% Storing Data
Geo_input_tier.b_can = b_can;
Geo_input_tier.b_can_e = b_can_e;
Geo_input_tier.b_can_fus = b_can_fus; 

%% HTP
% Entry Data, Surface LE sweep (in RAD)
% Leading Edge Sweept
Lambda_LE_HTP = Lambda_LE_HTP_e;
dihedral_HTP = dihedral_HTP_e;
% Relative distances between the inner and outter portion of HTP
x_y1_y2_HTP_CAD = x_loc_1R_y2_HTP_CAD - x_loc_1R_y1_HTP_CAD;
y_y1_y2_HTP_CAD = y_loc_1R_y2_HTP_CAD - y_loc_1R_y1_HTP_CAD;
z_y1_y2_HTP_CAD = z_loc_1R_y2_HTP_CAD - z_loc_1R_y1_HTP_CAD;
% Storing DATA
Geo_input_tier.Lambda_LE_HTP  = Lambda_LE_HTP;
Geo_input_tier.Lambda_LE_HTP_e  = Lambda_LE_HTP_e;
Geo_input_tier.Lambda_LE_HTP_k1_e  = Lambda_LE_HTP_k1_e;
Geo_input_tier.Lambda_LE_HTP_k2_e  = Lambda_LE_HTP_k2_e;
Geo_input_tier.dihedral_HTP = dihedral_HTP;
Geo_input_tier.dihedral_HTP_e = dihedral_HTP_e;
Geo_input_tier.dihedral_HTP_k1_e = dihedral_HTP_k1_e;
Geo_input_tier.dihedral_HTP_k2_e = dihedral_HTP_k2_e;

% Storing DATA
Geo_input_tier.x_loc_1R_y1_HTP_CAD = x_loc_1R_y1_HTP_CAD;
Geo_input_tier.y_loc_1R_y1_HTP_CAD = y_loc_1R_y1_HTP_CAD;
Geo_input_tier.z_loc_1R_y1_HTP_CAD = z_loc_1R_y1_HTP_CAD;
Geo_input_tier.x_loc_1R_y2_HTP_CAD = x_loc_1R_y2_HTP_CAD;
Geo_input_tier.y_loc_1R_y2_HTP_CAD = y_loc_1R_y2_HTP_CAD;
Geo_input_tier.z_loc_1R_y2_HTP_CAD = z_loc_1R_y2_HTP_CAD;
Geo_input_tier.x_y1_y2_HTP_CAD = x_y1_y2_HTP_CAD;
Geo_input_tier.y_y1_y2_HTP_CAD = y_y1_y2_HTP_CAD;
Geo_input_tier.z_y1_y2_HTP_CAD = z_y1_y2_HTP_CAD;
Geo_input_tier.y_loc_1R_yB1_HTP_CAD = y_loc_1R_yB1_HTP_CAD;
Geo_input_tier.y_loc_1R_yB2_HTP_CAD = y_loc_1R_yB2_HTP_CAD;

% Wingspan - HTP
% Central Section Section of fuselage with no wing
b_HTP_fus = b_HTP - b_HTP_e;
% Storing Data
Geo_input_tier.b_HTP = b_HTP;
Geo_input_tier.b_HTP_e = b_HTP_e;
Geo_input_tier.b_HTP_fus = b_HTP_fus; 

%% VTP
% Entry Data, Surface LE sweep (in RAD)
% Leading Edge Sweept
Lambda_LE_VTP = Lambda_LE_VTP_e;
dihedral_VTP = dihedral_VTP_e;
% Relative distances between the inner and outter portion of VTP
x_y1_y2_VTP_CAD = x_loc_1R_y2_VTP_CAD - x_loc_1R_y1_VTP_CAD;
y_y1_y2_VTP_CAD = y_loc_1R_y2_VTP_CAD - y_loc_1R_y1_VTP_CAD;
z_y1_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD - z_loc_1R_y1_VTP_CAD;
% Storing DATA
Geo_input_tier.Lambda_LE_VTP  = Lambda_LE_VTP;
Geo_input_tier.Lambda_LE_VTP_e  = Lambda_LE_VTP_e;
Geo_input_tier.Lambda_LE_VTP_k1_e  = Lambda_LE_VTP_k1_e;
Geo_input_tier.Lambda_LE_VTP_k2_e  = Lambda_LE_VTP_k2_e;
Geo_input_tier.dihedral_VTP = dihedral_VTP;
Geo_input_tier.dihedral_VTP_e = dihedral_VTP_e;
Geo_input_tier.dihedral_VTP_k1_e = dihedral_VTP_k1_e;
Geo_input_tier.dihedral_VTP_k2_e = dihedral_VTP_k2_e;

% Storing DATA
Geo_input_tier.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD;
Geo_input_tier.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD;
Geo_input_tier.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD;
Geo_input_tier.x_loc_1R_y2_VTP_CAD = x_loc_1R_y2_VTP_CAD;
Geo_input_tier.y_loc_1R_y2_VTP_CAD = y_loc_1R_y2_VTP_CAD;
Geo_input_tier.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD;
Geo_input_tier.x_y1_y2_VTP_CAD = x_y1_y2_VTP_CAD;
Geo_input_tier.y_y1_y2_VTP_CAD = y_y1_y2_VTP_CAD;
Geo_input_tier.z_y1_y2_VTP_CAD = z_y1_y2_VTP_CAD;
Geo_input_tier.z_loc_1R_yB1_VTP_CAD = z_loc_1R_yB1_VTP_CAD;
Geo_input_tier.z_loc_1R_yB2_VTP_CAD = z_loc_1R_yB2_VTP_CAD;

% Wingspan - VTP
% Central Section Section of fuselage with no wing
b_VTP_fus = b_VTP - b_VTP_e;
% Storing Data
Geo_input_tier.b_VTP = b_VTP;
Geo_input_tier.b_VTP_e = b_VTP_e;
Geo_input_tier.b_VTP_fus = b_VTP_fus; 

%% vee
% Entry Data, Surface LE sweep (in RAD)
% Leading Edge Sweept
Lambda_LE_vee = Lambda_LE_vee_e;
dihedral_vee = dihedral_vee_e;
% Relative distances between the inner and outter portion of vee
x_y1_y2_vee_CAD = x_loc_1R_y2_vee_CAD - x_loc_1R_y1_vee_CAD;
y_y1_y2_vee_CAD = y_loc_1R_y2_vee_CAD - y_loc_1R_y1_vee_CAD;
z_y1_y2_vee_CAD = z_loc_1R_y2_vee_CAD - z_loc_1R_y1_vee_CAD;
% Storing DATA
Geo_input_tier.Lambda_LE_vee  = Lambda_LE_vee;
Geo_input_tier.Lambda_LE_vee_e  = Lambda_LE_vee_e;
Geo_input_tier.Lambda_LE_vee_k1_e  = Lambda_LE_vee_k1_e;
Geo_input_tier.Lambda_LE_vee_k2_e  = Lambda_LE_vee_k2_e;
Geo_input_tier.dihedral_vee = dihedral_vee;
Geo_input_tier.dihedral_vee_e = dihedral_vee_e;
Geo_input_tier.dihedral_vee_k1_e = dihedral_vee_k1_e;
Geo_input_tier.dihedral_vee_k2_e = dihedral_vee_k2_e;

% Storing DATA
Geo_input_tier.x_loc_1R_y1_vee_CAD = x_loc_1R_y1_vee_CAD;
Geo_input_tier.y_loc_1R_y1_vee_CAD = y_loc_1R_y1_vee_CAD;
Geo_input_tier.z_loc_1R_y1_vee_CAD = z_loc_1R_y1_vee_CAD;
Geo_input_tier.x_loc_1R_y2_vee_CAD = x_loc_1R_y2_vee_CAD;
Geo_input_tier.y_loc_1R_y2_vee_CAD = y_loc_1R_y2_vee_CAD;
Geo_input_tier.z_loc_1R_y2_vee_CAD = z_loc_1R_y2_vee_CAD;
Geo_input_tier.x_y1_y2_vee_CAD = x_y1_y2_vee_CAD;
Geo_input_tier.y_y1_y2_vee_CAD = y_y1_y2_vee_CAD;
Geo_input_tier.z_y1_y2_vee_CAD = z_y1_y2_vee_CAD;
Geo_input_tier.y_loc_1R_yB1_vee_CAD = y_loc_1R_yB1_vee_CAD;
Geo_input_tier.y_loc_1R_yB2_vee_CAD = y_loc_1R_yB2_vee_CAD;

% Wingspan - vee
% Central Section Section of fuselage with no wing
b_vee_fus = b_vee - b_vee_e;
% Storing Data
Geo_input_tier.b_vee = b_vee;
Geo_input_tier.b_vee_e = b_vee_e;
Geo_input_tier.b_vee_fus = b_vee_fus; 

%% vee-2
% Entry Data, Surface LE sweep (in RAD)
% Leading Edge Sweept
Lambda_LE_vee2 = Lambda_LE_vee2_e;
dihedral_vee2 = dihedral_vee2_e;
% Relative distances between the inner and outter portion of vee
x_y1_y2_vee2_CAD = x_loc_1R_y2_vee2_CAD - x_loc_1R_y1_vee2_CAD;
y_y1_y2_vee2_CAD = y_loc_1R_y2_vee2_CAD - y_loc_1R_y1_vee2_CAD;
z_y1_y2_vee2_CAD = z_loc_1R_y2_vee2_CAD - z_loc_1R_y1_vee2_CAD;
% Storing DATA
Geo_input_tier.Lambda_LE_vee2  = Lambda_LE_vee2;
Geo_input_tier.Lambda_LE_vee2_e  = Lambda_LE_vee2_e;
Geo_input_tier.Lambda_LE_vee2_k1_e  = Lambda_LE_vee2_k1_e;
Geo_input_tier.Lambda_LE_vee2_k2_e  = Lambda_LE_vee2_k2_e;
Geo_input_tier.dihedral_vee2 = dihedral_vee2;
Geo_input_tier.dihedral_vee2_e = dihedral_vee2_e;
Geo_input_tier.dihedral_vee2_k1_e = dihedral_vee2_k1_e;
Geo_input_tier.dihedral_vee2_k2_e = dihedral_vee2_k2_e;

% Storing DATA
Geo_input_tier.x_loc_1R_y1_vee2_CAD = x_loc_1R_y1_vee2_CAD;
Geo_input_tier.y_loc_1R_y1_vee2_CAD = y_loc_1R_y1_vee2_CAD;
Geo_input_tier.z_loc_1R_y1_vee2_CAD = z_loc_1R_y1_vee2_CAD;
Geo_input_tier.x_loc_1R_y2_vee2_CAD = x_loc_1R_y2_vee2_CAD;
Geo_input_tier.y_loc_1R_y2_vee2_CAD = y_loc_1R_y2_vee2_CAD;
Geo_input_tier.z_loc_1R_y2_vee2_CAD = z_loc_1R_y2_vee2_CAD;
Geo_input_tier.x_y1_y2_vee2_CAD = x_y1_y2_vee2_CAD;
Geo_input_tier.y_y1_y2_vee2_CAD = y_y1_y2_vee2_CAD;
Geo_input_tier.z_y1_y2_vee2_CAD = z_y1_y2_vee2_CAD;
Geo_input_tier.y_loc_1R_yB1_vee2_CAD = y_loc_1R_yB1_vee2_CAD;
Geo_input_tier.y_loc_1R_yB2_vee2_CAD = y_loc_1R_yB2_vee2_CAD;

% Wingspan - vee
% Central Section Section of fuselage with no wing
b_vee2_fus = b_vee2 - b_vee2_e;
% Storing Data
Geo_input_tier.b_vee2 = b_vee2;
Geo_input_tier.b_vee2_e = b_vee2_e;
Geo_input_tier.b_vee2_fus = b_vee2_fus; 

%% 2VTP
% Entry Data, Surface LE sweep (in RAD)
% Leading Edge Sweept
Lambda_LE_2VTP = Lambda_LE_2VTP_e;
dihedral_2VTP = dihedral_2VTP_e;
% Relative distances between the inner and outter portion of 2VTP
x_y1_y2_2VTP_CAD = x_loc_1R_y2_2VTP_CAD - x_loc_1R_y1_2VTP_CAD;
y_y1_y2_2VTP_CAD = y_loc_1R_y2_2VTP_CAD - y_loc_1R_y1_2VTP_CAD;
z_y1_y2_2VTP_CAD = z_loc_1R_y2_2VTP_CAD - z_loc_1R_y1_2VTP_CAD;
% Storing DATA
Geo_input_tier.Lambda_LE_2VTP  = Lambda_LE_2VTP;
Geo_input_tier.Lambda_LE_2VTP_e  = Lambda_LE_2VTP_e;
Geo_input_tier.Lambda_LE_2VTP_k1_e  = Lambda_LE_2VTP_k1_e;
Geo_input_tier.Lambda_LE_2VTP_k2_e  = Lambda_LE_2VTP_k2_e;
Geo_input_tier.dihedral_2VTP = dihedral_2VTP;
Geo_input_tier.dihedral_2VTP_e = dihedral_2VTP_e;
Geo_input_tier.dihedral_2VTP_k1_e = dihedral_2VTP_k1_e;
Geo_input_tier.dihedral_2VTP_k2_e = dihedral_2VTP_k2_e;

% Storing DATA
Geo_input_tier.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_2VTP_CAD;
Geo_input_tier.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_2VTP_CAD;
Geo_input_tier.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_2VTP_CAD;
Geo_input_tier.x_loc_1R_y2_2VTP_CAD = x_loc_1R_y2_2VTP_CAD;
Geo_input_tier.y_loc_1R_y2_2VTP_CAD = y_loc_1R_y2_2VTP_CAD;
Geo_input_tier.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_2VTP_CAD;
Geo_input_tier.x_y1_y2_2VTP_CAD = x_y1_y2_2VTP_CAD;
Geo_input_tier.y_y1_y2_2VTP_CAD = y_y1_y2_2VTP_CAD;
Geo_input_tier.z_y1_y2_2VTP_CAD = z_y1_y2_2VTP_CAD;
Geo_input_tier.z_loc_1R_yB1_2VTP_CAD = z_loc_1R_yB1_2VTP_CAD;
Geo_input_tier.z_loc_1R_yB2_2VTP_CAD = z_loc_1R_yB2_2VTP_CAD;

% Wingspan - 2VTP
% Central Section Section of fuselage with no wing
b_2VTP_fus = b_2VTP - b_2VTP_e;
% Storing Data
Geo_input_tier.b_2VTP = b_2VTP;
Geo_input_tier.b_2VTP_e = b_2VTP_e;
Geo_input_tier.b_2VTP_fus = b_2VTP_fus; 

%% Chord CAD meassurements
%% w1
% cR_w1 = SF*0.226321;
% cT_w1 = SF*0.147462;
cR_w1 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;
cB_k1_w1 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_w1;
cB_k2_w1 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_w1;
cT_w1 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1;
% Storing Data
Geo_input_tier.cR_w1 = cR_w1;
Geo_input_tier.cB_k1_w1 = cB_k1_w1;
Geo_input_tier.cB_k2_w1 = cB_k2_w1;
Geo_input_tier.cT_w1 = cT_w1;

%% w2
% cR_w2 = SF*0.20453; 	
% cT_w2 = SF*0.102565;
cR_w2 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2;
cB_k1_w2 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_w2;
cB_k2_w2 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_w2;
cT_w2 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2;
% Storing Data
Geo_input_tier.cR_w2 = cR_w2;
Geo_input_tier.cB_k1_w2 = cB_k1_w2;
Geo_input_tier.cB_k2_w2 = cB_k2_w2;
Geo_input_tier.cT_w2 = cT_w2;

%% can
cR_can = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_can;
cB_k1_can = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_can;
cB_k2_can = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_can;
cT_can = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_can;
% Storing Data
Geo_input_tier.cR_can = cR_can;
Geo_input_tier.cB_k1_can = cB_k1_can;
Geo_input_tier.cB_k2_can = cB_k2_can;
Geo_input_tier.cT_can = cT_can;

%% HTP
cR_HTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_HTP;
cB_k1_HTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_HTP;
cB_k2_HTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_HTP;
cT_HTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_HTP;
% Storing Data
Geo_input_tier.cR_HTP = cR_HTP;
Geo_input_tier.cB_k1_HTP = cB_k1_HTP;
Geo_input_tier.cB_k2_HTP = cB_k2_HTP;
Geo_input_tier.cT_HTP = cT_HTP;

%% VTP
cR_VTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_VTP;
cB_k1_VTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_VTP;
cB_k2_VTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_VTP;
cT_VTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_VTP;
% Storing Data
Geo_input_tier.cR_VTP = cR_VTP;
Geo_input_tier.cB_k1_VTP = cB_k1_VTP;
Geo_input_tier.cB_k2_VTP = cB_k2_VTP;
Geo_input_tier.cT_VTP = cT_VTP;

%% V tail
cR_vee = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_vee;
cB_k1_vee = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_vee;
cB_k2_vee = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_vee;
cT_vee = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_vee;
% % Storing Data
% Geo_input_tier.cR_VTP = cR_VTP;
% Geo_input_tier.cT_VTP = cT_VTP;
% Storing Data
Geo_input_tier.cR_vee = cR_vee;
Geo_input_tier.cB_k1_vee = cB_k1_vee;
Geo_input_tier.cB_k2_vee = cB_k2_vee;
Geo_input_tier.cT_vee = cT_vee;

%% V tail2
cR_vee2 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_vee2;
cB_k1_vee2 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_vee2;
cB_k2_vee2 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_vee2;
cT_vee2 = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_vee2;
% Storing Data
% Geo_input_tier.cR_VTP = cR_VTP;
% Geo_input_tier.cT_VTP = cT_VTP;
Geo_input_tier.cR_vee2 = cR_vee2;
Geo_input_tier.cB_k1_vee2 = cB_k1_vee2;
Geo_input_tier.cB_k2_vee2 = cB_k2_vee2;
Geo_input_tier.cT_vee2 = cT_vee2;

%% 2VTP
cR_2VTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_2VTP;
cB_k1_2VTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_2VTP;
cB_k2_2VTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_2VTP;
cT_2VTP = SF*OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_2VTP;
% Storing Data
Geo_input_tier.cR_2VTP = cR_2VTP;
Geo_input_tier.cB_k1_2VTP = cB_k1_2VTP;
Geo_input_tier.cB_k2_2VTP = cB_k2_2VTP;
Geo_input_tier.cT_2VTP = cT_2VTP;

%% Control surfaces
% Determines the ammount of control surface in the w2 aileron & elevator
% For Elevons
% Control_surface 0 = ELEVON;
% Control_surface 1 = RUDDERVATOR;
% Control_surface 2 = 2-SURFACE;
% Control_surface 3 = 3-SURFACE;
%% Control surfaces
% Determines the ammount of control surface in the w2 aileron & elevator
% For Elevons
% Control_surface 0 = ELEVON;
% Control_surface 1 = RUDDERVATOR;
% Control_surface 2 = 2-SURFACE;
% Control_surface 3 = 3-SURFACE;
% if Control_surface == 0
%     K_y1_ail_w1 = 0.62; % Percentage of aileron from effective wing surface
%     K_y2_ail_w1 = 0.98; % Percentage of aileron from effective wing surface
%     K_y1_ele_w2 = 0; % complementary control surface
%     K_y2_ele_w2 = 0; % complementary control surface
%     K_y1_elevon_w1 = 0.62; % Percentage of elevon from effective wing surface
%     K_y2_elevon_w1 = 0.98; % Percentage of elevon from effective wing surface
%     K_y1_flap_w1 = 0.02; % Percentage of aileron from effective wing surface
%     K_y2_flap_w1 = 0.58; % Percentage of aileron from effective wing surface
%     K_y1_rudder_VTP = 0; % Percentage of rudder from effective VTP surface
%     K_y2_rudder_VTP = 0; % Percentage of rudder from effective VTP surface
%     K_y1_rudvtr_w2 = 0.02; % complementary control surface
%     K_y2_rudvtr_w2 = 0.98; % complementary control surface
%     K_y1_can_w2 = 0; % complementary control surface
%     K_y2_can_w2 = 0; % complementary control surface
%     
%     % Maximum deflections
%     delta_ail_min = -25*D2R;
%     delta_ail_max = 25*D2R;
%     delta_ele_min = -0*D2R;
%     delta_ele_max = 0*D2R;
%     delta_elevon_min = -25*D2R;
%     delta_elevon_max = 25*D2R;
%     delta_flap_min = -30*D2R;
%     delta_flap_max = 30*D2R;
%     delta_rudder_min = -0*D2R;
%     delta_rudder_max = 0*D2R;
%     delta_rudvtr_min = -25*D2R;
%     delta_rudvtr_max = 25*D2R;
%     delta_can_min = -0*D2R;
%     delta_can_max = 0*D2R;
% 
% elseif Control_surface == 1
%     K_y1_ail_w1 = 0.62; % Percentage of aileron from effective wing span
%     K_y2_ail_w1 = 0.98; % Percentage of aileron from effective wing span
%     K_y1_ele_w2 = 0; % complementary control surface
%     K_y2_ele_w2 = 0; % complementary control surface
%     K_y1_elevon_w1 = 0.0; % Percentage of elevon from effective wing span
%     K_y2_elevon_w1 = 0.0; % Percentage of elevon from effective wing span
%     K_y1_flap_w1 = 0.02; % Percentage of aileron from effective wing span
%     K_y2_flap_w1 = 0.58; % Percentage of aileron from effective wing span
%     K_y1_rudder_VTP = 0; % Percentage of rudder from effective VTP span
%     K_y2_rudder_VTP = 0; % Percentage of rudder from effective VTP span
%     K_y1_rudvtr_w2 = 0.02; % complementary control span
%     K_y2_rudvtr_w2 = 0.98; % complementary control span
%     K_y1_can_w2 = 0; % complementary control span
%     K_y2_can_w2 = 0; % complementary control span
%     
%     % Maximum deflections
%     delta_ail_min = -25*D2R;
%     delta_ail_max = 25*D2R;
%     delta_ele_min = -0*D2R;
%     delta_ele_max = 0*D2R;
%     delta_elevon_min = -0*D2R;
%     delta_elevon_max = 0*D2R;
%     delta_flap_min = -30*D2R;
%     delta_flap_max = 30*D2R;
%     delta_rudder_min = -0*D2R;
%     delta_rudder_max = 0*D2R;
%     delta_rudvtr_min = -25*D2R;
%     delta_rudvtr_max = 25*D2R;
%     delta_can_min = -0*D2R;
%     delta_can_max = 0*D2R;
%     
% elseif Control_surface == 2
%     K_y1_ail_w1 = 0.62; % Percentage of aileron from effective wing surface
%     K_y2_ail_w1 = 0.98; % Percentage of aileron from effective wing surface
%     K_y1_ele_w2 = 0; % complementary control surface
%     K_y2_ele_w2 = 0.98; % complementary control surface
%     K_y1_elevon_w1 = 0.0; % Percentage of elevon from effective wing surface
%     K_y2_elevon_w1 = 0.0; % Percentage of elevon from effective wing surface
%     K_y1_flap_w1 = 0.02; % Percentage of aileron from effective wing surface
%     K_y2_flap_w1 = 0.58; % Percentage of aileron from effective wing surface
%     K_y1_rudder_VTP = 0.02; % complementary control surface
%     K_y2_rudder_VTP = 0.98; % complementary control surface
%     K_y1_rudvtr_w2 = 0;
%     K_y2_rudvtr_w2 = 0;
%     K_y1_can_w2 = 0; 
%     K_y2_can_w2 = 0; 
%     
%     % Maximum deflections
%     delta_ail_min = -25*D2R;
%     delta_ail_max = 25*D2R;
%     delta_ele_min = -25*D2R;
%     delta_ele_max = 25*D2R;
%     delta_elevon_min = -0*D2R;
%     delta_elevon_max = 0*D2R;
%     delta_flap_min = -30*D2R;
%     delta_flap_max = 30*D2R;
%     delta_rudder_min = -25*D2R;
%     delta_rudder_max = 25*D2R;
%     delta_rudvtr_min = -0*D2R;
%     delta_rudvtr_max = 0*D2R;
%     delta_can_min = -0*D2R;
%     delta_can_max = 0*D2R;
%  
% elseif Control_surface == 3
%     K_y1_ail_w1 = 0.62; % Percentage of aileron from effective wing surface
%     K_y2_ail_w1 = 0.98; % Percentage of aileron from effective wing surface
%     K_y1_ele_w2 = 0; % complementary control surface
%     K_y2_ele_w2 = 0.98; % complementary control surface
%     K_y1_elevon_w1 = 0.0; % Percentage of elevon from effective wing surface
%     K_y2_elevon_w1 = 0.0; % Percentage of elevon from effective wing surface
%     K_y1_flap_w1 = 0.02; % Percentage of aileron from effective wing surface
%     K_y2_flap_w1 = 0.58; % Percentage of aileron from effective wing surface
%     K_y1_rudder_VTP = 0.02; % complementary control surface
%     K_y2_rudder_VTP = 0.98; % complementary control surface
%     K_y1_rudvtr_w2 = 0.0; % complementary control surface
%     K_y2_rudvtr_w2 = 0.0; % complementary control surface
%     K_y1_can_w2 = 0; % complementary control surface
%     K_y2_can_w2 = 0.98; % complementary control surface
%  
%     % Maximum deflections
%     delta_ail_min = -25*D2R;
%     delta_ail_max = 25*D2R;
%     delta_ele_min = -25*D2R;
%     delta_ele_max = 25*D2R;
%     delta_elevon_min = -0*D2R;
%     delta_elevon_max = 0*D2R;
%     delta_flap_min = -30*D2R;
%     delta_flap_max = 30*D2R;
%     delta_rudder_min = -25*D2R;
%     delta_rudder_max = 25*D2R;
%     delta_rudvtr_min = -0*D2R;
%     delta_rudvtr_max = 0*D2R;
%     delta_can_min = -25*D2R;
%     delta_can_max = 25*D2R;
%  
% end

K_y1_ail_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_ail_w1;
K_y2_ail_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_ail_w1;
K_y1_ele_HTP = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_ele_HTP;
K_y2_ele_HTP = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_ele_HTP;
K_y1_elevon_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_elevon_w1;
K_y2_elevon_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_elevon_w1;
K_y1_flap_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_flap_w1;
K_y2_flap_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_flap_w1;
K_y1_rudder_VTP = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_rudder_VTP;
K_y2_rudder_VTP = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_rudder_VTP;
K_y1_rudvtr_vee = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_rudvtr_vee;
K_y2_rudvtr_vee = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_rudvtr_vee;
K_y1_rudvtr_vee2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_rudvtr_vee2;
K_y2_rudvtr_vee2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_rudvtr_vee2;
K_y1_canard_can = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_canard_can;
K_y2_canard_can = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_canard_can;
delta_ail_min = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_ail_min;
delta_ail_max = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_ail_max;
delta_ele_min = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_ele_min;
delta_ele_max = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_ele_max;
delta_elevon_min = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_elevon_min;
delta_elevon_max = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_elevon_max;
delta_flap_min = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_flap_min;
delta_flap_max = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_flap_max;
delta_rudder_min = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_rudder_min;
delta_rudder_max = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_rudder_max;
delta_rudvtr_min = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_rudvtr_min;
delta_rudvtr_max = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_rudvtr_max;
delta_can_min = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_can_min;
delta_can_max = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_can_max;

% Storing DATA
Geo_input_tier.K_y1_ail_w1 = K_y1_ail_w1;
Geo_input_tier.K_y2_ail_w1 = K_y2_ail_w1;
Geo_input_tier.K_y1_ele_HTP = K_y1_ele_HTP;
Geo_input_tier.K_y2_ele_HTP = K_y2_ele_HTP;
Geo_input_tier.K_y1_elevon_w1 = K_y1_elevon_w1;
Geo_input_tier.K_y2_elevon_w1 = K_y2_elevon_w1;
Geo_input_tier.K_y1_flap_w1 = K_y1_flap_w1;
Geo_input_tier.K_y2_flap_w1 = K_y2_flap_w1;
Geo_input_tier.K_y1_rudder_VTP = K_y1_rudder_VTP;
Geo_input_tier.K_y2_rudder_VTP = K_y2_rudder_VTP;
Geo_input_tier.K_y1_rudvtr_vee = K_y1_rudvtr_vee;
Geo_input_tier.K_y2_rudvtr_vee = K_y2_rudvtr_vee;
Geo_input_tier.K_y1_rudvtr_vee2 = K_y1_rudvtr_vee2;
Geo_input_tier.K_y2_rudvtr_vee2 = K_y2_rudvtr_vee2;
Geo_input_tier.K_y1_canard_can = K_y1_canard_can;
Geo_input_tier.K_y2_canard_can = K_y2_canard_can;

K_ail = K_y2_ail_w1 - K_y1_ail_w1; % Percentage of aileron
K_ele = K_y2_ele_HTP - K_y1_ele_HTP; % complementary control surface
K_elevon = K_y2_elevon_w1 - K_y1_elevon_w1; % Percentage of elevon
K_flap = K_y2_flap_w1 - K_y1_flap_w1; % Percentage of flap
K_rudvtr = K_y2_rudvtr_vee - K_y1_rudvtr_vee; % Percentage of rudder-vator
K_rudvtr2 = K_y2_rudvtr_vee2 - K_y1_rudvtr_vee2; % Percentage of rudder-vator
K_rudder = K_y2_rudder_VTP - K_y1_rudder_VTP; % Percentage of rudder
K_can = K_y2_canard_can - K_y1_canard_can; % Percentage of canard

Geo_input_tier.K_ail = K_ail;
Geo_input_tier.K_ele = K_ele;
Geo_input_tier.K_elevon = K_elevon;
Geo_input_tier.K_flap = K_flap;
Geo_input_tier.K_rudder = K_rudder;
Geo_input_tier.K_rudvtr = K_rudvtr;
Geo_input_tier.K_rudvtr2 = K_rudvtr2;
Geo_input_tier.K_can = K_can;

Geo_input_tier.delta_ail_min = delta_ail_min;
Geo_input_tier.delta_ail_max = delta_ail_max;
Geo_input_tier.delta_ele_min = delta_ele_min;
Geo_input_tier.delta_ele_max = delta_ele_max;
Geo_input_tier.delta_elevon_min = delta_elevon_min;
Geo_input_tier.delta_elevon_max = delta_elevon_max;
Geo_input_tier.delta_flap_min = delta_flap_min;
Geo_input_tier.delta_flap_max = delta_flap_max;
Geo_input_tier.delta_rudder_min = delta_rudder_min;
Geo_input_tier.delta_rudder_max = delta_rudder_max;
Geo_input_tier.delta_rudvtr_min = delta_rudvtr_min;
Geo_input_tier.delta_rudvtr_max = delta_rudvtr_max;
Geo_input_tier.delta_can_min = delta_can_min;
Geo_input_tier.delta_can_max = delta_can_max;

%% AILERON
% cf_ail = 0.25; % percentage of control surface
% t_c_ail = 0.15; % Thinckness 2 chord ratio associated to the airfoil
cf_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_ail;
t_c_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_ail;
% Store DATA
Geo_input_tier.cf_ail = cf_ail;
Geo_input_tier.t_c_ail = t_c_ail; %

%% ELEVATOR
% cf_ele = 0.25; % percentage of control surface
% t_c_ele = 0.15; % Thinckness 2 chord ratio associated to the airfoil
cf_ele = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_ele;
t_c_ele = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_ele;
% Store DATA
Geo_input_tier.cf_ele = cf_ele;
Geo_input_tier.t_c_ele = t_c_ele; %

%% ELEVON
% cf_elevon = 0.25; % percentage of control surface
% t_c_elevon = 0.15; % Thinckness 2 chord ratio associated to the airfoil
cf_elevon = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_elevon;
t_c_elevon = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_elevon;
% Store DATA
Geo_input_tier.cf_elevon = cf_elevon;
Geo_input_tier.t_c_elevon = t_c_elevon; %

%% FLAP
% cf_flap = 0.25; % percentage of control surface
% t_c_flap = 0.15; % Thinckness 2 chord ratio associated to the airfoil
cf_flap = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_flap;
t_c_flap = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_flap;
% Store DATA
Geo_input_tier.cf_flap = cf_flap;
Geo_input_tier.t_c_flap = t_c_flap; %

%% RUDDER
% cf_rudder = 0.45; % percentage of control surface
% t_c_rudder = 0.12; % Thinckness 2 chord ratio associated to the airfoil
cf_rudder = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_rudder;
t_c_rudder = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_rudder;
% Store DATA
Geo_input_tier.cf_rudder = cf_rudder;
Geo_input_tier.t_c_rudder = t_c_rudder; %

%% RUDDERVATOR
% cf_rudvtr = 0.45; % percentage of control surface
% t_c_rudvtr = 0.12; % Thinckness 2 chord ratio associated to the airfoil
cf_rudvtr = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_rudvtr;
t_c_rudvtr = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_rudvtr;
% Store DATA
Geo_input_tier.cf_rudvtr = cf_rudvtr;
Geo_input_tier.t_c_rudvtr = t_c_rudvtr; %

%% RUDDERVATOR
% cf_rudvtr = 0.45; % percentage of control surface
% t_c_rudvtr = 0.12; % Thinckness 2 chord ratio associated to the airfoil
cf_rudvtr2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_rudvtr2;
t_c_rudvtr2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_rudvtr2;
% Store DATA
Geo_input_tier.cf_rudvtr2 = cf_rudvtr;
Geo_input_tier.t_c_rudvtr2 = t_c_rudvtr2; %

%% CANNARD
% cf_can = 0.45; % percentage of control surface
% t_c_can = 0.12; % Thinckness 2 chord ratio associated to the airfoil
cf_canard = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_canard;
t_c_canard = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_canard;
% Store DATA
Geo_input_tier.cf_canard = cf_canard;
Geo_input_tier.t_c_canard = t_c_canard; %

% %% Engine geometry
% % Geometry of nacelle
% l_nc = Prop_data.l_nc; % Length nacelle
% d_nc = Prop_data.d_nc; % depth nacelle

%% Location of the Xcg
% Center of gravity
x_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
y_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG;
z_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;

% Vertical distance from origin to minimu distance and maximum distance
% z_min_Xorigin_CAD = 81.179/1000;
% z_max_Xorigin_CAD = 186.215/1000;
% Store DATA
% Geo_input_tier.z_min_Xorigin_CAD = z_min_Xorigin_CAD;
% Geo_input_tier.z_max_Xorigin_CAD = z_max_Xorigin_CAD;

% save Geo_tier.mat