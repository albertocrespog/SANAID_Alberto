function Geo_input_tier = Generation_Input_Geometric_Data(Prop_data,conv_UNITS,CASE,AC_CONFIGURATION,SF)
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
D_prop = Prop_data.D_prop;

% Distance from origin in CATIa to Fuselage (offset value)
x_offset_CAD = SF*000/1000; %%No hay offset en x
z_offset_CAD = SF*0/1000; 
% Soting DATA
Geo_input_tier.x_offset_CAD = x_offset_CAD;
Geo_input_tier.z_offset_CAD = z_offset_CAD;

% Geometry
% Fuselage geomtry (Units in m)
w_fus = SF*0.2; % width
h_fus = SF*0.2; % height
l_fus = SF*2.096; % length
d_fus = (w_fus + h_fus)/2; % diameter
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
y_loc_1R_y1_w1_CAD = SF*0/1000; % distance from CAD refference point
y_loc_1R_y2_w1_CAD = SF*1200/1000; % distance from CAD refference point
z_loc_1R_y1_w1_CAD = SF*40/1000;% distance from CAD refference point

% Location of LE w2 
% referenced with rtespect ot the inner point (y1) and outter point (y2)
y_loc_1R_y1_w2_CAD = SF*80/1000; % distance from CAD refference point
y_loc_1R_y2_w2_CAD = SF*(353.5+80)/1000; % distance from CAD refference point
z_loc_1R_y1_w2_CAD = SF*-60/1000;% distance from CAD refference point

switch CASE
    case 1
        % SOLO USAMOS CASO 1 POR AHORA. Usaremos mas casos para las
        % modificaciones previstas (24/6/2019)
        % CASE1 Not making Any changes in the geometry, the airple is really stable in forward flight, trim conditions are really good but:
        % -  XCG is at x_XCG: 1.481 from the nose.
        % -  Aerodynamic center of airplane at 1,513 m from the nose
        % W1
        x_loc_1R_y1_w1_CAD = SF*700/1000; % distance from CAD refference point
        Lambda_LE_w1_e = 26.57*D2R;
        dihedral_w1_e = 0*D2R;
        b_w1 = y_loc_1R_y2_w1_CAD*2;
        b_w1_e = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CAD)*2;
        x_loc_1R_y2_w1_CAD = x_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(Lambda_LE_w1_e); % distance from CAD refference point
        z_loc_1R_y2_w1_CAD = z_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(dihedral_w1_e); % distance from CAD refference point
        
        % W2
        x_loc_1R_y1_w2_CAD = SF*1700/1000; % distance from CAD refference point
        Lambda_LE_w2_e = 24.38*D2R;
        dihedral_w2_e = 45*D2R;
        b_w2 = y_loc_1R_y2_w2_CAD*2;
        b_w2_e = (y_loc_1R_y2_w2_CAD - y_loc_1R_y1_w2_CAD)*2;
        x_loc_1R_y2_w2_CAD = x_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(Lambda_LE_w2_e); % distance from CAD refference point
        z_loc_1R_y2_w2_CAD = z_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(dihedral_w2_e); % distance from CAD refference point
    case 2
        % Case2) Removing sweep in wing
        % -  XCG is at x_XCG: 0.8997 from the nose: 0.1660 m behind the location of the engine thrust line
        % -  Aerodynamic center of airplane at 0.9935 m from the nose
        x_loc_1R_y1_w1_CAD = SF*715.629/1000; % distance from CAD refference point
        Lambda_LE_w1_e = 0*D2R;
        dihedral_w1_e = 4.9155*D2R;
        b_w1 = y_loc_1R_y2_w1_CAD*2;
        b_w1_e = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CAD)*2;
        x_loc_1R_y2_w1_CAD = x_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(Lambda_LE_w1_e); % distance from CAD refference point
        z_loc_1R_y2_w1_CAD = z_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(dihedral_w1_e); % distance from CAD refference point
        
        % W2
        x_loc_1R_y1_w2_CAD = SF*1494/1000; % distance from CAD refference point
        Lambda_LE_w2_e = 17.9289*D2R;
        dihedral_w2_e = 29.3514*D2R;
        b_w2 = y_loc_1R_y2_w2_CAD*2;
        b_w2_e = (y_loc_1R_y2_w2_CAD - y_loc_1R_y1_w2_CAD)*2;
        x_loc_1R_y2_w2_CAD = x_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(Lambda_LE_w2_e); % distance from CAD refference point
        z_loc_1R_y2_w2_CAD = z_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(dihedral_w2_e); % distance from CAD refference point        
    case 3
        % Case3) Putting wing backwards almos 20 cm and removin sweep
        % -  XCG is at x_XCG: 1.0080 from the nose: ONLY 0.0900 m behind the location of the engine thrust line
        % -  Aerodynamic center of airplane at 1.1189 m from the nose
        x_loc_1R_y1_w1_CAD = SF*900/1000; % distance from CAD refference point
        Lambda_LE_w1_e = 0*D2R;
        dihedral_w1_e = 4.9155*D2R;
        b_w1 = y_loc_1R_y2_w1_CAD*2;
        b_w1_e = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CAD)*2;
        x_loc_1R_y2_w1_CAD = x_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(Lambda_LE_w1_e); % distance from CAD refference point
        z_loc_1R_y2_w1_CAD = z_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(dihedral_w1_e); % distance from CAD refference point
        
        % W2
        x_loc_1R_y1_w2_CAD = SF*1494/1000; % distance from CAD refference point
        Lambda_LE_w2_e = 17.9289*D2R;
        dihedral_w2_e = 29.3514*D2R;
        b_w2 = y_loc_1R_y2_w2_CAD*2;
        b_w2_e = (y_loc_1R_y2_w2_CAD - y_loc_1R_y1_w2_CAD)*2;
        x_loc_1R_y2_w2_CAD = x_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(Lambda_LE_w2_e); % distance from CAD refference point
        z_loc_1R_y2_w2_CAD = z_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(dihedral_w2_e); % distance from CAD refference point        
        
    case 4
        % Case4) Putting wing backwards almos 20 cm and removinG sweep from both wing and Vtail
        % -  XCG is at x_XCG: 0.9833 from the nose: ONLY 0.0652 m behind the location of the engine thrust line
        % -  Aerodynamic center of airplane at 1.0942 m from the nose
        x_loc_1R_y1_w1_CAD = SF*900/1000; % distance from CAD refference point
        Lambda_LE_w1_e = 0*D2R;
        dihedral_w1_e = 4.9155*D2R;
        b_w1 = y_loc_1R_y2_w1_CAD*2;
        b_w1_e = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CAD)*2;
        x_loc_1R_y2_w1_CAD = x_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(Lambda_LE_w1_e); % distance from CAD refference point
        z_loc_1R_y2_w1_CAD = z_loc_1R_y1_w1_CAD + (b_w1_e/2)*tan(dihedral_w1_e); % distance from CAD refference point
        
        % W2
        x_loc_1R_y1_w2_CAD = SF*1494/1000; % distance from CAD refference point
        Lambda_LE_w2_e = 0*D2R;
        dihedral_w2_e = 29.3514*D2R;
        b_w2 = y_loc_1R_y2_w2_CAD*2;
        b_w2_e = (y_loc_1R_y2_w2_CAD - y_loc_1R_y1_w2_CAD)*2;
        x_loc_1R_y2_w2_CAD = x_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(Lambda_LE_w2_e); % distance from CAD refference point
        z_loc_1R_y2_w2_CAD = z_loc_1R_y1_w2_CAD + (b_w2_e/2)*tan(dihedral_w2_e); % distance from CAD refference point        
end

% Angles w1
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
Geo_input_tier.dihedral_w1 = dihedral_w1;
Geo_input_tier.dihedral_w1_e = dihedral_w1_e;

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

% Wingspan - w1
% Central Section Section of fuselage with no wing
b_w1_fus = b_w1 - b_w1_e;
% Storing Data
Geo_input_tier.b_w1 = b_w1;
Geo_input_tier.b_w1_e = b_w1_e;
Geo_input_tier.b_w1_fus = b_w1_fus; 

% Angles w2
% Entry Data, Surface LE sweep (in RAD)
% Leading Edge Sweept
Lambda_LE_w2 = Lambda_LE_w2_e;
dihedral_w2 = dihedral_w2_e;

% Angles v
dihedral_v = 0*D2R;

% Storing DATA
Geo_input_tier.Lambda_LE_w2  = Lambda_LE_w2;
Geo_input_tier.Lambda_LE_w2_e  = Lambda_LE_w2_e;
Geo_input_tier.dihedral_w2 = dihedral_w2;
Geo_input_tier.dihedral_w2_e = dihedral_w2_e;
Geo_input_tier.dihedral_v = dihedral_v;

% Relative distances between the inner and outter portion of w1
x_y1_y2_w2_CAD = x_loc_1R_y2_w2_CAD - x_loc_1R_y1_w2_CAD;
y_y1_y2_w2_CAD = y_loc_1R_y2_w2_CAD - y_loc_1R_y1_w2_CAD;
z_y1_y2_w2_CAD = z_loc_1R_y2_w2_CAD - z_loc_1R_y1_w2_CAD;

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

% Wingspan - w2
% Central Section Section of fuselage with no wing
b_w2_fus = b_w2 - b_w2_e; % distance from CAD refference point;
% Storing Data
Geo_input_tier.b_w2 = b_w2;
Geo_input_tier.b_w2_e = b_w2_e;
Geo_input_tier.b_w2_fus = b_w2_fus;

% CAD meassurements
% Chord
cR_w1 = SF*0.2;
cT_w1 = SF*0.2;
cR_w2 = SF*0.32; 	
cT_w2 = SF*0.16;

% Storing Data
Geo_input_tier.cR_w1 = cR_w1;
Geo_input_tier.cT_w1 = cT_w1;
Geo_input_tier.cR_w2 = cR_w2;
Geo_input_tier.cT_w2 = cT_w2;

%% Control surfaces
% Determines the ammount of control surface in the w2 aileron & elevator
% For Elevons
% Control_surface 0 = ELEVON;
% Control_surface 1 = RUDDERVATOR;
% Control_surface 2 = 2-SURFACE;
% Control_surface 3 = 3-SURFACE;
Control_surface = 1; %% Elegimos la de ruddervator
if Control_surface == 0
    K_y1_ail_w1 = 0.5; % Percentage of aileron from effective wing surface
    K_y2_ail_w1 = 1; % Percentage of aileron from effective wing surface
    K_y1_elevon_w1 = 0; % Percentage of elevon from effective wing surface
    K_y2_elevon_w1 = 1; % Percentage of elevon from effective wing surface
    K_y1_flap_w1 = 0.0; % Percentage of aileron from effective wing surface
    K_y2_flap_w1 = 0.0; % Percentage of aileron from effective wing surface
    K_y1_rudvtr_w2 = 0; % complementary control surface
    K_y2_rudvtr_w2 = 1; % complementary control surface
    K_y1_ele_w2 = 0; % complementary control surface
    K_y2_ele_w2 = 0; % complementary control surface
    K_y1_can_w2 = 0; % complementary control surface
    K_y2_can_w2 = 0; % complementary control surface
    
    K_ail = K_y2_ail_w1 - K_y1_ail_w1; % Percentage of aileron
    K_ele = K_y2_ele_w2 - K_y1_ele_w2; % complementary control surface
    K_can = K_y2_can_w2 - K_y1_can_w2; % Percentage of canard
    K_flap = K_y2_flap_w1 - K_y1_flap_w1; % Percentage of flap
    K_rudvtr = K_y2_rudvtr_w2 - K_y1_rudvtr_w2; % Percentage of rudder-vator
    K_elevon = K_y2_elevon_w1 - K_y1_elevon_w1; % Percentage of elevon

elseif Control_surface == 1
    K_y1_ail_w1 = 0.41; % Percentage of aileron from effective wing surface
    K_y2_ail_w1 = 0.98; % Percentage of aileron from effective wing surface
    K_y1_flap_w1 = 0.0; % Percentage of aileron from effective wing surface
    K_y2_flap_w1 = 0.0; % Percentage of aileron from effective wing surface
    K_y1_rudvtr_w2 = 0.02; % complementary control surface
    K_y2_rudvtr_w2 = 0.98; % complementary control surface
    K_y1_ele_w2 = 0; % complementary control surface
    K_y2_ele_w2 = 0; % complementary control surface
    K_y1_can_w2 = 0; % complementary control surface
    K_y2_can_w2 = 0; % complementary control surface
    
    K_ail = K_y2_ail_w1 - K_y1_ail_w1; % Percentage of aileron
    K_ele = K_y2_ele_w2 - K_y1_ele_w2; % complementary control surface
    K_can = K_y2_can_w2 - K_y1_can_w2; % Percentage of canard
    K_flap = K_y2_flap_w1 - K_y1_flap_w1; % Percentage of flap
    K_rudvtr = K_y2_rudvtr_w2 - K_y1_rudvtr_w2; % Percentage of rudder-vator
    
    % Maximum deflections
    delta_ail_min = -25*D2R;
    delta_ail_max = 25*D2R;
    delta_flap_min = -30*D2R;
    delta_flap_max = 30*D2R;
    delta_rudvtr_min = -25*D2R;
    delta_rudvtr_max = 25*D2R;
    
    % Storing DATA
    Geo_input_tier.K_y1_ail_w1 = K_y1_ail_w1;
    Geo_input_tier.K_y2_ail_w1 = K_y2_ail_w1;
    Geo_input_tier.K_y1_ele_w1 = K_y1_ele_w2;
    Geo_input_tier.K_y2_ele_w1 = K_y2_ele_w2;
    Geo_input_tier.K_y1_flap_w1 = K_y1_flap_w1;
    Geo_input_tier.K_y2_flap_w1 = K_y2_flap_w1;
    Geo_input_tier.K_y1_can_w1 = K_y1_can_w2;
    Geo_input_tier.K_y2_can_w1 = K_y2_can_w2;
    Geo_input_tier.K_y1_rudvtr_w1 = K_y1_rudvtr_w2;
    Geo_input_tier.K_y2_rudvtr_w1 = K_y2_rudvtr_w2;
    
    Geo_input_tier.K_ail = K_ail;
    Geo_input_tier.K_ele = K_ele;
    Geo_input_tier.K_can = K_can;
    Geo_input_tier.K_flap = K_flap;
    Geo_input_tier.K_rudvtr = K_rudvtr;
    
    Geo_input_tier.delta_ail_min = delta_ail_min;
    Geo_input_tier.delta_ail_max = delta_ail_max;
    Geo_input_tier.delta_flap_min = delta_flap_min;
    Geo_input_tier.delta_flap_max = delta_flap_max;
    Geo_input_tier.delta_rudvtr_min = delta_rudvtr_min;
    Geo_input_tier.delta_rudvtr_max = delta_rudvtr_max;
    
elseif Control_surface == 2
    K_y1_ail_w1 = 0.62; % Percentage of aileron from effective wing surface
    K_y2_ail_w1 = 0.98; % Percentage of aileron from effective wing surface
    K_y1_flap_w1 = 0.02; % Percentage of aileron from effective wing surface
    K_y2_flap_w1 = 0.58; % Percentage of aileron from effective wing surface
    K_y1_rudder_w2 = 0.02; % complementary control surface
    K_y2_rudder_w2 = 0.98; % complementary control surface
    K_y1_ele_w2 = 0; % complementary control surface
    K_y2_ele_w2 = 0; % complementary control surface
    
    K_ail = K_y2_ail_w1 - K_y1_ail_w1; % Percentage of aileron
    K_ele = K_y2_ele_w2 - K_y1_ele_w2; % complementary control surface
    K_flap = K_y2_flap_w1 - K_y1_flap_w1; % Percentage of flap
    K_rudder = K_y2_rudder_w2 - K_y1_rudder_w2; % Percentage of rudder-vator
    
    % Storing DATA
    Geo_input_tier.K_y1_ail_w1 = K_y1_ail_w1;
    Geo_input_tier.K_y2_ail_w1 = K_y2_ail_w1;
    Geo_input_tier.K_y1_ele_w1 = K_y1_ele_w2;
    Geo_input_tier.K_y2_ele_w1 = K_y2_ele_w2;
    Geo_input_tier.K_y1_flap_w1 = K_y1_flap_w1;
    Geo_input_tier.K_y2_flap_w1 = K_y2_flap_w1;
    Geo_input_tier.K_y1_can_w1 = K_y1_can_w2;
    Geo_input_tier.K_y2_can_w1 = K_y2_can_w2;
    Geo_input_tier.K_y1_rudvtr_w1 = K_y1_rudvtr_w2;
    Geo_input_tier.K_y2_rudvtr_w1 = K_y2_rudvtr_w2;
    
    Geo_input_tier.K_ail = K_ail;
    Geo_input_tier.K_ele = K_ele;
    Geo_input_tier.K_can = K_can;
    Geo_input_tier.K_flap = K_flap;
    Geo_input_tier.K_rudvtr = K_rudvtr;
   
    
elseif Control_surface == 3
    K_y1_ail_w1 = 0.62; % Percentage of aileron from effective wing surface
    K_y2_ail_w1 = 0.98; % Percentage of aileron from effective wing surface
    K_y1_flap_w1 = 0.02; % Percentage of aileron from effective wing surface
    K_y2_flap_w1 = 0.58; % Percentage of aileron from effective wing surface
    K_y1_rudvtr_w2 = 0.02; % complementary control surface
    K_y2_rudvtr_w2 = 0.98; % complementary control surface
    K_y1_ele_w2 = 0; % complementary control surface
    K_y2_ele_w2 = 0; % complementary control surface
    K_y1_can_w2 = 0; % complementary control surface
    K_y2_can_w2 = 0; % complementary control surface
    
    K_ail = K_y2_ail_w1 - K_y1_ail_w1; % Percentage of aileron
    K_ele = K_y2_ele_w2 - K_y1_ele_w2; % complementary control surface
    K_can = K_y2_can_w2 - K_y1_can_w2; % Percentage of canard
    K_flap = K_y2_flap_w1 - K_y1_flap_w1; % Percentage of flap
    K_rudvtr = K_y2_rudvtr_w2 - K_y1_rudvtr_w2; % Percentage of rudder-vator
    
    % Storing DATA
    Geo_input_tier.K_y1_ail_w1 = K_y1_ail_w1;
    Geo_input_tier.K_y2_ail_w1 = K_y2_ail_w1;
    Geo_input_tier.K_y1_ele_w1 = K_y1_ele_w2;
    Geo_input_tier.K_y2_ele_w1 = K_y2_ele_w2;
    Geo_input_tier.K_y1_flap_w1 = K_y1_flap_w1;
    Geo_input_tier.K_y2_flap_w1 = K_y2_flap_w1;
    Geo_input_tier.K_y1_can_w1 = K_y1_can_w2;
    Geo_input_tier.K_y2_can_w1 = K_y2_can_w2;
    Geo_input_tier.K_y1_rudvtr_w1 = K_y1_rudvtr_w2;
    Geo_input_tier.K_y2_rudvtr_w1 = K_y2_rudvtr_w2;
    
    Geo_input_tier.K_ail = K_ail;
    Geo_input_tier.K_ele = K_ele;
    Geo_input_tier.K_can = K_can;
    Geo_input_tier.K_flap = K_flap;
    Geo_input_tier.K_rudvtr = K_rudvtr;
end

%% AILERON
cf_ail = 0.25; % percentage of control surface
t_c_ail = 0.087; % Thinckness 2 chord ratio associated to the airfoil
% Store DATA
Geo_input_tier.cf_ail = cf_ail;
Geo_input_tier.t_c_ail = t_c_ail; %

%% FLAP
cf_flap = 0.25; % percentage of control surface
t_c_flap = 0.087; % Thinckness 2 chord ratio associated to the airfoil
% Store DATA
Geo_input_tier.cf_flap = cf_flap;
Geo_input_tier.t_c_flap = t_c_flap; %

%% RUDDERVATOR
cf_rudvtr = 0.375; % percentage of control surface
t_c_rudvtr = 0.07; % Thinckness 2 chord ratio associated to the airfoil
% Store DATA
Geo_input_tier.cf_rudvtr = cf_rudvtr;
Geo_input_tier.t_c_rudvtr = t_c_rudvtr; %

%% Engine geometry
% Geometry of nacelle
l_nc = Prop_data.l_nc; % Length nacelle
d_nc = Prop_data.d_nc; % depth nacelle

%% Location of the Xcg
% Vertical distance from origin to minimu distance and maximum distance
z_min_Xorigin_CAD = h_fus/2; % height
z_max_Xorigin_CAD = h_fus/2;
% Store DATA
Geo_input_tier.z_min_Xorigin_CAD = z_min_Xorigin_CAD;
Geo_input_tier.z_max_Xorigin_CAD = z_max_Xorigin_CAD;

save Geo_tier.mat