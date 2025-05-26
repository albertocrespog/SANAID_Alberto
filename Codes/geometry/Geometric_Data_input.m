function Geo_tier = Geometric_Data_input
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

% Conversion
D2R = pi/180;
R2D = 180/pi;

% Distance from origin in CATIa to Fuselage (offset value)
x_offset_CAD = -18.829/1000;
z_offset_CAD = +8.425/1000;
Geo_tier.x_offset_CAD = x_offset_CAD;
Geo_tier.z_offset_CAD = z_offset_CAD;

% Geometry
% Fuselage geomtry (Units in m)
w_fus = 0.26; % width
h_fus = 0.25; % height
l_fus = 1.679595; % length
d_fus = (w_fus + h_fus)/2; % diameter
% Storing Data
Geo_tier.w_fus = w_fus;
Geo_tier.h_fus = h_fus;
Geo_tier.l_fus = l_fus;
Geo_tier.d_fus = d_fus;

%% Geometric variations
% Entry Data, Dihedral (in RAD)
dihedral_w1_e = 5*D2R;
dihedral_w2_e = 29*D2R;
dihedral_w1 = dihedral_w1_e;
dihedral_w2 = dihedral_w2_e;
dihedral_v = 0*D2R;
% Storing Data
Geo_tier.dihedral_w1 = dihedral_w1;
Geo_tier.dihedral_w2 = dihedral_w2;
Geo_tier.dihedral_v = dihedral_v;

% CAD meassurements
% Chord
cR_w1 = 0.226321;
cT_w1 = 0.147462; 
cR_w2 = 0.20453; 	
cT_w2 = 0.102565; 
% Storing Data
Geo_tier.cR_w1 = cR_w1;
Geo_tier.cT_w1 = cT_w1;
Geo_tier.cR_w2 = cR_w2;
Geo_tier.cT_w2 = cT_w2;

% Wingspan - w1
b_w1 = 2.509573;
b_w1_fus = 0.250;
b_w1_e = b_w1 - b_w1_fus; % effective 
b_w1_s = (b_w1_e/cos(dihedral_w1_e)); % effective wingspan in y direction
% Wingspan - w2
b_w2 = 1.295646;
b_w2_fus = 0.212;
b_w2_e = b_w2 - b_w2_fus;
b_w2_s = (b_w2_e/cos(dihedral_w2_e)); % effective wingspan in y direction
% Storing Data
Geo_tier.b_w1 = b_w1;
Geo_tier.b_w2 = b_w2;
Geo_tier.b_w1_e = b_w1_e;
Geo_tier.b_w1_s = b_w1_s;
Geo_tier.b_w1_fus = b_w1_fus;
Geo_tier.b_w2_e = b_w2_e;
Geo_tier.b_w2_s = b_w2_s;
Geo_tier.b_w2_fus = b_w2_fus;

% Location of LE w1 
x_loc_LE_w1_CAD = 715.629/1000; % distance from CAD refference point
x_loc_LE_w1 = x_loc_LE_w1_CAD + x_offset_CAD; % corrected to the nose of the aircraft
y_loc_LE_w1 = b_w1_fus/2;
z_loc_LE_w1_CAD = 152.999/1000; % distance from CAD refference point
z_loc_LE_w1 = z_loc_LE_w1_CAD + z_offset_CAD; % corrected to the nose of the aircraft

Geo_tier.x_loc_LE_w1_CAD = x_loc_LE_w1_CAD;
Geo_tier.x_loc_LE_w1 = x_loc_LE_w1;
Geo_tier.y_loc_LE_w1 = y_loc_LE_w1;
Geo_tier.z_loc_LE_w1_CAD = z_loc_LE_w1_CAD;
Geo_tier.z_loc_LE_w1 = z_loc_LE_w1;

% Location of LE w2 
x_loc_LE_w2_CAD = 1494/1000; % distance from CAD refference point
x_loc_LE_w2 = x_loc_LE_w2_CAD + x_offset_CAD;
y_loc_LE_w2 = b_w2_fus/2;
z_loc_LE_w2_CAD = 162.796/1000;
z_loc_LE_w2 = z_loc_LE_w2_CAD + z_offset_CAD;

Geo_tier.x_loc_LE_w2_CAD = x_loc_LE_w2_CAD;
Geo_tier.x_loc_LE_w2 = x_loc_LE_w2;
Geo_tier.y_loc_LE_w2 = y_loc_LE_w2;
Geo_tier.z_loc_LE_w2_CAD = z_loc_LE_w2_CAD;
Geo_tier.z_loc_LE_w2 = z_loc_LE_w2;

% Entry Data, Surface LE sweep (in RAD)
% From trigonmetry
% Leading Edge Sweept
Lambda_LE_w1 = -7*D2R; 
Lambda_LE_w2 = 10.65*D2R;
Lambda_LE_w1_e = Lambda_LE_w1;
Lambda_LE_w2_e = Lambda_LE_w2;
% Storing DATA
Geo_tier.Lambda_LE_w1  = Lambda_LE_w1;
Geo_tier.Lambda_LE_w2  = Lambda_LE_w2;

%% Control surfaces
% Determines the ammount of control surface in the w2 aileron & elevator
% For Elevons
% Control_surface 0 = ELEVON;
% Control_surface 1 = RUDDERVATOR;
% Control_surface 2 = 2-SURFACE;
% Control_surface 3 = 3-SURFACE;
Control_surface = 1;

if Control_surface == 0
%     K_ail_w2 = 1; % Percentage of aileron
%     K_ele_w2 = 1; % complementary control surface
%     K_can_w1 = 0.0; % Percentage of canard
%     K_flap_w1 = 0.0;
%     K_rudvtr_w2 = 0.0;
%     
%     K_ail = K_ail_w2; % Percentage of aileron
%     K_ele = K_ele_w2; % complementary control surface
%     K_can = K_can_w1; % Percentage of canard
%     K_flap = K_flap_w1;
%     K_rudvtr = K_rudvtr_w2;
    
elseif Control_surface == 1
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
    
    
        
elseif Control_surface == 2
%     % For split surface Elevator and Ailero
%     K_ail_w1 = 0.50; % Percentage of aileron
%     K_flap_w1 =  1 - K_ail_w1; % Percentage of flaps
%     K_ele_w2 = 1; % complementary control surface
%     K_can_w1 = 0; % Percentage of canard
%     K_rudvtr_w2 = 0;
% 
%     K_ail = K_ail_w1; % Percentage of aileron
%     K_ele = K_ele_w2; % complementary control surface
%     K_can = K_can_w1; % Percentage of canard
%     K_flap = K_flap_w1;
%     K_rudvtr = K_rudvtr_w2;

elseif Control_surface == 3
%     % For split surface Elevator and Ailero
%     K_can_w1 = 1.0; % Percentage of canard
%     K_ail_w1 = 0.50; % Percentage of aileron
%     K_flap_w1 =  1 - K_ail_w1; % Percentage of flaps
%     K_ele_w2 = 1; % complementary control surface
%     K_rudvtr_w2 = 0;
% 
%     K_ail = K_ail_w1; % Percentage of aileron
%     K_ele = K_ele_w2; % complementary control surface
%     K_can = K_can_w1; % Percentage of canard
%     K_flap = K_flap_w1;
%     K_rudvtr = K_rudvtr_w2;
end

% Storage DATA

Geo_tier.K_y1_ail_w1 = K_y1_ail_w1;
Geo_tier.K_y2_ail_w1 = K_y2_ail_w1;
Geo_tier.K_y1_ele_w1 = K_y1_ele_w2;
Geo_tier.K_y2_ele_w1 = K_y2_ele_w2;
Geo_tier.K_y1_flap_w1 = K_y1_flap_w1;
Geo_tier.K_y2_flap_w1 = K_y2_flap_w1;
Geo_tier.K_y1_can_w1 = K_y1_can_w2;
Geo_tier.K_y2_can_w1 = K_y2_can_w2;
Geo_tier.K_y1_rudvtr_w1 = K_y1_rudvtr_w2;
Geo_tier.K_y2_rudvtr_w1 = K_y2_rudvtr_w2;

Geo_tier.K_ail = K_ail;
Geo_tier.K_ele = K_ele;
Geo_tier.K_can = K_can;
Geo_tier.K_flap = K_flap;
Geo_tier.K_rudvtr = K_rudvtr;

%% AILERON
cf_ail = 0.25; % percentage of control surface
t_c_ail = 0.15; % Thinckness 2 chord ratio associated to the airfoil

%% FLAP
cf_flap = 0.25; % percentage of control surface
t_c_flap = 0.15; % Thinckness 2 chord ratio associated to the airfoil

%% RUDDERVATOR
cf_rudvtr = 0.25; % percentage of control surface
t_c_rudvtr = 0.12; % Thinckness 2 chord ratio associated to the airfoil

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Engine properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propulsive model
% AXI 5360/24HD V2 GOLD LINE
l_eng = 0.104; %length
d_eng = 0.063; %diameter
n_eng=2; % Engine number.

Geo_tier.l_eng = l_eng;
Geo_tier.d_eng = d_eng;
Geo_tier.n_eng = n_eng;

% Approximate dimenstions of the nacelles
l_nc = l_eng*3;
d_nc= d_eng*2; 

%e= Engine number concerning to excel master table. From 1 to 5.
e=6;
load('Data.mat')
load('Data_Prop.mat')

%---------------------------ENGINE DATA----------------------------------%
W_eng = Data.engine(e).Weight;     % Engine Weight.
D_prop=0.6043;        % Propeller Diameter [m].
D_propin=D_prop*100/2.54; % Propeller Diameter [in].
A_prop=pi*(D_prop/2)^2;          % Swept area by propeller.

% Propeller data
% Datos genéricos Hélice para 22x12W - They are used to scale the Prop
b_p = (22*2.54/100);
c_p = 3/100;
b_p_c_p = b_p/c_p;
c_prop = D_prop/b_p_c_p;
S_prop = D_prop*c_prop;
AR_prop = (D_prop^2)/S_prop;
RPM_max = 145000/D_propin; % Max RPM by engine builder.

Geo_tier.W_eng = W_eng;
Geo_tier.D_prop = D_prop;
Geo_tier.D_propin = D_propin;
Geo_tier.A_prop = A_prop;
Geo_tier.RPM_max = RPM_max;
Geo_tier.c_prop = c_prop;
Geo_tier.S_prop = S_prop;
Geo_tier.AR_prop = AR_prop;

Geo_tier.l_eng = l_eng;
Geo_tier.d_eng = d_eng;
Geo_tier.l_nc = l_nc;
Geo_tier.d_nc = d_nc;

% fuel
x_fuel = 1; % 
z_fuel = 0; % 
Geo_tier.x_fuel = x_fuel;
Geo_tier.z_fuel = z_fuel;

% Payload
x_payload = 0.25;
z_payload = 0.25;
Geo_tier.x_payload = x_payload;
Geo_tier.z_payload = z_payload;

m_T0 = 18.8; % Kg
Geo_tier.m_T0 = m_T0;

save Geo_tier.mat