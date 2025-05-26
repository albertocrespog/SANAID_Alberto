function Geo_tier = Generate_Geometric_Data(Geo_tier)
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

% % Distance from origin in CATIa to Fuselage (offset value)
% x_offset_CAD = -18.829/1000;
% z_offset_CAD = +8.425/1000;

x_offset_CAD = Geo_tier.x_offset_CAD;
z_offset_CAD = Geo_tier.z_offset_CAD;

% Geometry
% Fuselage geomtry (Units in m)
% w_fus = 0.26; % width
% h_fus = 0.25; % height
% l_fus = 1.679595; % length
% d_fus = (w_fus + h_fus)/2; % diameter
% Storing Data
w_fus = Geo_tier.w_fus;
h_fus = Geo_tier.h_fus;
l_fus = Geo_tier.l_fus;
d_fus = Geo_tier.d_fus;

%% Geometric variations
% GEO_VAR = 1 then ProVant4 (Emergentia)
% GEO_VAR = 2 then Pepiño XXL
% GEO_VAR = 1;

% Entry Data, Dihedral (in RAD)
% dihedral_w1_e = 5*D2R;
% dihedral_w2_e = 29*D2R;
% dihedral_w1 = dihedral_w1_e;
% dihedral_w2 = dihedral_w2_e;
% dihedral_v = 0*D2R;
% Storing Data
dihedral_w1 = Geo_tier.dihedral_w1;
dihedral_w2 = Geo_tier.dihedral_w2;
dihedral_v = Geo_tier.dihedral_v;

% CAD meassurements
% Chord
% cR_w1 = 0.226321;
% cT_w1 = 0.147462; 
% cR_w2 = 0.20453; 	
% cT_w2 = 0.102565; 
% Storing Data
cR_w1 = Geo_tier.cR_w1;
cT_w1 = Geo_tier.cT_w1;
cR_w2 = Geo_tier.cR_w2;
cT_w2= Geo_tier.cT_w2;

% Wingspan - w1
% b_w1 = 2.509573;
% b_w1_fus = 0.250;
% b_w1_e = b_w1 - b_w1_fus; % effective 
% b_w1_s = (b_w1_e/cos(dihedral_w1_e)); % effective wingspan in y direction
% % Wingspan - w2
% b_w2 = 1.295646;
% b_w2_fus = 0.212;
% b_w2_e = b_w2 - b_w2_fus;
% b_w2_s = (b_w2_e/cos(dihedral_w2_e)); % effective wingspan in y direction
% Storing Data
b_w1 = Geo_tier.b_w1;
b_w2 = Geo_tier.b_w2;
b_w1_e = Geo_tier.b_w1_e;
b_w1_s = Geo_tier.b_w1_s;
b_w1_fus = Geo_tier.b_w1_s;
b_w2_e = Geo_tier.b_w2_e;
b_w2_s = Geo_tier.b_w2_s;
b_w2_fus = Geo_tier.b_w1_s;

% Wing area 1
S_w1_e = b_w1_e*(cR_w1+cT_w1)/2; % effective wing area: proyected in teh y-plane
S_w1_s = b_w1_s*(cR_w1+cT_w1)/2; % real wing area: proyected along the surface of wing
S_w1_fus = cR_w1*b_w1_fus; % w1 area within the fuselage
S_w1 = S_w1_e + S_w1_fus; % w1 Reference area
% Wing area 2
S_w2_e = b_w2_e*(cR_w2+cT_w2)/2; % effective wing area: proyected in teh y-plane
S_w2_s = b_w2_s*(cR_w2+cT_w2); % real wing area: proyected along the surface of wing
S_w2_fus = cR_w2*b_w2_fus; % w2 area within the fuselage
S_w2 = S_w2_e + S_w2_fus; % w2 Reference area
% Storing DATA
Geo_tier.S_w1 = S_w1;
Geo_tier.S_w1_e  = S_w1_e; 
Geo_tier.S_w1_s  = S_w1_s; 
Geo_tier.S_w1_fus  = S_w1_fus; 

Geo_tier.S_w2 = S_w2;
Geo_tier.S_w2_e  = S_w2_e; 
Geo_tier.S_w2_s  = S_w2_s; 
Geo_tier.S_w2_fus  = S_w2_fus; 

%%%%% IMPORTANT %%%%
% The user select the reference wing area that will be used 
% Reference area
S_ref = S_w1_e;
Geo_tier.S_ref = S_ref;

% Location of LE w1 
% x_loc_LE_w1_CAD = 715.629/1000; % distance from CAD refference point
% x_loc_LE_w1 = x_loc_LE_w1_CAD + x_offset_CAD; % corrected to the nose of the aircraft
% y_loc_LE_w1 = b_w1_fus/2;
% z_loc_LE_w1_CAD = 152.999/1000; % distance from CAD refference point
% z_loc_LE_w1 = z_loc_LE_w1_CAD + z_offset_CAD; % corrected to the nose of the aircraft

x_loc_LE_w1_CAD = Geo_tier.x_loc_LE_w1_CAD;
x_loc_LE_w1 = Geo_tier.x_loc_LE_w1;
y_loc_LE_w1 = Geo_tier.y_loc_LE_w1;
z_loc_LE_w1_CAD = Geo_tier.z_loc_LE_w1_CAD;
z_loc_LE_w1 = Geo_tier.z_loc_LE_w1;

% Location of LE w2 
% x_loc_LE_w2_CAD = 1494/1000; % distance from CAD refference point
% x_loc_LE_w2 = x_loc_LE_w2_CAD + x_offset_CAD;
% y_loc_LE_w2 = b_w2_fus/2;
% z_loc_LE_w2_CAD = 162.796/1000;
% z_loc_LE_w2 = z_loc_LE_w2_CAD + z_offset_CAD;

x_loc_LE_w2_CAD = Geo_tier.x_loc_LE_w2_CAD;
x_loc_LE_w2 = Geo_tier.x_loc_LE_w2;
y_loc_LE_w2 = Geo_tier.y_loc_LE_w2;
z_loc_LE_w2_CAD = Geo_tier.z_loc_LE_w2_CAD;
z_loc_LE_w2 = Geo_tier.z_loc_LE_w2;

x_w1_LE = x_loc_LE_w1;
x_w2_LE = x_loc_LE_w2;
y_w1_LE = y_loc_LE_w1;
y_w2_LE = y_loc_LE_w2;
z_w1_LE = z_loc_LE_w1;
z_w2_LE = z_loc_LE_w2;
% Storing DATA
Geo_tier.x_w1_LE = x_w1_LE; 	
Geo_tier.x_w2_LE = x_w2_LE;
Geo_tier.y_w1_LE = y_w1_LE; 	
Geo_tier.y_w2_LE = y_w2_LE;
Geo_tier.z_w1_LE = z_w1_LE; 	
Geo_tier.z_w2_LE = z_w2_LE;

% Aspect Ratio
AR_w1 = b_w1^2/S_w1;
AR_w1_e = b_w1_e^2/S_w1_e;
AR_w2 = b_w2^2/S_w2;
AR_w2_e = b_w2_e^2/S_w2_e;
% Storing DATA
Geo_tier.AR_w1 = AR_w1;
Geo_tier.AR_w1_e = AR_w1_e;
Geo_tier.AR_w2 = AR_w2;
Geo_tier.AR_w2_e = AR_w2_e;

% Taper Ratio
lambda_w1 = cT_w1/cR_w1;
lambda_w1_e = cT_w1/cR_w1;
lambda_w2 = cT_w2/cR_w2;
lambda_w2_e = cT_w2/cR_w2;
% Storing DATA
Geo_tier.lambda_w1 = lambda_w1;
Geo_tier.lambda_w1_e = lambda_w1_e;
Geo_tier.lambda_w2 = lambda_w2;
Geo_tier.lambda_w2_e = lambda_w2_e;

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

% Quarter Chord
Lambda_c4_w1 = Get_Nth_Lambda(1/4,AR_w1_e,Lambda_LE_w1,lambda_w1_e);
Lambda_c4_w2 = Get_Nth_Lambda(1/4,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
% Half chord
Lambda_c2_w1 = Get_Nth_Lambda(1/2,AR_w1_e,Lambda_LE_w1,lambda_w1_e);
Lambda_c2_w2 = Get_Nth_Lambda(1/2,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
% Trailing Edge
Lambda_TE_w1 = Get_Nth_Lambda(1,AR_w1_e,Lambda_LE_w1,lambda_w1_e);
Lambda_TE_w2 = Get_Nth_Lambda(1,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
% Storing DATA
Geo_tier.Lambda_TE_w1  = Lambda_TE_w1;
Geo_tier.Lambda_TE_w2  = Lambda_TE_w2;
Geo_tier.Lambda_TE_w1  = Lambda_TE_w1;
Geo_tier.Lambda_TE_w2  = Lambda_TE_w2;
Geo_tier.Lambda_c4_w1  = Lambda_c4_w1;
Geo_tier.Lambda_c4_w2  = Lambda_c4_w2;
Geo_tier.Lambda_c2_w1  = Lambda_c2_w1;
Geo_tier.Lambda_c2_w2  = Lambda_c2_w2;

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

% Geometric position of w1 MAC with respect to wing 1R
XYZ_MAC = Get_MAC_Coordinates(b_w1_e,lambda_w1_e,cR_w1,dihedral_w1,Lambda_c4_w1);
Geo_tier.cmac_w1 = XYZ_MAC.cbar;
Geo_tier.xbar_w1 = XYZ_MAC.xbar_w;
Geo_tier.ybar_w1 = XYZ_MAC.ybar_w;
Geo_tier.zbar_w1 = XYZ_MAC.zbar_w;

% Distances relative to the origin
x_xbar_w1 = x_loc_LE_w1 + Geo_tier.xbar_w1;
y_ybar_w1 = y_loc_LE_w1 + Geo_tier.ybar_w1;
z_zbar_w1 = z_loc_LE_w1 + Geo_tier.zbar_w1;

% Geometric position of w2 MAC
XYZ_MAC = Get_MAC_Coordinates(b_w2_e,lambda_w2_e,cR_w2,dihedral_w2,Lambda_c4_w2);
Geo_tier.cmac_w2 = XYZ_MAC.cbar;
Geo_tier.xbar_w2 = XYZ_MAC.xbar_w;
Geo_tier.ybar_w2 = XYZ_MAC.ybar_w;
Geo_tier.zbar_w2 = XYZ_MAC.zbar_w;

% Distances relative to the origin
x_xbar_w2 = x_loc_LE_w2 + Geo_tier.xbar_w2;
y_ybar_w2 = y_loc_LE_w2 + Geo_tier.ybar_w2;
z_zbar_w2 = z_loc_LE_w2 + Geo_tier.zbar_w2;

% Storing DATA
Geo_tier.x_xbar_w1 = x_xbar_w1;
Geo_tier.y_ybar_w1 = y_ybar_w1;
Geo_tier.z_zbar_w1 = z_zbar_w1;
Geo_tier.x_xbar_w2 = x_xbar_w2;
Geo_tier.y_ybar_w2 = y_ybar_w2;
Geo_tier.z_zbar_w2 = z_zbar_w2;

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
t_c_ail = 0.15; % Thinckness 2 chord ratio associated to the airfoil

% Calculates geometric information regarding the control surfaces
CS_geo = Get_ControlSurface_Coordinates(y_1R_y1_ail,y_2R_y1_ail,y_1R_y2_ail,y_2R_y2_ail,...
    x_1R_y1_w1,x_2R_y1_w1,y_offset_w1,Lambda_LE_w1,Lambda_TE_w1,Lambda_c4_w1,dihedral_w1,cf_ail);

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

%% FLAP
cf_flap = 0.25; % percentage of control surface
b_flap = b_w1_e*K_flap; % length of flap (both surfaces)
% inner and outter location of the control surface
% 1R identifies LE
% 2R identifies TE
% y1 identifies inner position
% y2 identifies outer position
y_1R_y1_flap = y_offset_w1 + (b_w1_e/2)*K_y1_flap_w1; % inner position from the center line
y_2R_y1_flap = y_1R_y1_ail; % inner position from the center line
y_1R_y2_flap = y_offset_w1 + (b_w1_e/2)*K_y2_flap_w1; % outter position from the center line
y_2R_y2_flap = y_1R_y2_flap; % outter position from the center line
t_c_flap = 0.15; % Thinckness 2 chord ratio associated to the airfoil

% Calculates geometric information regarding the control surfaces
CS_geo = Get_ControlSurface_Coordinates(y_1R_y1_ail,y_2R_y1_ail,y_1R_y2_ail,y_2R_y2_ail,...
    x_1R_y1_w1,x_2R_y1_w1,y_offset_w1,Lambda_LE_w1,Lambda_TE_w1,Lambda_c4_w1,dihedral_w1,cf_ail);

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


%% RUDDERVATOR
cf_rudvtr = 0.25; % percentage of control surface
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
t_c_rudvtr = 0.12; % Thinckness 2 chord ratio associated to the airfoil

% Calculates geometric information regarding the control surfaces
CS_geo = Get_ControlSurface_Coordinates(y_1R_y1_rudvtr,y_2R_y1_rudvtr,y_1R_y2_rudvtr,y_2R_y2_rudvtr,...
    x_1R_y1_w2,x_2R_y1_w2,y_offset_w2,Lambda_LE_w2,Lambda_TE_w2,Lambda_c4_w2,dihedral_w2,cf_rudvtr);

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

% Tail Volume Coefficients
l_xac_w1w2 = x_xbar_w2 - x_xbar_w1; % from xac_wing1 to xac_wing2
% fuselage length
lfus_b_w1 = l_fus/b_w1;
lfus_b_w2 = l_fus/b_w2;
% tail volume coefficient
Cw2 = S_w2*l_xac_w1w2/(Geo_tier.cmac_w2*S_w1);
% Storing DATA
Geo_tier.l_xac_w1w2 = l_xac_w1w2;
Geo_tier.lfus_b_w1 = lfus_b_w1;
Geo_tier.lfus_b_w2 = lfus_b_w2;
Geo_tier.Cw2 = Cw2;

% Assumes engines located at 1/4 of cT_w1
x_eng_xbar = x_loc_LE_w1 + x_1R_y2_w1 + cT_w1/4;
y_eng_ybar = b_w1/2;
z_eng_zbar = z_loc_LE_w1 + (y_eng_ybar - y_offset_w1)*tan(dihedral_w1);

Geo_tier.x_eng_xbar = x_eng_xbar;
Geo_tier.y_eng_xbar = y_eng_ybar;
Geo_tier.z_eng_xbar = z_eng_zbar;

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