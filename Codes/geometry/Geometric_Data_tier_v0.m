function Geo_tier = Geometric_Data_tier
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

% Geometry
% Fuselage geomtry (Units in m)
w_fus = 0.26; % width
h_fus = 0.25; % height
l_fus = 1.6797; % length
d_fus = (w_fus + h_fus)/2; % diameter

Geo_tier.w_fus = w_fus;
Geo_tier.h_fus = h_fus;
Geo_tier.l_fus = l_fus;
Geo_tier.d_fus = d_fus;

%% 
% Geometric variations
% GEO_VAR = 1 then ProVant4 (Emergentia)
% GEO_VAR = 2 then Pepiño XXL
 GEO_VAR = 1;

% if GEO_VAR == 1 
%     geo_var = geometry_variations;
%     dihedral_vtail = geo_var.dihedral_vtail;
%     S_VTP = geo_var.S_VTP;
%     S_HTP = geo_var.S_HTP;
% end

% Data from analysis (Jorge Carreño)
% S_VTP = 0.0228*6;
% S_HTP = 0.0391*4;
% dihedral_vtail = atan((S_VTP/2)/S_HTP/2);
% dihedral_vtail_deg = dihedral_vtail*R2D;
% S_Vtail_1 = (S_VTP/2)/sin(dihedral_vtail);

% Entry Data, Dihedral (in RAD)
dihedral_w1_e = 5*D2R;
dihedral_w2_e = 29*D2R;
dihedral_w1 = dihedral_w1_e;
dihedral_w2 = dihedral_w2_e;
dihedral_v = 0*D2R;
Geo_tier.dihedral_w1 = dihedral_w1;
Geo_tier.dihedral_w2 = dihedral_w2;
Geo_tier.dihedral_v = dihedral_v;

% CAD meassurements
% Chord
cR_w1 = 0.226321;
cT_w1 = 0.147462; 
cR_w2 = 0.20453; 	
cT_w2 = 0.102565; 

% % Equivalent Horizontal Tail Surface
% b_HTP_eqv = (S_HTP/2)/((cR_w2+cT_w2)/2);
% b_VTP_eqv = (S_VTP/2)/((cR_w2+cT_w2)/2);

% % Real surface area and span
% b_VTAIL_eqv = sqrt(b_HTP_eqv^2 + b_VTP_eqv^2);
% S_VTAIL_v1 = b_VTAIL_eqv*(cR_w1+cT_w1)/2;
% S_VTAIL = S_VTAIL_v1*2;

% Wingspan
b_w1 = 2.509573;
b_w1_fus = 0.250;
b_w1_e = b_w1 - b_w1_fus; % effective 
b_w2 = 1.295646;
b_w2_fus = 0.212;
b_w2_e = b_w2 - b_w2_fus;

% Wing area 1
S_w1_e = b_w1_e*(cR_w1+cT_w1)/2;
S_w1_fus = cR_w1*b_w1_fus;
S_w1 = S_w1_e + S_w1_fus;

% Wing area 2
S_w2_e = b_w2_e*(cR_w2+cT_w2)/2;
S_w2_fus = cR_w2*b_w2_fus;
S_w2 = S_w2_e + S_w2_fus;

% Location of LE w1 
x_loc_LE_w1 = 0.695289;
y_loc_LE_w1 = b_w1_fus/2;
z_loc_LE_w1 = 0.155658;

% Location of LE w2 
x_loc_LE_w2 = 1.476053;
y_loc_LE_w2 = b_w2_fus/2;
z_loc_LE_w2 = 0.161743;

%%
%%%%% IMPORTANT %%%%
% The user select the reference wing area that will be used 
% 
% Reference area
S_ref = S_w1_e;

% Storing DATA
Geo_tier.cR_w1 = cR_w1; 	
Geo_tier.cT_w1 = cT_w1;	
Geo_tier.cR_w2 = cR_w2; 	
Geo_tier.cT_w2 = cT_w2;	
Geo_tier.b_w1 = b_w1;
Geo_tier.b_w2 = b_w2;
Geo_tier.b_w1_e = b_w1_e;
Geo_tier.b_w2_e = b_w2_e;
Geo_tier.S_w1 = S_w1;
Geo_tier.S_w2 = S_w2;
Geo_tier.S_w1_e  = S_w1_e; 
Geo_tier.S_w2_e  = S_w2_e; 
Geo_tier.S_ref = S_ref;

% Aspect Ratio
AR_w1 = b_w1^2/S_w1;
AR_w1_e = b_w1_e^2/S_w1_e;
AR_w2 = b_w2^2/S_w2;
AR_w2_e = b_w2_e^2/S_w2_e;

Geo_tier.AR_w1 = AR_w1;
Geo_tier.AR_w1_e = AR_w1_e;
Geo_tier.AR_w2 = AR_w2;
Geo_tier.AR_w2_e = AR_w2_e;

% Taper Ratio
lambda_w1 = cT_w1/cR_w1;
lambda_w1_e = cT_w1/cR_w1;
lambda_w2 = cT_w2/cR_w2;
lambda_w2_e = cT_w2/cR_w2;

Geo_tier.lambda_w1 = lambda_w1;
Geo_tier.lambda_w1_e = lambda_w1_e;
Geo_tier.lambda_w2 = lambda_w2;
Geo_tier.lambda_w2_e = lambda_w2_e;

% Entry Data, Surface LE sweep (in RAD)
% From trigonmetry
Lambda_LE_w1 = -7*D2R;
Lambda_LE_w2 = 10.65*D2R;
Lambda_LE_w1_e = Lambda_LE_w1;
Lambda_LE_w2_e = Lambda_LE_w2;

% % Surface TE sweep (in RAD)
% Lambda_TE_w1_e = atan((cT_w1)/(b_w1_e/2));
% Lambda_TE_w1 = Lambda_TE_w1_e;
% Lambda_TE_w2_e = atan((cT_w2)/(b_w2_e/2));
% Lambda_TE_w2 = Lambda_TE_w2_e;

% Equivalent
% Nth = 1/4;
% Lambda_c4_w1_e = atan((1/AR_w1_e)*(AR_w1_e*tan(Lambda_LE_w1_e) - 4*(Nth)*((1-lambda_w1_e)/(1+lambda_w1_e))))
% Lambda_c4_w2_e = atan((1/AR_w2_e)*(AR_w2_e*tan(Lambda_LE_w2_e) - 4*(Nth)*((1-lambda_w2_e)/(1+lambda_w2_e))))
% Lambda_c4_w1 = atan((1/AR_w1)*(4*(Nth)*((1-lambda_w1)/(1+lambda_w1)) + AR_w1*tan(Lambda_LE_w1)));
% Lambda_c4_w2 = atan((1/AR_w2)*(4*(Nth)*((1-lambda_w2)/(1+lambda_w2)) + AR_w2*tan(Lambda_LE_w2)));
% Lambda_c4_w1 = Lambda_c4_w1_e;
% Lambda_c4_w2 = Lambda_c4_w2_e;
Lambda_c4_w1 = Get_Nth_Lambda(1/4,AR_w1_e,Lambda_LE_w1,lambda_w1_e);
Lambda_c4_w2 = Get_Nth_Lambda(1/4,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
Lambda_c2_w1 = Get_Nth_Lambda(1/2,AR_w1_e,Lambda_LE_w1,lambda_w1_e);
Lambda_c2_w2 = Get_Nth_Lambda(1/2,AR_w2_e,Lambda_LE_w2,lambda_w2_e);
Lambda_TE_w1 = Get_Nth_Lambda(1,AR_w1_e,Lambda_LE_w1,lambda_w1_e);
Lambda_TE_w2 = Get_Nth_Lambda(1,AR_w2_e,Lambda_LE_w2,lambda_w2_e);

% Geo_tier.Lambda_LE_w1_e  = Lambda_LE_w1_e;
% pause
% 
% Nth = 1/2;
% Lambda_c2_w1_e = atan((1/AR_w1_e)*(AR_w1_e*tan(Lambda_LE_w1_e) - 4*(Nth)*((1-lambda_w1_e)/(1+lambda_w1_e))))
% Lambda_c2_w2_e = atan((1/AR_w2_e)*(AR_w2_e*tan(Lambda_LE_w2_e) - 4*(Nth)*((1-lambda_w2_e)/(1+lambda_w2_e))))
% % Lambda_c2_w1 = atan((1/AR_w1)*(4*(Nth)*((1-lambda_w1)/(1+lambda_w1)) + AR_w1*tan(Lambda_LE_w1)));
% % Lambda_c2_w2 = atan((1/AR_w2)*(4*(Nth)*((1-lambda_w2)/(1+lambda_w2)) + AR_w2*tan(Lambda_LE_w2)));
% Lambda_c2_w1 = Lambda_c2_w1_e;
% Lambda_c2_w2 = Lambda_c2_w2_e;
% 
% Nth = 1;
% Lambda_TE_w1_e = atan((1/AR_w1_e)*(AR_w1_e*tan(Lambda_LE_w1_e) - 4*(Nth)*((1-lambda_w1_e)/(1+lambda_w1_e))));
% Lambda_TE_w2_e = atan((1/AR_w2_e)*(AR_w2_e*tan(Lambda_LE_w2_e) - 4*(Nth)*((1-lambda_w2_e)/(1+lambda_w2_e))));
% Lambda_TE_w1 = Lambda_TE_w1_e;
% Lambda_TE_w2 = Lambda_TE_w2_e;

% Equivalent 
% Geo_tier.Lambda_LE_w1_e  = Lambda_LE_w1_e;
% Geo_tier.Lambda_LE_w2_e  = Lambda_LE_w2_e;
% Geo_tier.Lambda_TE_w1_e  = Lambda_TE_w1_e;
% Geo_tier.Lambda_TE_w2_e  = Lambda_TE_w2_e;
% Geo_tier.Lambda_TE_w1_e  = Lambda_TE_w1_e;
% Geo_tier.Lambda_TE_w2_e  = Lambda_TE_w2_e;
% Geo_tier.Lambda_c4_w1_e  = Lambda_c4_w1_e;
% Geo_tier.Lambda_c4_w2_e  = Lambda_c4_w2_e;
% Geo_tier.Lambda_c2_w1_e  = Lambda_c2_w1_e;
% Geo_tier.Lambda_c2_w2_e  = Lambda_c2_w2_e;

Geo_tier.Lambda_LE_w1  = Lambda_LE_w1;
Geo_tier.Lambda_LE_w2  = Lambda_LE_w2;
Geo_tier.Lambda_TE_w1  = Lambda_TE_w1;
Geo_tier.Lambda_TE_w2  = Lambda_TE_w2;
Geo_tier.Lambda_TE_w1  = Lambda_TE_w1;
Geo_tier.Lambda_TE_w2  = Lambda_TE_w2;
Geo_tier.Lambda_c4_w1  = Lambda_c4_w1;
Geo_tier.Lambda_c4_w2  = Lambda_c4_w2;
Geo_tier.Lambda_c2_w1  = Lambda_c2_w1;
Geo_tier.Lambda_c2_w2  = Lambda_c2_w2;

% 
% 
% 
% cmac_w1_e = (2/3)*cR_w1*(1 + lambda_w1_e + lambda_w1_e^2)/(1+lambda_w1_e);
% cmac_w2_e = (2/3)*cR_w2*(1 + lambda_w2_e + lambda_w2_e^2)/(1+lambda_w2_e);
% % cmac_w1   = (2/3)*cR_w1*(1 + lambda_w1 + lambda_w1^2)/(1+lambda_w1);
% % cmac_w2   = (2/3)*cR_w2*(1 + lambda_w2 + lambda_w2^2)/(1+lambda_w2);
% cmac_w1   = cmac_w1_e;
% cmac_w2   = cmac_w2_e;
% 
% Geo_tier.cmac_w1_e = cmac_w1_e;
% Geo_tier.cmac_w2_e = cmac_w2_e;
% Geo_tier.cmac_w1 = cmac_w1;
% Geo_tier.cmac_w2 = cmac_w2;

% % Location of LE w2 
% x_loc_LE_w1 = 0.695289;
% y_loc_LE_w1 = b_w1_fus/2;
% z_loc_LE_w1 = 0.155658;
% 
% % Location of LE w2 
% x_loc_LE_w2 = 1.476053;
% y_loc_LE_w2 = b_w2_fus/2;
% z_loc_LE_w2 = 0.161743;

% Obtain:
% Mean Aerodynamic Chord
% x,y,z locations of the MAC

% Geometric position of w1 MAC
XYZ_MAC = Get_MAC_Coordinates(b_w1_e,lambda_w1_e,y_loc_LE_w1,cR_w1,dihedral_w1,Lambda_c4_w1)
Geo_tier.cmac_w1 = XYZ_MAC.cbar;
Geo_tier.xbar_w1 = XYZ_MAC.xbar_w;
Geo_tier.ybar_w1 = XYZ_MAC.ybar_w;
Geo_tier.zbar_w1 = XYZ_MAC.zbar_w;

% Geometric position of w1 MAC
XYZ_MAC = Get_MAC_Coordinates(b_w2_e,lambda_w2_e,y_loc_LE_w2,cR_w2,dihedral_w2,Lambda_c4_w2)
Geo_tier.cmac_w2 = XYZ_MAC.cbar;
Geo_tier.xbar_w2 = XYZ_MAC.xbar_w;
Geo_tier.ybar_w2 = XYZ_MAC.ybar_w;
Geo_tier.zbar_w2 = XYZ_MAC.zbar_w;

% 
% XYZ_MAC = Get_MAC_Coordinates(b_w1_e,lambda_w1_e,y_loc_LE_w1,cR_w1,dihedral_w1,Lambda_c4_w1)
% 
% ybar_w1_e = (b_w1_e/2)*((1 + 2*lambda_w1_e)/(3*(1+lambda_w1_e))) + y_loc_LE_w1
% ybar_w2_e = (b_w2_e/2)*((1 + 2*lambda_w2_e)/(3*(1+lambda_w2_e))) + y_loc_LE_w2;
% ybar_w1   = ybar_w1_e;
% ybar_w2   = ybar_w2_e;
% 
% Geo_tier.ybar_w1 = ybar_w1;
% Geo_tier.ybar_w1_e = ybar_w1_e;
% Geo_tier.ybar_w2 = ybar_w2;
% Geo_tier.ybar_w2_e = ybar_w2_e;
% 
% % X-locatrion of MAC
% n = 1/4; % assume is at 25% chord
% xbar_w1_e = n*cR_w1 + ybar_w1_e*tan(Lambda_c4_w1_e)
% xbar_w2_e = n*cR_w2 + ybar_w2_e*tan(Lambda_c4_w2_e);
% xbar_w1   = xbar_w1_e;
% xbar_w2   = xbar_w2_e;
% 
% Geo_tier.xbar_w1 = xbar_w1;
% Geo_tier.xbar_w1_e = xbar_w1_e;
% Geo_tier.xbar_w2 = xbar_w2;
% Geo_tier.xbar_w2_e = xbar_w2_e;
% 
% % Z-locatrion of MAC
% zbar_w1_e   = ybar_w1_e*tan(dihedral_w1_e)
% zbar_w2_e   = ybar_w2_e*tan(dihedral_w2_e);
% zbar_w1   = zbar_w1_e;
% zbar_w2   = zbar_w2_e;
% 
% Geo_tier.zbar_w1 = zbar_w1;
% Geo_tier.zbar_w1_e = zbar_w1_e;
% Geo_tier.zbar_w2 = zbar_w2;
% Geo_tier.zbar_w2_e = zbar_w2_e;


% % Vertical Stabilizer
% cR_v = 0.150;
% cT_v = 0.05;	
% b_v = 0.3;
% 
% S_v1 = b_v*(cR_v + cT_v)/2;
% S_v2 = b_v*(cR_v + cT_v)/2;
% S_v = S_v1 + S_v2;
% AR_v = b_v^2/S_v1; % Alargamiento de una superficie
% 
% lambda_v = cT_v/cR_v;
% Lambda_LE_v = 19.3*D2R;
% % Lambda_TE_v = atan((0)/(b_v/2))
% Lambda_TE_v = atan((cT_v)/(b_v));
% Lambda_c4_v = atan((1/AR_v)*(4*(1/4)*((1-lambda_v)/(1+lambda_v)) + AR_v*tan(Lambda_LE_v)));
% Lambda_c2_v = atan((1/AR_v)*(4*(1/2)*((1-lambda_v)/(1+lambda_v)) + AR_v*tan(Lambda_LE_v)));
% cmac_v = (2/3)*cR_v*(1+lambda_v+lambda_v^2)/(1+lambda_v);
% ybar_v = (b_v)*((1 + 2*lambda_v)/(3*(1 + lambda_v)));
% n = 1/4;
% xbar_v = n*cR_v + ybar_v*tan(Lambda_c4_v);

% Storing DATA
% Geo_tier.cR_v = cR_v; 	
% Geo_tier.cT_v = cT_v;	
% Geo_tier.b_v	 = b_v;
% Geo_tier.S_v1 = S_v1;
% Geo_tier.S_v2 = S_v2;
% Geo_tier.S_v  = S_v; 
% Geo_tier.AR_v = AR_v;
% Geo_tier.lambda_v  = lambda_v;
% Geo_tier.Lambda_LE_v  = Lambda_LE_v;
% Geo_tier.Lambda_TE_v  = Lambda_TE_v;
% Geo_tier.Lambda_c4_v  = Lambda_c4_v;
% Geo_tier.Lambda_c2_v  = Lambda_c2_v;
% Geo_tier.cmac_v = cmac_v;
% Geo_tier.ybar_v = ybar_v;
% Geo_tier.xbar_v = xbar_v;

% Control surfaces
% Determines the ammount of control surface in the w2 aileron & elevator
% For Elevons
% Control_surface 0 = ELEVON;
% Control_surface 1 = RUDDERVATOR;
% Control_surface 2 = 2-SURFACE;
% Control_surface 3 = 3-SURFACE;

Control_surface = 1;

if Control_surface == 0
    K_ail_w2 = 1; % Percentage of aileron
    K_ele_w2 = 1; % complementary control surface
    K_can_w1 = 0.0; % Percentage of canard
    K_flap_w1 = 0.0;
    K_rudvtr_w2 = 0.0;
    
    K_ail = K_ail_w2; % Percentage of aileron
    K_ele = K_ele_w2; % complementary control surface
    K_can = K_can_w1; % Percentage of canard
    K_flap = K_flap_w1;
    K_rudvtr = K_rudvtr_w2;
    
elseif Control_surface == 1
    K_ail_w1 = 0.4; % Percentage of aileron
    K_flap_w1 =  1 - K_ail_w1; % Percentage of flaps
    K_ele_w2 = 0; % complementary control surface
    K_can_w1 = 0;
    K_rudvtr_w2 = 1; % Percentage of ruddervator
    
    K_ail = K_ail_w1; % Percentage of aileron
    K_ele = K_ele_w2; % complementary control surface
    K_can = K_can_w1; % Percentage of canard
    K_flap = K_flap_w1;
    K_rudvtr = K_rudvtr_w2;
        
elseif Control_surface == 2
    % For split surface Elevator and Ailero
    K_ail_w1 = 0.50; % Percentage of aileron
    K_flap_w1 =  1 - K_ail_w1; % Percentage of flaps
    K_ele_w2 = 1; % complementary control surface
    K_can_w1 = 0; % Percentage of canard
    K_rudvtr_w2 = 0;

    K_ail = K_ail_w1; % Percentage of aileron
    K_ele = K_ele_w2; % complementary control surface
    K_can = K_can_w1; % Percentage of canard
    K_flap = K_flap_w1;
    K_rudvtr = K_rudvtr_w2;

elseif Control_surface == 3
    % For split surface Elevator and Ailero
    K_can_w1 = 1.0; % Percentage of canard
    K_ail_w1 = 0.50; % Percentage of aileron
    K_flap_w1 =  1 - K_ail_w1; % Percentage of flaps
    K_ele_w2 = 1; % complementary control surface
    K_rudvtr_w2 = 0;

    K_ail = K_ail_w1; % Percentage of aileron
    K_ele = K_ele_w2; % complementary control surface
    K_can = K_can_w1; % Percentage of canard
    K_flap = K_flap_w1;
    K_rudvtr = K_rudvtr_w2;
end

% Storage
Geo_tier.K_ail = K_ail;
Geo_tier.K_ele = K_ele;
Geo_tier.K_can = K_can;
Geo_tier.K_flap = K_flap;
Geo_tier.K_rudvtr = K_rudvtr;

% AILERON
cf_ail = 0.25; % percentage of control surface
b_ail = b_w1_e*K_ail_w1; % length of aileron's (both surfaces)
% inner and outter location of the control surface
y2_b2_ail = b_w1/2; % outter position from the center line
y1_b2_ail = y2_b2_ail - b_ail/2; % inner position  from the center 
t_c_ail = 0.15; % Thinckness 2 chord ratio associated to the airfoil
y_offset_w1 = b_w1_fus/2;
% Calculates geometric information regarding the control surfaces
CS_geo = Get_ControlSurface_Coordinates(b_w1,K_ail,cR_w1,cf_ail,y2_b2_ail,y1_b2_ail,y_offset_w1,Lambda_LE_w1,Lambda_TE_w1,lambda_w1,Lambda_c4_w1,AR_w1,dihedral_w1);
Geo_tier.cmac_ail = CS_geo.cmac_cs;
Geo_tier.xbar_ail = CS_geo.xbar_cs;
Geo_tier.ybar_ail = CS_geo.ybar_cs;
Geo_tier.zbar_ail = CS_geo.zbar_cs;
Geo_tier.b_ail = CS_geo.b_cs;
Geo_tier.cR_ail = CS_geo.cR_cs;
Geo_tier.cT_ail = CS_geo.cT_cs;
Geo_tier.c_ail = CS_geo.c_cs;
Geo_tier.S_ail = CS_geo.S_cs;
Geo_tier.lambda_ail = CS_geo.lambda_cs;
Geo_tier.Lambda_LE_ail = CS_geo.Lambda_LE_cs;
Geo_tier.x_1R_y1_ail = CS_geo.x_1R_y1_cs;
Geo_tier.x_2R_y1_ail = CS_geo.x_2R_y1_cs;
Geo_tier.x_1R_y2_ail = CS_geo.x_1R_y2_cs;
Geo_tier.x_2R_y2_ail = CS_geo.x_2R_y2_cs;

% FLAP
cf_flap = 0.40; % percentage of control surface
b_flap = b_w1_e*K_flap_w1; % length of control surface's (both surfaces)
% inner and outter location of the control surface
y2_b2_flap = y1_b2_ail; % outter position from the center line
y1_b2_flap = y2_b2_flap - b_flap/2; % inner position  from the center 
t_c_ail = 0.15; % Thinckness 2 chord ratio associated to the airfoil
y_offset_w1 = b_w1_fus/2;

% Calculates geometric information regarding the control surfaces
CS_geo = Get_ControlSurface_Coordinates(b_w1,K_flap,cR_w1,cf_flap,y2_b2_flap,y1_b2_flap,y_offset_w1,Lambda_LE_w1,Lambda_TE_w1,lambda_w1,Lambda_c4_w1,AR_w1,dihedral_w1);
Geo_tier.cmac_ail = CS_geo.cmac_cs;
Geo_tier.xbar_ail = CS_geo.xbar_cs;
Geo_tier.ybar_ail = CS_geo.ybar_cs;
Geo_tier.zbar_ail = CS_geo.zbar_cs;
Geo_tier.b_ail = CS_geo.b_cs;
Geo_tier.cR_ail = CS_geo.cR_cs;
Geo_tier.cT_ail = CS_geo.cT_cs;
Geo_tier.c_ail = CS_geo.c_cs;
Geo_tier.S_ail = CS_geo.S_cs;
Geo_tier.lambda_ail = CS_geo.lambda_cs;
Geo_tier.Lambda_LE_ail = CS_geo.Lambda_LE_cs;
Geo_tier.x_1R_y1_ail = CS_geo.x_1R_y1_cs;
Geo_tier.x_2R_y1_ail = CS_geo.x_2R_y1_cs;
Geo_tier.x_1R_y2_ail = CS_geo.x_1R_y2_cs;
Geo_tier.x_2R_y2_ail = CS_geo.x_2R_y2_cs;

% RUDDERVATOR
cf_rudvtr = 0.25; % percentage of control surface
b_rudvtr = b_w2_e*K_rudvtr; % length of control surface's (both surfaces)
% inner and outter location of the control surface
y2_b2_rudvtr = b_w2/2; % outter position from the center line
y1_b2_rudvtr = y2_b2_rudvtr - b_rudvtr/2; % inner position  from the center 
t_c_rudvtr = 0.15; % Thinckness 2 chord ratio associated to teh airfoil
y_offset_w2 = b_w2_fus/2;

% Calculates geometric information regarding the control surfaces
CS_geo = Get_ControlSurface_Coordinates(b_w2,K_rudvtr,cR_w2,cf_rudvtr,y2_b2_rudvtr,y1_b2_rudvtr,y_offset_w2,Lambda_LE_w2,Lambda_TE_w2,lambda_w2,Lambda_c4_w2,AR_w2,dihedral_w2);
Geo_tier.cmac_rudvtr = CS_geo.cmac_cs;
Geo_tier.xbar_rudvtr = CS_geo.xbar_cs;
Geo_tier.ybar_rudvtr = CS_geo.ybar_cs;
Geo_tier.zbar_rudvtr = CS_geo.zbar_cs;
Geo_tier.b_rudvtr = CS_geo.b_cs;
Geo_tier.cR_rudvtr = CS_geo.cR_cs;
Geo_tier.cT_rudvtr = CS_geo.cT_cs;
Geo_tier.c_rudvtr = CS_geo.c_cs;
Geo_tier.S_rudvtr = CS_geo.S_cs;
Geo_tier.lambda_rudvtr = CS_geo.lambda_cs;
Geo_tier.Lambda_LE_rudvtr = CS_geo.Lambda_LE_cs;
Geo_tier.x_1R_y1_rudvtr = CS_geo.x_1R_y1_cs;
Geo_tier.x_2R_y1_rudvtr = CS_geo.x_2R_y1_cs;
Geo_tier.x_1R_y2_rudvtr = CS_geo.x_1R_y2_cs;
Geo_tier.x_2R_y2_rudvtr = CS_geo.x_2R_y2_cs;

% % 
% % 
% % 
% % 
% % 
% % 
% % cf_rudvtr = 0.25; % percentage of control surface
% % cR_rudvtr = cR_w2*cf_rudvtr;
% % cT_rudvtr = cT_w2*cf_rudvtr;
% % c_rudvtr = (2/3)*cR_rudvtr*(1+lambda_w2 + lambda_w2^2)/(1+lambda_w2);
% % b_rudvtr = b_w2_e*K_rudvtr_w2;
% % b_rudvtr_b_w2 = b_rudvtr/b_w2;
% % % inner and outter location of the control surface
% % y2_b2_rudvtr = b_w2/2;
% % y1_b2_rudvtr = b_w2/2 - b_rudvtr/2;
% % 
% % 
% % % FLAP
% % cf_flap = 0.25; % percentage of control surface
% % cR_flap = cR_w1*cf_flap;
% % cT_flap = cT_w1*cf_flap;
% % c_flap = (2/3)*cR_flap*(1+lambda_w1+lambda_w1^2)/(1+lambda_w1);
% % b_flap = b_w1_e*K_flap_w1;
% % b_flap_b_w1 = b_flap/b_w1;
% % % inner and outter location of the control surface
% % y2_b2_flap = y1_b2_ail;
% % y1_b2_flap = b_w1/2 - b_ail/2;
% % t_c_flap = 0.15; % Thinckness 2 chord ratio associated to teh airfoil
% 
% % % RUDDER
% % cf_rud = 0.25; % percentage of control surface
% % cR_rud = cR_v*cf_rud;
% % cT_rud = cT_v*cf_rud;
% % c_rud = (2/3)*cR_rud*(1+lambda_v+lambda_v^2)/(1+lambda_v);
% % b_rud = b_v;
% % b_rud_b_v = b_rud/b_v;
% % % inner and outter location of the control surface
% % y2_b2_rud = b_v;
% % y1_b2_rud = b_v - b_rud;
% % t_c_rud = 0.12; % Thinckness 2 chord ratio associated to teh airfoil
% 
% % AILERON
% cf_ele = 0.25; % percentage of control surface
% cR_ail = cR_w1*cf_ail;
% cT_ail = cT_w1*cf_ail;
% c_ail = (2/3)*cR_ail*(1+lambda_w1+lambda_w1^2)/(1+lambda_w1);
% b_ail = b_w1_e*K_ail_w1;
% b_ail_b_w1 = b_ail/b_w1;
% % inner and outter location of the control surface
% y2_b2_ail = b_w2/2;
% y1_b2_ail = b_w2/2 - b_ail/2;
% t_c_ail = 0.15; % Thinckness 2 chord ratio associated to teh airfoil
% 
% 
% 
% % Storing DATA
% Geo_tier.cf_ail = cf_ail;
% Geo_tier.cf_rudvtr = cf_rudvtr;
% Geo_tier.cf_flap = cf_flap;
% 
% Geo_tier.c_rudvtr = c_rudvtr;
% Geo_tier.b_rudvtr = b_rudvtr;
% Geo_tier.y1_b2_rudvtr = y1_b2_rudvtr;
% Geo_tier.y2_b2_rudvtr = y2_b2_rudvtr;
% Geo_tier.cf_rudvtr = cf_rudvtr;
% Geo_tier.t_c_rudvtr = t_c_rudvtr;
% 
% Geo_tier.c_ail = c_ail;	
% Geo_tier.b_ail = b_ail;
% Geo_tier.y1_b2_ail = y1_b2_ail;
% Geo_tier.y2_b2_ail = y2_b2_ail;
% Geo_tier.cf_ail = cf_ail;
% Geo_tier.t_c_ail = t_c_ail;
% 
% Geo_tier.c_flap = c_flap;	
% Geo_tier.b_flap = b_flap;
% Geo_tier.y1_b2_flap = y1_b2_flap;
% Geo_tier.y2_b2_flap = y2_b2_flap;
% Geo_tier.cf_flap = cf_flap;
% Geo_tier.t_c_flap = t_c_flap;
% 
% Geo_tier.b_flap_b_w1 = b_flap_b_w1;
% Geo_tier.b_rudvtr_b_w2 = b_rudvtr_b_w2;
% Geo_tier.b_ail_b_w1 = b_ail_b_w1;

% Distances
x_w1_LE = 0.6968;
x_w1_xbar = x_w1_LE + xbar_w1;
x_w1_xbar_e = x_w1_LE + xbar_w1_e;

x_w2_LE = 1.4751;
x_w2_xbar = x_w2_LE + xbar_w2;
x_w2_xbar_e = x_w2_LE + xbar_w2_e;

% Positive up
% z_w1_LE = 0.108042396441479;
z_w1_LE = 0.1771;
% z_w2_LE = 0.093715591696784;
z_w2_LE = 0.1171;
z_w1_xbar = z_w1_LE + zbar_w1;
z_w1_xbar_e = z_w1_LE + zbar_w1_e;

% Storing DATA
Geo_tier.x_w1_LE = x_w1_LE;
Geo_tier.x_w1_xbar = x_w1_xbar;
Geo_tier.x_w1_xbar_e = x_w1_xbar_e;
Geo_tier.x_w2_LE = x_w2_LE;
Geo_tier.x_w2_xbar = x_w2_xbar;
Geo_tier.x_w2_xbar_e = x_w2_xbar_e;

Geo_tier.z_w1_LE = z_w1_LE;
Geo_tier.z_w2_LE = z_w2_LE;

Geo_tier.z_w1_xbar = z_w1_xbar;
Geo_tier.z_w1_xbar_e = z_w1_xbar_e;

% Geo_tier.x_v_LE = x_v_LE;
% Geo_tier.x_v_xbar = x_v_xbar;
% Geo_tier.z_v_LE = z_v_LE;
% Geo_tier.y_v_LE = y_v_LE;

% Tail Volume Coefficients
l_w1w2 = x_w2_xbar - x_w1_xbar; % from xac_wing1 to xac_wing2
l_w1w2_e = x_w2_xbar_e - x_w1_xbar_e; % from xac_wing to xac_wing2

% l_vt1 = x_v_xbar - x_w1_xbar; % from xac_wing to xac_wing1
% l_vt2 = x_v_xbar - x_w2_xbar; % from xac_wing to xac_wing2
% 
% l_vt1_e = x_v_xbar - x_w1_xbar_e; % from xac_wing to xac_wing1
% l_vt2_e = x_v_xbar - x_w2_xbar_e; % from xac_wing to xac_wing2

% fuselage length
lfus_b_w1 = l_fus/b_w1;
lfus_b_w2 = l_fus/b_w2;

Geo_tier.lfus_b_w1 = lfus_b_w2;
Geo_tier.lfus_b_w2 = lfus_b_w2;

% tail volume coefficient
% Cw1 = S_w1*l_w1w2/(cmac_w1*S_w1);
Cw2 = S_w2*l_w1w2/(cmac_w2*S_w2);
% CVT1 = S_v*l_vt1/(b_w1*S_w1);
% CVT2 = S_v*l_vt2/(b_w2*S_w2);

% Storing DATA
Geo_tier.l_w1w2 = l_w1w2;
Geo_tier.l_w1w2_e = l_w1w2_e;
% Geo_tier.l_vt1 = l_vt1;
% Geo_tier.l_vt2 = l_vt2;
% Geo_tier.l_vt1_e = l_vt1_e;
% Geo_tier.l_vt2_e = l_vt2_e;
Geo_tier.l_fus = l_fus;

% Geo_tier.Cw1 = Cw1;
Geo_tier.Cw2 = Cw2;
% Geo_tier.CVT1 = CVT1;
% Geo_tier.CVT2 = CVT2;

% engine
x_eng = 2.1;
z_eng = 0.02;
Geo_tier.x_eng = x_eng;
Geo_tier.z_eng = z_eng;

%% Propulsive model
% AXI 5360/24HD V2 GOLD LINE
l_eng = 0.104; %length
d_eng = 0.063; %diameter
n_eng=2; % Engine number.
n_fan=2; % Number of fans

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

Geo_tier.W_eng = W_eng;
Geo_tier.n_eng = n_eng;
Geo_tier.n_fan = n_fan;

% Propeller data
% Datos genéricos Hélice para 22x12W
b_p = (22*2.54/100);
c_p = 3/100;
b_p_c_p = b_p/c_p;
c_prop = D_prop/b_p_c_p;
S_prop = D_prop*c_prop;
AR_prop = (D_prop^2)/S_prop;
RPM_max = 145000/D_propin; % Max RPM by engine builder.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Engine properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Engine Data
% Geometry Engien Assume AXI 2826/12 GOLD LINE AND MULTIPLIES BY 2 DIMENSIONS
eng_dia = 0.0350;
eng_length = 0.0540;
Geo_tier.eng_dia = eng_dia;
Geo_tier.eng_length = eng_length;

save Geo_tier.mat