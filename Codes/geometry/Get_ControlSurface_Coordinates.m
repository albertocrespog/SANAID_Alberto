function CS_geo = Get_ControlSurface_Coordinates(y_1R_y1,y_2R_y1,y_1R_y2,y_2R_y2,...
    x_1R_y1_w,x_2R_y1_w,y_offset_w,Lambda_LE_w,Lambda_TE_w,Lambda_c4_w,dihedral_w,cf,AC_CONFIGURATION,VTP_cs)

% chord at each location
% _as area surface
% _cs control surface
% Location of the inner wing area with control surface (Area Surface)
% x-location
x_1R_y1_as = x_1R_y1_w  + (y_1R_y1 - y_offset_w)*tan(Lambda_LE_w); % defines inner position of wing chord LE
x_2R_y1_as = x_2R_y1_w  + (y_2R_y1 - y_offset_w)*tan(Lambda_TE_w); % defines inner position of wing chord TE
x_1R_y2_as = x_1R_y1_w  + (y_1R_y2 - y_offset_w)*tan(Lambda_LE_w); % defines inner position of wing chord LE
x_2R_y2_as = x_2R_y1_w  + (y_2R_y2 - y_offset_w)*tan(Lambda_TE_w); % defines inner position of wing chord TE
% Chord
c_y1_as = x_2R_y1_as - x_1R_y1_as; % defines chord of wing at beginning of control surface
c_y2_as = x_2R_y2_as - x_1R_y2_as; % defines chord of wing at beginning of control surface
% y-location
y_1R_y1_as = y_1R_y1;
y_2R_y1_as = y_2R_y1;
y_1R_y2_as = y_1R_y2;
y_2R_y2_as = y_2R_y2;
% z-location
z_1R_y1_as = (y_1R_y1 - y_offset_w)*tan(dihedral_w);
z_2R_y1_as = (y_2R_y1 - y_offset_w)*tan(dihedral_w);
z_1R_y2_as = (y_1R_y2 - y_offset_w)*tan(dihedral_w);
z_2R_y2_as = (y_2R_y2 - y_offset_w)*tan(dihedral_w);

% Storing DATA
CS_geo.x_1R_y1_as = x_1R_y1_as;
CS_geo.x_2R_y1_as = x_2R_y1_as;
CS_geo.x_1R_y2_as = x_1R_y2_as;
CS_geo.x_2R_y2_as = x_2R_y2_as;
CS_geo.y_1R_y1_as = y_1R_y1_as;
CS_geo.y_2R_y1_as = y_2R_y1_as;
CS_geo.y_1R_y2_as = y_1R_y2_as;
CS_geo.y_2R_y2_as = y_2R_y2_as;
CS_geo.z_1R_y1_as = z_1R_y1_as;
CS_geo.z_2R_y1_as = z_2R_y1_as;
CS_geo.z_1R_y2_as = z_1R_y2_as;
CS_geo.z_2R_y2_as = z_2R_y2_as;
CS_geo.c_y1_as = c_y1_as;
CS_geo.c_y2_as = c_y2_as;

% Geometry
lambda_as = c_y2_as/c_y1_as; % aileron taper ratio
b_as = y_1R_y2_as - y_1R_y1_as; % same span 
S_as = b_as*((c_y1_as + c_y2_as)/2); % total area effective surface
AR_as = b_as^2/S_as; % Aspect Ratio

% Storing DATA
CS_geo.lambda_as = lambda_as;
CS_geo.b_as = b_as;
CS_geo.S_as = S_as;
CS_geo.AR_as = AR_as;

%% Control Surface
% x-location
x_1R_y1_cs = x_1R_y1_as + c_y1_as*(1-cf); % Defines inner position of LE control Surface
x_2R_y1_cs = x_1R_y1_as + c_y1_as; % Defines inner position of TE control Surface
x_1R_y2_cs = x_1R_y2_as + c_y2_as*(1-cf); % Defines outter position of LE control Surface
x_2R_y2_cs = x_1R_y2_as + c_y2_as; % Defines outter position of TE control Surface
% Chord
c_y1_cs = x_2R_y1_cs - x_1R_y1_cs; % defines chord of wing at beginning of control surface
c_y2_cs = x_2R_y2_cs - x_1R_y2_cs; % defines chord of wing at beginning of control surface
% y-location
y_1R_y1_cs = y_1R_y1_as;
y_2R_y1_cs = y_2R_y1_as;
y_1R_y2_cs = y_1R_y2_as;
y_2R_y2_cs = y_2R_y2_as;
% z-location
z_1R_y1_cs = z_1R_y1_as;
z_2R_y1_cs = z_2R_y1_as;
z_1R_y2_cs = z_1R_y2_as;
z_2R_y2_cs = z_2R_y2_as;

% Storing DATA
CS_geo.x_1R_y1_cs = x_1R_y1_cs;
CS_geo.x_2R_y1_cs = x_2R_y1_cs;
CS_geo.x_1R_y2_cs = x_1R_y2_cs;
CS_geo.x_2R_y2_cs = x_2R_y2_cs;
CS_geo.y_1R_y1_cs = y_1R_y1_cs;
CS_geo.y_2R_y1_cs = y_2R_y1_cs;
CS_geo.y_1R_y2_cs = y_1R_y2_cs;
CS_geo.y_2R_y2_cs = y_2R_y2_cs;
CS_geo.z_1R_y1_cs = z_1R_y1_cs;
CS_geo.z_2R_y1_cs = z_2R_y1_cs;
CS_geo.z_1R_y2_cs = z_1R_y2_cs;
CS_geo.z_2R_y2_cs = z_2R_y2_cs;
CS_geo.c_y1_cs = c_y1_cs;
CS_geo.c_y2_cs = c_y2_cs;

% Geometry
lambda_cs = c_y2_cs/c_y1_cs; % aileron taper ratio
b_cs = y_1R_y2_cs - y_1R_y1_cs; % same span 
S_cs = b_as*((c_y1_cs + c_y2_cs)/2); % total area effective surface
AR_cs = b_cs^2/S_cs; % Aspect Ratio
% Storing DATA
CS_geo.lambda_cs = lambda_cs;
CS_geo.b_cs = b_cs;
CS_geo.S_cs = S_cs;
CS_geo.AR_cs = AR_cs;

% Calculate the Control Surface LE sweep (assuming sweep)
Lambda_LE_as = Lambda_LE_w;
Lambda_LE_as_c4 = Lambda_c4_w; % at the 1/4 of the chord of the contros surface
Lambda_TE_as = Lambda_TE_w;
% Calculate the Control Surface LE sweep (assuming sweep)
Lambda_LE_cs = Get_Nth_Lambda(1-cf,AR_cs,Lambda_LE_w,lambda_cs);
Lambda_LE_cs_c4 = Get_Nth_Lambda((1-(3/4)*cf),AR_cs,Lambda_LE_w,lambda_cs); % at the 1/4 of the chord of the contros surface
Lambda_TE_cs = Lambda_TE_w;
% Storing DATA
CS_geo.Lambda_LE_as = Lambda_LE_as;
CS_geo.Lambda_LE_as_c4 = Lambda_LE_as_c4;
CS_geo.Lambda_TE_as = Lambda_TE_as;
CS_geo.Lambda_LE_cs = Lambda_LE_cs;
CS_geo.Lambda_LE_cs_c4 = Lambda_LE_cs_c4;
CS_geo.Lambda_TE_cs = Lambda_TE_cs;

% Distances of the aerodynamic center of the CS and AS relative to the LE
% inner point
if VTP_cs == 1
    XYZ_MAC = Get_MAC_Coordinates_VTP(b_as,lambda_as,c_y1_as,dihedral_w,Lambda_LE_as_c4);
    CS_geo.cmac_as = XYZ_MAC.cbar;
    CS_geo.xbar_as = XYZ_MAC.xbar_w;
    CS_geo.ybar_as = XYZ_MAC.ybar_w;
    CS_geo.zbar_as = XYZ_MAC.zbar_w;

    XYZ_MAC = Get_MAC_Coordinates_VTP(b_cs,lambda_cs,c_y1_cs,dihedral_w,Lambda_LE_cs_c4);
    CS_geo.cmac_cs = XYZ_MAC.cbar;
    CS_geo.xbar_cs = XYZ_MAC.xbar_w;
    CS_geo.ybar_cs = XYZ_MAC.ybar_w;
    CS_geo.zbar_cs = XYZ_MAC.zbar_w;
else
    XYZ_MAC = Get_MAC_Coordinates(b_as,lambda_as,c_y1_as,dihedral_w,Lambda_LE_as_c4);
    CS_geo.cmac_as = XYZ_MAC.cbar;
    CS_geo.xbar_as = XYZ_MAC.xbar_w;
    CS_geo.ybar_as = XYZ_MAC.ybar_w;
    CS_geo.zbar_as = XYZ_MAC.zbar_w;

    XYZ_MAC = Get_MAC_Coordinates(b_cs,lambda_cs,c_y1_cs,dihedral_w,Lambda_LE_cs_c4);
    CS_geo.cmac_cs = XYZ_MAC.cbar;
    CS_geo.xbar_cs = XYZ_MAC.xbar_w;
    CS_geo.ybar_cs = XYZ_MAC.ybar_w;
    CS_geo.zbar_cs = XYZ_MAC.zbar_w;
end

% Distances relative to the origin
% x_xbar_cs = x_loc_LE_w2 + Geo_tier.xbar_w2;
% y_ybar_cs = y_loc_LE_w2 + Geo_tier.ybar_w2;
% z_zbar_cs = z_loc_LE_w2 + Geo_tier.zbar_w2;
