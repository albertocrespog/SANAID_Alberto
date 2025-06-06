function [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts] = ...
    Calculo_Stability_Derivatives_Abril2020_untouched(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
    Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim)

%% identifies the aerodynamic surfaces being used
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Nac = AC_CONFIGURATION.Nac;

twin_VTP = AC_CONFIGURATION.twin_VTP;
d_ail = AC_CONFIGURATION.d_ail;
d_ele = AC_CONFIGURATION.d_ele;
d_elevon = AC_CONFIGURATION.d_elevon;
d_flap = AC_CONFIGURATION.d_flap;
d_rudder = AC_CONFIGURATION.d_rudder;
d_rudvtr = AC_CONFIGURATION.d_rudvtr;
d_can = AC_CONFIGURATION.d_can;
AC_type = AC_CONFIGURATION.AC_type;

% Conversion units
g = conv_UNITS.g;
in2m = conv_UNITS.in2m;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;
ft2m = conv_UNITS.ft2m;
m2ft = conv_UNITS.m2ft;
W2hp = conv_UNITS.W2hp;
mps2ftps = conv_UNITS.mps2ftps;
kg2lb = conv_UNITS.kg2lb;
lb2kg = conv_UNITS.lb2kg;
ftpm2mps = conv_UNITS.ftpm2mps;
W2pftsec = conv_UNITS.W2pftsec;
m22ft2 = conv_UNITS.m22ft2;
rho_SI2rho_IMP = conv_UNITS.rho_SI2rho_IMP;
qmet2qimp = conv_UNITS.qmet2qimp;
N2lbf = conv_UNITS.N2lbf;


%% Mass and inertias
% m_TOW = Weight_tier.m_TOW;
m_TOW = conditions.m_TOW;
w_T0 = m_TOW*g;
alpha_f = conditions.alpha_f;
beta_f = conditions.beta_f;

% Number of engines
n_eng = Prop_data.n_eng;
case_AC = AC_CONFIGURATION.case_AC;

% Performance
alpha_max = Performance.alpha_max;
V_stall = Performance.V_stall;
V_min = Performance.V_min;
% Flight Conditions
h = conditions.h;
V = conditions.V;
rho = Performance.rho;
q_inf = 0.5*rho*V^2;                     %presi�n din�mica inicial
a = Performance.a;
% Mach Number
Mach = V/a;

CL_max_w1 = Performance.CL_max_w1;
CL_max_w1_ope = Performance.CL_max_w1_ope;
alpha_max_w1_ope = Performance.alpha_max_w1_ope;

% Xcg estimation
% x_XCG = XCG_data.x_XCG;
x_XCG = conditions.x_XCG;
y_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG;
z_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;

% Max values controil surfaces
delta_ail_min = Geo_tier.delta_ail_min;
delta_ail_max = Geo_tier.delta_ail_max;
delta_ele_min = Geo_tier.delta_ele_min;
delta_ele_max = Geo_tier.delta_ele_max;
delta_elevon_min = Geo_tier.delta_elevon_min;
delta_elevon_max = Geo_tier.delta_elevon_max;
delta_flap_min = Geo_tier.delta_flap_min;
delta_flap_max = Geo_tier.delta_flap_max;
delta_rudder_min = Geo_tier.delta_rudder_min;
delta_rudder_max = Geo_tier.delta_rudder_max;
delta_rudvtr_max = Geo_tier.delta_rudvtr_max;
delta_rudvtr_min = Geo_tier.delta_rudvtr_min;
delta_can_min = Geo_tier.delta_can_min;
delta_can_max = Geo_tier.delta_can_max;

% Fuselage geometry
Surf_TOT = Body_Geo.Surf_TOT;
Vol_TOT = Body_Geo.Vol_TOT;
z_Area_b_1_4 = Body_Geo.z_Area_b_1_4;
z_Area_b_3_4 = Body_Geo.z_Area_b_3_4;
x_Area_b_1_4 = Body_Geo.x_Area_b_1_4;
x_Area_b_3_4 = Body_Geo.x_Area_b_3_4;
x_Area_b_max = Body_Geo.x_Area_b_max;
Area_b_max = Body_Geo.Area_b_max;
w_Area_b_max = Body_Geo.w_Area_b_max;
h_Area_b_max = Body_Geo.h_Area_b_max;
Area_top = Body_Geo.Area_top;
Area_side = Body_Geo.Area_side;
length_fus = Body_Geo.l_fus;
Area_body = Body_Geo.Area_body;
x_Area_body = Body_Geo.x_Area_body;
CLa_fus = Body_Geo.CLa_fus;

Height_x_fus = Body_Geo.Height;
Bisectriz_z_x_fus = Body_Geo.Bisectriz_z;
Centroid_z_x_fus = Body_Geo.Centroid_z;
Angle_fus_x_fus = Body_Geo.Angle_fus;
Angle_fus_interp = Body_Geo.Angle_fus_interp;
x_interp = Body_Geo.x_interp;

length_x_position = Body_Geo.length_x_position;
width_x_position = Body_Geo.width_x_position;
height_x_position = Body_Geo.height_x_position;

% reference area
S_ref = Geo_tier.S_ref;
 
%% Geometry of Wing - W1
if W1 == 1
    % Design Selected incidences
    i_w1 = Design_criteria.i_w1;
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_w1 = Design_criteria.index_w1;
    %% Geometry
    % Wingspan
    b_w1 = Geo_tier.b_w1;
    b_w1_e = Geo_tier.b_w1_e;
    b_w1_pv = Geo_tier.b_w1_pv;
    % Chord
    cR_w1 = Geo_tier.cR_w1;
    cT_w1 = Geo_tier.cT_w1;
    cmac_w1 = Geo_tier.cmac_w1;
    cmac_w1_e = Geo_tier.cmac_w1;
    % Positions
    % Distances relative to the origin
    % Position (X,Y,Z) of the Xac for both CR and TOLD
    x_xbar_w1 = Geo_tier.x_xbar_w1;
    y_ybar_w1 = Geo_tier.y_ybar_w1;
    z_zbar_w1 = Geo_tier.z_zbar_w1;
    % Distances relative to the LE
    xbar_w1 = Geo_tier.xbar_w1;
    xbar_w1_e = Geo_tier.xbar_w1; % effective
    % Distances (X,Z) from origin to LE
    x_w1_LE = Geo_tier.x_w1_LE;
    z_w1_LE = Geo_tier.z_w1_LE;
    % Distances (X,Y,Z) from origin to LE
    x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
    y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
    z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
    % Distances (X,Y,Z) from origin to TE
    x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
    y_cR_w1_TE = Geo_tier.y_cR_w1_TE;
    z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
    % Areas
    S_w1 = Geo_tier.S_w1; 
    S_w1_e = Geo_tier.S_w1_e; % effective wing area: proyected in the y-plane
    S_w1_s = Geo_tier.S_w1_s; % real wing area: proyected along the surface of wing
    % prop wash area affecting surfaces
    S_w1_pw = Geo_tier.S_w1_pw;
    S_w1_afe = S_w1_pw;
    % Aspect Ratio
    AR_w1 = Geo_tier.AR_w1;
    AR_w1_e = Geo_tier.AR_w1_e;
    % Taper ratio
    lambda_w1 = Geo_tier.lambda_w1;
    lambda_w1_e = Geo_tier.lambda_w1_e;
    % Angles
    % Sweep LE
    Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
    Lambda_LE_w1_e = Geo_tier.Lambda_LE_w1_e;
    % Sweept c/4
    Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
    Lambda_c2_w1 = Geo_tier.Lambda_c2_w1;
    % dihedral
    dihedral_w1 = Geo_tier.dihedral_w1;
    dihedral_w1_e = Geo_tier.dihedral_w1_e;
    %% Aerodynamic properties
    % Condiciones de vuelo
    % Lift coefficient
    CL0_w1 = Aero.CL_0_w1_CR;
    CL0_w1_e = CL0_w1;
C
    alpha_CL_0_w1 = Aero.alpha_CL_0_w1_CR;
    % Drag
    CD0_w1 = Aero_TH.CD0_w1;
    % Moment coefficient (CM)
    CM0_w1 = Aero.CM_0_w1_CR;
    CM0_w1_e = CM0_w1;
end

%% W2
if HTP == 1
    % Design Selected incidences
    i_w2 = Design_criteria.i_w2;
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_w2 = Design_criteria.index_w2;
    %% Geometry
    % Wingspan
    b_w2 = Geo_tier.b_w2;
    b_w2_e = Geo_tier.b_w2_e;
    b_w2_pv = Geo_tier.b_w2_pv;
    % Chord
    cR_w2 = Geo_tier.cR_w2;
    cT_w2 = Geo_tier.cT_w2;
    cmac_w2 = Geo_tier.cmac_w2;
    cmac_w2_e = Geo_tier.cmac_w2;
    % Arm from Xac w1 to w2
    l_xac_w1w2 = Geo_tier.l_xac_w1w2;
    % Positions
    % Distances relative to the origin
    % Position (X,Y,Z) of the Xac for both CR and TOLD
    x_xbar_w2 = Geo_tier.x_xbar_w2;
    y_ybar_w2 = Geo_tier.y_ybar_w2;
    z_zbar_w2 = Geo_tier.z_zbar_w2;
    % Distances relative to the LE
    xbar_w2 = Geo_tier.xbar_w2;
    xbar_w2_e = Geo_tier.xbar_w2; % effective
    % Distances (X,Z) from origin to LE
    x_w2_LE = Geo_tier.x_w2_LE;
    z_w2_LE = Geo_tier.z_w2_LE;
    % Distances (X,Y,Z) from origin to LE
    x_cR_w2_LE = Geo_tier.x_cR_w2_LE;
    y_cR_w2_LE = Geo_tier.y_cR_w2_LE;
    z_cR_w2_LE = Geo_tier.z_cR_w2_LE;
    % Distances (X,Y,Z) from origin to TE
    x_cR_w2_TE = Geo_tier.x_cR_w2_TE;
    y_cR_w2_TE = Geo_tier.y_cR_w2_TE;
    z_cR_w2_TE = Geo_tier.z_cR_w2_TE;
    % Areas
    S_w2 = Geo_tier.S_w2;
    S_w2_e = Geo_tier.S_w2_e; % effective wing area: proyected in the y-plane
    S_w2_s = Geo_tier.S_w2_s; % real wing area: proyected along the surface of wing
    % prop wash area affecting surfaces
    S_w2_pw = Geo_tier.S_w2_pw;
    S_w2_afe = S_w2_pw;
    % Aspect Ratio
    AR_w2 = Geo_tier.AR_w2;
    AR_w2_e = Geo_tier.AR_w2_e;
    % Taper Ratio
    lambda_w2 = Geo_tier.lambda_w2;
    lambda_w2_e = Geo_tier.lambda_w2_e;
    % Angles
    % Sweep LE
    Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
    Lambda_LE_w2_e = Geo_tier.Lambda_LE_w2_e;
    % Sweep c/4
    Lambda_c4_w2 = Geo_tier.Lambda_c4_w2;
    Lambda_c2_w2 = Geo_tier.Lambda_c2_w2;
    % dihedral
    dihedral_w2 = Geo_tier.dihedral_w2;
    dihedral_w2_e = Geo_tier.dihedral_w2_e;
    % Tail Volume Coefficient
    Vbar_h = ((x_xbar_w2 - x_XCG)/cmac_w1)*(S_w2_e/S_ref);
    %% Aerodynamic properties
    % Condiciones de vuelo
    % Lift coefficient
    CL0_w2 = Aero.CL_0_w2_CR;
    CL0_w2_e = CL0_w2;
    CL_alpha_w2 = Aero.CL_alpha_w2_CR;
    CLalpha_w2_e = CL_alpha_w2;
    alpha_CL_0_w2 = Aero.alpha_CL_0_w2_CR;
    % Drag
    CD0_w2 = Aero_TH.CD0_w2;
    % Moment coefficient (CM)
    CM0_w2 = Aero.CM_0_w2_CR;
    CM0_w2_e = CM0_w2;
end

%% Vee tail
if Vee == 1
    % Design Selected incidences
    i_w2 = Design_criteria.i_w2;
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_w2 = Design_criteria.index_w2;
    %% Geometry
    % Wingspan
    b_w2 = Geo_tier.b_w2;
    b_w2_e = Geo_tier.b_w2_e;
    b_w2_pv = Geo_tier.b_w2_pv;
    % Chord
    cR_w2 = Geo_tier.cR_w2;
    cT_w2 = Geo_tier.cT_w2;
    cmac_w2 = Geo_tier.cmac_w2;
    cmac_w2_e = Geo_tier.cmac_w2;
    % Arm from Xac w1 to w2
    l_xac_w1w2 = Geo_tier.l_xac_w1w2;
    % Distances relative to the origin
    % Position (X,Y,Z) of the Xac for both CR and TOLD
    x_xbar_w2 = Geo_tier.x_xbar_w2;
    y_ybar_w2 = Geo_tier.y_ybar_w2;
    z_zbar_w2 = Geo_tier.z_zbar_w2;
    % Distances relative to the LE
    xbar_w2 = Geo_tier.xbar_w2;
    xbar_w2_e = Geo_tier.xbar_w2; % effective
    % Distances (X,Z) from origin to LE
    x_w2_LE = Geo_tier.x_w2_LE;
    z_w2_LE = Geo_tier.z_w2_LE;
    % Distances (X,Y,Z) from origin to LE
    x_cR_w2_LE = Geo_tier.x_cR_w2_LE;
    y_cR_w2_LE = Geo_tier.y_cR_w2_LE;
    z_cR_w2_LE = Geo_tier.z_cR_w2_LE;
    % Distances (X,Y,Z) from origin to TE
    x_cR_w2_TE = Geo_tier.x_cR_w2_TE;
    y_cR_w2_TE = Geo_tier.y_cR_w2_TE;
    z_cR_w2_TE = Geo_tier.z_cR_w2_TE;
    % Areas
    S_w2 = Geo_tier.S_w2;
    S_w2_e = Geo_tier.S_w2_e; % effective wing area: proyected in the y-plane
    S_w2_s = Geo_tier.S_w2_s; % real wing area: proyected along the surface of wing
    % prop wash area affecting surfaces
    S_w2_pw = Geo_tier.S_w2_pw;
    S_w2_afe = S_w2_pw;
    % Aspect Ratio
    AR_w2 = Geo_tier.AR_w2;
    AR_w2_e = Geo_tier.AR_w2_e;
    % Taper ratio
    lambda_w2 = Geo_tier.lambda_w2;
    lambda_w2_e = Geo_tier.lambda_w2_e;
    % Angles
    % Sweep LE
    Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
    Lambda_LE_w2_e = Geo_tier.Lambda_LE_w2_e;
    % Sweep c/4
    Lambda_c4_w2 = Geo_tier.Lambda_c4_w2;
    Lambda_c2_w2 = Geo_tier.Lambda_c2_w2;
    % dihedral
    dihedral_w2 = Geo_tier.dihedral_w2;
    dihedral_w2_e = Geo_tier.dihedral_w2_e;
    
    %% Tail Volume Coefficient
    Vbar_vee = ((x_xbar_w2 - x_XCG)/cmac_w1)*(S_w2_s/S_ref);
    
    %% Aerodynamic properties
    % Condiciones de vuelo
    % Lift coefficient
    CL0_w2 = Aero.CL_0_w2_CR;
    CL0_w2_e = CL0_w2;
    CL_alpha_wb_Vee = Aero.CL_alpha_w2_CR;
    CLalpha_w2_e = CL_alpha_wb_Vee;
    alpha_CL_0_w2 = Aero.alpha_CL_0_w2_CR;
    % Drag
    CD0_w2 = Aero_TH.CD0_w2;
    % Moment coefficient (CM)
    CM0_w2 = Aero.CM_0_w2_CR;
    CM0_w2_e = CM0_w2;
end

%% Canard
if Can == 1
    % Design Selected incidences
    i_can = Design_criteria.i_can;
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_can = Design_criteria.index_can;
    %% Geometry
    % Wingspan
    b_can = Geo_tier.b_can;
    b_can_e = Geo_tier.b_can_e;
    b_can_pv = Geo_tier.b_can_pv;
    % Chord
    cR_can = Geo_tier.cR_can;
    cT_can = Geo_tier.cT_can;
    cmac_can = Geo_tier.cmac_can;
    cmac_can_e = Geo_tier.cmac_can;
    % Arm from Xac w1 to w2
    l_xac_w1w2 = Geo_tier.l_xac_w1w2;
    % Distances relative to the origin
    % Position (X,Y,Z) of the Xac for both CR and TOLD
    x_xbar_can = Geo_tier.x_xbar_can;
    y_ybar_can = Geo_tier.y_ybar_can;
    z_zbar_can = Geo_tier.z_zbar_can;
    % Distances relative to the LE
    xbar_can = Geo_tier.xbar_can;
    xbar_can_e = Geo_tier.xbar_can; % effective
    % Distances (X,Z) from origin to LE
    x_can_LE = Geo_tier.x_can_LE;
    z_can_LE = Geo_tier.z_can_LE;
    % Distances (X,Y,Z) from origin to LE
    x_cR_can_LE = Geo_tier.x_cR_can_LE;
    y_cR_can_LE = Geo_tier.y_cR_can_LE;
    z_cR_can_LE = Geo_tier.z_cR_can_LE;
    % Distances (X,Y,Z) from origin to TE
    x_cR_can_TE = Geo_tier.x_cR_can_TE;
    y_cR_can_TE = Geo_tier.y_cR_can_TE;
    z_cR_can_TE = Geo_tier.z_cR_can_TE;
    % Areas
    S_can = Geo_tier.S_can;
    S_can_e = Geo_tier.S_can_e; % effective wing area: proyected in the y-plane
    S_can_s = Geo_tier.S_can_s; % real wing area: proyected along the surface of wing
    % prop wash area affecting surfaces
    S_can_pw = Geo_tier.S_can_pw;
    S_can_afe = S_can_pw;
    % Aspect Ratio
    AR_can = Geo_tier.AR_can;
    AR_can_e = Geo_tier.AR_can_e;
    % Taper ratio
    lambda_can = Geo_tier.lambda_can;
    lambda_can_e = Geo_tier.lambda_can_e;
    % Angles
    % Sweep LE
    Lambda_LE_can = Geo_tier.Lambda_LE_can;
    Lambda_LE_can_e = Geo_tier.Lambda_LE_can_e;
    % Sweep c/4
    Lambda_c4_can = Geo_tier.Lambda_c4_can;
    Lambda_c2_can = Geo_tier.Lambda_c2_can;
    % dihedral
    dihedral_can = Geo_tier.dihedral_can;
    dihedral_can_e = Geo_tier.dihedral_can_e;
    
    %% Tail Volume Coefficient
    Vbar_can = ((x_xbar_can + x_XCG)/cmac_w1)*(S_can_s/S_ref);
    
    %% Aerodynamic properties
    % Condiciones de vuelo
    % Lift coefficient
    CL0_can = Aero.CL_0_can_CR;
    CL0_can_e = CL0_can;
    CL_alpha_wb_Vee = Aero.CL_alpha_can_CR;
    CLalpha_can_e = CL_alpha_wb_Vee;
    alpha_CL_0_can = Aero.alpha_CL_0_can_CR;
    % Drag
    CD0_can = Aero_TH.CD0_can;
    % Moment coefficient (CM)
    CM0_can = Aero.CM_0_can_CR;
    CM0_can_e = CM0_can;
end
% 

%% VTP
if VTP == 1
    % Design Selected incidences
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_VTP = Design_criteria.index_VTP;
    %% Geometry
    % Wingspan
    b_VTP = Geo_tier.b_VTP;
    b_VTP_e = Geo_tier.b_VTP_e;
    b_VTP_s = Geo_tier.b_VTP_s;
    % Chord
    cR_VTP = Geo_tier.cR_VTP;
    cT_VTP = Geo_tier.cT_VTP;
    cmac_VTP = Geo_tier.cmac_VTP;
    cmac_VTP_e = Geo_tier.cmac_VTP;
    % Arm from Xac w1 to w2
    l_xac_w1VTP = Geo_tier.l_xac_w1VTP;
    % Distances relative to the origin
    % Position (X,Y,Z) of the Xac for both CR and TOLD
    x_xbar_VTP = Geo_tier.x_xbar_VTP;
    y_ybar_VTP = Geo_tier.y_ybar_VTP;
    z_zbar_VTP = Geo_tier.z_zbar_VTP;
    % Distances relative to the LE
    xbar_VTP = Geo_tier.xbar_VTP;
    xbar_VTP_e = Geo_tier.xbar_VTP;
    % Distances (X,Z) from origin to LE
    x_VTP_LE = Geo_tier.x_VTP_LE;
    z_VTP_LE = Geo_tier.z_VTP_LE;
    % Distances (X,Y,Z) from origin to LE
    x_cR_VTP_LE = Geo_tier.x_cR_w1_LE;
    y_cR_VTP_LE = Geo_tier.y_cR_w1_LE;
    z_cR_VTP_LE = Geo_tier.z_cR_w1_LE;
    % Distances (X,Y,Z) from origin to TE
    x_cR_VTP_TE = Geo_tier.x_cR_VTP_TE;
    y_cR_VTP_TE = Geo_tier.y_cR_VTP_TE;
    z_cR_VTP_TE = Geo_tier.z_cR_VTP_TE;
    % Areas
    S_VTP = Geo_tier.S_VTP;
    S_VTP_e = Geo_tier.S_VTP_e; % effective wing area: proyected in the y-plane
    S_VTP_s = Geo_tier.S_VTP_s; % real wing area: proyected along the surface of wing
    % prop wash area affecting surfaces
    S_VTP_pw = Geo_tier.S_VTP_pw;
    S_VTP_afe = S_VTP_pw;
    % Aspect Ratio
    AR_VTP = Geo_tier.AR_VTP;
    AR_VTP_e = Geo_tier.AR_VTP_e;
    % Taper ratio
    lambda_VTP = Geo_tier.lambda_VTP;
    lambda_VTP_e = Geo_tier.lambda_VTP_e;
    % Angles
    % Sweep LE
    Lambda_LE_VTP = Geo_tier.Lambda_LE_VTP;
    Lambda_LE_VTP_e = Geo_tier.Lambda_LE_VTP_e;
    % Sweep c/4
    Lambda_c4_VTP = Geo_tier.Lambda_c4_VTP;
    Lambda_c2_VTP = Geo_tier.Lambda_c2_VTP;
    % dihedral
    dihedral_VTP = Geo_tier.dihedral_VTP;
    dihedral_VTP_e = Geo_tier.dihedral_VTP_e;
    % Tail Volume Coefficient
    V_t = ((x_xbar_VTP - x_XCG)/cmac_w1)*(S_VTP_e/S_ref);
    
    %% Aerodynamic properties
    % Condiciones de vuelo
    CLalpha_VTP = Aero.CL_alpha_VTP_CR;
    CLalpha_VTP_e = CLalpha_VTP;
end

%% Properties of control surfaces
% aileron
if d_ail ==1
    % Calculates geometric information regarding the control surfaces:
    % chord and location of the aerodynamic center (x,y,z)
    cmac_ail = Geo_tier.cmac_ail;
    xbar_ail = Geo_tier.xbar_ail;
    ybar_ail = Geo_tier.ybar_ail;
    zbar_ail = Geo_tier.zbar_ail;
    
    % Distances relative to the origin
    % Store DATA
    x_xbar_ail = Geo_tier.x_xbar_ail;
    y_ybar_ail = Geo_tier.y_ybar_ail;
    z_zbar_ail = Geo_tier.z_zbar_ail;
    x_1R_y1_ail = Geo_tier.x_1R_y1_ail;
    x_2R_y1_ail = Geo_tier.x_2R_y1_ail;
    x_1R_y2_ail = Geo_tier.x_1R_y2_ail;
    x_2R_y2_ail = Geo_tier.x_2R_y2_ail;
    y_1R_y1_ail = Geo_tier.y_1R_y1_ail;
    y_2R_y1_ail = Geo_tier.y_2R_y1_ail;
    y_1R_y2_ail = Geo_tier.y_1R_y2_ail;
    y_2R_y2_ail = Geo_tier.y_2R_y2_ail;
    z_1R_y1_ail = Geo_tier.z_1R_y1_ail;
    z_2R_y1_ail = Geo_tier.z_2R_y1_ail;
    z_1R_y2_ail = Geo_tier.z_1R_y2_ail;
    z_2R_y2_ail = Geo_tier.z_2R_y2_ail;
    cR_ail = Geo_tier.cR_ail;
    cT_ail = Geo_tier.cT_ail;
    b_ail = Geo_tier.b_ail;
    S_ail = Geo_tier.S_ail;
end

% elevator
if d_ele == 1
    % Calculates geometric information regarding the control surfaces:
    % chord and location of the aerodynamic center (x,y,z)
    cmac_ele = Geo_tier.cmac_ele;
    xbar_ele = Geo_tier.xbar_ele;
    ybar_ele = Geo_tier.ybar_ele;
    zbar_ele = Geo_tier.zbar_ele;
    
    % Distances relative to the origin
    % Store DATA
    x_xbar_ele = Geo_tier.x_xbar_ele;
    y_ybar_ele = Geo_tier.y_ybar_ele;
    z_zbar_ele = Geo_tier.z_zbar_ele;
    x_1R_y1_ele = Geo_tier.x_1R_y1_ele;
    x_2R_y1_ele = Geo_tier.x_2R_y1_ele;
    x_1R_y2_ele = Geo_tier.x_1R_y2_ele;
    x_2R_y2_ele = Geo_tier.x_2R_y2_ele;
    y_1R_y1_ele = Geo_tier.y_1R_y1_ele;
    y_2R_y1_ele = Geo_tier.y_2R_y1_ele;
    y_1R_y2_ele = Geo_tier.y_1R_y2_ele;
    y_2R_y2_ele = Geo_tier.y_2R_y2_ele;
    z_1R_y1_ele = Geo_tier.z_1R_y1_ele;
    z_2R_y1_ele = Geo_tier.z_2R_y1_ele;
    z_1R_y2_ele = Geo_tier.z_1R_y2_ele;
    z_2R_y2_ele = Geo_tier.z_2R_y2_ele;
    cR_ele = Geo_tier.cR_ele;
    cT_ele = Geo_tier.cT_ele;
    b_ele = Geo_tier.b_ele;
    S_ele = Geo_tier.S_ele;
end

% elevon
if d_elevon == 1
    % Calculates geometric information regarding the control surfaces:
    % chord and location of the aerodynamic center (x,y,z)
    cmac_elevon = Geo_tier.cmac_elevon;
    xbar_elevon = Geo_tier.xbar_elevon;
    ybar_elevon = Geo_tier.ybar_elevon;
    zbar_elevon = Geo_tier.zbar_elevon;
    
    % Distances relative to the origin
    % Store DATA
    x_xbar_elevon = Geo_tier.x_xbar_elevon;
    y_ybar_elevon = Geo_tier.y_ybar_elevon;
    z_zbar_elevon = Geo_tier.z_zbar_elevon;
    x_1R_y1_elevon = Geo_tier.x_1R_y1_elevon;
    x_2R_y1_elevon = Geo_tier.x_2R_y1_elevon;
    x_1R_y2_elevon = Geo_tier.x_1R_y2_elevon;
    x_2R_y2_elevon = Geo_tier.x_2R_y2_elevon;
    y_1R_y1_elevon = Geo_tier.y_1R_y1_elevon;
    y_2R_y1_elevon = Geo_tier.y_2R_y1_elevon;
    y_1R_y2_elevon = Geo_tier.y_1R_y2_elevon;
    y_2R_y2_elevon = Geo_tier.y_2R_y2_elevon;
    z_1R_y1_elevon = Geo_tier.z_1R_y1_elevon;
    z_2R_y1_elevon = Geo_tier.z_2R_y1_elevon;
    z_1R_y2_elevon = Geo_tier.z_1R_y2_elevon;
    z_2R_y2_elevon = Geo_tier.z_2R_y2_elevon;
    cR_elevon = Geo_tier.cR_elevon;
    cT_elevon = Geo_tier.cT_elevon;
    b_elevon = Geo_tier.b_elevon;
    S_elevon = Geo_tier.S_elevon;
end

% flap
if d_flap == 1
    % Calculates geometric information regarding the control surfaces:
    % chord and location of the aerodynamic center (x,y,z)
    cmac_flap = Geo_tier.cmac_flap;
    xbar_flap = Geo_tier.xbar_flap;
    ybar_flap = Geo_tier.ybar_flap;
    zbar_flap = Geo_tier.zbar_flap;
    
    % Distances relative to the origin
    % Store DATA
    x_xbar_flap = Geo_tier.x_xbar_flap;
    y_ybar_flap = Geo_tier.y_ybar_flap;
    z_zbar_flap = Geo_tier.z_zbar_flap;
    x_1R_y1_flap = Geo_tier.x_1R_y1_flap;
    x_2R_y1_flap = Geo_tier.x_2R_y1_flap;
    x_1R_y2_flap = Geo_tier.x_1R_y2_flap;
    x_2R_y2_flap = Geo_tier.x_2R_y2_flap;
    y_1R_y1_flap = Geo_tier.y_1R_y1_flap;
    y_2R_y1_flap = Geo_tier.y_2R_y1_flap;
    y_1R_y2_flap = Geo_tier.y_1R_y2_flap;
    y_2R_y2_flap = Geo_tier.y_2R_y2_flap;
    z_1R_y1_flap = Geo_tier.z_1R_y1_flap;
    z_2R_y1_flap = Geo_tier.z_2R_y1_flap;
    z_1R_y2_flap = Geo_tier.z_1R_y2_flap;
    z_2R_y2_flap = Geo_tier.z_2R_y2_flap;
    cR_flap = Geo_tier.cR_flap;
    cT_flap = Geo_tier.cT_flap;
    b_flap = Geo_tier.b_flap;
    S_flap = Geo_tier.S_flap;
end

% rudder
if d_rudder == 1
    % Calculates geometric information regarding the control surfaces:
    % chord and location of the aerodynamic center (x,y,z)
    cmac_rudder = Geo_tier.cmac_rudder;
    xbar_rudder = Geo_tier.xbar_rudder;
    ybar_rudder = Geo_tier.ybar_rudder;
    zbar_rudder = Geo_tier.zbar_rudder;
    
    % Distances relative to the origin
    % Store DATA
    x_xbar_rudder = Geo_tier.x_xbar_rudder;
    y_ybar_rudder = Geo_tier.y_ybar_rudder;
    z_zbar_rudder = Geo_tier.z_zbar_rudder;
    x_1R_y1_rudder = Geo_tier.x_1R_y1_rudder;
    x_2R_y1_rudder = Geo_tier.x_2R_y1_rudder;
    x_1R_y2_rudder = Geo_tier.x_1R_y2_rudder;
    x_2R_y2_rudder = Geo_tier.x_2R_y2_rudder;
    y_1R_y1_rudder = Geo_tier.y_1R_y1_rudder;
    y_2R_y1_rudder = Geo_tier.y_2R_y1_rudder;
    y_1R_y2_rudder = Geo_tier.y_1R_y2_rudder;
    y_2R_y2_rudder = Geo_tier.y_2R_y2_rudder;
    z_1R_y1_rudder = Geo_tier.z_1R_y1_rudder;
    z_2R_y1_rudder = Geo_tier.z_2R_y1_rudder;
    z_1R_y2_rudder = Geo_tier.z_1R_y2_rudder;
    z_2R_y2_rudder = Geo_tier.z_2R_y2_rudder;
    cR_rudder = Geo_tier.cR_rudder;
    cT_rudder = Geo_tier.cT_rudder;
    b_rudder = Geo_tier.b_rudder;
    S_rudder = Geo_tier.S_rudder;
    
end

% ruddervator
if d_rudvtr == 1
    % Calculates geometric information regarding the control surfaces:
    % chord and location of the aerodynamic center (x,y,z)
    cmac_rudvtr = Geo_tier.cmac_rudvtr;
    xbar_rudvtr = Geo_tier.xbar_rudvtr;
    ybar_rudvtr = Geo_tier.ybar_rudvtr;
    zbar_rudvtr = Geo_tier.zbar_rudvtr;
    
    % Distances relative to the origin
    % Store DATA
    x_xbar_rudvtr = Geo_tier.x_xbar_rudvtr;
    y_ybar_rudvtr = Geo_tier.y_ybar_rudvtr;
    z_zbar_rudvtr = Geo_tier.z_zbar_rudvtr;
    x_1R_y1_rudvtr = Geo_tier.x_1R_y1_rudvtr;
    x_2R_y1_rudvtr = Geo_tier.x_2R_y1_rudvtr;
    x_1R_y2_rudvtr = Geo_tier.x_1R_y2_rudvtr;
    x_2R_y2_rudvtr = Geo_tier.x_2R_y2_rudvtr;
    y_1R_y1_rudvtr = Geo_tier.y_1R_y1_rudvtr;
    y_2R_y1_rudvtr = Geo_tier.y_2R_y1_rudvtr;
    y_1R_y2_rudvtr = Geo_tier.y_1R_y2_rudvtr;
    y_2R_y2_rudvtr = Geo_tier.y_2R_y2_rudvtr;
    z_1R_y1_rudvtr = Geo_tier.z_1R_y1_rudvtr;
    z_2R_y1_rudvtr = Geo_tier.z_2R_y1_rudvtr;
    z_1R_y2_rudvtr = Geo_tier.z_1R_y2_rudvtr;
    z_2R_y2_rudvtr = Geo_tier.z_2R_y2_rudvtr;
    cR_rudvtr = Geo_tier.cR_rudvtr;
    cT_rudvtr = Geo_tier.cT_rudvtr;
    b_rudvtr = Geo_tier.b_rudvtr;
    S_rudvtr = Geo_tier.S_rudvtr;
    
end

% cannard
if d_can == 1
    % Calculates geometric information regarding the control surfaces:
    % chord and location of the aerodynamic center (x,y,z)
    cmac_can = Geo_tier.cmac_can;
    xbar_can = Geo_tier.xbar_can;
    ybar_can = Geo_tier.ybar_can;
    zbar_can = Geo_tier.zbar_can;
    
    % Distances relative to the origin
    % Store DATA
    x_xbar_can = Geo_tier.x_xbar_can;
    y_ybar_can = Geo_tier.y_ybar_can;
    z_zbar_can = Geo_tier.z_zbar_can;
    x_1R_y1_can = Geo_tier.x_1R_y1_can;
    x_2R_y1_can = Geo_tier.x_2R_y1_can;
    x_1R_y2_can = Geo_tier.x_1R_y2_can;
    x_2R_y2_can = Geo_tier.x_2R_y2_can;
    y_1R_y1_can = Geo_tier.y_1R_y1_can;
    y_2R_y1_can = Geo_tier.y_2R_y1_can;
    y_1R_y2_can = Geo_tier.y_1R_y2_can;
    y_2R_y2_can = Geo_tier.y_2R_y2_can;
    z_1R_y1_can = Geo_tier.z_1R_y1_can;
    z_2R_y1_can = Geo_tier.z_2R_y1_can;
    z_1R_y2_can = Geo_tier.z_1R_y2_can;
    z_2R_y2_can = Geo_tier.z_2R_y2_can;
    cR_can = Geo_tier.cR_can;
    cT_can = Geo_tier.cT_can;
    b_can = Geo_tier.b_can;
    S_can = Geo_tier.S_can;
    
end

% Propeller/Engine distances
z_d_T = Geo_tier.z_d_T;
x_d_T = Geo_tier.x_d_T; % Positive for engine behind Xcg : x_d_T = x_eng_xbar - x_XCG; 
y_d_T = Geo_tier.y_d_T;

% Thrust line inclination angle (for futuire version with tilting rotors)
phi_T = 0;
% The perpendicular distance between the thrustline and the airplane center of gravity is given by:
d_T = x_d_T*sin(phi_T) + z_d_T*cos(phi_T);

% if HTP == 1
%     %% Tail Volume Coefficient
%     V_h = ((x_xbar_w2 - x_XCG)/cmac_w1)*(S_w2_e/S_ref);
% end
% if VTP == 1
%     %% Tail Volume Coefficient
%     V_t = ((x_xbar_VTP - x_XCG)/cmac_w1)*(S_VTP_e/S_ref);
% end
% 
% if Vee == 1
%     %% Tail Volume Coefficient
%     V_vee = ((x_xbar_w2 - x_XCG)/cmac_w1)*(S_w2_s/S_ref);
% end

%% Aerodynamic properties
% Condiciones de vuelo
C_D0 = Aero_TH.CD0;
C_D1 = Aero_TH.CD1;
C_D2 = Aero_TH.CD2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CONFIGURACI�N PUSHER. MODELO DE HELIC�PTERO%%%%%%%%%%%%%%%
D_prop = Prop_data.D_prop;
R_prop = D_prop/2;
S_heli = (pi*R_prop^2);              %superficie de la h�lice

% Estimation of Desired Thrust
CL = w_T0/(q_inf*S_ref); % assume equilibry steady state flight
CD = C_D0 + C_D1*CL + C_D2*CL^2;
D = CD*q_inf*S_ref;

%% Propulsive Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%DERIVADA PROPULSIVA LONGITUDINAL%%%%%%%%%%%%%%%%
CD_Total_prel = C_D0 + C_D1*CL + C_D2*CL^2;  %polar aeronave
Drag = q_inf*S_ref*CD_Total_prel;
Fdes = Drag; % assume T = D

%% Propulsion information
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine
% 4: CONSUMO ESPECIFICO: % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
% 5: DERIVACION(TURBOFANES) By-pass
% 6: EFICIENCIA DE LA HELICE (ETA_P)
% 7: AVION CIVIL = 1/MILITAR = 2 - Normativa
% 8: CAPACIDAD CALORIFICA DEL COMBUSTIBLE
% 9: DIAMETRO DEL MOTOR(TURBOHELICES)
propul = AC_CONFIGURATION.propulsion;
type_engine = propul(1);
% Type of prop: fixed pitch or variable pitch
type_prop = AC_CONFIGURATION.type_prop; 
bypass_ratio = propul(5); % By-pass Ratio

% Prop Wash Effect
prop_wash_effect = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;
%% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
switch type_engine
    case 1 % TIPO DE MOTOR --> 1_TURBOFAN 
%         n = Selection_J_CT_F_des(V,Prop_data,Drag,rho);
%         [Propulsion] = get_EngineProperties_v3(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,n_eng)
        [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION);
%         d_CT_d_J = Propulsion.d_CT_d_J;
        d_CT_d_V = Propulsion.d_CT_d_V;
%         d_CP_d_J = Propulsion.d_CP_d_J;
%         d_CP_d_V = Propulsion.d_CP_d_V;
%         d_T_d_J = Propulsion.d_T_d_J;
        d_T_d_V = Propulsion.d_T_d_V;
%         d_P_d_J = Propulsion.d_P_d_J;
%         d_P_d_V = Propulsion.d_P_d_V;
%         n = Propulsion.RPS;
%         d_etap_d_V = Propulsion.d_etap_d_V;
%         d_etap_d_J = Propulsion.d_etap_d_J;
%         J = Propulsion.J;
%         etha_emp = Propulsion.etha_emp;
        T_tot = Propulsion.Ti;
        T_eng = Propulsion.Ti_eng;
        %% Estimation of data for induced velocity associated to engine configuration
        if prop_wash_effect == 1
            v_i = Propulsion.v_i; % Induced velocity at prop disk
        else
            v_i = 0;
        end
        R_inf = Propulsion.R_inf; % Radius of prop wash at infinity
        v_inf = Propulsion.v_inf; % induced velocity at propwash at infinity
        % v_i = -(1/2)*V + sqrt(1/4*V^2 + T_eng/(2*rho*S_heli));
        % R_inf = R_heli*sqrt((V + v_i)/(V + 2*v_i));
        % Propeller variables
    case 2 % TIPO DE MOTOR --> 2_TURBOPROP
        [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION);
%         n = Selection_J_CT_F_des(V,Prop_data,Drag,rho);
%         [Propulsion] = get_EngineProperties_v3(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,n_eng);
%         d_CT_d_J = Propulsion.d_CT_d_J;
        d_CT_d_V = Propulsion.d_CT_d_V;
%         d_CP_d_J = Propulsion.d_CP_d_J;
        d_CP_d_V = Propulsion.d_CP_d_V;
%         d_T_d_J = Propulsion.d_Ti_d_J;
        d_T_d_V = Propulsion.d_T_d_V;
%         d_P_d_J = Propulsion.d_Pi_d_J;
        d_P_d_V = Propulsion.d_P_d_V;
        
%         d_etap_d_V = Propulsion.d_etap_d_V;
%         d_etap_d_J = Propulsion.d_etap_d_J;
%         J = Propulsion.J;
%         etha_emp = Propulsion.etha_emp;
        
        T_tot = Propulsion.Ti;
        T_eng = Propulsion.Ti_eng;
        %% Estimation of data for induced velocity associated to engine configuration
        if prop_wash_effect == 1
            v_i = Propulsion.v_i; % Induced velocity at prop disk
        else
            v_i = 0;
        end
        R_inf = Propulsion.R_inf; % Radius of prop wash at infinity
        v_inf = Propulsion.v_inf; % induced velocity at propwash at infinity
        % v_i = -(1/2)*V + sqrt(1/4*V^2 + T_eng/(2*rho*S_heli));
        % R_inf = R_heli*sqrt((V + v_i)/(V + 2*v_i));
        % Propeller variables
    case 3 % TIPO DE MOTOR --> 3_PISTON 
%          n = Selection_J_CT_F_des(V,Prop_data,Drag,rho);
%         [Propulsion] = get_EngineProperties_v3(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,n_eng);
        [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION);
        d_CT_d_J = Propulsion.d_CT_d_J;
        d_CT_d_V = Propulsion.d_CT_d_V;
        d_CP_d_J = Propulsion.d_CP_d_J;
        d_CP_d_V = Propulsion.d_CP_d_V;
        d_T_d_J = Propulsion.d_T_d_J;
        d_T_d_V = Propulsion.d_T_d_V;
        d_P_d_J = Propulsion.d_P_d_J;
        d_P_d_V = Propulsion.d_P_d_V;
        
        d_etap_d_V = Propulsion.d_etap_d_V;
        d_etap_d_J = Propulsion.d_etap_d_J;
        J = Propulsion.J;
        etha_emp = Propulsion.etha_emp;

        T_tot = Propulsion.Ti;
        T_eng = Propulsion.Ti_eng;
        %% Estimation of data for induced velocity associated to engine configuration
        if prop_wash_effect == 1
            v_i = Propulsion.v_i; % Induced velocity at prop disk
        else
            v_i = 0;
        end
        R_inf = Propulsion.R_inf; % Radius of prop wash at infinity
        v_inf = Propulsion.v_inf; % induced velocity at propwash at infinity
        % v_i = -(1/2)*V + sqrt(1/4*V^2 + T_eng/(2*rho*S_heli));
        % R_inf = R_heli*sqrt((V + v_i)/(V + 2*v_i));
        % Propeller variables
    case 4 % TIPO DE MOTOR --> 4_ELECTRIC_PROP
%         n = Selection_J_CT_F_des(V,Prop_data,Drag,rho);
%         [Propulsion] = get_EngineProperties_v3(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,n_eng)
        [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION);
        d_CT_d_J = Propulsion.d_CT_d_J;
        d_CT_d_V = Propulsion.d_CT_d_V;
        d_CP_d_J = Propulsion.d_CP_d_J;
        d_CP_d_V = Propulsion.d_CP_d_V;
        d_T_d_J = Propulsion.d_T_d_J;
        d_T_d_V = Propulsion.d_T_d_V;
        d_P_d_J = Propulsion.d_P_d_J;
        d_P_d_V = Propulsion.d_P_d_V;
        n = Propulsion.RPS;
        d_etap_d_V = Propulsion.d_etap_d_V;
        d_etap_d_J = Propulsion.d_etap_d_J;
        J = Propulsion.J;
        etha_emp = Propulsion.etha_emp;
        T_tot = Propulsion.Ti;
        T_eng = Propulsion.Ti_eng;
        
        %% Estimation of data for induced velocity associated to engine configuration
        if prop_wash_effect == 1
            v_i = Propulsion.v_i; % Induced velocity at prop disk
        else
            v_i = 0;
        end
        R_inf = Propulsion.R_inf; % Radius of prop wash at infinity
        v_inf = Propulsion.v_inf; % induced velocity at propwash at infinity
        % v_i = -(1/2)*V + sqrt(1/4*V^2 + T_eng/(2*rho*S_heli));
        % R_inf = R_heli*sqrt((V + v_i)/(V + 2*v_i));
        % Propeller variables
end

% v_i = 0;

if Vee == 1
    %%%%%%%%%%%%%%%%%%%%%%C�LCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kh = (1-(z_w2_LE - z_w1_LE)/b_w1)/(2*l_xac_w1w2/b_w1)^(1/3);          %eq 3.43 PAMADI
    kl = (10-3*lambda_w1)/7;                                     %eq 3.43 PAMADI
    ka = 1/AR_w1_e -1/(1 + AR_w1_e^(1.7));                        %eq 3.43 PAMADI
    deps_dalpha_clean = 4.44*(kh*kl*ka*sqrt(cos(Lambda_c4_w1)))^1.19;  %eq 3.43 PAMADI
    deps_dalpha_poff = 0; % Assumes not affecting
    deps_dalpha_Vee = deps_dalpha_clean + deps_dalpha_poff;
    
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash = 1 - deps_dalpha_Vee;
    % eps_vee_0 = deps_dalpha*(alpha_CL_0_w1*D2R - i_w1);
    eps_w2_0 = deps_dalpha_Vee*(alpha_CL_0_w1*D2R - i_w1);
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_w2 = eps_w2_0 + deps_dalpha_Vee*alpha_f;
    
    %% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection
    
    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_w2 - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to w2 MAC
    l = x_xbar_w2 - x_cR_w1_TE; % x distance from the trailing edge of the wing to the w2 MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_w1 + atan(z_TE/l);
    gamma_w2 = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_w1 - eps_vee);
    x_w2_wake = a*cos(gamma_w2 - alpha_f - i_w1 - eps_w2);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    Deltaz_w2_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_w2_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    z_w2_wake = a*sin(gamma_w2 - alpha_f - i_w1 + eps_w2);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_w2_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_w2_wake)/(2*Deltaz_w2_wake)))^2)/((x_w2_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_w2_power = 0;
    % NOTE- need to take into account power effects
    % eta_vee = eta_vee_power_off + eta_vee_power
    eta_w2 = eta_w2_power_off + eta_w2_power;
else
    deps_dalpha_Vee = 0;
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash = 1 - deps_dalpha_Vee;
    eps_w2_0 = 0;
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_w2 = eps_w2_0 + deps_dalpha_Vee*alpha_f;
end

if HTP == 1
    %%%%%%%%%%%%%%%%%%%%%%C�LCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kh = (1-(z_w2_LE - z_w1_LE)/b_w1)/(2*l_xac_w1w2/b_w1)^(1/3);          %eq 3.43 PAMADI
    kl = (10-3*lambda_w1)/7;                                     %eq 3.43 PAMADI
    ka = 1/AR_w1_e -1/(1 + AR_w1_e^(1.7));                        %eq 3.43 PAMADI
    deps_dalpha_clean = 4.44*(kh*kl*ka*sqrt(cos(Lambda_c4_w1)))^1.19;  %eq 3.43 PAMADI
    deps_dalpha_poff = 0; % Assumes not affecting
    deps_dalpha_h = deps_dalpha_clean + deps_dalpha_poff;
    
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash = 1 - deps_dalpha_h;
    % eps_vee_0 = deps_dalpha*(alpha_CL_0_w1*D2R - i_w1);
    eps_w2_0 = deps_dalpha_h*(alpha_CL_0_w1*D2R - i_w1);
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_w2 = eps_w2_0 + deps_dalpha_h*alpha_f;
    
    %% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection
    
    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_w2 - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to w2 MAC
    l = x_xbar_w2 - x_cR_w1_TE; % x distance from the trailing edge of the wing to the w2 MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_w1 + atan(z_TE/l);
    gamma_w2 = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_w1 - eps_vee);
    x_w2_wake = a*cos(gamma_w2 - alpha_f - i_w1 - eps_w2);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    Deltaz_w2_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_w2_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    z_w2_wake = a*sin(gamma_w2 - alpha_f - i_w1 + eps_w2);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_w2_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_w2_wake)/(2*Deltaz_w2_wake)))^2)/((x_w2_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_w2_power = 0;
    
    % NOTE- need to take into account power effects
    % eta_vee = eta_vee_power_off + eta_vee_power
    eta_w2 = eta_w2_power_off + eta_w2_power;
else
    deps_dalpha_h = 0;
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash = 1 - deps_dalpha_h;
    eps_w2_0 = 0;
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_w2 = eps_w2_0 + deps_dalpha_h*alpha_f;
end

if Can == 1
    %%%%%%%%%%%%%%%%%%%%%%C�LCULO DEL UPWASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deps_dalpha_clean = upwash_calc(AR, X_w, X_c, c_m);
    deps_dalpha_poff = 0; % Assumes not affecting
    deps_dalpha_can = deps_dalpha_clean + deps_dalpha_poff;
    
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    upwash = 1 + deps_dalpha_can;
    eps_can_0 = deps_dalpha_can*(alpha_CL_0_w1*D2R - i_w1);
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_can = eps_can_0 + deps_dalpha_can*alpha_f;
    
    %% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection
    
    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_can - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to w2 MAC
    l = x_xbar_can - x_cR_w1_TE; % x distance from the trailing edge of the wing to the w2 MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the Canard aerodynamic center
    gamma_can = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    x_can_wake = a*cos(gamma_can - alpha_f - i_w1 - eps_can);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    Deltaz_can_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_can_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    z_can_wake = a*sin(gamma_can - alpha_f - i_w1 + eps_can);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_can_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_can_wake)/(2*Deltaz_can_wake)))^2)/((x_can_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_can_power = 0;
    
    % NOTE- need to take into account power effects
    eta_can = eta_can_power_off + eta_can_power;
else
    deps_dalpha_can = 0;
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    upwash = 1 + deps_dalpha_can;
    eps_can_0 = 0;
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_can = eps_can_0 + deps_dalpha_can*alpha_f;
end

if VTP == 1
    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_VTP - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to w2 MAC
    l = x_xbar_VTP - x_cR_w1_TE; % x distance from the trailing edge of the wing to the w2 MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_w1 + atan(z_TE/l);
    gamma_VTP = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_w1 - eps_vee);
    x_VTP_wake = acos(gamma_VTP - alpha_f - i_w1 - eps_w2);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    Deltaz_VTP_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_VTP_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    z_VTP_wake = a*sin(gamma_VTP - alpha_f - i_w1 + eps_w2);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_VTP_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_VTP_wake)/(2*Deltaz_VTP_wake)))^2)/((x_VTP_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_VTP_power = 0;
    
    % NOTE- need to take into account power effects
    eta_VTP = eta_VTP_power_off + eta_VTP_power;
end

% Upwash influencing in prop
% Location of the prop disk (source of thrust)
x_prop_cF = Geo_tier.x_prop_cF; % location of prop
AR_w1 = Geo_tier.AR_w1;
depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);

if W1 == 1
    S_w1_afe = S_w1_pw;
    S_w1_no_afe = S_w1 - S_w1_afe;
    q_w1_no_afe = q_inf;
    if Posicion_Palanca == 0
        V_w1 = V;
        q_w1_afe    = q_inf;
    else
        V_w1 = V + v_i;
        q_w1_afe       = 0.5*rho*V_w1^2;
    end
    eta_w1_afe = (q_w1_afe/q_inf);
    eta_w1_no_afe = (q_w1_no_afe/q_inf);
    eta_w1_afe_S_w1_afe_S_ref = eta_w1_afe*(S_w1_afe/S_ref);
    eta_w1_no_afe_S_w1_no_afe_S_ref = eta_w1_no_afe*(S_w1_no_afe/S_ref);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Correction of Lift curve slope of w1 associated to prop wash
    % eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
    % eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
    CLalpha_w1_e_pw = eta_w1_afe_S_w1_afe_S_ref*CLalpha_w1_e + eta_w1_no_afe_S_w1_no_afe_S_ref*CLalpha_w1_e;
    
    CL0_w1_e_corrected = eta_w1_afe_S_w1_afe_S_ref*CL0_w1_e + eta_w1_no_afe_S_w1_no_afe_S_ref*CL0_w1_e;
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CLalpha_w1 = CLalpha_w1_e_pw;
end

if Can == 1
    S_can_afe = S_can_pw;
    S_can_no_afe = S_can - S_can_afe;
    q_can_no_afe = q_inf;
    if Posicion_Palanca == 0
        V_can = V;
        q_can_afe    = q_inf;
    else
        V_can = V + v_i;
        q_can_afe       = 0.5*rho*V_can^2;
    end
    
    eta_can_afe = (q_can_afe/q_inf);
    eta_can_no_afe = (q_w1_no_afe/q_inf);
    eta_can_afe_S_can_afe_S_ref = eta_can_afe*(S_can_afe/S_ref);
    eta_can_no_afe_S_can_no_afe_S_ref = eta_can_no_afe*(S_can_no_afe/S_ref);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Correction of Lift curve slope of w1 associated to prop wash
    % eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
    % eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
    CLalpha_can_e_pw = eta_can_afe_S_can_afe_S_ref*CLalpha_can_e + eta_can_no_afe_S_can_no_afe_S_ref*CLalpha_can_e;
    CL0_can_e_corrected = eta_can_afe_S_can_afe_S_ref*CL0_can_e + eta_can_no_afe_S_can_no_afe_S_ref*CL0_can_e;
    Stab_Der_parts.CL0_can = CL0_can_e_corrected;
    Stab_Der_parts.CLalpha_can = CLalpha_can_e_pw;
end

if HTP == 1
    S_w2_afe = S_w2_pw;
    S_w2_no_afe = S_w2_s - S_w2_afe;
    q_w2_no_afe    = eta_w2*q_inf;
    %presion dinamica no afectada
    
    if Posicion_Palanca == 0
        V_w2 = V;
        q_w2_afe    = q_inf;
    else
        V_w2 = sqrt(eta_w2)*V + 2*v_i;
        q_w2_afe       = 0.5*rho*V_w2^2;
    end
    
    if Posicion_Palanca == 0
        V_w1 = V;
        q_w1_afe    = q_inf;
    else
        V_w1 = V + v_i;
        q_w1_afe       = 0.5*rho*V_w1^2;
    end
   
    eta_w2_afe = (q_w2_afe/q_inf);
    eta_w2_no_afe = (q_w2_no_afe/q_inf);
    eta_w2_afe_S_w2_afe_S_ref = eta_w2_afe*(S_w2_afe/S_ref);
    eta_w2_no_afe_S_w2_no_afe_S_ref = eta_w2_no_afe*(S_w2_no_afe/S_ref);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Correction of Lift curve slope of w1 associated to prop wash
    % eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
    % eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
    CLalpha_HTP_e_pw = eta_w2_afe_S_w2_afe_S_ref*CLalpha_w2_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CLalpha_w2_e;
    CL0_HTP_e_corrected = eta_w2_afe_S_w2_afe_S_ref*CL0_w2_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CL0_w2_e;
    Stab_Der_parts.CL0_HTP = CL0_HTP_e_corrected;
    Stab_Der_parts.CLalpha_HTP = CLalpha_HTP_e_pw;
end

if Vee == 1
    S_w2_afe = S_w2_pw;
    S_w2_no_afe = S_w2_s - S_w2_afe;
    q_w2_no_afe    = eta_w2*q_inf;
    %presion dinamica no afectada
    
    if Posicion_Palanca == 0
        V_w2 = V;
        q_w2_afe    = q_inf;
    else
        V_w2 = sqrt(eta_w2)*V + 2*v_i;
        q_w2_afe       = 0.5*rho*V_w2^2;
    end
    
    eta_w2_afe = (q_w2_afe/q_inf);
    eta_w2_no_afe = (q_w2_no_afe/q_inf);
    eta_w2_afe_S_w2_afe_S_ref = eta_w2_afe*(S_w2_afe/S_ref);
    eta_w2_no_afe_S_w2_no_afe_S_ref = eta_w2_no_afe*(S_w2_no_afe/S_ref);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Correction of Lift curve slope of w1 associated to prop wash
    % eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
    % eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
    CLalpha_Vee_e_pw = eta_w2_afe_S_w2_afe_S_ref*CLalpha_w2_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CLalpha_w2_e;    
    CL0_Vee_e_corrected = eta_w2_afe_S_w2_afe_S_ref*CL0_w2_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CL0_w2_e;
    Stab_Der_parts.CL0_Vee = CL0_Vee_e_corrected;
    Stab_Der_parts.CLalpha_Vee = CLalpha_Vee_e_pw;
end
if VTP == 1
    if twin_VTP == 1
        S_VTP_afe = 2*S_VTP_pw;
        S_VTP_no_afe = S_VTP_s - S_VTP_afe;
        q_VTP_no_afe    = eta_VTP*q_inf;                    %presion dinamica no afectada
        if Posicion_Palanca == 0
            V_VTP = V;
            q_VTP_afe    = q_inf;
        else
            V_VTP = sqrt(eta_VTP)*V + 2*v_i;
            q_VTP_afe       = 0.5*rho*V_VTP^2;
        end
        eta_VTP_afe = (q_VTP_afe/q_inf);
        eta_VTP_no_afe = (q_VTP_no_afe/q_inf);
        eta_VTP_afe_S_VTP_afe_S_ref = eta_VTP_afe*(S_VTP_afe/S_ref);
        eta_VTP_no_afe_S_VTP_no_afe_S_ref = eta_VTP_no_afe*(S_VTP_no_afe/S_ref);
    else
        S_VTP_afe = S_VTP_pw;
        S_VTP_no_afe = S_VTP_s - S_VTP_afe;
        q_VTP_no_afe    = eta_VTP*q_inf;                    %presion dinamica no afectada
        if Posicion_Palanca == 0
            V_VTP = V;
            q_VTP_afe    = q_inf;
        else
            V_VTP = sqrt(eta_VTP)*V + 2*v_i;
            q_VTP_afe       = 0.5*rho*V_VTP^2;
        end
        eta_VTP_afe = (q_VTP_afe/q_inf);
        eta_VTP_no_afe = (q_VTP_no_afe/q_inf);
        eta_VTP_afe_S_VTP_afe_S_ref = eta_VTP_afe*(S_VTP_afe/S_ref);
        eta_VTP_no_afe_S_VTP_no_afe_S_ref = eta_VTP_no_afe*(S_VTP_no_afe/S_ref);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Correction of Lift curve slope of w1 associated to prop wash
    CLalpha_VTP_e_pw = eta_VTP_afe_S_VTP_afe_S_ref*CLalpha_VTP_e + eta_VTP_no_afe_S_VTP_no_afe_S_ref*CLalpha_VTP_e;
    CL0_VTP_e_corrected = eta_VTP_afe_S_VTP_afe_S_ref*CL0_w2_e + eta_VTP_no_afe_S_VTP_no_afe_S_ref*CL0_w2_e;
    Stab_Der_parts.CL0_VTP = CL0_VTP_e_corrected;
    Stab_Der_parts.CLalpha_VTP = CLalpha_VTP_e_pw;
end

%% Correction of Lift curve slope of w1 associated to fuselage interaction
% Pendiente de sustentaci�n del morro.
Fineness_Ratio = length_fus/w_Area_b_max;
%digitaliazacion figura PAMADI CAP3
x_f_k2_k1 = [4.,5.,6.,8.,10.,12.,14.,16.,18.,20.];
y_f_k2_k1 = [.77,.825,.865,.91,.94,.955,.965,.97,.973,.975];
f_k2_k1  = interp1(x_f_k2_k1,y_f_k2_k1,Fineness_Ratio,'spline');
x1 = x_Area_b_max;

if W1 == 1
    %% W1 body interference
    aN_w1 = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
        
%     f_k2_k1         = k2_k1_calc(Body_Geo.l_fus, Geo_tier.d_fus);      % Fuselage apparent mass coefficient. Pamadi, Figure 3.6
%     Body_Geo.CLa_fus    = 2*�*(Body_Geo.S_front/Sref);         % Pendiente de sustentaci�n del morro aislado

    aN_w1 = CLa_fus; %  corrigiendo con la integrtaci�n del fuselaje a partir de modelo XFLR5
    KN_w1 = (aN_w1/CLalpha_w1_e_pw)*(S_w1/S_w1_e);  %eq 3.25 PAMADI
    CLalpha_fus = aN_w1;
    KWB_w1 = 0.1714*(w_Area_b_max/b_w1)^2 + 0.8326*(w_Area_b_max/b_w1) + 0.9974; %eq 3.27 PAMADI
    KBW_w1 = 0.7810*(w_Area_b_max/b_w1)^2 + 1.1976*(w_Area_b_max/b_w1) + 0.0088; %eq 3.28 PAMADI
    a_bw_w1 = KBW_w1*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.34 PAMADI
    a_wb_w1 = KWB_w1*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.33 PAMADI
    a_WB_w1 = (KN_w1 + KWB_w1 + KBW_w1)*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.24 PAMADI
    CLalpha_WB_w1 = a_WB_w1; % contribution of w1 to fuselage and fuselage to w1
    CL_alpha_wb_w1 = CLalpha_WB_w1;
    CL_alpha_w1_corrected = CL_alpha_wb_w1;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_w1 = CL_alpha_wb_w1;
end

%% Canard body interference
if Can == 1
    aN_can = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    aN_can = CLa_fus; %  corrigiendo con la integrtaci�n del fuselaje a partir de modelo XFLR5
    KN_can = (aN_Can/CLalpha_can_e_pw)*(S_can/S_can_e);  %eq 3.25 PAMADI
    CLalpha_fus = aN_can;
    KWB_can = 0.1714*(w_Area_b_max/b_can)^2 + 0.8326*(w_Area_b_max/b_can) + 0.9974; %eq 3.27 PAMADI
    KBW_can = 0.7810*(w_Area_b_max/b_can)^2 + 1.1976*(w_Area_b_max/b_can) + 0.0088; %eq 3.28 PAMADI
    a_bw_can = KBW_can*CLalpha_can_e_pw*(S_can_e/S_can); %eq 3.34 PAMADI
    a_wb_can = KWB_can*CLalpha_can_e_pw*(S_can_e/S_can); %eq 3.33 PAMADI
    a_WB_can = (KN_can + KWB_can + KBW_can)*CLalpha_can_e_pw*(S_can_e/S_can); %eq 3.24 PAMADI
    CLalpha_WB_can = a_WB_can; % contribution of ww1 to fuselage and fuselage to w1
    CL_alpha_wb_can = CLalpha_WB_can;
    %% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
    CL_alpha_wb_can = CLalpha_can_e_pw;
    CL_alpha_can_corrected = CL_alpha_wb_can;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_can = CL_alpha_wb_can;
end

%% W2 body interference (HTP or Vee tail)
if HTP == 1
    aN_w2 = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    aN_w2 = CLa_fus; %  corrigiendo con la integrtaci�n del fuselaje a partir de modelo XFLR5
    KN_w2 = (aN_w2/CLalpha_HTP_e_pw)*(S_w2/S_w2_s);  %eq 3.25 PAMADI
    CLalpha_fus = aN_w2;
    KWB_w2 = 0.1714*(w_Area_b_max/b_w2)^2 + 0.8326*(w_Area_b_max/b_w2) + 0.9974; %eq 3.27 PAMADI
    KBW_w2 = 0.7810*(w_Area_b_max/b_w2)^2 + 1.1976*(w_Area_b_max/b_w2) + 0.0088; %eq 3.28 PAMADI
    a_bw_w2 = KBW_w1*CLalpha_HTP_e_pw*(S_w2_s/S_w2); %eq 3.34 PAMADI
    a_wb_w2 = KWB_w1*CLalpha_HTP_e_pw*(S_w2_s/S_w2); %eq 3.33 PAMADI
    a_WB_w2 = (KN_w2 + KWB_w2 + KBW_w1)*CLalpha_HTP_e_pw*(S_w2_s/S_w2); %eq 3.24 PAMADI
    CLalpha_WB_w2 = a_WB_w2; % contribution of w2 to fuselage and fuselage to w2
    CLalpha_wb_w2 = CLalpha_WB_w2;
    %% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
    CL_alpha_wb_HTP = CLalpha_HTP_e_pw;
    CL_alpha_HTP_corrected = CLalpha_wb_w2;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_HTP = CL_alpha_wb_HTP;
end

%% W2 body interference (HTP or Vee tail)
if Vee==1
    aN_w2 = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    aN_w2 = CLa_fus; %  corrigiendo con la integrtaci�n del fuselaje a partir de modelo XFLR5
    KN_w2 = (aN_w2/CLalpha_Vee_e_pw)*(S_w2/S_w2_s);  %eq 3.25 PAMADI
    CLalpha_fus = aN_w2;
    KWB_w2 = 0.1714*(w_Area_b_max/b_w2)^2 + 0.8326*(w_Area_b_max/b_w2) + 0.9974; %eq 3.27 PAMADI
    KBW_w2 = 0.7810*(w_Area_b_max/b_w2)^2 + 1.1976*(w_Area_b_max/b_w2) + 0.0088; %eq 3.28 PAMADI
    a_bw_w2 = KBW_w1*CLalpha_Vee_e_pw*(S_w2_s/S_w2); %eq 3.34 PAMADI
    a_wb_w2 = KWB_w1*CLalpha_Vee_e_pw*(S_w2_s/S_w2); %eq 3.33 PAMADI
    a_WB_w2 = (KN_w2 + KWB_w2 + KBW_w1)*CLalpha_Vee_e_pw*(S_w2_s/S_w2); %eq 3.24 PAMADI
    CLalpha_WB_w2 = a_WB_w2; % contribution of w2 to fuselage and fuselage to w2
    CLalpha_wb_w2 = CLalpha_WB_w2;
    %% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
    CL_alpha_wb_Vee = CLalpha_Vee_e_pw;
    CL_alpha_Vee_corrected = CLalpha_wb_w2;
    % Correction of CLalphais for wing with no dihedral
    %     CL_alpha_wb_Vee = CLalpha_Vee_e_pw*(cos(dihedral_w2_e))^2;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_Vee = CL_alpha_wb_Vee;

end

CL_alpha_wb_nac = 0;
% Storing Values
Stab_Der_parts.CL_alpha_fus = CLalpha_fus;
Stab_Der_parts.CL_alpha_nac = CL_alpha_wb_nac;

%%%%%%%%%%%%%%%%%CALCULO CENTRO AERODINAMICO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_f = sqrt(1 - Mach^2);
Check_method = beta_f*AR_w1_e;

% Digitalizacion figura 4.3.2.2-35 Datcom
x_f_xi = [0.0,.025,.050,.075,.10,.15,.20,.30,.40,.50,.60,.70,.80];
y_f_xi = [0.0,.056,.101,.130,.152,.190,.22,.266,.301,.33,.348,1.365,.375];
f_xi  = interp1(x_f_xi,y_f_xi,w_Area_b_max/b_w1,'spline');
% Revised version �lvaro
chi1    = chi_calc(b_w1, w_Area_b_max);
% Aerodynamic center of fuselage
% area derivative
dSdX = diff(Area_body)./diff(x_Area_body);
% Distance from the leading edge of the fuselage to the leading edge of the
% exposed wing root
l_N = x_cR_w1_LE;
%  Calculates where fuselage ceases to be potential
i=1;
while dSdX(i) > 0
    x_0_v = i;
    i=i+1;
end
x_0 = x_Area_body(x_0_v+1);

% Aerodynimic centers based in the leading edge
int_Sb = @(x) interp1(x_Area_body(1:end-1),dSdX,x).*(l_N - x);     %eq 3.36 PAMADI
Xac_cre_N = -(1/(cmac_w1_e*Area_b_max))*quad(@(x) int_Sb(x),0,x_0); %eq 3.36 PAMADI

if Check_method > 4
    Xac_cre_BW = 1/4 + (b_w1 - w_Area_b_max)/(2*cmac_w1_e)*f_xi*tan(Lambda_c4_w1); %eq 3.39 PAMADI
else
    Warning = 'WARNING!!! Revise method for estimating Xac with WB - Code in PAUSE';
    %     disp(Warning)
    %     pause
    % From Fig 3.20 obtain Xac_cre_BW for beta*Ae = 4
    Xac_cre_BW_4 = 1/4 + (b_w1 - w_Area_b_max)/(2*cmac_w1_e)*f_xi*tan(Lambda_c4_w1); %eq 3.39 PAMADI
    beta_Ae_4 = 4;
    % From Fig 3.21 obtain Xac_cre_BW for beta*Ae = 0
    % x-axis of Fig 3.21 (TFM fig 3.6)
    x_axis_F_3_21 = 1/4*(AR_w1_e*(1 + lambda_w1_e)*tan(Lambda_LE_w1));
    Xac_cre_BW_0 = 0.0; % Fig 3.21 with above results
    beta_Ae_0 = 0;
    % Linear interpolation
    slope_Xac_cre_BW_0 = (Xac_cre_BW_4 - Xac_cre_BW_0)/(beta_Ae_4 - beta_Ae_0); % Fig 3.21
    Xac_cre_BW = Xac_cre_BW_0 + slope_Xac_cre_BW_0*beta_f*AR_w1_e;
    %% Revisado �lvaro
    Xca_cre_BW4     = 1/4 + (b_w1 - w_Area_b_max)/(2*cmac_w1_e)*chi1*tan(Lambda_c4_w1);
    Xca_cre_BW0     = xac_cre_BW0_calc(AR_w1_e, lambda_w1, Lambda_LE_w1);
    Xca_cre_BW      = interp1([0 4], [Xca_cre_BW0 Xca_cre_BW4], Check_method, 'linear');
end

Xac_cre_WB = xbar_w1/cmac_w1_e; %eq 3.38 PAMADI % assumes the influence of the body on the location of the wing aerodynamic center is small
XAC_WB_cre = ((Xac_cre_N*aN_w1 + Xac_cre_WB*a_wb_w1 + Xac_cre_BW*a_bw_w1)/a_WB_w1)*cmac_w1_e; %eq 3.32 PAMADI
% Xace at the appex (w1 croot LE)
XAC_WB_LE = XAC_WB_cre*(cmac_w1_e/cmac_w1);                                        %eq 3.36 PAMADI
% Convertin Xac appex to reference origin
XAC_WB = XAC_WB_LE + x_w1_LE;                                                     %eq 3.36 PAMADI
x_xbar_wb_w1 = XAC_WB;

%% CL_alpha of ac
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail

% Identifies which CLalpha is used, the prop-wash correction or the body
% influence
body_interference = 0;

if W1 == 1
    switch body_interference
        case 0
            CL_alpha_w1 = CLalpha_w1_e_pw;
        case 1
            CL_alpha_w1 = CL_alpha_wb_w1;
    end
end

if Can == 1
    switch body_interference
        case 0
            CL_alpha_can = CLalpha_can_e_pw;
        case 1
            CL_alpha_can = CL_alpha_wb_can;
    end
end

if HTP == 1
    switch body_interference
        case 0
            CL_alpha_HTP = CLalpha_HTP_e_pw;
        case 1
            CL_alpha_HTP = CL_alpha_wb_HTP;
    end
end

if Vee == 1
    switch body_interference
        case 0
            CL_alpha_Vee = CLalpha_Vee_e_pw;
        case 1
            CL_alpha_Vee = CL_alpha_wb_Vee;
    end
end

if VTP == 1
    switch body_interference
        case 0
            CLalpha_VTP = CLalpha_VTP_e_pw;
        case 1
            CLalpha_VTP = CL_alpha_wb_VTP;
    end
end

if AC_type == 1
    CL_alpha_ac = CL_alpha_w1;
    X_NP = (CL_alpha_w1*XAC_WB)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1; 
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;    
elseif AC_type == 2
    CL_alpha_ac = CL_alpha_w1 + CL_alpha_HTP*(downwash);
    X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_HTP*(downwash)*x_xbar_w2)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_HTP_e_corrected + CL_alpha_HTP*(i_w2 - eps_w2);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_HTP_e_corrected = CL0_HTP_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_HTP = CL0_HTP_e_corrected;
    Stab_Der_parts.CL_alpha_htp = CL_alpha_HTP;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
elseif AC_type == 3
    CL_alpha_ac = CL_alpha_w1 + CL_alpha_can*(upwash) + CL_alpha_HTP*(downwash);
    X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can + CL_alpha_HTP*(downwash)*x_xbar_w2)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_HTP_e_corrected + CL_alpha_HTP*(i_w2 - eps_w2) + ...
        CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_HTP_e_corrected = CL0_HTP_e_corrected;
    Trim_ITER.CL0_can_e_corrected = CL0_can_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_HTP = CL0_HTP_e_corrected;
    Stab_Der_parts.CL_alpha_htp = CL_alpha_htp;
    Stab_Der_parts.CL0_can = CL0_can_e_corrected;
    Stab_Der_parts.CL_alpha_can = CL_alpha_can;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
elseif AC_type == 4
    CL_alpha_ac = CL_alpha_w1 + CL_alpha_Vee*(downwash);
    X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_Vee*(downwash)*x_xbar_w2)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_Vee_e_corrected + CL_alpha_Vee*(i_w2 - eps_w2);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_Vee_e_corrected = CL0_Vee_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_Vee = CL0_Vee_e_corrected;
    Stab_Der_parts.CL_alpha_Vee = CL_alpha_Vee;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
elseif AC_type == 5
    CL_alpha_ac = CL_alpha_w1 + CL_alpha_can*(upwash) + CL_alpha_Vee*(downwash);
    X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can + CL_alpha_Vee*(downwash)*x_xbar_w2)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_Vee_e_corrected + CL_alpha_Vee*(i_w2 - eps_w2) + ...
        CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_Vee_e_corrected = CL0_Vee_e_corrected;
    Trim_ITER.CL0_can_e_corrected = CL0_can_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_Vee = CL0_Vee_e_corrected;
    Stab_Der_parts.CL_alpha_Vee = CL_alpha_Vee;
    Stab_Der_parts.CL0_can = CL0_can_e_corrected;
    Stab_Der_parts.CL_alpha_can = CL_alpha_can;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
end

%% To replace by wind tunnel models
% Fuselage and nacelle lift coefficient at zero angle of attack
CL0_fus = 0;
CL0_nac = 0;
Stab_Der_parts.CL0_fus = CL0_fus;
Stab_Der_parts.CL0_nac = CL0_nac;

%% NOTE IMPORTANT
%% Estimation of SM with only CLalpha contribution
SM_Excel_No_fus = (X_NP - x_XCG)/cmac_w1;
x_XCG_Excel = x_XCG;
Trim_ITER.x_XCG_Excel = x_XCG_Excel;
Trim_ITER.X_NP = X_NP;
Trim_ITER.SM_Excel_No_fus = SM_Excel_No_fus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cm_0_fus%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nose section of Fuselage
% Digitalizacion Fig 3.9(b) Pamadi
%CALCULO Cmo eq 3.8 PAMADI
int_low_anglef = 0;
% wing_body_CLalpha = 0.0785;
int_bf = @(x) interp1(length_x_position,width_x_position,x);
alpha_CL_0_we_fus = (alpha_CL_0_w1*D2R + i_w1);
int_up_anglef = length_fus;
f_anglef = @(x) interp1(x_interp,Angle_fus_interp,x);
int_Cm0 = quad(@(x) ((int_bf(x)).^2.*(alpha_CL_0_we_fus + f_anglef(x))),int_low_anglef,int_up_anglef);
CM0_fus = (f_k2_k1/(36.5*S_w1*cmac_w1))*int_Cm0;
Stab_Der_parts.CM0_fus = CM0_fus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cm_alpha%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nose section of Fuselage
% Digitalizacion Fig 3.9(b) Pamadi
% x_cre_1 = 0.8;
y_pdf_1 = 14.43;
y_real_1 = 1;
x_d_epsu_d_alpha_1 = [0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4];
y_d_epsu_d_alpha_1 = (y_real_1/y_pdf_1)*[4.56,3.39,2.62,2.16,1.62,1.0,0.98,0.87,0.66,0.61];
% f_d_epsu_d_alpha_1  = interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,x_cre_1,'spline');

% Fuselage section right before wing
% Digitalizacion Fig 3.9(a) Pamadi
% x1bar_cre = 0.2;
y_pdf_2 = 57.88;
y_real_2 = 4;
x_d_epsu_d_alpha_2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
y_d_epsu_d_alpha_2 = (y_real_2/y_pdf_2)*[52.03,36.43,28.43,22.99,18.88,16.22,14.37,13.03,11.70,10.98];
% f_d_epsu_d_alpha_2  = interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,x1bar_cre,'spline');

% selection of segments
% x1bar_sec = x_w1_LE - cmac_w_e;
x1bar_sec = x_w1_LE - cR_w1;

wing_body_CLalpha = 0.0785;
% int_bf = @(x) interp1(length_x_position,width_x_position,x);
% M�todo III eq 3.7 PAMADI
% Section from tip-nose to front (from graphic)
f_d_epsu_d_alpha_f1  = @(x) interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,x,'linear','extrap').*((x_w1_LE) - x);
int_low_1 = 0;
int_up_1 = (x_w1_LE - cR_w1);
% Determines if there are sections farther away from the LE of the wing.
% The limit is 1 wing root chord - cR_w1
if int_up_1 > 0
    % if there are sections computes the integral
    int_Cmalpha_1 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f1(x)),int_low_1,int_up_1);
else
    % No sections
    int_Cmalpha_1 = 0;
end
%   CMalpha_f1 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_1;
CMalpha_f1 = (CL_alpha_wb_w1/(wing_body_CLalpha*R2D))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_1;

% Section from right before wing (from graphic)
f_d_epsu_d_alpha_f2  = @(x) interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,x,'linear','extrap').*((x_w1_LE) - x);

% Determines if there are sections right aft the LE of the wing.
% The limit is 1 wing root chord - cR_w1
if int_up_1 > 0
    % if there are sections computes the integral with the limits
    int_low_2 = (x_w1_LE - cR_w1);
    int_up_2 = x_w1_LE;
    int_Cmalpha_2 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f2(x)),int_low_2,int_up_2);
else
    % if there are sections no sections farther away, computes if there are any portions right aft the wing less that 1 wing root chord from LE
    int_low_2 = 0;
    int_up_2 = x_w1_LE;
    if (int_up_2- int_low_2)>0
        int_Cmalpha_2 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f2(x)),int_low_2,int_up_2);
    else
        int_Cmalpha_2 = 0;
    end
end

%   CMalpha_f2 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_2;
CMalpha_f2 = (a_WB_w1/(wing_body_CLalpha*R2D))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_2;

% Section from behind the wing (from eq 3.10 PAMADI)
int_low_3 = x_w1_LE + cmac_w1_e;
int_up_3 = length_fus;
% Determines if there are sections behind the wing.
% The limit is 1 wing root chord - cR_w1
if (int_up_3-int_low_3) > 0
    % if there are sections computes the integral
    x_3 = (x_w2_LE  + cmac_w2/4) - (x_w1_LE + cmac_w1_e);
    f_d_epsu_d_alpha_f3 = @(x) ((x-(x_w1_LE + cmac_w1_e))/x_3).*downwash;
    int_Cmalpha_3 = quad(@(x) ((int_bf(x)).^2).*f_d_epsu_d_alpha_f3(x),int_low_3,int_up_3);
else
    % No sections
    int_Cmalpha_3 = 0;
end

%   CMalpha_f3 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_3;
CMalpha_f3 = (a_WB_w1/(wing_body_CLalpha*R2D))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_3;
C_Malpha_fus = CMalpha_f1 + CMalpha_f2 + CMalpha_f3;
Stab_Der_parts.CM_alpha_fus = C_Malpha_fus;

Angle_fus_x_fus = -Angle_fus_x_fus;
Angle_fus_interp = -Angle_fus_interp;

if W1 == 1
    alpha_CL_0_w1 = Aero.alpha_CL_0_w1_CR;
end
if HTP == 1 | Vee == 1
    alpha_CL_0_w2 = Aero.alpha_CL_0_w2_CR;
end
if Can == 1
    alpha_CL_0_can = Aero.alpha_CL_0_can_CR;
end

%% Propulsive Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%DERIVADA PROPULSIVA LONGITUDINAL%%%%%%%%%%%%%%%%
CD_Total_prel = C_D0 + C_D1*CL + C_D2*CL^2;  %polar aeronave
Drag = q_inf*S_ref*CD_Total_prel;

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CTx1 %%%%%%%%%%%%%%%%%%%%%%%%%
phiT = 0; %	is the thrust line inclination angle.
T_set = Propulsion.Ti;
CTx1 = CD_Total_prel;
% CTx1 = T_set*cos(phiT+trim_alpha)/q_inf*S_ref;
Stab_Der.T = T_set;

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CMT1 %%%%%%%%%%%%%%%%%%%%%%%%%
%% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
switch type_engine
    case 1 % TIPO DE MOTOR --> 1_TURBOFAN 
        % Increment of Pitching Moment Coefficient due to the Lift Component of the Propeller Normal Force Aproximada a 0
        Delta_CM_Nprop = 0;
        Delta_CM_Tprop = -(d_T/cmac_w1)*CTx1;
        CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;
    case 2 % TIPO DE MOTOR --> 2_TURBOPROP
        % Increment of Pitching Moment Coefficient due to the Lift Component of the Propeller Normal Force Aproximada a 0
        Delta_CM_Nprop = 0;
        Delta_CM_Tprop = -(d_T/cmac_w1)*CTx1;
        CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;
    case 3 % TIPO DE MOTOR --> 3_PISTON 
        CMT1 = - (d_T/cmac_w1)*CTx1;
    case 4 % TIPO DE MOTOR --> 4_ELECTRIC_PROP
        % Increment of Pitching Moment Coefficient due to the Lift Component of the Propeller Normal Force Aproximada a 0
        Delta_CM_Nprop = 0;
        Delta_CM_Tprop = -(d_T/cmac_w1)*CTx1;
        CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;
end
%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CMTalpha %%%%%%%%%%%%%%%%%%%%%%%%%
%% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
switch type_engine
    case 1 % TIPO DE MOTOR --> 1_TURBOFAN
        n_j = n_eng; % number of jet engines
        switch bypass_ratio
            case 1 % BPR from 0 to 1
                k_gas = 0.0003;
            case 2 % BPR from 1 to 2
                k_gas = 0.0007;
            case 3 % BPR from 2 to 4
                k_gas = 0.0009;
            case 4 % BPR from 4 to 6
                k_gas = 0.0011;
        end
        slugs_s2kg_s = 14.593903; % Convertsfrom slugs/s to kg/s
        m_dot_gass = T_eng*k_gas*slugs_s2kg_s;
        m_dot_cool = 0.006*m_dot_gass;
        m_dot = m_dot_gass + m_dot_cool;
        % is the inlet cross sectional area for Max Thrust
        % effect of thrustline offset on longitudinal stability
        dCM_dCL_TL = 0 ; % Negligible for Jet engines
        % effect of propeller or inlet normal force on lomngitudinal stability
        dCM_dCL_N = 0.035*m_dot*(-x_d_T)*depsu_dalpha/(S_w1*cmac_w1*rho*V*CL_alpha_ac);
        Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;        
        CMTalpha = Delta_CM_CL_T*CL_alpha_ac;
%         CMTalpha = 0;
        % Data from NACA REport 640 - properties of airfoil
        w_R_30 = 0.0525*2;
        w_R_60 = 0.073*2;
        w_R_90 = 0.045*2;        

    case 2 % TIPO DE MOTOR --> 2_TURBOPROP
        Nprop = 1; % number of props
        % Data from NACA REport 640 - properties of airfoil
        w_R_30 = 0.0525*2;
        w_R_60 = 0.073*2;
        w_R_90 = 0.045*2;        
        % Relation between propeller geometry - pitch and beta
        pitch =12*D2R; % Pitch of propeller 22x12W
        beta_3_4R = atan(pitch/(2*pi*(2/3)*R_prop));
        % DATCOM - PROPELLER INFLOW FACTOR - FIGURE 4.6.1-25B
        % Nominal blade angle at 0.75R Radius in deg
        Beta_blade =[15.,20.,25.,30.,35.,40.,60.];
        % Data for CNa
        beta_blade = beta_3_4R*R2D;
        % N_blades = [2.,3.,4.,6.];
        n_blades=2;
        if n_blades==2
            CNa_n = [.08,.10,.115,.126,.140,.150,.192];
        elseif n_blades==3
            CNa_n = [.11,.139,.160,.182,.20,.216,.275];
        elseif n_blades==5
            CNa_n = [.136,.172,.20,.226,.25,.272,.35];
        elseif n_blades==6
            CNa_n = [.196,.237,.275,.315,.35,.382,.5];
        end
        C_Na_807  = interp1(Beta_blade,CNa_n,beta_blade,'spline');
        % Airplane Design; Roskam, Jan; Part VI, pag 342 (2032 PDF)
        % depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop, cr_we, XLE_w, l_fus, downw, CLa_WB*S_w/Sref); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
        % n_blades = 2;
        % beta_blade = 20; % beta_blade in degrees!!
        % CNa_807         = CNa_807_calc(n_blades, beta_blade)
        KN              = 262*(2*w_R_30/D_prop) + 262*(2*w_R_60/D_prop) + 135*(2*w_R_90/D_prop);
        dCN_dalpha  = C_Na_807*(1 + 0.8*(KN/80.7 - 1));
        
        % Distance from the end of wing to the end of fuselage
        % x_1 = x_m_propeller - (x_w_LE + cmac_w1);
        % l_h = l_ht - (3*cmac_w1/4);
        % depsu_dalpha = (1-deps_dalpha)*(x1/l_h)-1;
        
        % % Location of the prop disk (source of thrust)
        % x_prop_cF = Geo_tier.x_prop_cF;
        % AR_w1 = Geo_tier.AR_w1;
        % depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);
        
        % rho_met2imp = conv_UNITS.rho_met2imp;
        % W2hp = conv_UNITS.W2hp;
        
        T_set_eng = Propulsion.Ti_eng;
        CT_p_eng = T_set_eng/(q_inf*S_ref);
        % P_SHP = P_SET*W2hp;
        Pi_eng = Propulsion.Pi_eng; % power per engine
        P_SHP = Pi_eng*W2hp;
        % rho_met2imp = 0.0017645;
        W_current_lb = m_TOW*2.20462;
        % The intermediate calculation parameter is given by:
        dTc_dCL = (3/2)*550*P_SHP*sqrt(rho*rho_SI2rho_IMP)*sqrt(CL)/...
            (sqrt((2*W_current_lb/(S_w1*m22ft2))^3)*(D_prop*m2ft)^2);
        % effect of thrustline offset on longitudinal stability
        dCM_dCL_TL = (Nprop*(2*((D_prop)^2)*d_T)/(S_ref*cmac_w1))*dTc_dCL;
        % DATCOM -
        % ----PROPELLER INFLOW FACTOR
        %      ----FIGURE 4.6.1-25B
        f_inflow_input_vec = [0.,1.,2.,3.,4.,6.,8.,10.,14.,19.,22.];
        f_inflow_input = S_ref*(CT_p_eng)/(8*(R_prop^2));
        f_inflow_vec = [1.0,1.55,1.94,2.20,2.40,2.75,3.05,3.30,3.75,4.25,4.54];
        finflow  = interp1(f_inflow_input_vec,f_inflow_vec,f_inflow_input,'spline');
        % effect of propeller or inlet normal force on longitudinal stability is given by:
        l_prop = (-1)*x_d_T*cos(phi_T) - (-1)*y_d_T*sin(phi_T);
        dCM_dCL_N = (pi/4)*finflow*Nprop*l_prop*((D_prop)^2)...
            *dCN_dalpha*(1 - depsu_dalpha)/(S_w1*cmac_w1*CL_alpha_ac);
        Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;
        CMTalpha = Delta_CM_CL_T*CL_alpha_ac;
    case 3 % TIPO DE MOTOR --> 3_PISTON
        Nprop = 1; % number of props
        % Data from NACA REport 640 - properties of airfoil
        w_R_30 = 0.0525*2;
        w_R_60 = 0.073*2;
        w_R_90 = 0.045*2;       
        % Relation between propeller geometry - pitch and beta
        pitch =12*D2R; % Pitch of propeller 22x12W
        beta_3_4R = atan(pitch/(2*pi*(2/3)*R_prop));
        % DATCOM - PROPELLER INFLOW FACTOR - FIGURE 4.6.1-25B
        % Nominal blade angle at 0.75R Radius in deg
        Beta_blade =[15.,20.,25.,30.,35.,40.,60.];
        % Data for CNa
        beta_blade = beta_3_4R*R2D;
        % N_blades = [2.,3.,4.,6.];
        n_blades=2;
        if n_blades==2
            CNa_n = [.08,.10,.115,.126,.140,.150,.192];
        elseif n_blades==3
            CNa_n = [.11,.139,.160,.182,.20,.216,.275];
        elseif n_blades==5
            CNa_n = [.136,.172,.20,.226,.25,.272,.35];
        elseif n_blades==6
            CNa_n = [.196,.237,.275,.315,.35,.382,.5];
        end
        C_Na_807  = interp1(Beta_blade,CNa_n,beta_blade,'spline');
        % Airplane Design; Roskam, Jan; Part VI, pag 342 (2032 PDF)
        % depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop, cr_we, XLE_w, l_fus, downw, CLa_WB*S_w/Sref); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
        % n_blades = 2;
        % beta_blade = 20; % beta_blade in degrees!!
        % CNa_807         = CNa_807_calc(n_blades, beta_blade)
        KN              = 262*(2*w_R_30/D_prop) + 262*(2*w_R_60/D_prop) + 135*(2*w_R_90/D_prop);
        dCN_dalpha  = C_Na_807*(1 + 0.8*(KN/80.7 - 1));
        
        % Distance from the end of wing to the end of fuselage
        % x_1 = x_m_propeller - (x_w_LE + cmac_w1);
        % l_h = l_ht - (3*cmac_w1/4);
        % depsu_dalpha = (1-deps_dalpha)*(x1/l_h)-1;
        
        % % Location of the prop disk (source of thrust)
        % x_prop_cF = Geo_tier.x_prop_cF;
        % AR_w1 = Geo_tier.AR_w1;
        % depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);
        
        % rho_met2imp = conv_UNITS.rho_met2imp;
        % W2hp = conv_UNITS.W2hp;
        
        T_set_eng = Propulsion.Ti_eng;
        CT_p_eng = T_set_eng/(q_inf*S_ref);
        % P_SHP = P_SET*W2hp;
        Pi_eng = Propulsion.Pi_eng; % power per engine
        P_SHP = Pi_eng*W2hp;
        % rho_met2imp = 0.0017645;
        W_current_lb = m_TOW*2.20462;
        % The intermediate calculation parameter is given by:
        dTc_dCL = (3/2)*550*P_SHP*sqrt(rho*rho_SI2rho_IMP)*sqrt(CL)/...
            (sqrt((2*W_current_lb/(S_w1*m22ft2))^3)*(D_prop*m2ft)^2);
        % effect of thrustline offset on longitudinal stability
        dCM_dCL_TL = (Nprop*(2*((D_prop)^2)*d_T)/(S_ref*cmac_w1))*dTc_dCL;
        % DATCOM -
        % ----PROPELLER INFLOW FACTOR
        %      ----FIGURE 4.6.1-25B
        f_inflow_input_vec = [0.,1.,2.,3.,4.,6.,8.,10.,14.,19.,22.];
        f_inflow_input = S_ref*(CT_p_eng)/(8*(R_prop^2));
        f_inflow_vec = [1.0,1.55,1.94,2.20,2.40,2.75,3.05,3.30,3.75,4.25,4.54];
        finflow  = interp1(f_inflow_input_vec,f_inflow_vec,f_inflow_input,'spline');
        % effect of propeller or inlet normal force on longitudinal stability is given by:
        l_prop = (-1)*x_d_T*cos(phi_T) - (-1)*y_d_T*sin(phi_T);
        dCM_dCL_N = (pi/4)*finflow*Nprop*l_prop*((D_prop)^2)...
            *dCN_dalpha*(1 - depsu_dalpha)/(S_w1*cmac_w1*CL_alpha_ac);
        Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;
        CMTalpha = Delta_CM_CL_T*CL_alpha_ac;
    case 4 % TIPO DE MOTOR --> 4_ELECTRIC_PROP
        Nprop = 1; % number of props
        % Data from NACA REport 640 - properties of airfoil
        w_R_30 = 0.0525*2;
        w_R_60 = 0.073*2;
        w_R_90 = 0.045*2;        
        % Relation between propeller geometry - pitch and beta
        pitch =12*D2R; % Pitch of propeller 22x12W
        beta_3_4R = atan(pitch/(2*pi*(2/3)*R_prop));
        % DATCOM - PROPELLER INFLOW FACTOR - FIGURE 4.6.1-25B
        % Nominal blade angle at 0.75R Radius in deg
        Beta_blade =[15.,20.,25.,30.,35.,40.,60.];
        % Data for CNa
        beta_blade = beta_3_4R*R2D;
        % N_blades = [2.,3.,4.,6.];
        n_blades=2;
        if n_blades==2
            CNa_n = [.08,.10,.115,.126,.140,.150,.192];
        elseif n_blades==3
            CNa_n = [.11,.139,.160,.182,.20,.216,.275];
        elseif n_blades==5
            CNa_n = [.136,.172,.20,.226,.25,.272,.35];
        elseif n_blades==6
            CNa_n = [.196,.237,.275,.315,.35,.382,.5];
        end
        C_Na_807  = interp1(Beta_blade,CNa_n,beta_blade,'spline');
        
        % Airplane Design; Roskam, Jan; Part VI, pag 342 (2032 PDF)
        % depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop, cr_we, XLE_w, l_fus, downw, CLa_WB*S_w/Sref); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
        % n_blades = 2;
        % beta_blade = 20; % beta_blade in degrees!!
        % CNa_807         = CNa_807_calc(n_blades, beta_blade)
        KN              = 262*(2*w_R_30/D_prop) + 262*(2*w_R_60/D_prop) + 135*(2*w_R_90/D_prop);
        dCN_dalpha  = C_Na_807*(1 + 0.8*(KN/80.7 - 1));
        
        % Distance from the end of wing to the end of fuselage
        % x_1 = x_m_propeller - (x_w_LE + cmac_w1);
        % l_h = l_ht - (3*cmac_w1/4);
        % depsu_dalpha = (1-deps_dalpha)*(x1/l_h)-1;
        
        % % Location of the prop disk (source of thrust)
        % x_prop_cF = Geo_tier.x_prop_cF;
        % AR_w1 = Geo_tier.AR_w1;
        % depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);
        
        % rho_met2imp = conv_UNITS.rho_met2imp;
        % W2hp = conv_UNITS.W2hp;
        
        T_set_eng = Propulsion.Ti_eng;
        CT_p_eng = T_set_eng/(q_inf*S_ref);
        % P_SHP = P_SET*W2hp;
        Pi_eng = Propulsion.Pi_eng; % power per engine
        P_SHP = Pi_eng*W2hp;
        % rho_met2imp = 0.0017645;
        W_current_lb = m_TOW*2.20462;
        % The intermediate calculation parameter is given by:
        dTc_dCL = (3/2)*550*P_SHP*sqrt(rho*rho_SI2rho_IMP)*sqrt(CL)/...
            (sqrt((2*W_current_lb/(S_w1*m22ft2))^3)*(D_prop*m2ft)^2);
        % effect of thrustline offset on longitudinal stability
        dCM_dCL_TL = (Nprop*(2*((D_prop)^2)*d_T)/(S_ref*cmac_w1))*dTc_dCL;
        % DATCOM -
        % ----PROPELLER INFLOW FACTOR
        %      ----FIGURE 4.6.1-25B
        f_inflow_input_vec = [0.,1.,2.,3.,4.,6.,8.,10.,14.,19.,22.];
        f_inflow_input = S_ref*(CT_p_eng)/(8*(R_prop^2));
        f_inflow_vec = [1.0,1.55,1.94,2.20,2.40,2.75,3.05,3.30,3.75,4.25,4.54];
        finflow  = interp1(f_inflow_input_vec,f_inflow_vec,f_inflow_input,'spline');
        % effect of propeller or inlet normal force on longitudinal stability is given by:
        l_prop = (-1)*x_d_T*cos(phi_T) - (-1)*y_d_T*sin(phi_T);
        dCM_dCL_N = (pi/4)*finflow*Nprop*l_prop*((D_prop)^2)...
            *dCN_dalpha*(1 - depsu_dalpha)/(S_w1*cmac_w1*CL_alpha_ac);
        Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;
        CMTalpha = Delta_CM_CL_T*CL_alpha_ac;
end
% CMTalpha = 0;
% Stab_Der.CMTalpha = CMTalpha;

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas%%%%%%%%%%%%%%%%%%%%%%%%%
% Coeffcicients of second order Power model
% P = A_power*V^2 + B_power*V + C_power

%conversiones de unidades imperiales
% qmet2qimp = 10.76391*(3.28084^2)*(1.2250/0.0023772);
% W2pftsec = 0.7375621;
% m22ft2 = 10.76391;
% rho_SI2rho_IMP = (0.0023772/1.2250);
% qmet2qimp = (3.28084^2)*(1.2250/0.0023772);
% qmet2qimp = m22ft2*rho_SI2rho_IMP;
% N2lbf = 0.2248089;
% m2ft = 3.28083;

% For future versions
alpha = 0;
beta = 0;
% RPS = Propulsion.RPS;
% adimensional_P_Heli =rho*(RPS^3)*(D_prop^5);
% adimensional_P_Aero =q_inf*S_ref;
% adimensional_conv = adimensional_P_Aero/adimensional_P_Heli;
% dcP_du = d_CP_d_V*adimensional_conv;
% ctxu = (1/(q_inf*qmet2qimp*S_w1*m22ft2))*(2*A_power*V + B_power)*W2pftsec;

%% Type engine
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CTxu %%%%%%%%%%%%%%%%%%%%%%%%%
switch type_engine
    case 1 % TIPO DE MOTOR --> 1_TURBOFAN 
                CTxu = (V/(q_inf*S_ref))*d_T_d_V - 2*CTx1;
    case 2 % TIPO DE MOTOR --> 2_TURBOPROP 
        switch type_prop
            case 1
                %         CTxu = (1/(q_inf*S_ref))*d_P_d_V - 2*CTx1  %por helice a paso fijo
                CTxu = d_CP_d_V - 2*CTx1;  %por helice a paso fijo
                CTxu = - 3*CTx1 + (CTx1*J*d_etap_d_J/etha_emp);   %por helice a paso fijo
            case 2
                CTxu = - 3*CTx1;  %por helice a paso variable
        end
    case 3 % TIPO DE MOTOR --> 3_PISTON
        switch type_prop
            case 1
                %         CTxu = (1/(q_inf*S_ref))*d_P_d_V - 2*CTx1  %por helice a paso fijo
                CTxu = d_CP_d_V - 2*CTx1;  %por helice a paso fijo
                CTxu = - 3*CTx1 + (CTx1*J*d_etap_d_J/etha_emp);   %por helice a paso fijo
            case 2
                CTxu = - 3*CTx1;  %por helice a paso variable
        end
    case 4% TIPO DE MOTOR --> 4_ELECTRIC_PROP
        switch type_prop
            case 1
%                 CTxu = (1/(q_inf*S_ref))*d_P_d_V - 2*CTx1;  %por helice a paso fijo
                CTxu = d_CP_d_V - 2*CTx1;  %por helice a paso fijo
                CTxu = - 3*CTx1 + (CTx1*J*d_etap_d_J/etha_emp);   %por helice a paso fijo
            case 2
                CTxu = - 3*CTx1;  %por helice a paso variable
        end
end

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CTxalpha %%%%%%%%%%%%%%%%%%%%%%%%%
CTxalpha = 0;

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CMtu %%%%%%%%%%%%%%%%%%%%%%%%%
CMtu = - (d_T/cmac_w1)*CTxu;

Stab_Der.CTx1 = CTx1;
Stab_Der.CTxu = CTxu;
Stab_Der.CTxalpha = CTxalpha;
Stab_Der.CMt1 = CMT1;
Stab_Der.CMtu = CMtu;
Stab_Der.CMtalpha = CMTalpha;
Stab_Der.CL = CL;
Stab_Der.CD = CD;

%% NOTE Corrections of x_XCG taking into avvount fuselage CMalpha contribution
% CM_alpha CONTRIBUTION ONLY TAKING INTO ACCOUNT CLalpha
SM_Excel_No_fus;
CM_alpha_ac_No_fus = -CL_alpha_ac*SM_Excel_No_fus;
CM_alpha_ac_fus = CM_alpha_ac_No_fus + CMTalpha + C_Malpha_fus; % Total CMalpha with fuselage and propulsive contribution - Uses Munk's contribution
SM_Excel_fus = -CM_alpha_ac_fus/CL_alpha_ac; % Real SM with fuselage contribution

% Shift in wing+fuselage aerodynamic center from wing aerodynamic center
% caused by the so called Munk effect
X_NP_fus = (x_XCG - cmac_w1*CM_alpha_ac_fus/CL_alpha_ac);
Delta_X_NP_ac_fus = X_NP_fus - X_NP;
Percentage_inc_X_NP = ((Delta_X_NP_ac_fus)/X_NP)*100;

Trim_ITER.CM_alpha_ac_No_fus = CM_alpha_ac_No_fus;
Trim_ITER.CM_alpha_ac_fus = CM_alpha_ac_fus;
Trim_ITER.SM_Excel_fus = SM_Excel_fus;
Trim_ITER.X_NP_fus = X_NP_fus;
Trim_ITER.Delta_X_NP_ac_fus = Delta_X_NP_ac_fus;
Trim_ITER.Percentage_inc_X_NP = Percentage_inc_X_NP;

% Adjusting x_XCG to ensure longitudinal stability due to the fuselage
% contribution

% Gathers all the flags
Munk_fuselage_constribution = OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution;
if Munk_fuselage_constribution == 1;
    x_XCG_fus = X_NP_fus - SM_des*cmac_w1;
    CM_alpha_ac = CM_alpha_ac_fus;
    x_XCG = x_XCG_fus;
    X_NP = X_NP + Delta_X_NP_ac_fus;
    SM = (X_NP - x_XCG)/cmac_w1;
else
    x_XCG_No_fus = X_NP - SM_des*cmac_w1;
    CM_alpha_ac = CM_alpha_ac_No_fus;
    x_XCG = x_XCG_No_fus;
    X_NP = X_NP;
    SM = (X_NP - x_XCG)/cmac_w1;
end

TRIM_RESULTS.CM_alpha_ac = CM_alpha_ac;
TRIM_RESULTS.x_XCG = x_XCG;
TRIM_RESULTS.X_NP = X_NP;
TRIM_RESULTS.SM = SM;
SM
% % x_XCG with no fuselage contribution
% % Revised values
% study_var_xcg = conditions.study_var_xcg;
% if study_var_xcg == 1
%     CM_alpha_ac = -CL_alpha_ac*((X_NP-x_XCG_No_fus)/cmac_w1) + CMTalpha + C_Malpha_fus;
%     x_XCG = x_XCG_No_fus;
% else
%     CM_alpha_ac = -CL_alpha_ac*((X_NP-x_XCG_fus)/cmac_w1) + CMTalpha + C_Malpha_fus;
%     x_XCG = x_XCG_fus;
% end

Stab_Der.CL_alpha_ac = CL_alpha_ac;
Stab_Der.CM_alpha_ac = CM_alpha_ac;

% SM_actual = -CM_alpha_ac/CL_alpha_ac
% SM_real = (X_NP - x_XCG_fus)/cmac_w1
% pause

% Store DATA
% Trim_ITER.CM_alpha_ac_CLalpha = CM_alpha_ac_CLalpha;
% Trim_ITER.CM_alpha_ac_fus = CM_alpha_ac_fus;
Trim_ITER.CMTalpha = CMTalpha;
Trim_ITER.CMalpha_f = C_Malpha_fus;
% Trim_ITER.SM_Excel_fus = SM_Excel_fus;
% Trim_ITER.CM_alpha_ac = CM_alpha_ac;
% Trim_ITER.SM_actual = SM_actual;
% Trim_ITER.SM_real = SM_real;

%%%%%%%%%%%%%%%%%%%%%%%CALCULO CM0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if W1 == 1
    CM_0_w1_ae = Aero.CM_0_w1_CR;
    CM0_w1_e_corrected = (eta_w1_afe_S_w1_afe_S_ref*CM_0_w1_ae + eta_w1_no_afe_S_w1_no_afe_S_ref*CM_0_w1_ae)*(cmac_w1/cmac_w1);
    CM0_w1_CL0_w1 = ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CL_alpha_w1*i_w1);
    CM_0_w1_wrt_XCG = CM_0_w1_ae + CM0_w1_CL0_w1;
    Trim_ITER.CM0_w1_e_corrected = CM0_w1_e_corrected;
    Trim_ITER.CM0_w1_CL0_w1 = CM0_w1_CL0_w1;
    Stab_Der_parts.CM0_w1_ae = CM0_w1_e_corrected;
    Stab_Der_parts.CM0_w1_CL0 = CM0_w1_CL0_w1;
    Stab_Der_parts.CM_0_w1 = CM_0_w1_wrt_XCG;
end
if HTP == 1
    CM_0_w2_ae = Aero.CM_0_w1_CR;
    CM0_w2_e_corrected = (eta_w2_afe_S_w2_afe_S_ref*CM_0_w2_ae + eta_w2_no_afe_S_w2_no_afe_S_ref*CM_0_w2_ae)*(cmac_w2/cmac_w1);
    CM0_w2_CL0_w2 = ((x_XCG - x_xbar_w2)/cmac_w1)*(CL0_HTP_e_corrected + CL_alpha_HTP*(i_w2 - eps_w2));
    CM_0_w2_wrt_XCG = CM_0_w2_ae + CM0_w2_CL0_w2;
    Trim_ITER.CM0_w2_CL0_w2 = CM0_w2_CL0_w2;
    Trim_ITER.CM0_w2_e_corrected = CM0_w2_e_corrected;
    Trim_ITER.CM0_w2_CL0_w2 = CM0_w2_CL0_w2;
    Stab_Der_parts.CM0_w2_ae = CM0_w2_e_corrected;
    Stab_Der_parts.CM0_w2_CL0 = CM0_w2_CL0_w2;
    Stab_Der_parts.CM_0_w2 = CM_0_w2_wrt_XCG;
end
    
if Vee == 1
    CM_0_w2_ae = Aero.CM_0_w1_CR;
    CM0_w2_e_corrected = (eta_w2_afe_S_w2_afe_S_ref*CM_0_w2_ae + eta_w2_no_afe_S_w2_no_afe_S_ref*CM_0_w2_ae)*(cmac_w2/cmac_w1);
    CM0_w2_CL0_w2 = ((x_XCG - x_xbar_w2)/cmac_w1)*(CL0_Vee_e_corrected + CL_alpha_Vee*(i_w2 - eps_w2));
    CM_0_w2_wrt_XCG = CM_0_w2_ae + CM0_w2_CL0_w2;
    Trim_ITER.CM0_w2_CL0_w2 = CM0_w2_CL0_w2;
    Trim_ITER.CM0_w2_e_corrected = CM0_w2_e_corrected;
    Trim_ITER.CM0_w2_CL0_w2 = CM0_w2_CL0_w2;
    Stab_Der_parts.CM0_w2_ae = CM0_w2_e_corrected;
    Stab_Der_parts.CM0_w2_CL0 = CM0_w2_CL0_w2;
    Stab_Der_parts.CM_0_Vee = CM_0_w2_wrt_XCG;
end
if Can == 1
    CM_0_can_ae = Aero.CM_0_w1_can;
    CM0_can_e_corrected = (eta_can_afe_S_can_afe_S_ref*CM_0_can_ae + eta_can_no_afe_S_can_no_afe_S_ref*CM_0_can_ae)*(cmac_can/cmac_w1);
    CM0_can_CL0_can = ((x_XCG - x_xbar_can)/cmac_w1)*(CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can));
    CM_0_can_wrt_XCG = CM_0_can_ae + CM0_can_CL0_w1;
    Trim_ITER.CM0_can_CL0_can = CM0_can_CL0_can;
    Trim_ITER.CM0_can_e_corrected = CM0_can_e_corrected;
    Stab_Der_parts.CM0_can_ae = CM0_can_e_corrected;
    Stab_Der_parts.CM0_can_CL0 = CM0_can_CL0_can;
    Stab_Der_parts.CM_0_can = CM_0_can_wrt_XCG;
end

Trim_ITER.CM0_f = CM0_fus;
Trim_ITER.CMT1 = CMT1;

% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
if AC_type == 1
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CMT1;
elseif AC_type == 2
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_w2_wrt_XCG + CMT1;
elseif AC_type == 3
%     CM0_w1_wrt_XCG = CM0_w1_e_corrected + ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CL_alpha_w1_corrected*i_w1);
%     CM0_can_wrt_XCG = CM0_can_e_corrected*(cmac_can/cmac_w1) + ((x_XCG - x_xbar_can)/cmac_w1)*(CL0_can_e_corrected + CLalpha_can*(i_can + eps_can));
%     CM0_w2_wrt_XCG = CM0_w2_e_corrected*(cmac_w2/cmac_w1) + ((x_XCG - x_xbar_w2)/cmac_w1)*(CL0_HTP_e_corrected + CL_alpha_wb_HTP*(i_w2 - eps_w2));
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_can_wrt_XCG + CM0_w2_wrt_XCG + CMT1
elseif AC_type == 4
%     CM0_w1_wrt_XCG = CM0_w1_e_corrected + ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CL_alpha_w1_corrected*i_w1);
%     CM0_w2_wrt_XCG = CM0_w2_e_corrected*(cmac_w2/cmac_w1) + ((x_XCG - x_xbar_w2)/cmac_w1)*(CL0_Vee_e_corrected + CL_alpha_wb_Vee*(i_w2 - eps_w2));
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_w2_wrt_XCG + CMT1;
elseif AC_type == 5
%     CM0_w1_wrt_XCG = CM0_w1_e_corrected + ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CL_alpha_w1_corrected*i_w1);
%     CM0_can_wrt_XCG = CM0_can_e_corrected*(cmac_can/cmac_w1) + ((x_XCG - x_xbar_can)/cmac_w1)*(CL0_can_e_corrected + CLalpha_can*(i_can + eps_can));
%     CM0_w2_wrt_XCG = CM0_w2_e_corrected*(cmac_w2/cmac_w1) + ((x_XCG - x_xbar_w2)/cmac_w1)*(CL0_Vee_e_corrected + CL_alpha_wb_Vee*(i_w2 - eps_w2));
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_can_wrt_XCG + CM_0_w2_wrt_XCG + CMT1;
end

Trim_ITER.CM0_ac = CM0_ac;
Stab_Der.CM0_ac = CM0_ac;

%% Control allocation for longitudinal


% ASpro CLdeltae, CDdeltae, CMdeltae
%% Estimaci�n con proyecciones en el plano horizontal
if d_rudvtr == 1
    cf_rudvtr = Geo_tier.cf_rudvtr;
    cf_d_rudvtr = cf_rudvtr;
    t_c_rudvtr = Geo_tier.t_c_rudvtr;
    t_c_d_rudvtr = t_c_rudvtr;
    K_y1_rudvtr_w2 = Geo_tier.K_y1_rudvtr_w2;
    K_y2_rudvtr_w2 = Geo_tier.K_y2_rudvtr_w2;
    K_y1 = K_y1_rudvtr_w2;
    K_y2 = K_y2_rudvtr_w2;
end

if d_ele == 1
    cf_ele = Geo_tier.cf_ele;
    cf_d_ele = cf_ele;
    t_c_ele = Geo_tier.t_c_ele;
    t_c_d_ele = t_c_ele;
    K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
    K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;
    K_y1 = K_y1_ele_w2;
    K_y2 = K_y2_ele_w2;
end

if d_rudder == 1
    cf_rudder = Geo_tier.cf_rudder;
    cf_d_rudder = cf_rudder;
    t_c_rudder = Geo_tier.t_c_rudder;
    t_c_d_rudder = t_c_rudder;
    K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
    K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
    K_y1 = K_y1_rudder_VTP;
    K_y2 = K_y2_rudder_VTP;
end

if d_ail == 1
    cf_ail = Geo_tier.cf_ail;
    cf_d_ail = cf_ail;
    t_c_ail = Geo_tier.t_c_ail;
    t_c_d_ail = t_c_ail;
    K_y1_ail_w1 = Geo_tier.K_y1_ail_w1;
    K_y2_ail_w1 = Geo_tier.K_y2_ail_w1;
    K_y1 = K_y1_ail_w1;
    K_y2 = K_y2_ail_w1;
end

if d_flap == 1
    cf_flap = Geo_tier.cf_flap;
    cf_d_flap = cf_flap;
    t_c_flap = Geo_tier.t_c_flap;
    t_c_d_flap = t_c_flap;
    K_y1_flap_w1 = Geo_tier.K_y1_flap_w1;
    K_y2_flap_w1 = Geo_tier.K_y2_flap_w1;
    K_y1 = K_y1_flap_w1;
    K_y2 = K_y2_flap_w1;
end

lambda_w2 = Geo_tier.lambda_w2;
AR_w2_e = Geo_tier.AR_w2_e;

% % theoretical V-Tail sectional lift curve slope at Mach equal to zero
% k_prima              = 1;
% Cla_M0               = 2*pi + 5.0525*t_c_d;
% Kb_w2                = Kb_calc(K_y1, K_y2, lambda_w2);              % Fig 3.36 Pamadi
% Clde_Cldetheo_w2     = Cldelta_Cldeltatheory_calc(cf_d, Cla_M0);% Fig 3.37b Pamadi
% Cldetheo_w2          = Cldeltatheory_calc(cf_d, t_c_d);               % Fig 3.37a Pamadi
% alphaCL_alphaCl_w2   = alphaCL_alphaCl_calc(cf_d,AR_w2_e);              % Fig 3.35 Pamadi
% alpha_delta_e       = Kb_w2*Clde_Cldetheo_w2*Cldetheo_w2/Cla_M0*alphaCL_alphaCl_w2;

% Revision to determine CLalpha of the horizontal control surface
if d_rudvtr == 1
    CL_alpha_wb_long = CL_alpha_wb_Vee;
elseif d_ele == 1
    CL_alpha_wb_long = CL_alpha_wb_HTP;
end
% 
% CL_delta_e          = alpha_delta_e*CL_alpha_wb_long*(S_w2_e/S_ref);
% CD_delta_e          = 2*CL0_ac*C_D2*(eta_w2_afe_S_w2_afe_S_ref*CL_delta_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CL_delta_e);
% CM_delta_e            = - CL_delta_e*((x_xbar_w2 - x_XCG)/cmac_w1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPLET FOR DIFFERENT CONTROL SURFACES!
%% CL_delta_rv
if d_rudvtr == 1
    % theoretical V-Tail sectional lift curve slope at Mach equal to zero
    k_prima              = 1;
    Cla_M0               = 2*pi + 5.0525*t_c_d_rudvtr;
    Kb_w2                = Kb_calc(K_y1_rudvtr_w2, K_y2_rudvtr_w2, lambda_w2);              % Fig 3.36 Pamadi
    Clde_Cldetheo_w2     = Cldelta_Cldeltatheory_calc(cf_d_rudvtr, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_w2          = Cldeltatheory_calc(cf_d_rudvtr, t_c_d_rudvtr);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_w2   = alphaCL_alphaCl_calc(cf_d_rudvtr,AR_w2_e);              % Fig 3.35 Pamadi
    alpha_delta_e       = Kb_w2*Clde_Cldetheo_w2*Cldetheo_w2/Cla_M0*alphaCL_alphaCl_w2;

    %% CL_delta_rv
    Balance_rv = cf_rudvtr;
    %% CL_delta_rv
    Balance_rv = cf_rudvtr;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_rv = 0.83 + 0.30714*Balance_rv; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    CL_i_vee = (eta_w2_afe_S_w2_afe_S_ref*CL_alpha_wb_Vee + eta_w2_no_afe_S_w2_no_afe_S_ref*CL_alpha_wb_Vee);
    
    CL_delta_rv_0 = f_val_rv*CL_i_vee*alpha_delta_e;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    
    CL_delta_rv = CL_delta_rv_0*K_prima2;
    CL_delta_e = CL_delta_rv;
    
    %% CD_delta_rv
    CD_i_vee          = 2*CL0_ac*C_D2*(eta_w2_afe_S_w2_afe_S_ref*CL_alpha_wb_Vee + eta_w2_no_afe_S_w2_no_afe_S_ref*CL_alpha_wb_Vee);
    CD_delta_rv       = alpha_delta_e*CD_i_vee;
    CD_delta_e       = CD_delta_rv;
    
    %% CMdelta_rv
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_vee = -CL_i_vee*((x_xbar_w2 - x_XCG)/cmac_w1);
    CM_delta_rv_0 = f_val_rv*CM_i_vee*alpha_delta_e;
    CM_delta_rv = CM_delta_rv_0*K_prima2;
    CM_delta_e = CM_delta_rv;
    
elseif d_ele == 1
    % theoretical V-Tail sectional lift curve slope at Mach equal to zero
    k_prima              = 1;
    Cla_M0               = 2*pi + 5.0525*t_c_d_ele;
    Kb_w2                = Kb_calc(K_y1_ele_w2, K_y2_ele_w2, lambda_w2);              % Fig 3.36 Pamadi
    Clde_Cldetheo_w2     = Cldelta_Cldeltatheory_calc(cf_d_ele, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_w2          = Cldeltatheory_calc(cf_d_ele, t_c_d_ele);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_w2   = alphaCL_alphaCl_calc(cf_d_ele,AR_w2_e);              % Fig 3.35 Pamadi
    alpha_delta_e       = Kb_w2*Clde_Cldetheo_w2*Cldetheo_w2/Cla_M0*alphaCL_alphaCl_w2;

    Balance_e = cf_ele;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_e = 0.83 + 0.30714*Balance_e; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    CL_i_w2 = (eta_w2_afe_S_w2_afe_S_ref*CL_alpha_wb_HTP + eta_w2_no_afe_S_w2_no_afe_S_ref*CL_alpha_wb_HTP);
    
    CL_delta_e_0 = f_val_e*CL_i_w2*alpha_delta_e;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    CL_delta_e = CL_delta_e_0*K_prima2;
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
        %% CD_delta
    CD_i_h          = 2*CL0_ac*C_D2*(eta_w2_afe_S_w2_afe_S_ref*CL_alpha_wb_HTP + eta_w2_no_afe_S_w2_no_afe_S_ref*CL_alpha_wb_HTP);
    CD_delta_e       = alpha_delta_e*CD_i_h;
    
    %% CMdelta_rv
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_w2 = -CL_i_w2*((x_xbar_w2 - x_XCG)/cmac_w1);
    CM_delta_e_0 = f_val_e*CM_i_w2*alpha_delta_e;
    CM_delta_e = CM_delta_e_0*K_prima2;
end

if d_ail ==1
    % theoretical wing sectional lift curve slope at Mach equal to zero
    k_prima              = 1;
    Cla_M0               = 2*pi + 5.0525*t_c_d_ail;
    Kb_w1                = Kb_calc(K_y1_ail_w1, K_y2_ail_w1, lambda_w1);              % Fig 3.36 Pamadi
    Clde_Cldetheo_w1     = Cldelta_Cldeltatheory_calc(cf_d_ail, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_w1          = Cldeltatheory_calc(cf_d_ail, t_c_d_ail);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_w1   = alphaCL_alphaCl_calc(cf_d_ail,AR_w1_e);              % Fig 3.35 Pamadi
    alpha_delta_ail       = Kb_w1*Clde_Cldetheo_w1*Cldetheo_w1/Cla_M0*alphaCL_alphaCl_w1;

    Balance_ail = cf_ail;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_ail = 0.83 + 0.30714*Balance_ail; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    CL_i_w1 = (eta_w1_afe_S_w1_afe_S_ref*CL_alpha_w1 + eta_w1_no_afe_S_w1_no_afe_S_ref*CL_alpha_w1);
    
    CL_delta_ail_0 = f_val_ail*CL_i_w1*alpha_delta_ail;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    CL_delta_ail = CL_delta_ail_0*K_prima2;
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    
    %% CD_delta
    CD_i_h          = 2*CL0_ac*C_D2*(eta_w1_afe_S_w1_afe_S_ref*CL_alpha_w1 + eta_w1_no_afe_S_w1_no_afe_S_ref*CL_alpha_w1);
    CD_delta_ail       = alpha_delta_ail*CD_i_h;
    
    %% CMdelta_rv
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_w1 = -CL_i_w1*((x_xbar_w2 - x_XCG)/cmac_w1);
    CM_delta_ail_0 = f_val_ail*CM_i_w1*alpha_delta_ail;
    CM_delta_ail = CM_delta_ail_0*K_prima2;
end

if d_flap ==1
    % theoretical wing sectional lift curve slope at Mach equal to zero
    k_prima              = 1;
    Cla_M0               = 2*pi + 5.0525*t_c_d_flap;
    Kb_w1                = Kb_calc(K_y1_flap_w1, K_y2_flap_w1, lambda_w1);              % Fig 3.36 Pamadi
    Clde_Cldetheo_w1     = Cldelta_Cldeltatheory_calc(cf_d_flap, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_w1          = Cldeltatheory_calc(cf_d_flap, t_c_d_flap);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_w1   = alphaCL_alphaCl_calc(cf_d_flap,AR_w1_e);              % Fig 3.35 Pamadi
    alpha_delta_flap       = Kb_w1*Clde_Cldetheo_w1*Cldetheo_w1/Cla_M0*alphaCL_alphaCl_w1;

    Balance_flap = cf_flap;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_flap = 0.83 + 0.30714*Balance_flap; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    CL_i_w1 = (eta_w1_afe_S_w1_afe_S_ref*CL_alpha_w1 + eta_w1_no_afe_S_w1_no_afe_S_ref*CL_alpha_w1);
    
    CL_delta_flap_0 = f_val_flap*CL_i_w1*alpha_delta_flap;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    CL_delta_flap = CL_delta_flap_0*K_prima2;
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    
    %% CD_delta
    CD_i_h          = 2*CL0_ac*C_D2*(eta_w1_afe_S_w1_afe_S_ref*CL_alpha_w1 + eta_w1_no_afe_S_w1_no_afe_S_ref*CL_alpha_w1);
    CD_delta_flap       = alpha_delta_flap*CD_i_h;
    
    %% CMdelta_rv
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_w1 = -CL_i_w1*((x_xbar_w2 - x_XCG)/cmac_w1);
    CM_delta_flap_0 = f_val_ail*CM_i_w1*alpha_delta_flap;
    CM_delta_flap = CM_delta_flap_0*K_prima2;
end

Trim_ITER.CL_delta_e= CL_delta_e;
Trim_ITER.CD_delta_e= CD_delta_e;
Trim_ITER.CM_delta_e= CM_delta_e;

Stab_Der.CL_delta_e= CL_delta_e;
Stab_Der.CD_delta_e= CD_delta_e;
Stab_Der.CM_delta_e= CM_delta_e;

Stab_Der.CL_delta_ail= CL_delta_ail;
Stab_Der.CD_delta_ail= CD_delta_ail;
Stab_Der.CM_delta_ail= CM_delta_ail;

Stab_Der.CL_delta_flap= CL_delta_flap;
Stab_Der.CD_delta_flap= CD_delta_flap;
Stab_Der.CM_delta_flap= CM_delta_flap;

% Control derivatives
CX_delta_e = -CD_delta_e;
CZ_delta_e = -CL_delta_e;

Stab_Der.CX_delta_e = CX_delta_e;
Stab_Der.CZ_delta_e = CZ_delta_e;

%% Aircraft type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
switch AC_type
    case 1 % AC_type = 1 - flying wing
        CL_delta = CL_delta_e;
        CD_delta = CD_delta_e;
        CM_delta = CM_delta_e;
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        CL_delta = CL_delta_e;
        CD_delta = CD_delta_e;
        CM_delta = CM_delta_e;
    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        CL_delta = CL_delta_can + CL_delta_e;
        CD_delta = CD_delta_can + CD_delta_e;
        CM_delta = CM_delta_can + CM_delta_e;
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        CL_delta = CL_delta_rv;
        CD_delta = CD_delta_rv;
        CM_delta = CM_delta_rv;
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        CL_delta = CL_delta_can + CL_delta_rv;
        CD_delta = CD_delta_can + CD_delta_rv;
        CM_delta = CM_delta_can + CM_delta_rv;
end

% A_mat = [CL_alpha_ac, CL_delta; CM_alpha_ac, CM_delta];
% B_mat = [w_T0/(q_inf*S_ref) - CL0_ac; - CM0_ac]
% trim = B_mat\A_mat
% trim_deg = trim*R2D
% pause
%% Resolution of Trim Conditions
% num_delta_e = w_T0/(q_inf*S_ref) - CL0_ac + CL_alpha_ac*CM0_ac/CM_alpha_ac;
% den_delta_e = -CL_alpha_ac*CM_delta/CM_alpha_ac + CL_delta;
% trim_delta_e = num_delta_e/den_delta_e;
% trim_delta_e_deg = trim_delta_e*R2D;
% trim_alpha = - (CM_delta_e*trim_delta_e + CM0_ac)/CM_alpha_ac;
% trim_alpha_deg = trim_alpha*R2D;

CL_needed = w_T0/(q_inf*S_ref);
trim_alpha = (CL_needed*CM_delta - CL0_ac*CM_delta + CL_delta*CM0_ac)/(CL_alpha_ac*CM_delta - CL_delta*CM_alpha_ac);
trim_delta_e = -(CL_needed*CM_alpha_ac - CL0_ac*CM_alpha_ac + CL_alpha_ac*CM0_ac)/(CL_alpha_ac*CM_delta - CL_delta*CM_alpha_ac);
trim_alpha_deg = trim_alpha*R2D;
trim_delta_e_deg = trim_delta_e*R2D;

TRIM_RESULTS.X_NP = X_NP;
TRIM_RESULTS.CL0_ac = CL0_ac;
TRIM_RESULTS.CL_alpha_ac = CL_alpha_ac;
TRIM_RESULTS.CL_delta = CL_delta;
TRIM_RESULTS.CM0_ac = CM0_ac;
TRIM_RESULTS.CM_alpha_ac = CM_alpha_ac;
TRIM_RESULTS.CM_delta = CM_delta;
TRIM_RESULTS.trim_alpha = trim_alpha;
TRIM_RESULTS.trim_delta_e = trim_delta_e;
TRIM_RESULTS.trim_alpha_deg = trim_alpha_deg;
TRIM_RESULTS.trim_delta_e_deg = trim_delta_e_deg;
% TRIM_RESULTS.SM_actual = SM_actual;
TRIM_RESULTS.SM = SM;
TRIM_RESULTS.delta_T = Propulsion.delta_T;

% Euler angle
theta1 = trim_alpha;

% Calculo de valores totales
% Total Lift coefficient
CL_needed = w_T0/(q_inf*S_ref);
CL = CL0_ac + CL_alpha_ac*trim_alpha + CL_delta*trim_delta_e;
CD_trim_delta = abs(CD_delta*trim_delta_e);
CD_Total = C_D0 + C_D1*CL + C_D2*CL^2 + CD_trim_delta;  %polar aeronave

%% Aircraft Type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
switch AC_type
    case 1
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL_alpha_w1*trim_alpha;
        CL_Total = CL_w1;
        Trim_ITER.CL_w1 = CL_w1;
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
        
    case 2
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL_alpha_w1*trim_alpha;
        CL_HTP = CL0_HTP_e_corrected + CL_alpha_HTP*(i_w2 - eps_w2) + CL_alpha_HTP*(downwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_HTP;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_HTP = CL_HTP;
        
        % angle of attack of HTP
        alpha_HTP = ((i_w2 - eps_w2) + (downwash)*trim_alpha);
        alpha_HTP_deg = alpha_HTP*R2D;
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
    case 3
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL_alpha_w1*trim_alpha;
        CL_HTP = CL0_HTP_e_corrected + CL_alpha_HTP*(i_w2 - eps_w2) + CL_alpha_HTP*(downwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_can = CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can) + CL_alpha_can(upwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_HTP + CL_can;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_HTP = CL_HTP;
        Trim_ITER.CL_can = CL_can;
        
        % angle of attack of HTP
        alpha_HTP = ((i_w2 - eps_w2) + (downwash)*trim_alpha);
        alpha_HTP_deg = alpha_HTP*R2D;
        % angle of attack of Can
        alpha_can = ((i_can + eps_can) + (upwash)*trim_alpha);
        alpha_can_deg = alpha_HTP*R2D;
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
        
    case 4
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL_alpha_w1*trim_alpha;
        CL_Vee = CL0_Vee_e_corrected + CL_alpha_Vee*(i_w2 - eps_w2) + CL_alpha_Vee*(downwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_Vee;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_Vee = CL_Vee;
        % angle of attack of HTP
        alpha_Vee = ((i_w2 - eps_w2) + (downwash)*trim_alpha);
        alpha_Vee_deg = alpha_Vee*R2D;
        
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_rudvtr_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
    case 5
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_corrected*i_w1 + CL_alpha_wb_w1*trim_alpha;
        CL_HTP = CL0_Vee_e_corrected + CL_alpha_wb_Vee*(i_w2 - eps_w2) + CLalpha_Vee_e_pw*(downwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_can = CL0_can_e_corrected + CLalpha_can*(i_can + eps_can) + CLalpha_can_e_pw*(upwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_HTP + CL_can;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_Vee = CL_Vee;
        Trim_ITER.CL_can = CL_can;
        
        % angle of attack of HTP
        alpha_Vee = ((i_w2 - eps_w2) + (downwash)*trim_alpha);
        alpha_Vee_deg = alpha_Vee*R2D;
        % angle of attack of Can
        alpha_can = ((i_can + eps_can) + (upwash)*trim_alpha);
        alpha_can_deg = alpha_Vee*R2D;
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_rudvtr_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
end

% Sotore DATA
Trim_ITER.CL = CL;
Trim_ITER.CL_needed = CL_needed;
Trim_ITER.CD_trim_delta = CD_trim_delta;
Trim_ITER.CL_Total = CL_Total;
Trim_ITER.CD_Total = CD_Total;

x_offset_engine = Geo_tier.x_eng_xbar - TRIM_RESULTS.x_XCG ;
Trim_ITER.x_offset_engine = x_offset_engine;

% Total CM - Checking that is zero

CM = CM0_ac + CM_alpha_ac*trim_alpha + CM_delta*trim_delta_e;
TRIM_RESULTS.CM = CM;
Stab_Der.CM = CM;

%% Most Forward XCG
SM_min = 0.10;
x_XCG_rwd = X_NP - (-CL_alpha_ac*SM_min - (CMTalpha + C_Malpha_fus))*(cmac_w1/(-CL_alpha_ac)); % x_XCG with fuselage contribution
% Storing DATA
TRIM_RESULTS.x_XCG_fwd = x_XCG_fwd;
TRIM_RESULTS.x_XCG_rwd = x_XCG_rwd;

Stab_Der.CL_Total = CL_Total;

%% Conversion of variables
if VTP == 1
    % Arm from Zcg to Zac_w2
    l_zcg_VTP = z_zbar_VTP - z_XCG;
    l_xcg_VTP = x_xbar_VTP - x_XCG;
    z_v = (l_zcg_VTP*cos(trim_alpha) - l_xcg_VTP*sin(trim_alpha));
end

if HTP == 1 | Vee == 1
    % Arm from Zcg to Zac_w2
    l_xcg_w2 = x_xbar_w2 - x_XCG;
    l_zcg_w2 = z_zbar_w2 - z_XCG;
    z_v = (l_zcg_w2*cos(trim_alpha) - l_xcg_w2*sin(trim_alpha));
end

if Can == 1
    % Arm from Zcg to Zac_can
    l_xcg_can = x_XCG - x_xbar_can;
    l_zcg_can = z_XCG - z_zbar_can;
    z_v = (l_zcg_can*cos(trim_alpha) - l_xcg_can*sin(trim_alpha));
end

model_conversion = 1;
if model_conversion == 1
    % oswald factor
    taper_e = 1;
    a1 = 0.0004;
    a2 = -0.0080;
    a3 = 0.0501;
    a4 =0.8642;
    lambda1 = AR_w1_e*taper_e/(cos(Lambda_LE_w1_e));
    R = a1*lambda1^3 + a2*lambda1^2 + a3*lambda1 + a4;
    e_e = 1.1*CLalpha_w1_e/(R*CLalpha_w1_e + (1-R)*pi*AR_w1_e);
    
    modelo.general.qinf = q_inf;
    modelo.general.CL = CL;
    modelo.general.CL_w = CL_w1;
    %     modelo.general.CL_h = CL_HTP;
    modelo.general.mtow = m_TOW;
    modelo.general.w_w0 = 1;
    modelo.general.Sref = S_ref;
    modelo.general.Minf = Mach;
    modelo.general.Vinf = V;
    modelo.general.rhoinf = rho;
    modelo.general.Xcg = x_XCG;
    modelo.general.Zcg = z_XCG;
    modelo.general.CLa_WB = CL_alpha_wb_w1;
    modelo.general.downwash = deps_dalpha_Vee;
    modelo.general.upwash = deps_dalpha_can;
    modelo.general.l_fus = Geo_tier.l_fus;
    
    modelo.conversion.m22ft2 = m22ft2;
    modelo.conversion.m2ft = m2ft;
    modelo.propulsion.engines = n_eng; %numero de motores
    modelo.propulsion.blades = 2; %numero de palas
    modelo.propulsion.beta = 0.25; %[deg]
    modelo.propulsion.X = x_d_T;
    modelo.propulsion.wR30 = w_R_30;
    modelo.propulsion.wR60 = w_R_60;
    modelo.propulsion.wR90 = w_R_90;
    modelo.propulsion.D = D_prop;
    
    modelo.ala.S = S_w1;
    modelo.ala.AR = AR_w1;
    modelo.ala.ARwe = AR_w1_e;
    modelo.ala.XLE = x_w1_LE;
    modelo.ala.Xca = x_xbar_w1;
    %vertical distance between the fuselage line and wing root (positive if wing above center line)
    % modelo.ala.Zca1 = z_w1_LE + ybar_w*tan(Sigma_w*D2R);
    modelo.ala.Zca1 = z_zbar_w1;
    %vertical distance between the fuselage line and wing root (positive if wing below center line)
    modelo.ala.Zca = - modelo.ala.Zca1;
    modelo.ala.LAMc2 = Lambda_c2_w1;
    modelo.ala.LAMc4 = Lambda_c4_w1;
    modelo.ala.diedro = dihedral_w1;
    modelo.ala.b = b_w1;
    modelo.ala.TR = lambda_w1;
    modelo.ala.TR_we = lambda_w1_e;
    modelo.ala.xca = xbar_w1;
    modelo.ala.le_y = (b_w1/2)*tan(Lambda_LE_w1);
    modelo.ala.ct = cT_w1;
    modelo.ala.MAC = cmac_w1;
    modelo.ala.MAC_e = cmac_w1_e;
    
    y_1R_y1_ail = Geo_tier.y_1R_y1_ail;
    y_1R_y2_ail = Geo_tier.y_1R_y2_ail;
    cf_ail = Geo_tier.cf_ail;
    
    modelo.ala.y1_b2 = y_1R_y2_ail/(b_w1_e/2); %outer aleron coordinate
    modelo.ala.y0_b2 = y_1R_y1_ail/(b_w1_e/2); %inner aleron coordinate
    
    % modelo.ala.cm_c = c_ail/cmac_w1;
    modelo.ala.cm_c = cf_ail;
    modelo.ala.t_c = 0.15;
    %% OJO assumes is the same as 3D
    modelo.ala.Cla = CL_alpha_wb_w1;
    modelo.ala.CLa = CL_alpha_wb_w1;
    modelo.ala.CLa_we = CL_alpha_wb_w1;
    modelo.ala.oswald = e_e;
    modelo.ala.CD0 = C_D0;
    modelo.ala.CD1 = C_D1;
    modelo.ala.CD2 = C_D2;
    
    modelo.fuselaje.x = x_Area_body;
    modelo.fuselaje.D_x = Area_body;
    modelo.fuselaje.vol = Vol_TOT;
    modelo.fuselaje.CLa = CLalpha_fus;
    modelo.fuselaje.D = h_Area_b_max;
    modelo.fuselaje.l = length_fus;
    modelo.fuselaje.Sside = Area_side;
    modelo.fuselaje.W = Area_b_max;
    
    % Real area, span and AR for the VTAIL
    S_w1_s = Geo_tier.S_w1_s;
    S_w2_s = Geo_tier.S_w2_s;
    b_w1_s = Geo_tier.b_w1_s;
    b_w2_s = Geo_tier.b_w2_s;
    AR_w1_s = Geo_tier.AR_w1_s;
    AR_w2_s = Geo_tier.AR_w2_s;
    
    S_w2_pv = Geo_tier.S_w2_pv;
    S_w2_ph = Geo_tier.S_w2_ph;
    
    % Defines equivalent parameters for VTAIL as projection of HTP and VTP
    
    switch AC_type
        case 1 % AC_type = 1 - flying wing
            
        case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
            if VTP == 1
                % VTP
                modelo.vertical.S = Geo_tier.S_VTP;
                modelo.vertical.eta = eta_VTP_no_afe;
                modelo.vertical.Cla = CLalpha_VTP;
                modelo.vertical.CLa = CLalpha_VTP;
                modelo.vertical.b = (b_VTP_s);
                modelo.vertical.Xca = x_xbar_VTP;
                modelo.vertical.Zca = l_zcg_VTP;
                modelo.vertical.AR = AR_VTP;
                modelo.vertical.TR = lambda_VTP;
                modelo.vertical.cm_c = 0.25*cmac_VTP;
                modelo.vertical.t_c = 0.12; %naca 0012
                modelo.vertical.LAMc2 = Lambda_c2_VTP;
                modelo.vertical.LAMc4 = Lambda_c4_VTP;
                modelo.vertical.diedro = dihedral_VTP;
                modelo.vertical.xca = xbar_VTP;
                modelo.vertical.le_y = (b_VTP)*tan(Lambda_LE_VTP);
                modelo.vertical.ct = cT_VTP;
                modelo.vertical.MAC = cmac_VTP;
                modelo.vertical.MAC_e = cmac_VTP_e;
                modelo.vertical.y0_b2 = K_y1_rudder_VTP;
                modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            end
            
            if HTP == 1
                modelo.horizontal.S = Geo_tier.S_w2;
                modelo.horizontal.eta = eta_w2_no_afe;
                modelo.horizontal.Cla = CL_alpha_HTP;
                modelo.horizontal.CLa = CL_alpha_HTP;
                modelo.horizontal.b = (b_w2_s);
                modelo.horizontal.Xca = x_xbar_w2;
                modelo.horizontal.Zca = l_zcg_w2;
                modelo.horizontal.AR = AR_w2;
                modelo.horizontal.TR = lambda_w2;
                modelo.horizontal.cm_c = 0.25*cmac_w2;
                modelo.horizontal.t_c = 0.12; %naca 0012
                modelo.horizontal.LAMc2 = Lambda_c2_w2;
                modelo.horizontal.LAMc4 = Lambda_c4_w2;
                modelo.horizontal.diedro = dihedral_w2;
                modelo.horizontal.xca = xbar_w2;
                modelo.horizontal.le_y = (b_w2/2)*tan(Lambda_LE_w2);
                modelo.horizontal.ct = cT_w2;
                modelo.horizontal.MAC = cmac_w2;
                modelo.horizontal.MAC_e = cmac_w2_e;
                modelo.horizontal.y0_b2 = K_y1_ele_w2;
                modelo.horizontal.y1_b2 = K_y2_ele_w2;
                %                 modelo.horizontal.y0_b2 = z_1R_y1_ele;
                %                 modelo.horizontal.y1_b2 = z_1R_y2_ele;
                
            end
            
            
        case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
            
            if VTP == 1
                % VTP
                modelo.vertical.S = Geo_tier.S_VTP;
                modelo.vertical.eta = eta_VTP_no_afe;
                modelo.vertical.Cla = CLalpha_VTP;
                modelo.vertical.CLa = CLalpha_VTP;
                modelo.vertical.b = (b_VTP_s);
                modelo.vertical.Xca = x_xbar_VTP;
                modelo.vertical.Zca = l_zcg_VTP;
                modelo.vertical.AR = AR_VTP;
                modelo.vertical.TR = lambda_VTP;
                modelo.vertical.cm_c = 0.25*cmac_VTP;
                modelo.vertical.t_c = 0.12; %naca 0012
                modelo.vertical.LAMc2 = Lambda_c2_VTP;
                modelo.vertical.LAMc4 = Lambda_c4_VTP;
                modelo.vertical.diedro = dihedral_VTP;
                modelo.vertical.xca = xbar_VTP;
                modelo.vertical.le_y = (b_VTP)*tan(Lambda_LE_VTP);
                modelo.vertical.ct = cT_VTP;
                modelo.vertical.MAC = cmac_VTP;
                modelo.vertical.MAC_e = cmac_VTP_e;
                
                modelo.vertical.y0_b2 = K_y1_rudder_VTP;
                modelo.vertical.y1_b2 = K_y2_rudder_VTP;
                %                 modelo.vertical.y0_b2 = z_1R_y1_rudder;
                %                 modelo.vertical.y1_b2 = z_1R_y2_rudder;
                
                
                
            end
            
            if HTP == 1 | Vee == 1
                modelo.horizontal.S = Geo_tier.S_w2;
                modelo.horizontal.eta = eta_w2_no_afe;
                modelo.horizontal.Cla = CLalpha_wb_HTP;
                modelo.horizontal.CLa = CLalpha_wb_HTP;
                modelo.horizontal.b = (b_w2_s);
                modelo.horizontal.Xca = x_xbar_w2;
                modelo.horizontal.Zca = l_zcg_w2;
                modelo.horizontal.AR = AR_w2;
                modelo.horizontal.TR = lambda_w2;
                modelo.horizontal.cm_c = 0.25*cmac_w2;
                modelo.horizontal.t_c = 0.12; %naca 0012
                modelo.horizontal.LAMc2 = Lambda_c2_w2;
                modelo.horizontal.LAMc4 = Lambda_c4_w2;
                modelo.horizontal.diedro = dihedral_w2;
                modelo.horizontal.xca = xbar_w2;
                modelo.horizontal.le_y = (b_w2/2)*tan(Lambda_LE_w2);
                modelo.horizontal.ct = cT_w2;
                modelo.horizontal.MAC = cmac_w2;
                modelo.horizontal.MAC_e = cmac_w2_e;
                %                 modelo.horizontal.y0_b2 = z_1R_y1_ele;
                %                 modelo.horizontal.y1_b2 = z_1R_y2_ele;
                modelo.horizontal.y0_b2 = K_y1_ele_w2;
                modelo.horizontal.y1_b2 = K_y2_ele_w2;
                
            end
            
        case 4 % AC_type = 4 - 2 surface: wing + V-tail
            % Estimates properties as proyectos from virtual HTP and virtual
            % VTP
            modelo.vertical.S = Geo_tier.S_w2_pv;
            modelo.vertical.eta = eta_w2_no_afe;
            modelo.vertical.Cla = CL_alpha_wb_Vee*(sin(dihedral_w2))^2;
            modelo.vertical.Cla = CL_alpha_wb_Vee;
            modelo.vertical.CLa = CL_alpha_wb_Vee*(sin(dihedral_w2))^2;
            modelo.vertical.CLa = CL_alpha_wb_Vee;
            modelo.vertical.b = (b_w2_s/2)*sin(dihedral_w2);
            modelo.vertical.Xca = x_xbar_w2;
            modelo.vertical.Zca = l_zcg_w2;
            modelo.vertical.AR = (((b_w2_s/2)*sin(dihedral_w2))^2)/(S_w2_pv/2);
            modelo.vertical.TR = lambda_w2;
            modelo.vertical.cm_c = 0.25*cmac_w2;
            modelo.vertical.t_c = 0.12; %naca 0012
            modelo.vertical.LAMc2 = Lambda_c2_w2;
            modelo.vertical.LAMc4 = Lambda_c4_w2;
            modelo.vertical.diedro = dihedral_w2;
            modelo.vertical.xca = xbar_w2;
            modelo.vertical.le_y = (b_w2/2)*tan(Lambda_LE_w2);
            modelo.vertical.ct = cT_w2;
            modelo.vertical.MAC = cmac_w2;
            modelo.vertical.MAC_e = cmac_w2_e;
            %             modelo.vertical.y0_b2 = K_y1_rudder_VTP;
            %             modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            modelo.vertical.y0_b2 = z_1R_y1_rudvtr;
            modelo.vertical.y1_b2 = z_1R_y2_rudvtr;
            
            modelo.horizontal.S = Geo_tier.S_w2_ph;
            modelo.horizontal.eta = eta_w2_no_afe;
            modelo.horizontal.Cla = CL_alpha_wb_Vee*(cos(dihedral_w2))^2;
            modelo.horizontal.CLa = CL_alpha_wb_Vee*(cos(dihedral_w2))^2;
            modelo.horizontal.b = (b_w2_s)*cos(dihedral_w2);
            modelo.horizontal.Xca = x_xbar_w2;
            modelo.horizontal.Zca = l_zcg_w2;
            modelo.horizontal.AR = (((b_w2_s)*cos(dihedral_w2))^2)/(S_w2_ph);
            modelo.horizontal.TR = lambda_w2;
            modelo.horizontal.cm_c = 0.25*cmac_w2;
            modelo.horizontal.t_c = 0.12; %naca 0012
            modelo.horizontal.LAMc2 = Lambda_c2_w2;
            modelo.horizontal.LAMc4 = Lambda_c4_w2;
            modelo.horizontal.diedro = dihedral_w2;
            modelo.horizontal.xca = xbar_w2;
            modelo.horizontal.le_y = (b_w2/2)*tan(Lambda_LE_w2);
            modelo.horizontal.ct = cT_w2;
            modelo.horizontal.MAC = cmac_w2;
            modelo.horizontal.MAC_e = cmac_w2_e;
            modelo.horizontal.y0_b2 = z_1R_y1_rudvtr;
            modelo.horizontal.y1_b2 = z_1R_y2_rudvtr;
            modelo.horizontal.y0_b2 = K_y1_rudvtr_w2;
            modelo.horizontal.y1_b2 = K_y2_rudvtr_w2;
            
            modelo.vee.S = Geo_tier.S_w2_s;
            modelo.vee.Cla = CL_alpha_wb_Vee*(cos(dihedral_w2))^2;
            modelo.vee.CLa = CL_alpha_wb_Vee*(cos(dihedral_w2))^2;
            modelo.vee.b = b_w2_s;
            modelo.vee.Xca = x_xbar_w2;
            modelo.vee.Zca = l_zcg_w2;
            modelo.vee.dihedral_w2 = dihedral_w2;
            
        case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
            % Estimates properties as proyectos from virtual HTP and virtual
            % VTP
            modelo.vertical.S = Geo_tier.S_w2_pv;
            modelo.vertical.eta = eta_w2_no_afe;
            modelo.vertical.Cla = CL_alpha_wb_Vee*(sin(dihedral_w2))^2;
            modelo.vertical.Cla = CL_alpha_wb_Vee;
            modelo.vertical.CLa = CL_alpha_wb_Vee*(sin(dihedral_w2))^2;
            modelo.vertical.CLa = CL_alpha_wb_Vee;
            modelo.vertical.b = (b_w2_s/2)*sin(dihedral_w2);
            modelo.vertical.Xca = x_xbar_w2;
            modelo.vertical.Zca = l_zcg_w2;
            modelo.vertical.AR = (((b_w2_s/2)*sin(dihedral_w2))^2)/(S_w2_pv/2);
            modelo.vertical.TR = lambda_w2;
            modelo.vertical.cm_c = 0.25*cmac_w2;
            modelo.vertical.t_c = 0.12; %naca 0012
            modelo.vertical.LAMc2 = Lambda_c2_w2;
            modelo.vertical.LAMc4 = Lambda_c4_w2;
            modelo.vertical.diedro = dihedral_w2;
            modelo.vertical.xca = xbar_w2;
            modelo.vertical.le_y = (b_w2/2)*tan(Lambda_LE_w2);
            modelo.vertical.ct = cT_w2;
            modelo.vertical.MAC = cmac_w2;
            modelo.vertical.MAC_e = cmac_w2_e;
            modelo.vertical.y0_b2 = K_y1_rudder_VTP;
            modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            %             modelo.vertical.y0_b2 = z_1R_y1_rudvtr;
            %             modelo.vertical.y1_b2 = z_1R_y2_rudvtr;
            
            modelo.horizontal.S = Geo_tier.S_w2_ph;
            modelo.horizontal.eta = eta_w2_no_afe;
            modelo.horizontal.Cla = CL_alpha_wb_Vee*(cos(dihedral_w2))^2;
            modelo.horizontal.CLa = CL_alpha_wb_Vee*(cos(dihedral_w2))^2;
            modelo.horizontal.b = (b_w2_s)*cos(dihedral_w2);
            modelo.horizontal.Xca = x_xbar_w2;
            modelo.horizontal.Zca = l_zcg_w2;
            modelo.horizontal.AR = (((b_w2_s)*cos(dihedral_w2))^2)/(S_w2_ph);
            modelo.horizontal.TR = lambda_w2;
            modelo.horizontal.cm_c = 0.25*cmac_w2;
            modelo.horizontal.t_c = 0.12; %naca 0012
            modelo.horizontal.LAMc2 = Lambda_c2_w2;
            modelo.horizontal.LAMc4 = Lambda_c4_w2;
            modelo.horizontal.diedro = dihedral_w2;
            modelo.horizontal.xca = xbar_w2;
            modelo.horizontal.le_y = (b_w2/2)*tan(Lambda_LE_w2);
            modelo.horizontal.ct = cT_w2;
            modelo.horizontal.MAC = cmac_w2;
            modelo.horizontal.MAC_e = cmac_w2_e;
            %             modelo.horizontal.y0_b2 = z_1R_y1_ele;
            %             modelo.horizontal.y1_b2 = z_1R_y2_ele;
            modelo.horizontal.y0_b2 = K_y1_ele_w2;
            modelo.horizontal.y1_b2 = K_y2_ele_w2;
            
            modelo.vee.S = Geo_tier.S_w2_s;
            modelo.vee.Cla = CL_alpha_wb_Vee*(cos(dihedral_w2))^2;
            modelo.vee.CLa = CL_alpha_wb_Vee*(cos(dihedral_w2))^2;
            modelo.vee.b = b_w2_s;
            modelo.vee.Xca = x_xbar_w2;
            modelo.vee.Zca = l_zcg_w2;
            modelo.vee.dihedral_w2 = dihedral_w2;
    end
    
end


if only_trim == 0
    
    Trim_ITER.SM_des = SM_des;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%Pitch Derivatives%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%CLq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CLalpha_fus = 2*f_k2_k1*Area_b_max/(Vol_TOT^(2/3));  %eq 4.495 PAMADI
    CLalpha_Bp = CLalpha_fus*(Vol_TOT^(2/3))/Area_b_max; %eq 4.494 PAMADI
    CLq_B = 2*CLalpha_Bp*(1 + x_XCG/length_fus);           %eq 4.493 PAMADI
    CLq_fus = CLq_B*((length_fus*Area_b_max)/(S_w1*cmac_w1)); %eq 4.488 PAMADI
    
    if W1 == 1
        x_xbar_w1_e = Geo_tier.x_xbar_w1_e;
        B_prandalt = sqrt(1 - (Mach^2)*(cos(Lambda_c4_w1))^2);     %eq 4.510 PAMADI
        chi = (x_xbar_w1_e - x_XCG)/cmac_w1_e; %eq 4.490 PAMADI
        % CLalpha_wb_w1 already includes (KWB_w1 + KBW_w1)
        CLq_e = CL_alpha_wb_w1*(0.5 + 2*chi);   %eq 4.489 PAMADI
        CLq_w1 = ((S_w1_e*cmac_w1_e)/(S_ref*cmac_w1))*CLq_e;        
        % ROSKAM AAA
        % Curve lift slope with no flaps
        CL_alpha_clean = CL_alpha_wb_w1;
        % increment of wing-fuselage lift curve slope due to power
        DeltaCL_alpha_clean = 0;
        CLq_w_M0 = (CL_alpha_clean + DeltaCL_alpha_clean)*(0.5 + 2*chi);
        int_CLq_w1 = (AR_w1 + 2*cos(Lambda_c4_w1))/(AR_w1*B_prandalt + 2*cos(Lambda_c4_w1));
        CLq_w1 = ((S_w1_e*cmac_w1_e)/(S_ref*cmac_w1))*int_CLq_w1*CLq_w_M0;
        % W1 and w2 with no dynamic pressure correction
        CLq_w1_pw = eta_w1_afe_S_w1_afe_S_ref*CLq_w1 + eta_w1_no_afe_S_w1_no_afe_S_ref*CLq_w1;
    else
        CLq_w1_pw = 0;
    end
    
    if HTP == 1
        x_xbar_w2_e = Geo_tier.x_xbar_w2_e;
        CLq_w2 = 2*CLalpha_w2_e*Vbar_h; %eq 10.72 ROSKAM
        % W1 and w2 with no dynamic pressure correction
        CLq_HTP_pw = eta_w2_afe_S_w2_afe_S_ref*CLq_w2 + eta_w2_no_afe_S_w2_no_afe_S_ref*CLq_w2;
    else
        CLq_HTP_pw = 0;
    end
    
    if Vee == 1
        x_xbar_w2_e = Geo_tier.x_xbar_w2_e;
        CLq_w2 = 2*CLalpha_w2_e*Vbar_vee; %eq 10.72 ROSKAM
        % W1 and w2 with no dynamic pressure correction
        CLq_Vee_pw = eta_w2_afe_S_w2_afe_S_ref*CLq_w2 + eta_w2_no_afe_S_w2_no_afe_S_ref*CLq_w2;
    else
        CLq_Vee_pw = 0;
    end
    
    if Can == 1
        x_xbar_can_e = Geo_tier.x_xbar_can_e;
        CLq_can = -2*CLalpha_can_e*Vbar_can; %eq 10.72 ROSKAM
        % Can with no dynamic pressure correction
        CLq_can_pw = eta_can_afe_S_can_afe_S_ref*CLq_can + eta_can_no_afe_S_can_no_afe_S_ref*CLq_can;
    else
        CLq_can_pw = 0;
    end
    
    % NAcelle contribution
    if Nac== 1
        CLq_nac_pw = 0;
    else
        CLq_nac_pw = 0;
    end
    
    % Total Derivative CLq
    CLq = CLq_w1_pw + CLq_fus + CLq_HTP_pw + CLq_Vee_pw + CLq_can_pw + CLq_nac_pw;
    Stab_Der.CLq = CLq;
    % Stores Derivatives per parts
    Stab_Der_parts.CL_q = CLq;
    Stab_Der_parts.CL_q_w1 = CLq_w1_pw;
    Stab_Der_parts.CL_q_fus = CLq_fus;
    Stab_Der_parts.CL_q_HTP = CLq_HTP_pw;
    Stab_Der_parts.CL_q_Vee = CLq_Vee_pw;
    Stab_Der_parts.CL_q_can = CLq_can_pw;
    Stab_Der_parts.CL_q_nac = CLq_nac_pw;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%Cmq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fuselage contribution
    % Centroid fuselage
    int_Sb1 = @(x) interp1(x_Area_body,Area_body,x).*(x);               %eq 4.514 PAMADI
    x_c_fus = (1/Vol_TOT)*quad(@(x) int_Sb1(x),0,length_fus);           %eq 4.514 PAMADI
    x_m = x_XCG;
    xm_1 = x_m/length_fus;                      %eq 4.513 PAMADI
    xc_1 = x_c_fus/length_fus;                  %eq 4.513 PAMADI
    VB_1 = Vol_TOT/(Area_b_max*length_fus);     %eq 4.513 PAMADI
    
    % C�lculo del centro aerodinamico del fuselaje
    % Derivada del �rea:
    dSdX = diff(Area_body)./diff(x_Area_body);
    %  Calculates where fuselage ceases to be potential
    i=1;
    while dSdX(i) > 0
        x_0_v = i;
        i=i+1;
    end
    x_0 = x_Area_body(x_0_v+1);
    int_Sb_CMq = @(x) interp1(x_Area_body(1:end-1),dSdX,x).*(x_m - x);    %eq 4.516 PAMADI
    CMalpha_B = (2*f_k2_k1/Vol_TOT)*quad(@(x) int_Sb_CMq(x),0,x_0);       %eq 4.516 PAMADI
    CMalpha_p = CMalpha_B*VB_1;                                           %eq 4.515 PAMADI
    CMq_B = 2*CMalpha_p*((1-xm_1)^2 - VB_1*(xc_1-xm_1))/(1-xm_1-VB_1);    %eq 4.512 PAMADI
    
    % W1 and w2 with no dynamic pressure correction %eq 4.502 PAMADI
    % Ojo revisar por la contribuci�n de la presi�n din�mica
    CMq_fus = CMq_B*(Area_b_max/S_w1)*((length_fus/cmac_w1)^2);                   %eq 4.502 PAMADI
    
    if W1 == 1
        c1 = AR_w1_e^3*(tan(Lambda_c4_w1))^2;          %eq 4.505 PAMADI
        c2 = 3/B_prandalt;                           %eq 4.506 PAMADI
        c3 = AR_w1_e*B_prandalt + 6*cos(Lambda_c4_w1); %eq 4.507 PAMADI
        c4 = AR_w1_e + 6*cos(Lambda_c4_w1);            %eq 4.508 PAMADI
        c5 = AR_w1_e + 2*cos(Lambda_c4_w1);            %eq 5.509 PAMADI
        % Assumes that the CLalpha_w1 in 2D is the same as 3D
        CLalpha_w1_2D = CLalpha_w1;
        CMq_e_M02 = -0.7*CLalpha_w1_2D*cos(Lambda_c4_w1)*((AR_w1_e*(0.5*chi + 2*chi^2))/c5 +...
            (c1/(24*c4)) + 1/8);                    %eq 4.504 PAMADI
        CMq_e = CMq_e_M02*(c1/c3 + c2)/(c1/c4 + 3); %eq 4.503 PAMADI
        
        % W1 and w2 with no dynamic pressure correction %eq 4.502 PAMADI
        % Ojo revisar por la contribuci�n de la presi�n din�mica
        CMq_w1 = (KWB_w1 + KBW_w1)*((cmac_w1_e/cmac_w1)^2)*CMq_e; %eq 4.502 PAMADI
        CMq_w1_pw = eta_w1_afe_S_w1_afe_S_ref*CMq_w1 + eta_w1_no_afe_S_w1_no_afe_S_ref*CMq_w1;
    else
        CMq_w1_pw = 0;
    end
    
    if HTP == 1
        %% Revisar OJO - contribuici�n de presi�n din�mica en ala w1
        CMq_w2 = - 2*CLalpha_w2_e*Vbar_h*((x_xbar_w2 - x_XCG)/cmac_w1);              %eq 10.78 ROSKAM
        CMq_HTP_pw = eta_w2_afe_S_w2_afe_S_ref*CMq_w2 + eta_w2_no_afe_S_w2_no_afe_S_ref*CMq_w2;
    else
        CMq_HTP_pw = 0;
    end
    
    if Vee == 1
        %% Revisar OJO - contribuici�n de presi�pn din�mica en ala w1
        CMq_w2 = - 2*CLalpha_w2_e*Vbar_vee*((x_xbar_w2 - x_XCG)/cmac_w1);              %eq 10.78 ROSKAM
        CMq_Vee_pw = eta_w2_afe_S_w2_afe_S_ref*CMq_w2 + eta_w2_no_afe_S_w2_no_afe_S_ref*CMq_w2;
    else
        CMq_Vee_pw = 0;
    end
    
    if Can == 1
        %% Revisar OJO - contribuici�n de presi�pn din�mica en ala w1
        CMq_w2 = - 2*CLalpha_can_e*Vbar_can*((x_xbar_can - x_XCG)/cmac_w1);              %eq 10.78 ROSKAM
        CMq_can_pw = eta_can_afe_S_can_afe_S_ref*CMq_can + eta_can_no_afe_S_can_no_afe_S_ref*CMq_can;
    else
        CMq_can_pw = 0;
    end
    
    % NAcelle contribution
    if Nac== 1
        CMq_nac_pw = 0;
    else
        CMq_nac_pw = 0;
    end
    
    % Total Derivative CMq
    CMq = CMq_w1_pw + CMq_fus + CMq_HTP_pw + CMq_Vee_pw + CMq_can_pw + CMq_nac_pw;
    % Storing DATA
    Stab_Der.CMq = CMq;
    %     Stab_Der.CMq_w1_pw = CMq_w1_pw;
    %     Stab_Der.CMq_w2_pw = CMq_w2_pw;
    %     Stab_Der.CMq_can_pw = CMq_can_pw;
    %     Stab_Der.CMq_fus = CMq_fus;
    
    % Stores Derivatives per parts
    Stab_Der_parts.CM_q = CMq;
    Stab_Der_parts.CM_q_w1 = CMq_w1_pw;
    Stab_Der_parts.CM_q_fus = CMq_fus;
    Stab_Der_parts.CM_q_HTP = CMq_HTP_pw;
    Stab_Der_parts.CM_q_Vee = CMq_Vee_pw;
    Stab_Der_parts.CM_q_can = CMq_can_pw;
    Stab_Der_parts.CM_q_nac = CMq_nac_pw;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%Cdq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The airplane drag-coefficient-due-to-pitch-rate derivative is negligible:

    CDq_w1_pw = 0;
    CDq_fus = 0;
    CDq_HTP_pw = 0;
    CDq_Vee_pw = 0;
    CDq_can_pw = 0;
    CDq_nac_pw = 0;
    % Total Derivative CDq
    CDq = CDq_w1_pw + CDq_fus + CDq_HTP_pw + CDq_Vee_pw + CDq_can_pw + CDq_nac_pw;
    
    % Storing DATA
    Stab_Der.CDq = CDq;
    
    % Stores Derivatives per parts
    Stab_Der_parts.CD_q = CDq;
    Stab_Der_parts.CD_q_w1 = CDq_w1_pw;
    Stab_Der_parts.CD_q_fus = CDq_fus;
    Stab_Der_parts.CD_q_HTP = CDq_HTP_pw;
    Stab_Der_parts.CD_q_Vee = CDq_Vee_pw;
    Stab_Der_parts.CD_q_can = CDq_can_pw;
    Stab_Der_parts.CD_q_nac = CDq_nac_pw;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Alpha punto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% CL_alphapunto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fuselage contribution
    CLalpha_B = 2*f_k2_k1*(Area_b_max/Vol_TOT^(2/3));      %eq 4.532 PAMADI
    CLalpha_B_p = CLalpha_B*(Vol_TOT^(2/3)/(Area_b_max));    %eq 4.531 PAMADI
    CL_alphapunto_B = 2*CLalpha_B_p*(Vol_TOT/(length_fus*Area_b_max));  %eq 4.530 PAMADI
    CL_alphapunto_fus = CL_alphapunto_B*((length_fus*Area_b_max)/(S_ref*cmac_w1));      %eq 4.527 PAMADI
    
    if W1 == 1
        beta = sqrt(1 - Mach^2);
        tau = beta*AR_w1_e;
        % Only valid for tau < 4
        if tau<4
            CL_g = (-pi*AR_w1_e/(2*beta^2))*(0.0013*tau^4 - 0.0122*tau^3 + ...
                0.0317*tau^2 + 0.0186*tau - 0.0004);                %eq 4.529 PAMADI
        else
            CL_g = 0;
        end
        %% OJO! using he distance from apex of wing LE root to Xac
        CL_alphapunto_e1 = 1.5*(xbar_w1/cmac_w1_e)*CLalpha_w1_e; %eq 4.528 PAMADI
        CL_alphapunto_e2 = 3*CL_g;                              %eq 4.528 PAMADI
        % CL_alphapunto_e2 = 0; %correccion
        CL_alphapunto_e = CL_alphapunto_e1 + CL_alphapunto_e2;  %eq 4.528 PAMADI
        
        % W1 and w2 with no dynamic pressure correction
        % Ojo revisar por la contribuci�n de la presi�n din�mica
        CL_alphapunto_w1 = (KWB_w1 + KBW_w1)*(cmac_w1_e/cmac_w1)*CL_alphapunto_e;    %eq 4.527 PAMADI
        CL_alphapunto_w1_pw = eta_w1_afe_S_w1_afe_S_ref*CL_alphapunto_w1 +...
            eta_w1_no_afe_S_w1_no_afe_S_ref*CL_alphapunto_w1;
        
    else
        CL_alphapunto_w1_pw = 0;
    end
    
    if HTP == 1
        % Modified version Pamadi to account for different dynamic pressure
        % W1 and w2 with no dynamic pressure correction
        CL_alphapunto_w2 = 2*CLalpha_w2_e*Vbar_h*deps_dalpha_h; %eq 4.525 PAMADI
        % Revisar OJO - contribuici�n de presi�pn din�mica en ala w2
        CL_alphapunto_HTP_pw = eta_w2_afe_S_w2_afe_S_ref*CL_alphapunto_w2 + ...
            eta_w2_no_afe_S_w2_no_afe_S_ref*CL_alphapunto_w2;
    else
        CL_alphapunto_HTP_pw = 0;
    end
    
    if Vee == 1
        % W1 and w2 with no dynamic pressure correction
        CL_alphapunto_w2 = 2*CLalpha_w2_e*Vbar_vee*deps_dalpha_Vee; %eq 4.525 PAMADI
        % Revisar OJO - contribuici�n de presi�pn din�mica en ala w2
        CL_alphapunto_Vee_pw = eta_w2_afe_S_w2_afe_S_ref*CL_alphapunto_w2 + ...
            eta_w2_no_afe_S_w2_no_afe_S_ref*CL_alphapunto_w2;
    else
        CL_alphapunto_Vee_pw = 0;
    end
    
    if Can == 1
        CL_alphapunto_can = 2*CLalpha_can_e*Vbar_can*deps_dalpha_can; %eq 4.525 PAMADI
        % Revisar OJO - contribuici�n de presi�pn din�mica en ala can
        CL_alphapunto_can_pw = eta_can_afe_S_can_afe_S_ref*CL_alphapunto_can + ...
            eta_can_no_afe_S_can_no_afe_S_ref*CL_alphapunto_can;
    else
        CL_alphapunto_can_pw = 0;
    end
    
    % NAcelle contribution
    if Nac== 1
        CL_alphapunto_nac_pw = 0;
    else
        CL_alphapunto_nac_pw = 0;
    end
    
    % Total Derivative CLalphapunto
    CL_alphapunto = CL_alphapunto_w1_pw + CL_alphapunto_fus + CL_alphapunto_HTP_pw + ...
        CL_alphapunto_Vee_pw + CL_alphapunto_can_pw + CL_alphapunto_nac_pw;
    
    % Storing DATA
    Stab_Der.CL_alphapunto = CL_alphapunto;
    %     Stab_Der.CL_alphapunto_w1_pw = CL_alphapunto_w1_pw;
    %     Stab_Der.CL_alphapunto_w2_pw = CL_alphapunto_w2_pw;
    %     Stab_Der.CL_alphapunto_can_pw = CL_alphapunto_can_pw;
    %     Stab_Der.CL_alphapunto_fus = CL_alphapunto_fus;
    
    % Stores Derivatives per parts
    Stab_Der_parts.CL_alphapunto = CL_alphapunto;
    Stab_Der_parts.CL_alphapunto_w1 = CL_alphapunto_w1_pw;
    Stab_Der_parts.CL_alphapunto_fus = CL_alphapunto_fus;
    Stab_Der_parts.CL_alphapunto_HTP = CL_alphapunto_HTP_pw;
    Stab_Der_parts.CL_alphapunto_Vee = CL_alphapunto_Vee_pw;
    Stab_Der_parts.CL_alphapunto_can = CL_alphapunto_can_pw;
    Stab_Der_parts.CL_alphapunto_nac = CL_alphapunto_nac_pw;
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CM_alphapunto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CM_alphapunto_B = 2*CMalpha_p*((xc_1 - xm_1)/(1 - xm_1 - VB_1))*(Vol_TOT/(length_fus*Area_b_max)); %eq 4.543 PAMADI
    CM_alphapunto_fus = CM_alphapunto_B*((Area_b_max*length_fus^2)/(S_w1*cmac_w1^2));      %eq 4.539 PAMADI
    
    if W1 == 1
        % Only valid for tau < 4
        if tau<4
            CMo_g = (pi*AR_w1_e/(2*beta^2))*(0.0008*tau^4 - 0.0075*tau^3 +...
                0.0185*tau^2 + 0.0128*tau - 0.0003); %eq 4.542 PAMADI
        else
            CMo_g = 0;
        end
        
        %% OJO! using he distance from apex of wing LE root to Xac
        CM_alphapunto_e_pp_1 = -(81/32)*((xbar_w1/cmac_w1_e)^2)*CLalpha_w1_e;        %eq 4.541 PAMADI
        CM_alphapunto_e_pp_2 = (9/2)*CMo_g;                                         %eq 4.541 PAMADI
        CM_alphapunto_e_pp =  CM_alphapunto_e_pp_1 + CM_alphapunto_e_pp_2;          %eq 4.541 PAMADI
        CM_alphapunto_e = CM_alphapunto_e_pp + ((x_xbar_w1_e - x_XCG)/cmac_w1_e)*CL_alphapunto_e; %eq 4.540 PAMADI
        
        % W1 and w2 with no dynamic pressure correction %eq 4.539 PAMADI
        % Ojo revisar por la contribuci�n de la presi�n din�mica
        CM_alphapunto_w1 = (KWB_w1 + KBW_w1)*(cmac_w1_e^2/cmac_w1^2)*CM_alphapunto_e;    %eq 4.539 PAMADI
        CM_alphapunto_w1_pw = eta_w1_afe_S_w1_afe_S_ref*CM_alphapunto_w1 +...
            eta_w1_no_afe_S_w1_no_afe_S_ref*CM_alphapunto_w1;
        
    else
        CM_alphapunto_w1_pw = 0;
    end
    
    if HTP == 1
        CM_alphapunto_w2 = - CL_alphapunto_w2*((x_xbar_w2 - x_XCG)/cmac_w1);      %eq 4.537 PAMADI
        % W1 and w2 with no dynamic pressure correction
        %% Revisar OJO - contribuici�n de presi�n din�mica en ala w1
        CM_alphapunto_HTP_pw = eta_w2_afe_S_w2_afe_S_ref*CM_alphapunto_w2 + ...
            eta_w2_no_afe_S_w2_no_afe_S_ref*CM_alphapunto_w2;
    else
        CM_alphapunto_HTP_pw = 0;
    end
    
    if Vee == 1
        CM_alphapunto_w2 = - CL_alphapunto_w2*((x_xbar_w2 - x_XCG)/cmac_w1);      %eq 4.537 PAMADI
        % W1 and w2 with no dynamic pressure correction
        %% Revisar OJO - contribuici�n de presi�n din�mica en ala w1
        CM_alphapunto_Vee_pw = eta_w2_afe_S_w2_afe_S_ref*CM_alphapunto_w2 + ...
            eta_w2_no_afe_S_w2_no_afe_S_ref*CM_alphapunto_w2;
    else
        CM_alphapunto_Vee_pw = 0;
    end
    
    if Can == 1
        CM_alphapunto_can = - CL_alphapunto_can*((x_xbar_can - x_XCG)/cmac_w1);      %eq 4.537 PAMADI
        % W1 and w2 with no dynamic pressure correction
        %% Revisar OJO - contribuici�n de presi�n din�mica en ala w1
        CM_alphapunto_can_pw = eta_can_afe_S_can_afe_S_ref*CM_alphapunto_can + ...
            eta_can_no_afe_S_can_no_afe_S_ref*CM_alphapunto_can;
    else
        CM_alphapunto_can_pw = 0;
    end
    
    % NAcelle contribution
    if Nac== 1
        CM_alphapunto_nac_pw = 0;
    else
        CM_alphapunto_nac_pw = 0;
    end
    
    % Total Derivative CMalphapunto
    CM_alphapunto = CM_alphapunto_w1_pw + CM_alphapunto_fus + CM_alphapunto_HTP_pw + ...
        CM_alphapunto_Vee_pw + CM_alphapunto_can_pw + CM_alphapunto_nac_pw;
    % Storing DATA
    Stab_Der.CM_alphapunto = CM_alphapunto;
    %     Stab_Der.CM_alphapunto_w1_pw = CM_alphapunto_w1_pw;
    %     Stab_Der.CM_alphapunto_w2_pw = CM_alphapunto_w2_pw;
    %     Stab_Der.CM_alphapunto_can_pw = CM_alphapunto_can_pw;
    %     Stab_Der.CM_alphapunto_fus = CM_alphapunto_fus;
    
    % Stores Derivatives per parts
    Stab_Der_parts.CM_alphapunto = CM_alphapunto;
    Stab_Der_parts.CM_alphapunto_w1 = CM_alphapunto_w1_pw;
    Stab_Der_parts.CM_alphapunto_fus = CM_alphapunto_fus;
    Stab_Der_parts.CM_alphapunto_HTP = CM_alphapunto_HTP_pw;
    Stab_Der_parts.CM_alphapunto_Vee = CM_alphapunto_Vee_pw;
    Stab_Der_parts.CM_alphapunto_can = CM_alphapunto_can_pw;
    Stab_Der_parts.CM_alphapunto_nac = CM_alphapunto_nac_pw;
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Alpha punto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CD_alphapunto = 0;              % vale 0 para vuelo subs�nico
    CX_alphapunto = -CD_alphapunto;              % vale 0 para vuelo subs�nico
    CZ_alphapunto = -CL_alphapunto; % Roskam 3.14
    
    Stab_Der.CD_alphapunto = CD_alphapunto;
    Stab_Der.CX_alphapunto = CX_alphapunto;
    Stab_Der.CZ_alphapunto = CZ_alphapunto;
    Stab_Der.CL_alphapunto = CL_alphapunto;
    Stab_Der.CM_alphapunto = CM_alphapunto;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CDalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CD_alpha = (C_D1 + 2*C_D2*CL)*CL_alpha_ac;
    CX_alfa = - CD_alpha + CL_alpha_ac*trim_alpha + CL; % Roskam 3.128
    CX_alfa = - CD_alpha  + CL; % Steady State Flight condition Roskam 3.128
    CX_alpha = - CD_alpha + CL; % Pamadi 4.447
    % CD_trim_delta = abs(CD_delta*trim_delta_e);
    % CD = C_D0 + C_D1*CL + C_D2*CL^2 + CD_trim_delta;  %polar aeronave
    CZ_alpha = - CL_alpha_ac -CD_alpha*trim_alpha - CD; % Roskam 3.131
    CZ_alpha = - CL_alpha_ac - CD; % Steady State Flight condition Roskam 3.131
    CM_alpha = CM_alpha_ac;
    
    Stab_Der.CD_alpha = CD_alpha;
    Stab_Der.CX_alpha = CX_alpha;
    Stab_Der.CZ_alpha = CZ_alpha;
    Stab_Der.CM_alpha = CM_alpha;
    
    %Stab_Der.CXalfapunto = CXalphapunto;
    %Stab_Der.CZalfapunto = CZalphapunto;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PITCH ANGLE DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CZ_teta = -CL*sin(theta1); % Pamadi 4.440
    CX_teta = -CL*cos(theta1); % Pamadi 4.434
    CM_teta = 0;
    
    %% Assumption!!! NEED TO CHECK!!
    CD_teta = -CX_teta;
    CL_teta = -CZ_teta;
    
    Stab_Der.CX_teta=CX_teta;
    Stab_Der.CD_teta=CD_teta;
    Stab_Der.CZ_teta=CZ_teta;
    Stab_Der.CL_teta=CL_teta;
    Stab_Der.CM_teta=CM_teta;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPEED DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CLu = (Mach/(1 - Mach^2))*CL;                   %para vuelo subs�nico
    CZu = -2*CL - CLu ;
    % Czu = -2*CL - Clu - 2*m1*q0;
    CDu = 0;
    CXu = -2*CD - CDu;  % Pamadi 4.432
    CMu = 0;                      %0 para vuelo subs�nico
    % CMu = 2*CM + CMu;                      %0 para vuelo subs�nico
    
    Stab_Der.CDu = CDu;
    Stab_Der.CZu = CZu;
    Stab_Der.CXu = CXu;
    Stab_Der.CLu = CLu;
    Stab_Der.CMu = CMu;
    
    CZq = -CLq;
    CXq = 0;
        
    Stab_Der.CXq = CXq;
    Stab_Der.CZq = CZq;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%COEFICIENTES Stab_Der ESTABILIDAD EST�TICA Stab_DerLATERAL.DIRECCIONAL %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%derivadas en funcion de beta%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Cyb
    % ASPro Cyb
    Stab_Der = getCybeta(modelo,Stab_Der);
    Cyb_w = Stab_Der.Cyb_w;
    Cyb_fus = Stab_Der.Cyb_fus;
    Cyb_v = Stab_Der.Cyb_v;
    Cyb = Stab_Der.Cyb;
    
    qinf    = modelo.general.qinf;
    W       = modelo.general.mtow*modelo.general.w_w0*9.8065;
    Sref    = modelo.general.Sref;
    C_L      = modelo.general.CL;
    C_Lw   = modelo.general.CL_w;
    %     C_Lh   = modelo.general.CL_h;
    
    AR_w    = modelo.ala.AR;
    
    z_w     = modelo.ala.Zca;
    LAM_w   = modelo.ala.LAMc2;
    LAMc4_w = modelo.ala.LAMc4;
    diedro  = modelo.ala.diedro;
    
    S_v     = modelo.vertical.S;
    eta_v   = modelo.vertical.eta;
    CLa_v   = modelo.vertical.CLa;
    b_v     = modelo.vertical.b;
    Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca, 'pchip');
    
    vol_fus = modelo.fuselaje.vol;
    CLa_fus = modelo.fuselaje.CLa;
    
    
    %% Cy_beta = Cy_beta_fus + Cy_beta_wing + Cy_beta_vert
    
    if W1 ==1
        %% APORTE ALA
        % - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
        Cy_beta_w1       = -0.00573*abs(dihedral_w1*180/pi);
        
        % - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
        C_L = W/qinf/S_ref;
        Cy_beta_w1       = (CL_w1^2)*(6*tan(Lambda_c2_w1)*sin(Lambda_c2_w1))/(pi*AR_w1*(AR_w1 + 4*cos(Lambda_c2_w1)));
    else
        Cy_beta_w1 = 0;
    end
    
    %% APORTE FUSELAJE
    Dfus_w1  = interp1(x_Area_body, Area_body, xbar_w1, 'pchip');
    % vertical distance between the fuselage line and wing root (positive if wing below center line)
    Ki  = Ki_calc(-z_zbar_w1, Dfus_w1);
    % - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
    %Cy_beta_fus     = -2*Ki*S0/Sref;
    
    % - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
    A_refFus        = Vol_TOT^(2/3);
    Cy_beta_fus     = -Ki*CLalpha_fus*A_refFus/S_ref;
    
    if VTP == 1
        twin_VTP = AC_CONFIGURATION.twin_VTP;
        %% APORTE VERTICAL
        % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
        % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)
        
        % Calculo del side-wash: deflexion de la corriente debida a la presencia
        % del ala
        sidewash = 0.724 + 3.06*(S_VTP/S_ref)/(1+cos(Lambda_c4_w1)) + 0.4*(-z_zbar_w1)/Dfus_w1 + 0.009*AR_w1;
        % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
        % LAMc4: flecha del ala en c/4
        % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))
        
        % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
        % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
        Dfus_VTP  = interp1(x_Area_body, Area_body, x_xbar_VTP, 'pchip');
        
        k       = k_calc(b_VTP_s, Dfus_VTP);
        
        % Single Vertical tail
        Cy_beta_VTP    = -k*CLalpha_VTP*sidewash*S_VTP/S_ref;
        if twin_VTP ==1
            % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
            % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
            Dfus_VTP  = interp1(x_Area_body, Area_body, x_xbar_VTP, 'pchip');
            Cybeta_v_Cybeta_eff = Cybeta_v_Cybeta_eff_calc(b_w2,length_fus,Dfus_VTP/2,b_VTP);
            Cy_beta_VTP = 2*Cy_beta_VTP*Cybeta_v_Cybeta_eff*(S_VTP/S_ref);
        end
        
    else
        Cy_beta_VTP = 0;
    end
    
    if Vee == 1
        %% APORTE VERTICAL
        % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
        % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)
        
        % Calculo del side-wash: deflexion de la corriente debida a la presencia
        % del ala
        sidewash = 0.724 + 3.06*(S_w2_pv/S_ref)/(1 + cos(Lambda_c4_w1)) + 0.4*(-z_zbar_w1)/Dfus_w1 + 0.009*AR_w1;
        % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
        % LAMc4: flecha del ala en c/4
        % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))
        
        % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
        % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
        Dfus_Vee  = interp1(x_Area_body, Area_body, x_xbar_w2, 'pchip');
        k       = k_calc(b_w2_s, Dfus_Vee);
        % Single Vertical tail
        Cy_beta_Vee    = -k*(CLalpha_Vee_e_pw*(tan(dihedral_w2)))*sidewash*S_w2/S_ref;
    else
        Cy_beta_Vee = 0;
    end
    %% DERIVADA TOTAL
    Cy_beta         = Cy_beta_w1 + Cy_beta_fus + Cy_beta_VTP + Cy_beta_Vee;
    Stab_Der.Cyb_w = Cy_beta_w1;
    Stab_Der.Cyb_fus = Cy_beta_fus;
    Stab_Der.Cyb_VTP = Cy_beta_VTP;
    Stab_Der.Cyb_Vee = Cy_beta_Vee;
    Stab_Der.Cyb = Cy_beta;
    
    % Stores Derivatives per parts
    Stab_Der_parts.Cy_beta = Cy_beta;
    Stab_Der_parts.Cy_beta_w1 = Cy_beta_w1;
    Stab_Der_parts.Cy_beta_fus = Cy_beta_fus;
%     Stab_Der_parts.Cy_beta_HTP = Cy_beta_HTP;
    Stab_Der_parts.Cy_beta_Vee = Cy_beta_Vee;
%     Stab_Der_parts.Cy_beta_can = Cy_beta_can;
%     Stab_Der_parts.Cy_beta_nac = Cy_beta_nac;
    
    %% Clb
    % ASPro Clb
    % IMPORTANT MISSING TO INCLUDE
    % The wing-fuselage, horizontal tail-fuselage and/or canard-fuselage contributions to the dihedral effect are found from:
    Stab_Der = getClbeta(modelo,trim_alpha,Stab_Der);
    Clb_wb = Stab_Der.Clb_wb;
    Clb_v = Stab_Der.Clb_v;
    Clb = Stab_Der.Clb;
    
    Mach       = modelo.general.Minf;
    % qinf    = modelo.general.qinf;
    % W       = modelo.general.mtow*modelo.general.w_w0*9.8065;
    Sref    = modelo.general.Sref;
    % C_L      = modelo.general.CL;
    C_Lw   = modelo.general.CL_w;
    %     C_Lh   = modelo.general.CL_h;
    
    D_fus   = modelo.fuselaje.D;
    if modelo.vertical.Xca < modelo.fuselaje.l
        Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca);
    else
        Dfus_v  = 0;
    end
    
    AR_w    = modelo.ala.AR;
    b_w     = modelo.ala.b;
    Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.ala.Xca, 'pchip');
    z_w     = modelo.ala.Zca;
    TR_w    = modelo.ala.TR;
    LAMc2_w = modelo.ala.LAMc2;
    LAMc4_w = modelo.ala.LAMc4;
    lt_w    = modelo.ala.Xca - modelo.ala.xca + modelo.ala.le_y(end) + modelo.ala.ct/2;
    diedro  = modelo.ala.diedro;
    
    S_v     = modelo.vertical.S;
    b_v     = modelo.vertical.b;
    Zca_v     = modelo.vertical.Zca;
    CLa_v   = modelo.vertical.CLa;
    Xca_v    = modelo.vertical.Xca;
    eta_v   = modelo.vertical.eta;
    
    x_XCG = modelo.general.Xcg;
    z_XCG = modelo.general.Zcg;
    
    %% FUSELAJE
    Dfus_w1  = interp1(x_Area_body, Area_body, x_xbar_w1, 'pchip');
    Cl_beta_fus = 1.2*sqrt(AR_w1)*(-z_zbar_w1)/b_w1*2*Dfus_w1/b_w1; % 1/rad
    
    if W1 ==1
        %%  APORTE DEL ALA-FUSELAJE
        % Estimacion de Cl_beta_wf extraida de:
        %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
        %   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)
        Dfus_w1  = interp1(x_Area_body, Area_body, x_xbar_w1, 'pchip');
        
        DClbeta_tau         = -180/pi*0.005*sqrt(AR_w1)*(Dfus_w1/b_w1)^2; % 1/rad
        Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(lambda_w1, AR_w1, Lambda_c2_w1); % 1/rad
        K_MLAM              = K_MLAM_calc(Mach, AR_w1, Lambda_c2_w1);
        
        Kf                  = Kf_calc(lt_w, b_w1, AR_w1, Lambda_c2_w1); % Corrected with DATCOM
        Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_w1, lambda_w1); % 1/rad
        Clbeta_die          = 180/pi*Clbeta_die_calc(lambda_w1, AR_w1, Lambda_c2_w1); % 1/rad
        K_Mdie              = K_Mdie_calc(Mach, AR_w, Lambda_c2_w1); % Corrected with DATCOM
        
        % C_L = W/qinf/Sref;
        %% APORTE ALA
        Cl_beta_w1 = C_Lw*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR) + diedro*(Clbeta_die*K_Mdie + DClbeta_tau);
        % die: Angulo de diedro en radianes
    else
        Cl_beta_w1 = 0;
    end
    
    if HTP == 1
        %%  APORTE DEL ALA-FUSELAJE
        % Estimacion de Cl_beta_wf extraida de:
        %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
        %   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)
        DClbeta_tau         = -180/pi*0.005*sqrt(AR_w2)*(Dfus_w/b_w2)^2; % 1/rad
        Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(lambda_w1, AR_w2, Lambda_c2_w2); % 1/rad
        K_MLAM              = K_MLAM_calc(Mach, AR_w2, Lambda_c2_w2);
        lt_w2 = Geo_tier.x_cT_w2_TE - Geo_tier.cT_w2/2; 
        Kf                  = Kf_calc(lt_w2, b_w2, AR_w2, Lambda_c2_w2); % Corrected with DATCOM
        Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_w2, lambda_w2); % 1/rad
        Clbeta_die          = 180/pi*Clbeta_die_calc(lambda_w2, AR_w2, Lambda_c2_w2); % 1/rad
        K_Mdie              = K_Mdie_calc(Mach, AR_w2, Lambda_c2_w2); % Corrected with DATCOM
        
        % C_L = W/qinf/Sref;
        %% APORTE ALA
        Cl_beta_w2 = C_Lw*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR) + diedro*(Clbeta_die*K_Mdie + DClbeta_tau);
        % die: Angulo de diedro en radianes
    else
        Cl_beta_w2 = 0;
    end
    
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))
    
    % % Calculo de la relacion de aspecto efectiva del vertical
    % AvB_Av      = AvB_Av_calc(b_v, Wfus_v, TR_v);
    % KH          = KH_calc(S_h, S_v);
    % AvHB_AvB    = AvHB_AvB_calc(xac_h, b_v, zac_h, c_v);
    %
    % ARv_eff = AvB_Av*AR_v*(1+KH*(AvHB_AvB-1));
    %
    % % b_v:      Envergadura del vertical
    % %
    % % r1:       Radio medio del fuselaje en donde se encuentra situado el
    % %           vertical.
    % %
    % % lambda_v: Estrechamiento del vertical.
    % %
    % % ARv:      Relaci�n de aspecto del vertical.
    % %
    % % St:       Superficie del estabilizador horizontal.
    % %
    % % Sv:       Superficie del estabilizador vertical.
    % %
    % % xac_h:    Distancia relativa entre el centro aerodinamico del estabilizador y el borde de ataque del estabilizador vertical.
    % %
    % % zac_h:    Distancia del plano del estabilizador horizontal a la linea
    % %           central del fuselaje.
    % %
    % % AvB_Av:   Ratio of the vertical tail aspect ratio in the presence of the
    % %           fuselage to that of the isolated vertical tail.
    % %
    % % AvHB_AvB: Ratio of the vertical tail aspect ratio in the presence of the
    % %           fuselage and the horizontal tail to that in the presence of the fuselage alone.
    % %
    % % KH:       Factor that accounts for the relative size of the horizontal and the
    % %           vertical tail.
    %
    % % Calculo de la pendiente de sustentacion del vertical
    % beta    = sqrt(1-M^2);
    % Cla_M   = Cla_0/beta; % Cla_0 Pendiente de la curva de sustentacion del perfil a angulo de ataque nulo
    % k_cl    = Cla_M/2/pi;
    % Cla_AR  = Cla_AR_calc(LAMc2_v, k_cl, beta);
    % CLa_v   = ARv_eff*Cla_AR;     %CLa_v = 2*pi*ARv_eff/(2 + sqrt((ARv_eff*beta/k)^2*(1+(tan(LAMvc2)/beta)^2)+4));
    
    
    %% DERIVADA DE ESTABILIDAD COMPLETA
    if VTP == 1
        %% APORTE VERTICAL
        % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
        % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)
        
        % Calculo del side-wash: deflexion de la corriente debida a la presencia
        % del ala
        sidewash = 0.724 + 3.06*(S_VTP/S_ref)/(1+cos(Lambda_c4_w1)) + 0.4*(-z_zbar_w1)/Dfus_w1 + 0.009*AR_w1;
        % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
        % LAMc4: flecha del ala en c/4
        % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))
        
        % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
        % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
        Dfus_VTP  = interp1(x_Area_body, Area_body, x_xbar_VTP, 'pchip');
        k       = k_calc(b_VTP_s, Dfus_VTP);
        
%         % Single Vertical tail
%         Cy_beta_VTP    = -k*CLalpha_VTP*sidewash*S_VTP/S_ref;
%         if twin_VTP ==1
%             % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
%             % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
%             Dfus_VTP  = interp1(x_Area_body, Area_body, x_xbar_VTP, 'pchip');
%             Cybeta_v_Cybeta_eff = Cybeta_v_Cybeta_eff_calc(b_w2,length_fus,Dfus_VTP/2,b_VTP);
%             Cy_beta_VTP = 2*Cy_beta_VTP*Cybeta_v_Cybeta_eff*(S_VTP/S_ref);
%         end
        z_zbar_VTP = Geo_tier.z_zbar_VTP;
        Cl_beta_VTP    = Cy_beta_VTP*((z_zbar_VTP - z_XCG)*cos(alpha)-(xbar_VTP - x_XCG)*sin(alpha))/b_w1;
    else
        Cl_beta_VTP = 0;
    end
    
    if Vee == 1
        %% APORTE VERTICAL
        % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
        % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)
        
        % Calculo del side-wash: deflexion de la corriente debida a la presencia
        % del ala
        sidewash = 0.724 + 3.06*(S_w2_pv/S_ref)/(1 + cos(Lambda_c4_w1)) + 0.4*(-z_zbar_w1)/Dfus_w1 + 0.009*AR_w1;
        % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
        % LAMc4: flecha del ala en c/4
        % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))
        
        % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
        % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
        Dfus_Vee  = interp1(x_Area_body, Area_body, x_xbar_w2, 'pchip');
        k       = k_calc(b_w2_s, Dfus_Vee);
        
        % V tail
%         Cy_beta_Vee    = -k*(CLalpha_Vee_e_pw*(tan(dihedral_w2)))*sidewash*S_w2/S_ref;
        Cl_beta_Vee    = Cy_beta_Vee*((z_zbar_w2 - z_XCG)*cos(alpha)-(x_xbar_w2 - x_XCG)*sin(alpha))/b_w1;
        % Single Vertical tail
        
    else
        Cl_beta_Vee = 0;
    end
    
    %% DERIVADA TOTAL
    Cl_beta         = Cl_beta_w1 + Cl_beta_fus + Cl_beta_VTP + Cl_beta_Vee;
    
    % %Clbwb
    % in_K_f = (x_w_LE + mamparo_F_front + cmac_w1/2)/(b_w);
    % K_f = 0.95; % Approx from Fig 3.98 Pamadi
    % KM_A = 1;
    % KM_diedro = 1;
    % ClbCl_flecha = 0.000; % Approx from Fig 3.96 Pamadi
    % ClbCl_AR = -0.001; % Approx from Fig 3.99 Pamadi
    % Clb_diedro   = -0.0003;
    % Delta_Clb_diedro = -0.0005*sqrt(AR_w)*(w_Area_b_max/b_w)^2;
    % Delta_Clb_zw = (1.2*sqrt(AR_w)/57.3)*(-z_ala_centerline/b_w)*(2*w_Area_b_max/b_w);
    % Clb_wb = CL*(ClbCl_flecha*R2D*KM_A*K_f + ClbCl_AR*R2D) + Sigma_w*D2R*(Clb_diedro*(R2D^2)*KM_diedro + ...
    %     Delta_Clb_diedro*(R2D^2)) + Delta_Clb_zw;
    % %Clbv
    % % horizontal and vertical distance betwen the CG and aerocynamic center
    % l_vt1 = x_v_xbar - x_XCG;
    % z_v_xbar = z_v_LE + ybar_v;
    % z_vt1 = z_v_xbar - Zcg;
    % z_v = (z_vt1*cos(trim_alpha) - l_vt1*sin(trim_alpha));
    % Clb_v = Cyb_v*z_v/b_w;
    % Clb = Clb_wb + Clb_v;
    
    
    %% Cnb
    % ASPro Cnb
    Stab_Der = getCnbeta(modelo,Stab_Der);
    Cnb_b = Stab_Der.Cnb_b;
    Cnb_w = Stab_Der.Cnb_w;
    Cnb_wb = Stab_Der.Cnb_wb;
    Cnb_v = Stab_Der.Cnb_v;
    Cnb = Stab_Der.Cnb;
    
    % %Cnbdiedro
    % % Strip theory approximation
    % Cnb_diedro = -Sigma_w*D2R*(CL-CD_alpha)/4; % Rectangular wings with constant Chord - Pamadi Eq 3.270
    % % approximation that solves the problem of strip theory of ignoring the
    % % induced drag effects
    % Cnb_diedro = -0.075*Sigma_w*D2R*CL; % Rectangular wings with constant Chord - Pamadi Eq 3.271
    % Cnb_flecha = CL^2/(4*pi*AR_w); % Appproximation Pamadi 3.299 with all terms with Sweep ignored
    % % Input data to the subfigures in Fig 3.73 Pamadi
    % x_m = (x_XCG + mamparo_F_front)/length_fus;
    % inFig_3_73 = (length_fus^2)/Area_side;
    % h_1 = z_Area_b_1_4;
    % h_2 = z_Area_b_3_4;
    % in_h1_h2 = sqrt(h_1/h_2);
    % in_h_bf_max = h_Area_b_max/w_Area_b_max;
    % Knn=0.00125;                       %de grafica de altura/bfmax
    % % Figura 3.74 - Pamadi
    % Re_fus = rho*V*length_fus/(mu);
    % Re_fus_f_Re_fus = [1e6,7e6, 30e6, 50e6, 80e6 350e6];
    % K_Ri_f_Re_fus = [1,1.4, 1.7 1.8, 1.9,2.2];
    % Kri_f_Re_fus  = interp1(Re_fus_f_Re_fus,K_Ri_f_Re_fus,Re_fus,'spline');
    % Kri=Kri_f_Re_fus;                         %sale del numero de reynolds
    % Cnb_b = -Knn*Kri*(Area_side/S_w)*(length_fus/b_w)*R2D;
    % % Cnbw_b = -1.3*(Vol_TOT/(S_w*b_w))*(h_Area_b_max/w_Area_b_max);
    % Cnb_w = Cnb_diedro + Cnb_flecha;
    % Cnb_wb = Cnb_diedro + Cnb_flecha + Cnb_b;
    % Cnb_v = - Cyb_v*((x_v_xbar - x_XCG)*cos(trim_alpha) + (z_v_xbar - Zcg)*sin(trim_alpha))/b_w;
    % Cnb =  Cnb_wb + Cnb_v;
    
    %%%%%%%%%%%%%%%derivadas propulsivas laterales-direccionales%%%%%%%%%%%%%%%
    % % Data from NACA REport 640
    % w_R_30 = 0.0525*2;
    % w_R_60 = 0.073*2;
    % w_R_90 = 0.045*2;
    % Nprop = 1;
    % CNalpha_p_KN = 0.1;
    % KN = 262*w_R_30 + 262*w_R_60 + 135*w_R_90;
    % dCN_dalpha = CNalpha_p_KN*(1 + 0.8*(KN/80.7) - 1);
    % CyTb = (pi/4)*Nprop*(D_prop^2)*dCN_dalpha/S_w1
    % % Thrust lines inclination angle
    % psi_T = 0;
    % l_prop = (x_XCG-x_m_propeller)*cos(psi_T);
    % CNTb = (pi/4)*Nprop*l_prop*(D_prop^2)*dCN_dalpha/(S_w1*b_w1)
    % Stab_Der.CyTb = CyTb;
    % Stab_Der.CNTb = CNTb;
    % pause
    
    %% REVIEW
    %ASpro CyTb
    Stab_Der = getCyTbeta(modelo,Stab_Der);
    CyTb = Stab_Der.CyTb;
    
    %% REVIEW
    %Aspro CNTb
    Stab_Der = getCnTbeta(modelo,Stab_Der);
    CNTb = Stab_Der.CNTb;
    
    %%%%%%%%%%%%%%%%%%%%%%%control alerones%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Cydeltaa
    % Cydeltaa=0;
    %ASpro Cydeltaa
    Stab_Der = getCyda(modelo,Stab_Der);
    Cydeltaa = Stab_Der.Cydeltaa;
    
    % % % Cldeltaa
    % METODO_AIL = 1;
    % if METODO_AIL == 0;
    %     CLalpha_v_2D = 5.9656;
    %     Cldelta_theory = 4.1666;
    %     Cldelta_Cldelta_theory = 0.45;
    %     K1 = 1.04;
    %     K2 = 0.0714;
    %     Cldeltaa = Cldelta_theory*Cldelta_Cldelta_theory*(CLalpha_v_2D/CLalpha_v)*K1*K2;
    % elseif METODO_AIL == 1;
    % %    M�todo II
    %     tau_a = 0.45;
    % %   ail_int = b_w/2 - b_ail;
    %     ail_int = b_w1/2 - 0.60;
    %     ail_ext = b_w1/2;
    %     int_ca = @(x) cmac_w1.*x;
    %     Cldeltaa = (2*CLalpha_w1*tau_a/(S_w1*b_w1))*quad(@(x)int_ca(x),ail_int,ail_ext);
    % end
    
    %% Cldeltaa
    %ASpro Cldelta
    Stab_Der = getClda(modelo,Stab_Der);
    Cldeltaa = Stab_Der.Cldeltaa;
    
    %% Cndeltaa
    % eta = ((b_w/2)- b_ail)/(b_w/2);
    % K11=-0.15;
    % Cndeltaa = 2*K11*CL*Cldeltaa;
    %ASpro Cndeltaa
    Stab_Der = getCnda(modelo,Stab_Der);
    Cndeltaa = Stab_Der.Cndeltaa;
    
    %%%%%%%%%%%%%%%%%%%%%%%control rudder%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %derivadas en deltar
    if d_rudder ==1
        %% Cydeltar
        Stab_Der = getCydr(modelo,Stab_Der);
        Cydeltar = Stab_Der.Cydeltar;
        
        %% Cldr
        Stab_Der = getCldr(modelo,trim_alpha,Stab_Der);
        Cldeltar = Stab_Der.Cldeltar;
        
        %%  Cndeltar
        Stab_Der = getCndr(modelo,trim_alpha,Stab_Der);
        Cndeltar = Stab_Der.Cndeltar;
    end
    
    if d_rudvtr == 1
        %% Cydeltarv - ruddervator
        Stab_Der = getCydrv(modelo,Stab_Der);
        Cydeltarv = Stab_Der.Cydeltarv;
        Cydeltar = Cydeltarv;
        Stab_Der.Cydeltar = Cydeltar;
        
        %% Cydeltarv - ruddervator
        Stab_Der = getCldrv(modelo,trim_alpha,Stab_Der);
        Cldeltarv = Stab_Der.Cldeltarv;
        Cldeltar = Cldeltarv;
        Stab_Der.Cldeltar = Cldeltar;
        
        %% Cndeltarv - ruddervator
        Stab_Der = getCndrv(modelo,trim_alpha,Stab_Der);
        Cndeltarv = Stab_Der.Cndeltarv;
        Cndeltar = Cndeltarv;
        Stab_Der.Cndeltar = Cndeltar;
    end
    
    % pause
    % z_v     = modelo.vertical.Zca;
    % l_v     = modelo.vertical.Xca - modelo.general.Xcg;
    %
    % l_xcg_w2 = x_xbar_w2 - x_XCG;
    % % Arm from Zcg to Zac_w2
    % l_zcg_w2 = z_zbar_w2 - z_XCG;
    % z_v = (l_zcg_w2*cos(trim_alpha) - l_xcg_w2*sin(trim_alpha));
    %
    % z_vt1 = (z_v*cos(trim_alpha) - l_v*sin(trim_alpha))
    % Cldeltar = Cydeltar*(z_vt1/b_w1_e)
    % Cndeltar = -Cydeltar*((x_v_xbar - x_XCG)/b_w)*eta_h_no_afe;
    
    %%%%%%%%%%%%%%%%%%%Derivadas en funcion de p%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %Cyp
    % k = 1.0; % Fig 3.75 Pamadi
    % in_beta_Clp_S0_CL0 = beta*AR_w_e/(k);
    % beta_Clp_S0_CL0 = -0.175;
    % beta_Clp_S0_CL0 = -0.575 - 0.150;
    % k = CLalpha_w_2D/(2*pi);
    % Clp_Sigma0_CL0 = beta_Clp_S0_CL0*(k/beta);
    % % beta_Clp_S0_CL0 = -0.575;
    % Cyp_w_diedro = 3*sin(Sigma_w*D2R)*(1 - 4*(z_v/b_w)*sin(Sigma_w*D2R))*Clp_Sigma0_CL0;
    % %oswald factor
    % taper_e = 1;
    % a1 = 0.0004;
    % a2 = -0.0080;
    % a3 = 0.0501;
    % a4 =0.8642;
    % lambda1 = AR_w_e*taper_e/(cos(Lambda_LE_w_e));
    % R = a1*lambda1^3 + a2*lambda1^2 + a3*lambda1 + a4;
    % e_e = 1.1*CLalpha_we/(R*CLalpha_we + (1-R)*pi*AR_w_e);
    % % Pamadi method > 1, use approximation
    % e_e = 0.8;
    % aw1 = CLalpha_we/(pi*AR_w_e*e_e);
    % aw2 = e_e*aw1;
    % K = (1-aw1)/(1-aw2);
    % Cyp_CL_CL0_M0 = -0.4; % approx Pamadi Fig 4.24 (not enought data)
    % Cyp_CL_CL0 = ((AR_w_e + beta*cos(Lambda_c4_w1_e))*(AR_w_e*beta + beta*cos(Lambda_c4_w1_e)))/...
    %     ((AR_w_e*beta + 4*cos(Lambda_c4_w1_e))*(AR_w_e + cos(Lambda_c4_w1_e)));
    % Cyp_w1 = K*Cyp_CL_CL0*Cyp_CL_CL0_M0;
    % Cyp_w = Cyp_w1 + Cyp_w_diedro;
    % Cyp_v =  2*((z_v - z_vt1)/b_w)*Cyb_v; % Eq4.545 Pamadi
    % Cyp = Cyp_v + Cyp_w;
    
    %ASpro Cyp
    Stab_Der = getCyp(modelo,trim_alpha,Stab_Der);
    Cyp_w_diedro = Stab_Der.Cyp_w_diedro;
    Cyp_w = Stab_Der.Cyp_w;
    Cyp_v = Stab_Der.Cyp_v;
    Cyp = Stab_Der.Cyp;
    
    % %Clp
    % % z = ((z_v_LE - Zcg)*cos(trim_alpha) - (x_v_xbar - x_XCG)*sin(trim_alpha));
    % % aproxuimation for untwisted rectangular constant chord wing: Eq 4.575
    % % Pamadi
    % Clp_w = -(1/6)*(CLalpha_w_CR + CD);
    % zp = (z_w_LE - Zcg)/b_w;
    % Clp_Sigma_Clp_Sigma0 = 1 - zp*sin(Sigma_w*D2R) + 3*(zp^2)*(sin(Sigma_w*D2R))^2;
    % Clp_w = Clp_Sigma0_CL0*Clp_Sigma_Clp_Sigma0;
    % Clp_v = abs(2*(z_v/b_w)*((z_v - z_vt1)/b_w))*Cyb_v;
    % Clp = Clp_v + Clp_w;
    
    %ASpro Clp
    Stab_Der = getClp(modelo,trim_alpha,Stab_Der);
    Clp_w = Stab_Der.Clp_w;
    Clp_v = Stab_Der.Clp_v;
    Clp = Stab_Der.Clp;
    
    
    % %Cnp
    % chi1 = (chi*tan(Lambda_c4_w1_e)/AR_w_e) + (1/12)*(tan(Lambda_c4_w1_e))^2;
    % Cnp_CL_CL0_M0 = -(AR_w_e + 6*(AR_w_e + cos(Lambda_c4_w1_e))*chi1)/...
    %     (6*(AR_w_e + 4*cos(Lambda_c4_w1_e)));
    % Cnp_CL_CL0 = (((AR_w_e + 4*cos(Lambda_c4_w1_e))/((AR_w_e*beta + 4*cos(Lambda_c4_w1_e))))*...
    %     ((AR_w_e*beta + 0.5*((AR_w_e*beta + cos(Lambda_c4_w1_e)))*(tan(Lambda_c4_w1_e)^2))/...
    %     (AR_w_e + 0.5*((AR_w_e + cos(Lambda_c4_w1_e)))*(tan(Lambda_c4_w1_e)^2))))*...
    %     Cnp_CL_CL0_M0;
    % Cnp_w = Clp_w*tan(trim_alpha)*(K-1) + K*Cnp_CL_CL0*CL;
    % Cnp_v = -(2/b_w)*(x_v_xbar - x_XCG)*cos(trim_alpha) + z_vt1*sin(trim_alpha)*...
    %     ((z_v-z_vt1)/b_w)*Cyb_v;
    % Cnp = Cnp_v + Cnp_w;
    
    %ASpro Cnp
    Stab_Der = getCnp(modelo,trim_alpha,Stab_Der);
    Cnp_w = Stab_Der.Cnp_w;
    Cnp_v = Stab_Der.Cnp_v;
    Cnp = Stab_Der.Cnp;
    
    %%%%%%%%%%%%%Derivadas en respecto a r%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Cyr
    % % rectangular wing approx
    % Cyr_w = 0.143*CL - 0.05;
    % Cyr_v = -(2/b_w)*((x_v_xbar - x_XCG)*cos(trim_alpha) + z_vt1*sin(trim_alpha))*...
    %     Cyb_v;
    % Cyr=Cyr_w + Cyr_v;
    
    %ASpro Cyr
    Stab_Der = getCyr(modelo,trim_alpha,Stab_Der);
    Cyr_w = Stab_Der.Cyr_w;
    Cyr_v = Stab_Der.Cyr_v;
    Cyr = Stab_Der.Cyr;
    
    % % Clr
    % % rectangular wing approx Pamadi 4.612
    % % Clr_w = CL/3;
    % Clr_CL_CL0_M0 = 0.266;
    % Num_clr = 1 + (AR_w_e*(1-beta^2))/(2*beta*(AR_w_e*beta + 2*cos(Lambda_c4_w1_e))) + ...
    %     ((AR_w_e*beta + 2*cos(Lambda_c4_w1_e))/(AR_w_e*beta + 4*cos(Lambda_c4_w1_e)))*...
    %     ((tan(Lambda_c4_w1_e))^2)/8;
    % Den_clr = 1 + ((AR_w_e + 2*cos(Lambda_c4_w1_e))/(AR_w_e + 4*cos(Lambda_c4_w1_e)))*...
    %     ((tan(Lambda_c4_w1_e))^2)/8;
    % Clr_CL_CL0 = (Num_clr/Den_clr)*Clr_CL_CL0_M0;
    % DeltaClr_Sigma = (1/12)*((pi*AR_w_e*sin(Lambda_c4_w1_e))/(AR_w_e + 4*cos(Lambda_c4_w1_e)));
    % Clr_w = CL*Clr_CL_CL0 + DeltaClr_Sigma*Sigma_w*D2R;
    % Clr_v = -(2/b_w^2)*((x_v_xbar - x_XCG)*cos(trim_alpha) + z_vt1*sin(trim_alpha))*...
    %     (z_vt1*cos(trim_alpha) + (x_v_xbar - x_XCG)*sin(trim_alpha))*Cyb_v;
    % Clr=Clr_w + Clr_v;
    
    %ASpro Clr
    Stab_Der = getClr(modelo,trim_alpha,Stab_Der);
    Clr_w = Stab_Der.Clr_w;
    Clr_v = Stab_Der.Clr_v;
    Clr = Stab_Der.Clr;
    
    % % Cnr
    % Cnr_w = -(1/3)*(C_D0 + CD_alpha*trim_alpha)
    % Cnr_v = (2/b_w^2)*(((x_v_xbar - x_XCG)*cos(trim_alpha) + z_vt1*sin(trim_alpha))^2)*Cyb_v
    % Cnr = Cnr_v + Cnr_w
    
    %ASpro Cnr
    Stab_Der = getCnr(modelo,trim_alpha,Stab_Der);
    Cnr_w = Stab_Der.Cnr_w;
    Cnr_v = Stab_Der.Cnr_v;
    Cnr = Stab_Der.Cnr;
    
    %%%%%%%%%%%%%%%derivadas respecto a beta punto%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %Cybpunto
    % in_Sigma_cybpunto1 = z_vt1/(b_w/2);
    % in_Sigma_cybpunto2 = w_Area_b_max/(b_w/2);
    % sigma_beta_alpha = -0.0005;
    % sigma_beta_Sigma = -1; % Approx Fig 4.31 Pamadi
    % sigma_beta_WB = 0.15;
    % sigma_beta = sigma_beta_alpha*(alpha_wb_TRIM_LLT) + sigma_beta_Sigma*Sigma_w + sigma_beta_WB;
    % Cybpunto = 2*CLalpha_v*sigma_beta*(S_v/S_w)*...
    %     ((x_v_xbar - x_XCG)*cos(trim_alpha) + z_vt1*sin(trim_alpha))/b_w;
    
    % ASpro Cybpunto
    Stab_Der = getCybdot(modelo,trim_alpha,Stab_Der);
    Cybpunto = Stab_Der.Cybpunto;
    
    %Clbpunto
    % Clbpunto = Cybpunto*(z_vt1*cos(trim_alpha) - (x_v_xbar - x_XCG)*sin(trim_alpha))/b_w
    
    %ASpro Clbpunto
    Stab_Der = getClbdot(modelo,trim_alpha,Stab_Der);
    Clbpunto = Stab_Der.Clbpunto;
    
    %Cnbpunto
    % Cnbpunto = -Cybpunto*((x_v_xbar - x_XCG)*cos(trim_alpha) + z_vt1*sin(trim_alpha))/b_w;
    
    %ASpro nbpunto
    Stab_Der = getCnbdot(modelo,trim_alpha,Stab_Der);
    Cnbpunto = Stab_Der.Cnbpunto;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stab_Dyn_LatDir = lateral_directional_analysis(V,rho,S_w1,b_w1,cmac_w1,m_TOW,...
    %     Stab_Der,q_inf,theta1,conv_UNITS,Stability_Pamadi,CL,Weight_tier);
    %
    % %%%%%%%%%%%%%%%%%%%%TRIM LATERAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [Trim_ITER_LAT] = Calculo_Trim_ITER_LAT_v3(Geo_tier,conv_UNITS,V,CL,rho,...
    %     Stab_Der,Weight_tier,m_TOW)
    % pause
    
    % [Trim_ITER_LAT_Viraje,Fig] = Calculo_Trim_ITER_LAT_Viraje_v3(Geo_tier,conv_UNITS,...
    %     V,CL,rho,Stab_Der,Weight_tier,m_TOW);
end
% save Stability_Derivatives_v3.mat Stab_Der
%save Stablity_Derivatites_16.mat Stab_Der
