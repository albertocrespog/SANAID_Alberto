function [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts] = ...
    Calculo_Stability_Derivatives_Noviembre2021(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
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
q_inf = 0.5*rho*V^2;                     %presión dinámica inicial
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
    % Positions
    % Distances relative to the origin
    % Position (X,Y,Z) of the Xac for both CR and TOLD
    Xac_cre_W_B = get_xac_cre_wing(lambda_w1, AR_w1,Lambda_LE_w1,Mach);
    Geo_tier.x_xbar_w1 = Xac_cre_W_B*cR_w1 + x_w1_LE;
    Geo_tier.xbar_w1 = Xac_cre_W_B*cR_w1;
    

    x_xbar_w1 = Geo_tier.x_xbar_w1;
    y_ybar_w1 = Geo_tier.y_ybar_w1;
    z_zbar_w1 = Geo_tier.z_zbar_w1;
    % Distances relative to the LE
    xbar_w1 = Geo_tier.xbar_w1;
    xbar_w1_e = Geo_tier.xbar_w1; % effective
    %% Aerodynamic properties
    % Condiciones de vuelo
    % Lift coefficient
    CL0_w1 = Aero.CL_0_w1_CR;
    CL0_w1_e = CL0_w1;
    CLalpha_w1 = Aero.CL_alpha_w1_CR;
    CLalpha_w1_e = CLalpha_w1;
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
    % Distances relative to the LE
    xbar_w2 = Geo_tier.xbar_w2;
    xbar_w2_e = Geo_tier.xbar_w2; % effective
    %% Tail Volume Coefficient
    Vbar_vee = ((x_xbar_w2 - x_XCG)/cmac_w1)*(S_w2_s/S_ref);
    
    %% Aerodynamic properties
    % Condiciones de vuelo
    % Lift coefficient
    CL0_w2 = Aero.CL_0_w2_CR;
    CL0_w2_e = CL0_w2;
    CL_alpha_wb_Vee = Aero.CL_alpha_w2_CR;
    CLalpha_w2_e = CL_alpha_wb_Vee;
    CL_alpha_wb_Vee_0 = Aero.CL_alpha_w2_CR;
    CLalpha_vee = Aero.CLalpha_vee; %%multiplicado por cos(diedro)^2;
    CYbeta_vee = Aero.CYbeta_vee;    
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
    CLalpha_can_e = Aero.CL_alpha_can_CR; %% ??
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
%%%%%%%%%%%%%%%%%CONFIGURACIÓN PUSHER. MODELO DE HELICÓPTERO%%%%%%%%%%%%%%%
D_prop = Prop_data.D_prop;
R_prop = D_prop/2;
S_heli = (pi*R_prop^2);              %superficie de la hélice

% Estimation of Desired Thrust
CL = w_T0/(q_inf*S_ref); % assume equilibry steady state flight
CD = C_D0 + C_D1*CL + C_D2*CL^2;
D = CD*q_inf*S_ref;
Stab_Der.CD = CD;

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

%%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL DOWN-WASH y UPWASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Effects = effects(AC_CONFIGURATION,Geo_tier,Design_criteria,conditions,Aero_TH,Aero,conv_UNITS);

% % Upwash influencing in prop 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Stab_Der_parts, afe] = propwash_influence(AC_CONFIGURATION,Geo_tier,Posicion_Palanca,prop_wash_effect,conditions,Propulsion,Performance,Aero,Effects,OUTPUT_read_XLSX);
%Tanto CL0 como CLalpha en Stab_Der_parts ya están multiplicadas por eta*S_/S_ref; 
%Las derivadas en Aero no, son las del elemento individualmente.


body_interference = 0;
[Stab_Der_parts,  Trim_ITER] = getCLalpha(AC_CONFIGURATION,Body_Geo,Geo_tier,body_interference,AC_type,Stab_Der_parts, Effects,Design_criteria);
%Con body_interference=0, lo que hace es llamar a las
%obtenidas en propwash_influence, es decir, a las corregidas con
%eta*S_/S_w1 y calcula CLalpha y CL_0 de la aeronave completa. OK.


CL0_fus = 0;
CL0_nac = 0;
Stab_Der_parts.CL0_fus = CL0_fus;
Stab_Der_parts.CL0_nac = CL0_nac;

Angle_fus_x_fus = -Angle_fus_x_fus;
Angle_fus_interp = -Angle_fus_interp;


% %% Propulsive Model
[Stab_Der,Trim_ITER] = get_propulsive_derivatives(type_engine,type_prop,bypass_ratio, Propulsion, Aero_TH,Prop_data,Geo_tier,Stab_Der_parts,afe,conv_UNITS,conditions,Performance,Trim_ITER,Stab_Der);
%Se calculan las derivadas propulsivas, no he tocado nada.


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cm_0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%CALCULO CM0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Trim_ITER,Stab_Der_parts,Stab_Der] = getCM0(AC_type,AC_CONFIGURATION,Geo_tier,afe,Effects,Stab_Der_parts,Body_Geo,....
    Aero,conv_UNITS,Design_criteria,conditions,Stab_Der,Trim_ITER);
%Calcula el CM0 de la aeronave completa.

% [Cmalpha_wf, Cmalpha_WB,Cmalpha_wb_nelson] = get_cmalpha_wf_comparacion(AC_type,AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,Effects,afe)
[Trim_ITER,Stab_Der_parts,TRIM_RESULTS,Stab_Der] = get_cmalpha_wf_def(AC_type,AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,Effects,afe,SM_des,OUTPUT_read_XLSX,Stab_Der,Trim_ITER);
%Se calcula el punto neutro, el SM deseado y no deseado, el CMalpha del fus
%y de la aeronave completa, todos los cálculos con y sin tener en cuenta el fuselaje.


% %%%%%%%%%%%%%%%%%%%%%%%CALCULO CONTROL LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der,Trim_ITER] = get_long_control(AC_CONFIGURATION,Geo_tier,Stab_Der_parts,Aero_TH,conditions,afe,Performance,Stab_Der,Trim_ITER);

[TRIM_RESULTS,Trim_ITER,Stab_Der] = resolution_trim_conditions(AC_type,AC_CONFIGURATION,Stab_Der,Stab_Der_parts,Design_criteria,Effects,...
    Geo_tier,conv_UNITS,conditions,Performance,Aero_TH,TRIM_RESULTS,Propulsion,Trim_ITER);

%% Conversion of variables
% x_XCG = TRIM_RESULTS.x_XCG_des;
x_XCG = conditions.x_XCG;
trim_alpha = TRIM_RESULTS.trim_alpha;
CL_w1 = Trim_ITER.CL_w1;
CL_alpha_wb_w1 = Stab_Der_parts.CL_alpha_wb_w1;

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
    CL = Stab_Der.CL;
    deps_dalpha_Vee = Effects.deps_dalpha_Vee;
    deps_dalpha_can = Effects.deps_dalpha_can;
    deps_dalpha_h = Effects.deps_dalpha_h;
    w_R_30 = Trim_ITER.w_R_30;
    w_R_60 = Trim_ITER.w_R_60;
    w_R_90 = Trim_ITER.w_R_90;
    
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
    if HTP == 1
    modelo.general.downwash = deps_dalpha_h;
    elseif Vee == 1
        modelo.general.downwash = deps_dalpha_Vee;
    end
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
    modelo.fuselaje.CLa = Stab_Der_parts.CL_alpha_fus;
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
                eta_VTP_no_afe = afe.eta_VTP_no_afe;
                K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
                K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
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
                modelo.vertical.cm_c = 0.25;
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
                eta_w2_no_afe = afe.eta_w2_no_afe;
                K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
                K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;

                modelo.horizontal.S = Geo_tier.S_w2;
                modelo.horizontal.eta = eta_w2_no_afe;
                modelo.horizontal.Cla = CL_alpha_HTP;
                modelo.horizontal.CLa = CL_alpha_HTP;
                modelo.horizontal.b = (b_w2_s);
                modelo.horizontal.Xca = x_xbar_w2;
                modelo.horizontal.Zca = l_zcg_w2;
                modelo.horizontal.AR = AR_w2;
                modelo.horizontal.TR = lambda_w2;
                modelo.horizontal.cm_c = 0.25;
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
                eta_VTP_no_afe = afe.eta_VTP_no_afe;
                K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
                K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
                
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
                modelo.vertical.cm_c = 0.25;
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
                eta_w2_no_afe = afe.eta_w2_no_afe;
%                 CLalpha_wb_HTP 
                K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
                K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;
                
                modelo.horizontal.S = Geo_tier.S_w2;
                modelo.horizontal.eta = eta_w2_no_afe;
                modelo.horizontal.Cla = CLalpha_wb_HTP;
                modelo.horizontal.CLa = CLalpha_wb_HTP;
                modelo.horizontal.b = (b_w2_s);
                modelo.horizontal.Xca = x_xbar_w2;
                modelo.horizontal.Zca = l_zcg_w2;
                modelo.horizontal.AR = AR_w2;
                modelo.horizontal.TR = lambda_w2;
                modelo.horizontal.cm_c = 0.25;
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
            eta_w2_no_afe = afe.eta_w2_no_afe;
            K_y1_rudvtr_w2 = Geo_tier.K_y1_rudvtr_w2;
            K_y2_rudvtr_w2 = Geo_tier.K_y2_rudvtr_w2;
            
            modelo.vertical.S = Geo_tier.S_w2_pv;
            modelo.vertical.eta = eta_w2_no_afe;
            modelo.vertical.Cla = CL_alpha_wb_Vee;
            modelo.vertical.CLa = CL_alpha_wb_Vee;
            modelo.vertical.b = (b_w2_s/2)*sin(dihedral_w2);
            modelo.vertical.Xca = x_xbar_w2;
            modelo.vertical.Zca = l_zcg_w2;
            modelo.vertical.AR = (((b_w2_s/2)*sin(dihedral_w2))^2)/(S_w2_pv/2);
            modelo.vertical.TR = lambda_w2;
            modelo.vertical.cm_c = 0.25;
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
            modelo.horizontal.Cla = CL_alpha_wb_Vee;
            modelo.horizontal.CLa = CL_alpha_wb_Vee;
            modelo.horizontal.b = (b_w2_s)*cos(dihedral_w2);
            modelo.horizontal.Xca = x_xbar_w2;
            modelo.horizontal.Zca = l_zcg_w2;
            modelo.horizontal.AR = (((b_w2_s)*cos(dihedral_w2))^2)/(S_w2_ph);
            modelo.horizontal.TR = lambda_w2;
            modelo.horizontal.cm_c = 0.25;
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
            modelo.vee.b = b_w2_s;
            modelo.vee.Xca = x_xbar_w2;
            modelo.vee.Zca = l_zcg_w2;
            modelo.vee.dihedral_w2 = dihedral_w2;
            
        case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
            % Estimates properties as proyectos from virtual HTP and virtual
            % VTP
            eta_w2_no_afe = afe.eta_w2_no_afe;
            K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
            K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
            K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
            K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;
            
            modelo.vertical.S = Geo_tier.S_w2_pv;
            modelo.vertical.eta = eta_w2_no_afe;
            modelo.vertical.Cla = CL_alpha_wb_Vee;
            modelo.vertical.CLa = CL_alpha_wb_Vee;
            modelo.vertical.b = (b_w2_s/2)*sin(dihedral_w2);
            modelo.vertical.Xca = x_xbar_w2;
            modelo.vertical.Zca = l_zcg_w2;
            modelo.vertical.AR = (((b_w2_s/2)*sin(dihedral_w2))^2)/(S_w2_pv/2);
            modelo.vertical.TR = lambda_w2;
            modelo.vertical.cm_c = 0.25;
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
            modelo.horizontal.Cla = CL_alpha_wb_Vee;
            modelo.horizontal.CLa = CL_alpha_wb_Vee;
            modelo.horizontal.b = (b_w2_s)*cos(dihedral_w2);
            modelo.horizontal.Xca = x_xbar_w2;
            modelo.horizontal.Zca = l_zcg_w2;
            modelo.horizontal.AR = (((b_w2_s)*cos(dihedral_w2))^2)/(S_w2_ph);
            modelo.horizontal.TR = lambda_w2;
            modelo.horizontal.cm_c = 0.25;
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
            modelo.vee.Cla = CL_alpha_wb_Vee;
            modelo.vee.CLa = CL_alpha_wb_Vee;
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

 [Stab_Der, Stab_Der_parts] = get_pitch_derivatives(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,Body_Geo,Aero);
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Alpha punto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der,Stab_Der_parts] = get_alphadot_derivatives(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,Effects,Aero_TH);
 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PITCH ANGLE DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Stab_Der = get_pitch_angle_derivatives(Stab_Der, TRIM_RESULTS);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPEED DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der] = get_speed_derivatives(Stab_Der,conditions,Performance);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%COEFICIENTES Stab_Der ESTABILIDAD ESTÁTICA Stab_DerLATERAL.DIRECCIONAL %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%derivadas en funcion de beta%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der_parts,Stab_Der] = get_beta_derivatives(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,conditions);

%     
    %%%%%%%%%%%%%%%%%%%%%%%control alerones%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Stab_Der = get_aileron_latdir_der(modelo,Stab_Der);
    Cydeltaa = Stab_Der.Cydeltaa;
    Cldeltaa = Stab_Der.Cldeltaa;
    Cndeltaa = Stab_Der.Cndeltaa;
    
    %%%%%%%%%%%%%%%%%%%%%%%control rudder%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %derivadas en deltar
    if d_rudder ==1
        Stab_Der = get_deltar_latdir_deriv(modelo,Stab_Der,Geo_tier,afe,trim_alpha);
        Cydeltar = Stab_Der.Cydeltar;
        Cldeltar = Stab_Der.Cldeltar;
        Cndeltar = Stab_Der.Cndeltar;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%control ruddervator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if d_rudvtr == 1
        Stab_Der = get_delta_rv_latdir_deriv(modelo,Stab_Der,afe,trim_alpha,Aero,Geo_tier);
        
        Cydeltarv = Stab_Der.Cydeltarv;
        Cydeltar = Cydeltarv;
        Stab_Der.Cydeltar = Cydeltar;
        
        Cldeltarv = Stab_Der.Cldeltarv;
        Cldeltar = Cldeltarv;
        Stab_Der.Cldeltar = Cldeltar;
        
        Cndeltarv = Stab_Der.Cndeltarv;
        Cndeltar = Cndeltarv;
        Stab_Der.Cndeltar = Cndeltar;
        
    end
    
   
    %%%%%%%%%%%%%%%%%%%Derivadas en funcion de p%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %ASpro Cyp
    Stab_Der = getCyp(modelo,trim_alpha,Stab_Der,Stab_Der_parts,Trim_ITER);
    Cyp_w_diedro = Stab_Der.Cyp_w_diedro;
    Cyp_w = Stab_Der.Cyp_w;
    Cyp_v = Stab_Der.Cyp_v;
    Cyp = Stab_Der.Cyp;

    %ASpro Clp
    Stab_Der = getClp(modelo,trim_alpha,Stab_Der,Stab_Der_parts);
    Clp_w = Stab_Der.Clp_w;
    Clp_v = Stab_Der.Clp_v;
    Clp = Stab_Der.Clp;
    
    %ASpro Cnp
    Stab_Der = getCnp(modelo,trim_alpha,Stab_Der,Stab_Der_parts,Trim_ITER);
    Cnp_w = Stab_Der.Cnp_w;
    Cnp_v = Stab_Der.Cnp_v;
    Cnp = Stab_Der.Cnp;
  
    %%%%%%%%%%%%%Derivadas en respecto a r%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %ASpro Cyr
    Stab_Der = getCyr (modelo,trim_alpha,Stab_Der,Stab_Der_parts);
    Cyr_w = Stab_Der.Cyr_w;
    Cyr_v = Stab_Der.Cyr_v;
    Cyr = Stab_Der.Cyr;
        
    %ASpro Clr
    Stab_Der = getClr(modelo,trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER);
    Clr_w = Stab_Der.Clr_w;
    Clr_v = Stab_Der.Clr_v;
    Clr = Stab_Der.Clr;

    
    %ASpro Cnr
    Stab_Der = getCnr(modelo,trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER);
    Cnr_w = Stab_Der.Cnr_w;
    Cnr_v = Stab_Der.Cnr_v;
    Cnr = Stab_Der.Cnr;
    
    %%%%%%%%%%%%%%%derivadas respecto a beta punto%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Roskam page 149: Except for airplanes in the high subsonic speed
    %range, the beta_dot derivatives are frequently considered
    %negligible.
    
    % ASpro Cybpunto
    Stab_Der = getCybdot(modelo,trim_alpha,Stab_Der,AC_CONFIGURATION, Geo_tier, Aero);
    Cybpunto = Stab_Der.Cybpunto;

    
    %ASpro Clbpunto
    Stab_Der = getClbdot(modelo,trim_alpha,Stab_Der);
    Clbpunto = Stab_Der.Clbpunto;
        
    %ASpro nbpunto
    Stab_Der = getCnbdot(modelo,trim_alpha,Stab_Der);
    Cnbpunto = Stab_Der.Cnbpunto;

end

