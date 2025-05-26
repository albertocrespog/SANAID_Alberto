function [TRIM_RESULTS,Trim_ITER] = Calculo_Trimado_2019(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,...
    Prop_data,conv_UNITS,Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,...
    V_Performance,Data_ATM,XCG_data)

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

% % Trim results
% trim_alpha = TRIM_RESULTS.trim_alpha;
% theta1 = trim_alpha;

%% Mass and inertias
m_TOW = Weight_tier.m_TOW;
w_T0 = m_TOW*g;
alpha_f = conditions.alpha_f;
beta_f = conditions.beta_f;

% Number of engines
n_eng = Prop_data.n_eng;

% Performance
alpha_max = V_Performance.alpha_max;
CL_max_w1 = V_Performance.CL_max_w1;
CL_max_w1_ope = V_Performance.CL_max_w1_ope;
alpha_max_w1_ope = V_Performance.alpha_max_w1_ope;
V_stall = V_Performance.V_stall;
V_min = V_Performance.V_min;
% Flight Conditions
h = conditions.h;
V = conditions.V;
rho = Data_ATM.rho;
q_inf = 0.5*rho*V^2;                      %presión dinámica inicial
a = Data_ATM.a;
% Mach Number
M = V/a;

% % Determines XCG either for Stability in Forward Flight (From calculation of Neutra Point)
% if XCG_FF == 1
%     % XCG from design
%     x_XCG = NP.XCG_sm;
%     y_XCG = XCG_data.z_XCG;
%     z_XCG = XCG_data.z_XCG;
% else
    % XCG from design
    x_XCG = XCG_data.x_XCG;
    y_XCG = XCG_data.y_XCG;
    z_XCG = XCG_data.z_XCG;
% end

% Design Selected incidences
i_w1 = Design_criteria.i_w1;
i_w2 = Design_criteria.i_w2;
% Select txt files associated with each aerodynamic study
% Make sure to check with case asssigned in read_aero_files_Aug2018.m
index_w1 = Design_criteria.index_w1;
index_w2 = Design_criteria.index_w2;

% Max values controil surfaces
delta_rudvtr_min = Geo_tier.delta_rudvtr_min;
delta_rudvtr_max = Geo_tier.delta_rudvtr_max;

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

%% Geometry
% distances
b_w1 = Geo_tier.b_w1;
b_w2 = Geo_tier.b_w2;
cR_w1 = Geo_tier.cR_w1;
cT_w1 = Geo_tier.cT_w1;
cmac_w1 = Geo_tier.cmac_w1;
cmac_w1_e = Geo_tier.cmac_w1;
cmac_w2 = Geo_tier.cmac_w2;
cR_w2 = Geo_tier.cR_w2;
cT_w2 = Geo_tier.cT_w2;
cmac_w2_e = Geo_tier.cmac_w2;
% Arm from Xac w1 to w2
l_xac_w1w2 = Geo_tier.l_xac_w1w2;

% Positions
x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
x_cR_w2_LE = Geo_tier.x_cR_w2_LE;
% Position of the Xac for both CR and TOLD
x_xbar_w1 = Geo_tier.x_xbar_w1;
y_ybar_w1 = Geo_tier.y_ybar_w1;
z_zbar_w1 = Geo_tier.z_zbar_w1;
xbar_w1 = Geo_tier.xbar_w1;
xbar_w1_e = Geo_tier.xbar_w1;
x_w1_LE = Geo_tier.x_w1_LE;
z_w1_LE = Geo_tier.z_w1_LE;
x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
y_cR_w1_TE = Geo_tier.y_cR_w1_TE;
z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
% Position of the Xac for both CR and TOLD
x_xbar_w2 = Geo_tier.x_xbar_w2;
y_ybar_w2 = Geo_tier.y_ybar_w2;
z_zbar_w2 = Geo_tier.z_zbar_w2;
xbar_w2 = Geo_tier.xbar_w2;
xbar_w2_e = Geo_tier.xbar_w2;
x_w2_LE = Geo_tier.x_w2_LE;
z_w2_LE = Geo_tier.z_w2_LE;
x_cR_w2_TE = Geo_tier.x_cR_w2_TE;
y_cR_w2_TE = Geo_tier.y_cR_w2_TE;
z_cR_w2_TE = Geo_tier.z_cR_w2_TE;

% Propeller distances
z_d_T = Geo_tier.z_d_T;
x_d_T = Geo_tier.x_d_T;

% Thrust line inclination angle (for futuire version with tilting rotors)
phi_T =0;
% The perpendicular distance between the thrustline and the airplane center of gravity is given by:
d_T = x_d_T*sin(phi_T) + z_d_T*cos(phi_T);

% Areas
S_w1 = Geo_tier.S_w1;
S_w1_e = Geo_tier.S_w1_e;
% prop wash area affecting surfaces
S_w1_pw = Geo_tier.S_w1_pw;
S_w1_afe = S_w1_pw;
S_w2 = Geo_tier.S_w2; 
S_w2_e = Geo_tier.S_w2_e;
S_w2_e = Geo_tier.S_w2_e; 
% prop wash area affecting surfaces
S_w2_pw = Geo_tier.S_w2_pw;
S_w2_afe = S_w2_pw;
S_ref = Geo_tier.S_ref;

% Adimensional
AR_w1 = Geo_tier.AR_w1;
AR_w1_e = Geo_tier.AR_w1_e;
AR_w2 = Geo_tier.AR_w2;
AR_w2_e = Geo_tier.AR_w2_e;
lambda_w1 = Geo_tier.lambda_w1;
lambda_w1_e = Geo_tier.lambda_w1_e;
lambda_w2 = Geo_tier.lambda_w2;
lambda_w2_e = Geo_tier.lambda_w2_e;

% angles
Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
Lambda_c4_w2 = Geo_tier.Lambda_c4_w2;


%% Aerodynamic properties
% Condiciones de vuelo
C_D0 = Aero_TH.CD0;
C_D1 = Aero_TH.CD1;
C_D2 = Aero_TH.CD2;

CD0_w1 = Aero_TH.CD0_w1;
CD0_w2 = Aero_TH.CD0_w2;
alpha_CL_0_w1 = Aero.alpha_CL_0_w1_CR;
alpha_CL_0_w2 = Aero.alpha_CL_0_w2_CR;

% Cruise
CLalpha_w1 = Aero.CL_alpha_w1_CR;
CLalpha_w1_e = CLalpha_w1;
CLalpha_w2 = Aero.CL_alpha_w2_CR;
CLalpha_w2_e = CLalpha_w2;

% CL0_wb = Aero3D.CL0_wb_CR;
CL0_w2 = Aero.CL_0_w2_CR;
CL0_w1 = Aero.CL_0_w1_CR;
CL0_w1_e = CL0_w1;
CL0_w2_e = CL0_w2;

% CM
CM0_w1 = Aero.CM_0_w1_CR;
CM0_w1_e = CM0_w1;
CM0_w2 = Aero.CM_0_w2_CR;
CM0_w2_e = CM0_w2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CONFIGURACIÓN PUSHER. MODELO DE HELICÓPTERO%%%%%%%%%%%%%%%
D_prop = Prop_data.D_prop; 
R_prop = D_prop/2;
S_heli = (pi*R_prop^2);              %superficie de la hélice

% Estimation of Desired Thrust
CL = w_T0/(q_inf*S_ref);
CD = C_D0 + C_D1*CL + C_D2*CL^2;
D = CD*q_inf*S_ref;
Fdes = D;
[Propulsion] = get_EngineProperties_v3(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,n_eng);
T_tot = Propulsion.Ti;
T_eng = Propulsion.Ti_eng;

%% Estimation of data for induced velocity associated to engine configuration
% v_i = -(1/2)*V + sqrt(1/4*V^2 + T_eng/(2*rho*S_heli));
% R_inf = R_heli*sqrt((V + v_i)/(V + 2*v_i));
% Propeller variables
v_i = Propulsion.v_i; % Induced velocity at prop disk
R_inf = Propulsion.R_inf; % Radius of prop wash at infinity
v_inf = Propulsion.v_inf; % induced velocity at propwash at infinity

%%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kh = (1-(z_w2_LE - z_w1_LE)/b_w1)/(2*l_xac_w1w2/b_w1)^(1/3);          %eq 3.43 PAMADI
kl = (10-3*lambda_w1)/7;                                     %eq 3.43 PAMADI
ka = 1/AR_w1_e -1/(1 + AR_w1_e^(1.7));                        %eq 3.43 PAMADI
deps_dalpha_clean = 4.44*(kh*kl*ka*sqrt(cos(Lambda_c4_w1)))^1.19;  %eq 3.43 PAMADI
deps_dalpha_poff = 0; % Assumes not affecting
deps_dalpha = deps_dalpha_clean + deps_dalpha_poff;

% NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
downwash = 1 - deps_dalpha;
% eps_vee_0 = deps_dalpha*(alpha_CL_0_w1*D2R - i_w1); 
eps_w2_0 = deps_dalpha*(alpha_CL_0_w1*D2R - i_w1); 
% eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
eps_w2 = eps_w2_0 + deps_dalpha*alpha_f;

% Upwash influencing in prop
% Location of the prop disk (source of thrust)
x_prop_cF = Geo_tier.x_prop_cF; % location of prop
AR_w1 = Geo_tier.AR_w1;
depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);
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
x_w2_wake = acos(gamma_w2 - alpha_f - i_w1 - eps_w2);
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

if Posicion_Palanca == 0
    S_w1_afe = S_w1_pw;
    S_w2_afe = S_w2_pw;
    S_w1_no_afe = S_w1 - S_w1_afe;
    S_w2_no_afe = S_w2;
    V_w1 = V;
    V_w2 = V;
    q_w1_afe    = q_inf;
    q_w2_afe    = eta_w2*q_inf;                    %presion dinamica afectada
    q_w1_no_afe = q_inf;
    q_w2_no_afe = eta_w2*q_inf;                    %presion dinamica no afectada
else
    S_w1_afe = S_w1_pw;
    S_w2_afe = S_w2_pw;
    S_w1_no_afe = S_w1 - S_w1_afe;
    S_w2_no_afe = S_w2;
    V_w1 = V + v_i;
    V_w2 = sqrt(eta_w2)*V + 2*v_i;
    q_w1_afe       = 0.5*rho*V_w1^2;
    q_w2_afe       = 0.5*rho*V_w2^2;                 %presion dinamica afectada
    q_w1_no_afe   = q_inf;
    q_w2_no_afe    = eta_w2*q_inf;                    %presion dinamica no afectada
end

eta_w1_afe = (q_w1_afe/q_inf);
eta_w2_afe = (q_w2_afe/q_inf);
eta_w1_no_afe = (q_w1_no_afe/q_inf);
eta_w2_no_afe = (q_w2_no_afe/q_inf);
eta_w1_afe_S_w1_afe_S_ref = eta_w1_afe*(S_w1_afe/S_ref);
eta_w1_no_afe_S_w1_no_afe_S_ref = eta_w1_no_afe*(S_w1_no_afe/S_ref);
eta_w2_afe_S_w2_afe_S_ref = eta_w2_afe*(S_w2_afe/S_ref);
eta_w2_no_afe_S_w2_no_afe_S_ref = eta_w2_no_afe*(S_w2_no_afe/S_ref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correction of Lift curve slope of w1 associated to prop wash
% eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
% eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
CLalpha_w1_e_pw = eta_w1_afe_S_w1_afe_S_ref*CLalpha_w1_e + eta_w1_no_afe_S_w1_no_afe_S_ref*CLalpha_w1_e;
CLalpha_w2_e_pw = eta_w2_afe_S_w2_afe_S_ref*CLalpha_w2_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CLalpha_w2_e;

%% Correction of Lift curve slope of w1 associated to fuselage interaction
% Pendiente de sustentación del morro.
Fineness_Ratio = length_fus/w_Area_b_max;
%digitaliazacion figura PAMADI CAP3
x_f_k2_k1 = [4.,5.,6.,8.,10.,12.,14.,16.,18.,20.];
y_f_k2_k1 = [.77,.825,.865,.91,.94,.955,.965,.97,.973,.975];
f_k2_k1  = interp1(x_f_k2_k1,y_f_k2_k1,Fineness_Ratio,'spline');    
x1 = x_Area_b_max;

%% W1 body interference
aN_w1 = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
aN_w1 = CLa_fus; %  corrigiendo con la integrtación del fuselaje a partir de modelo XFLR5
KN_w1 = (aN_w1/CLalpha_w1_e_pw)*(S_w1/S_w1_e);  %eq 3.25 PAMADI    
CLalpha_fus = aN_w1;
KWB_w1 = 0.1714*(w_Area_b_max/b_w1)^2 + 0.8326*(w_Area_b_max/b_w1) + 0.9974; %eq 3.27 PAMADI
KBW_w1 = 0.7810*(w_Area_b_max/b_w1)^2 + 1.1976*(w_Area_b_max/b_w1) + 0.0088; %eq 3.28 PAMADI
a_bw_w1 = KBW_w1*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.34 PAMADI
a_wb_w1 = KWB_w1*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.33 PAMADI
a_WB_w1 = (KN_w1 + KWB_w1 + KBW_w1)*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.24 PAMADI
CLalpha_WB_w1 = a_WB_w1; % contribution of ww1 to fuselage and fuselage to w1
CLalpha_wb_w1 = CLalpha_WB_w1;

%% W2 body interference
aN_w2 = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
aN_w2 = CLa_fus; %  corrigiendo con la integrtación del fuselaje a partir de modelo XFLR5
KN_w2 = (aN_w2/CLalpha_w2_e_pw)*(S_w2/S_w2_e);  %eq 3.25 PAMADI    
CLalpha_fus = aN_w2;
KWB_w2 = 0.1714*(w_Area_b_max/b_w2)^2 + 0.8326*(w_Area_b_max/b_w2) + 0.9974; %eq 3.27 PAMADI
KBW_w2 = 0.7810*(w_Area_b_max/b_w2)^2 + 1.1976*(w_Area_b_max/b_w2) + 0.0088; %eq 3.28 PAMADI
a_bw_w2 = KBW_w1*CLalpha_w2_e_pw*(S_w2_e/S_w2); %eq 3.34 PAMADI
a_wb_w2 = KWB_w1*CLalpha_w2_e_pw*(S_w2_e/S_w2); %eq 3.33 PAMADI
a_WB_w2 = (KN_w2 + KWB_w2 + KBW_w1)*CLalpha_w2_e_pw*(S_w2_e/S_w2); %eq 3.24 PAMADI
CLalpha_WB_w2 = a_WB_w2; % contribution of ww1 to fuselage and fuselage to w1
CLalpha_wb_w2 = CLalpha_WB_w2;
%% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
CLalpha_wb_w2 = CLalpha_w2_e_pw;

%%%%%%%%%%%%%%%%%CALCULO CENTRO AERODINAMICO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_f = sqrt(1 - M^2);
Check_method = beta_f*AR_w1_e;

% Digitalizacion figura 4.3.2.2-35 Datcom 
x_f_xi = [0.0,.025,.050,.075,.10,.15,.20,.30,.40,.50,.60,.70,.80];
y_f_xi = [0.0,.056,.101,.130,.152,.190,.22,.266,.301,.33,.348,1.365,.375];
f_xi  = interp1(x_f_xi,y_f_xi,w_Area_b_max/b_w1,'spline');
% Revised version Álvaro
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

% Aerodynimic centers based in teh leading edge 
int_Sb = @(x) interp1(x_Area_body(1:end-1),dSdX,x).*(l_N - x);     %eq 3.36 PAMADI
Xac_cre_N = -(1/(cmac_w1_e*Area_b_max))*quad(@(x) int_Sb(x),0,x_0); %eq 3.36 PAMADI
  
if Check_method > 4
    Xac_cre_BW = 1/4 + (b_w1 - w_Area_b_max)/(2*cmac_w1_e)*f_xi*tan(Lambda_c4_w1); %eq 3.39 PAMADI
else
    Warning = 'WARNING!!! Revise method for estimating Xac with WB - Code in PAUSE';
    disp(Warning)
    pause
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
    %% Revisado Álvaro
    Xca_cre_BW4     = 1/4 + (b_w - w_Area_b_max)/(2*cmac_w1_e)*chi1*tan(Lambda_c4_w1);
    Xca_cre_BW0     = xac_cre_BW0_calc(AR_w1_e, lambda_w1, Lambda_LE_w1);
    Xca_cre_BW      = interp1([0 4], [Xca_cre_BW0 Xca_cre_BW4], Check_method, 'linear');
end
    
Xac_cre_WB = xbar_w1/cmac_w1_e; %eq 3.38 PAMADI % assumes the influence of the body on the location of the wing aerodynamic center is small
XAC_WB_cre = ((Xac_cre_N*aN_w1 + Xac_cre_WB*a_wb_w1 + Xac_cre_BW*a_bw_w1)/a_WB_w1)*cmac_w1_e; %eq 3.32 PAMADI
% Xace at teh appex (w1 croot LE)
XAC_WB_LE = XAC_WB_cre*(cmac_w1_e/cmac_w1);                                        %eq 3.36 PAMADI
% Convertin Xac appex to reference origin
XAC_WB = XAC_WB_LE + x_w1_LE;                                                     %eq 3.36 PAMADI  
x_xbar_wb_w1 = XAC_WB;

CL_alpha_w1w2b = CLalpha_wb_w1 + CLalpha_w2_e_pw*(downwash);
%eq PAMADI (sin canard)
X_NP = (CLalpha_wb_w1*XAC_WB + CLalpha_w2_e_pw*(downwash)*x_xbar_w2)/CL_alpha_w1w2b;

%% NOTE IMPORTANT
%% Estimation of SM with only CLalpha contribution
SM_CLalpha = (X_NP - x_XCG)/cmac_w1;
% Data
Trim_ITER.X_NP = X_NP;
Trim_ITER.SM_CLalpha = SM_CLalpha;

%%%%%%%%%%%%%%%%%%%%%%%%CALCULO CL0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CL0_w1_e_corrected = eta_w1_afe_S_w1_afe_S_ref*CL0_w1_e + eta_w1_no_afe_S_w1_no_afe_S_ref*CL0_w1_e;
CLalpha_w1_corrected = CLalpha_wb_w1;
CL0_w2_e_corrected = eta_w2_afe_S_w2_afe_S_ref*CL0_w2_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CL0_w2_e;
CLalpha_w2_corrected = CLalpha_wb_w2;

CL0_w1w2b = CL0_w1_e_corrected + CLalpha_w1_corrected*i_w1 + CL0_w2_e_corrected + ...
    CLalpha_w2_corrected*(i_w2 - eps_w2);

% Storing DATA
Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
Trim_ITER.CLalpha_w1_corrected = CLalpha_w1_corrected;
Trim_ITER.CL0_w2_e_corrected = CL0_w2_e_corrected;
Trim_ITER.CLalpha_w2_corrected = CLalpha_w2_corrected;
Trim_ITER.CL0_w1w2 = CL0_w1w2b;

%%%%%%%%%%%%%CALCULO CLdeltae %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METODO_HTP = 1;
% if METODO_HTP == 0;
%     d_alpha_d_delta_e = 0.46;
%     CL_delta_e_h = CLalpha_h*d_alpha_d_delta_e;
%     CL_delta_e = CL_delta_e_h*(eta_h_afe*S_h_afe/S_w + eta_h_no_afe*S_h_no_afe/S_w);

% ASpro CLdeltae, CDdeltae
%% Estimación con proyecciones en el plano horizontal
cf_rudvtr = Geo_tier.cf_rudvtr;
t_c_rudvtr = Geo_tier.t_c_rudvtr;
K_y1_rudvtr_w1 = Geo_tier.K_y1_rudvtr_w1;
K_y2_rudvtr_w1 = Geo_tier.K_y2_rudvtr_w1;
lambda_w2 = Geo_tier.lambda_w2;
AR_w2_e = Geo_tier.AR_w2_e;

% theoretical V-Tail sectional lift curve slope at Mach equal to zero
Cla_M0 = 2*pi + 5.0525*t_c_rudvtr;
Kb_w2                = Kb_calc(K_y1_rudvtr_w1, K_y2_rudvtr_w1, lambda_w2);              % Fig 3.36 Pamadi
Clde_Cldetheo_w2     = Cldelta_Cldeltatheory_calc(cf_rudvtr, Cla_M0);% Fig 3.37b Pamadi
Cldetheo_w2          = Cldeltatheory_calc(cf_rudvtr, t_c_rudvtr);               % Fig 3.37a Pamadi
alphaCL_alphaCl_w2   = alphaCL_alphaCl_calc(cf_rudvtr,AR_w2_e);              % Fig 3.35 Pamadi
alpha_delta_e       = Kb_w2*Clde_Cldetheo_w2*Cldetheo_w2/Cla_M0*alphaCL_alphaCl_w2;
CL_delta_e          = alpha_delta_e*CLalpha_w2_corrected;
CD_delta_e          = 2*CL0_w1w2b*C_D2*(eta_w2_afe_S_w2_afe_S_ref*CL_delta_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CL_delta_e);
% CD_delta_e          = (2*C_D2*W_TO_1/(q_inf*S_w) + C_D1)*CL_delta_e;

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
x1bar_sec = x_w1_LE - cmac_w1_e;

wing_body_CLalpha = 0.0785;
int_bf = @(x) interp1(length_x_position,width_x_position,x);
% Método III eq 3.7 PAMADI
% Section from tip-nose to front (from graphic)
f_d_epsu_d_alpha_f1  = @(x) interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,x,'linear','extrap').*((x_w1_LE) - x);
int_low_1 = 0;
int_up_1 = (x_w1_LE - cmac_w1_e);
int_Cmalpha_1 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f1(x)),int_low_1,int_up_1);
%   CMalpha_f1 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_1;
CMalpha_f1 = (CLalpha_wb_w1/(wing_body_CLalpha*R2D))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_1;

% Section from right before wing (from graphic)
f_d_epsu_d_alpha_f2  = @(x) interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,x,'linear','extrap').*((x_w1_LE) - x);
int_low_2 = (x_w1_LE - cmac_w1_e);
int_up_2 = x_w1_LE;
int_Cmalpha_2 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f2(x)),int_low_2,int_up_2);
%   CMalpha_f2 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_2;
CMalpha_f2 = (a_WB_w1/(wing_body_CLalpha*R2D))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_2;

% Section from behind the wing (from eq 3.10 PAMADI)
int_low_3 = x_w1_LE + cmac_w1_e;
int_up_3 = length_fus;
x_3 = (x_w2_LE  + cmac_w2/4) - (x_w1_LE + cmac_w1_e);
f_d_epsu_d_alpha_f3 = @(x) ((x-(x_w1_LE + cmac_w1_e))/x_3).*downwash;
int_Cmalpha_3 = quad(@(x) ((int_bf(x)).^2).*f_d_epsu_d_alpha_f3(x),int_low_3,int_up_3);
%   CMalpha_f3 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_3;
CMalpha_f3 = (a_WB_w1/(wing_body_CLalpha*R2D))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_3;
CMalpha_f = CMalpha_f1 + CMalpha_f2 + CMalpha_f3;

Angle_fus_x_fus = -Angle_fus_x_fus;
Angle_fus_interp = -Angle_fus_interp;

alpha_CL_0_w1 = Aero.alpha_CL_0_w1_CR;
alpha_CL_0_w2 = Aero.alpha_CL_0_w2_CR;

%CALCULO Cmo eq 3.8 PAMADI
int_low_anglef = 0;
alpha_CL_0_we_fus = (alpha_CL_0_w1*D2R + i_w1);
int_up_anglef = length_fus;
f_anglef = @(x) interp1(x_interp,Angle_fus_interp,x);
int_Cm0 = quad(@(x) ((int_bf(x)).^2.*(alpha_CL_0_we_fus + f_anglef(x))),int_low_anglef,int_up_anglef);
CM0_f = (f_k2_k1/(36.5*S_w1*cmac_w1))*int_Cm0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%DERIVADA PROPULSIVA LONGITUDINAL%%%%%%%%%%%%%%%%
% Increment of Pitching Moment Coefficient due to the Lift Component of the Propeller Normal Force Aproximada a 0
Delta_CM_Nprop = 0;
% change in airplane pitching moment coefficient due to the lift component of the propeller thrust force.
% T_avail = CT*q_inf*S_w;
T_avail = Propulsion.Ti;

% Storing DATA
x_d_T = Geo_tier.x_d_T;
y_d_T = Geo_tier.y_d_T;
z_d_T = Geo_tier.z_d_T;

Delta_CM_Tprop = - T_avail*z_d_T/(q_inf*S_w1*cmac_w1);
CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;

%%%%%%%%%%%%%%%%%%%%%%%%%%% CMTalpha%%%%%%%%%%%%%%%%
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

T_SET = Propulsion.Ti_eng;
CT_p_eng = T_SET/(q_inf*S_ref);
% P_SHP = P_SET*W2hp;
Pi_eng = Propulsion.Pi_eng; % power per engine
P_SHP = Pi_eng*W2hp;
% rho_met2imp = 0.0017645;

W_current_lb = m_TOW*2.20462;

% The intermediate calculation parameter is given by:
dTc_dCL = (3/2)*550*P_SHP*sqrt(rho*rho_SI2rho_IMP)*sqrt(CL)/...
    (sqrt((2*W_current_lb/(S_w1*m22ft2))^3)*(D_prop*m2ft)^2);
% dCM_dCL_TL = (Nprop*(2*((D_prop*m2ft)^2)*d_T*m2ft)/(S_w1*m22ft2*cmac_w1*m2ft))*dTc_dCL;
dCM_dCL_TL = (Nprop*(2*((D_prop)^2)*d_T)/(S_ref*cmac_w1))*dTc_dCL;
% DATCOM - 
% ----PROPELLER INFLOW FACTOR
%      ----FIGURE 4.6.1-25B
f_inflow_input_vec = [0.,1.,2.,3.,4.,6.,8.,10.,14.,19.,22.]; 
f_inflow_input = S_ref*(CT_p_eng)/(8*(R_prop^2));
f_inflow_vec = [1.0,1.55,1.94,2.20,2.40,2.75,3.05,3.30,3.75,4.25,4.54]; 
finflow  = interp1(f_inflow_input_vec,f_inflow_vec,f_inflow_input,'spline');


dCM_dCL_N = (pi/4)*finflow*Nprop*(x_d_T)*m2ft*((D_prop*m2ft)^2)...
    *dCN_dalpha*(1 - depsu_dalpha)/(S_w1*m22ft2*cmac_w1*m2ft*CL_alpha_w1w2b);
Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;
CMTalpha = Delta_CM_CL_T*CL_alpha_w1w2b;
% CMTalpha = 0;

%% NOTE Corrections of x_XCG taking into avvount fuselage CMalpha contribution
% CM_alpha CONTRIBUTION ONLY TAKING INTO ACCOUNT CLalpha
CM_alpha_w1w2b_CLalpha = -CL_alpha_w1w2b*SM_CLalpha; 
CM_alpha_w1w2b_fus = CM_alpha_w1w2b_CLalpha + CMTalpha + CMalpha_f; % Total CMañpha with fuselage contribution
SM_actual_fus = -CM_alpha_w1w2b_fus/CL_alpha_w1w2b; % Real SM fith fuselage contribution

% Adjusting x_XCG to ensure longitudinal stability due to the fuselage
% contribution
x_XCG_no_fus = x_XCG; % x_XCG with no fuselage contribution
x_XCG_fus = X_NP - (-CL_alpha_w1w2b*SM_des - (CMTalpha + CMalpha_f))*(cmac_w1/(-CL_alpha_w1w2b)); % x_XCG with fuselage contribution
TRIM_RESULTS.x_XCG = x_XCG;

% Revised values
CM_alpha_w1w2b = -CL_alpha_w1w2b*((X_NP-x_XCG_fus)/cmac_w1) + CMTalpha + CMalpha_f;
SM_actual = -CM_alpha_w1w2b/CL_alpha_w1w2b;
% update the x_XCG
x_XCG = x_XCG_fus;
TRIM_RESULTS.x_XCG = x_XCG;
    
% Store DATA
Trim_ITER.CM_alpha_w1w2b_CLalpha = CM_alpha_w1w2b_CLalpha;
Trim_ITER.CM_alpha_w1w2b_fus = CM_alpha_w1w2b_fus;
Trim_ITER.CMTalpha = CMTalpha;
Trim_ITER.CMalpha_f = CMalpha_f;
Trim_ITER.SM_actual_fus = SM_actual_fus;
Trim_ITER.CM_alpha_w1w2b = CM_alpha_w1w2b;
Trim_ITER.SM_actual = SM_actual;

%%%%%%%%%%%%%%%%%%%%%%%CALCULO CM0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CM_0_w1 = Aero.CM_0_w1_CR;
CM_0_w2 = Aero.CM_0_w1_CR;

CM0_w1_e_corrected = eta_w1_afe_S_w1_afe_S_ref*CM_0_w1 + eta_w1_no_afe_S_w1_no_afe_S_ref*CM_0_w1;
CM0_w2_e_corrected = eta_w2_afe_S_w2_afe_S_ref*CM_0_w2 + eta_w2_no_afe_S_w2_no_afe_S_ref*CM_0_w2;

Trim_ITER.CM0_f = CM0_f;
Trim_ITER.CM0_w1_e_corrected = CM0_w1_e_corrected;
Trim_ITER.CM0_w2_e_corrected = CM0_w2_e_corrected*(cmac_w2/cmac_w1);
Trim_ITER.CM0_w1_CL0_w1 = ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CLalpha_w1_corrected*i_w1);
Trim_ITER.CM0_w2_CL0_w2 = ((x_XCG - x_xbar_w2)/cmac_w1)*(CL0_w2_e_corrected + CLalpha_w2_corrected*(i_w2 - eps_w2));
Trim_ITER.CMT1 = CMT1;

CM0_w1w2b = CM0_f + CM0_w1_e_corrected + CM0_w2_e_corrected*(cmac_w2/cmac_w1) + ...
    ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CLalpha_w1_corrected*i_w1) + ...
    ((x_XCG - x_xbar_w2)/cmac_w1)*(CL0_w2_e_corrected + CLalpha_w2_corrected*(i_w2 - eps_w2)) + CMT1;

Trim_ITER.CM0_w1w2b = CM0_w1w2b;

% Estiamtion of rest of lingitudinal Stability derivatives and TRim
% conditions
CM_delta_e = - CL_delta_e*((x_xbar_w2 - x_XCG)/cmac_w1);

Trim_ITER.CL_delta_e= CL_delta_e;
Trim_ITER.CD_delta_e= CD_delta_e;
Trim_ITER.CM_delta_e= CM_delta_e;

%% Resolution of Trim Conditions
num_delta_e = w_T0/(q_inf*S_ref) - CL0_w1w2b + CL_alpha_w1w2b*CM0_w1w2b/CM_alpha_w1w2b;
den_delta_e = -CL_alpha_w1w2b*CM_delta_e/CM_alpha_w1w2b + CL_delta_e;
trim_delta_e = num_delta_e/den_delta_e;
trim_delta_e_deg = trim_delta_e*R2D;
trim_alpha = - (CM_delta_e*trim_delta_e + CM0_w1w2b)/CM_alpha_w1w2b;
trim_alpha_deg = trim_alpha*R2D;

TRIM_RESULTS.X_NP = X_NP;
TRIM_RESULTS.CL0_w1w2b = CL0_w1w2b;
TRIM_RESULTS.CL_alpha_w1w2b = CL_alpha_w1w2b;
TRIM_RESULTS.CL_delta_e = CL_delta_e;
TRIM_RESULTS.CM0_w1w2b = CM0_w1w2b;
TRIM_RESULTS.CM_alpha_w1w2b = CM_alpha_w1w2b;
TRIM_RESULTS.CM_delta_e = CM_delta_e;
TRIM_RESULTS.trim_alpha_deg = trim_alpha_deg;
TRIM_RESULTS.trim_delta_e_deg = trim_delta_e_deg;

% Calculo de valores totales
% Total Lift coefficient
CL_needed = w_T0/(q_inf*S_ref);
CL = CL0_w1w2b + CL_alpha_w1w2b*trim_alpha + CL_delta_e*trim_delta_e;

% CL_w1
CL_w1_b = CL0_w1_e_corrected + CLalpha_w1_corrected*i_w1 + CLalpha_wb_w1*trim_alpha;
% CL_w2
CL_w2 = CL0_w2_e_corrected + CLalpha_w2_corrected*(i_w2 - eps_w2) + CLalpha_w2_e_pw*(downwash)*trim_alpha + CL_delta_e*trim_delta_e;
% Total CL
CL_Total = CL_w1_b + CL_w2;

% Sotore DATA
Trim_ITER.CL = CL;
Trim_ITER.CL_needed = CL_needed;
Trim_ITER.CL_w1_b = CL_w1_b;
Trim_ITER.CL_w2 = CL_w2;
Trim_ITER.CL_Total = CL_Total;

% angle of attack of HTP
alpha_w2 = ((i_w2 - eps_w2) + (downwash)*trim_alpha);
alpha_w2_deg = alpha_w2*R2D;

x_offset_engine = Geo_tier.x_eng_xbar - TRIM_RESULTS.x_XCG ;
Trim_ITER.x_offset_engine = x_offset_engine;

% Total CM - Checking that is zero
CM = CM0_w1w2b + CM_alpha_w1w2b*trim_alpha + CM_delta_e*trim_delta_e;
TRIM_RESULTS.CM = CM;


%% Most Forward XCG
x_XCG_fwd =  -(delta_rudvtr_min + (CM0_w1w2b/CM_delta_e))/(alpha_max_w1_ope*D2R*CL_alpha_w1w2b/(CM_delta_e*cmac_w1)) + X_NP;
%% Most Forward XCG
SM_min = 0.10;
x_XCG_rwd = X_NP - (-CL_alpha_w1w2b*SM_min - (CMTalpha + CMalpha_f))*(cmac_w1/(-CL_alpha_w1w2b)); % x_XCG with fuselage contribution
% Storing DATA
TRIM_RESULTS.x_XCG_fwd = x_XCG_fwd;
TRIM_RESULTS.x_XCG_rwd = x_XCG_rwd;

Trim_ITER.SM_des = SM_des;

% if x_XCG_variation ==1
%     W_i =
%     W_f = 
%     from 


% % Plot conditions
% SAVE_FIGS=0; % if = 1saves figures
% LS = 2; % Lines size
% TS = 10; % Text size
% 
% % defines the variation of mass (fuel being burn)
% Delta_m = 0.050; %kg
% M_vec = (m_TO_1:-Delta_m:(m_TO_1-mF_trayecto));
% 
% % Solves the trim conditions for the entired range of weights
% for i=1:length(M_vec)
%     W_TO_1m = g*M_vec(i);
%     num_delta_e_vec(i) = W_TO_1m/(q_inf*S_w) - CL0_wbh + CL_alpha_wbh*CM0_wbh/CM_alpha_wbh;
%     den_delta_e_vec(i) = -CL_alpha_wbh*CM_delta_e/CM_alpha_wbh + CL_delta_e;
%     trim_delta_e_vec(i) = num_delta_e_vec(i)/den_delta_e_vec(i);
%     trim_delta_e_deg_vec(i) = trim_delta_e_vec(i)*R2D;
%     trim_alpha_vec(i) = - (CM_delta_e*trim_delta_e_vec(i) + CM0_wbh)/CM_alpha_wbh;
%     trim_alpha_deg_vec(i) = trim_alpha_vec(i)*R2D;
%     
%     % Cálculo de Xcg
%     if dead_weight_small == 1
%         Xcg_vec(i) = Calc_Xcg_vs_W_payload(M_vec(i));
%     else
%         Xcg_vec(i) = Calc_Xcg_vs_W(M_vec(i));
%     end
%     
%     % Staticv margin
%     SM_vec(i) = (X_NP - Xcg_vec(i))/cmac_w;
%     
%     % Calculo de valores totales
%     CL_vec(i) = CL0_wbh + CL_alpha_wbh*trim_alpha_vec(i) + CL_delta_e*trim_delta_e_vec(i);
%     
%     % CL del ala
%     CL_wb_vec(i) = CL0_we*conv_we + CLalpha_we*i_w*(D2R)*conv_we + CLalpha_wb*trim_alpha_vec(i);
%     
%     % CL del HTP
%     CL_h_vec(i) = CL0_h*(eta_h_afe*S_h_afe/S_w + eta_h_no_afe*S_h_no_afe/S_w) + ...
%         CLalpha_h*(i_t*D2R - epsilon_0)*(eta_h_afe*S_h_afe/S_w + eta_h_no_afe*S_h_no_afe/S_w) + ...
%         CL_delta_e*trim_delta_e_vec(i) + ...
%         CLalpha_h*(downwash)*(eta_h_afe*S_h_afe/S_w + eta_h_no_afe*S_h_no_afe/S_w)*trim_alpha_vec(i);
%     
%     alpha_h_deg_vec(i) = ((i_t*D2R - epsilon_0) + (downwash)*trim_alpha_vec(i))*R2D;
%     
%     % Solves the real angle of attack from XFLR5 results
%     if CRUISE == 0
%         alpha_we_TRIM_LLT_vec(i) = interp1(DATA(31).CL,DATA(31).alpha,CL_wb_vec(i)*S_w/S_w_e,'spline');
%         alpha_wb_TRIM_LLT_vec(i) = alpha_we_TRIM_LLT_vec(i) - i_w;
%     elseif CRUISE == 1
%         alpha_we_TRIM_LLT_vec(i)  = interp1(DATA(34).CL,DATA(34).alpha,CL_wb_vec(i)*S_w/S_w_e,'spline');
%         alpha_wb_TRIM_LLT_vec(i) = alpha_we_TRIM_LLT_vec(i) - i_w;
%     elseif CRUISE == 2
%         alpha_we_TRIM_LLT_vec(i)  = interp1(DATA(34).CL,DATA(34).alpha,CL_wb_vec(i)*S_w/S_w_e,'spline');
%         alpha_wb_TRIM_LLT_vec(i) = alpha_we_TRIM_LLT_vec(i) - i_w;
%     end
% end
% 
% % vStores results
% Trim_ITER.trim_alpha_m_max_deg = trim_alpha_deg_vec(1);
% Trim_ITER.trim_alpha_m_min_deg = trim_alpha_deg_vec(end);
% Trim_ITER.trim_delta_m_max_deg = trim_delta_e_deg_vec(1);
% Trim_ITER.trim_delta_m_min_deg = trim_delta_e_deg_vec(end);
% Trim_ITER.trim_alpha_wb_m_max_deg = alpha_we_TRIM_LLT_vec(1)*R2D;
% Trim_ITER.trim_alpha_wb_m_min_deg = alpha_we_TRIM_LLT_vec(end)*R2D;
% TRIM_RESULTS.CASO = CASE_iter_num;
% TRIM_RESULTS.m_TO_1 =  m_TO_1;
% TRIM_RESULTS.T = NP.T;
% TRIM_RESULTS.deltaP = NP.deltaP;
% TRIM_RESULTS.V = V;
% TRIM_RESULTS.V_h = V_h;
% TRIM_RESULTS.eta_h_afe = eta_h_afe;
% TRIM_RESULTS.Xcg = Xcg;
% TRIM_RESULTS.X_NP = X_NP;
% TRIM_RESULTS.SM = SM;
% TRIM_RESULTS.CL = CL;
% TRIM_RESULTS.CL_Total = CL_Total;
% TRIM_RESULTS.CL_wb = CL_wb ;
% TRIM_RESULTS.CL_h = CL_h;
% TRIM_RESULTS.trim_alpha_deg = trim_alpha_deg;
% TRIM_RESULTS.alpha_we_TRIM_LLT = alpha_we_TRIM_LLT;
% TRIM_RESULTS.alpha_wb_TRIM_LLT = alpha_wb_TRIM_LLT;
% TRIM_RESULTS.trim_delta_e_deg = trim_delta_e_deg;
% TRIM_RESULTS.alpha_h = alpha_h;
% TRIM_RESULTS_show = TRIM_RESULTS
% 
% % Figure of results of trim analysis
% Fig = 0;
% Fig = Fig + 1;
% figure(Fig)
% plot(M_vec,trim_alpha_deg_vec,'b','LineWidth', LS)
% hold on
% plot(M_vec,trim_delta_e_deg_vec,'g','LineWidth', LS)
% hold off
% title('W vs Trim angles (solution)')
% xlabel('M (kg)')
% ylabel('\alpha_{trim} & \delta_{e_{trim}} (deg)')
% h_legend=legend('\alpha_{trim}','\delta_{e_{trim}}');
% set(h_legend, 'Location','Best','FontSize',TS)
% grid on
% if SAVE_FIGS==1
%     if dead_weight_small == 1
%         prefix = strcat('Cefiro3_Trim_longitudinal_PI_alpha_delta_e_calc_deadweight_0_4');
%     else
%         prefix = strcat('Cefiro3_Trim_longitudinal_PI_alpha_delta_e_calc_deadweight_3_2');
%     end
%     name   = strcat(prefix,sufix);
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end
% 
% Fig = Fig + 1;
% figure(Fig)
% plot(M_vec,alpha_we_TRIM_LLT_vec,'b','LineWidth', LS)
% hold on
% plot(M_vec,trim_delta_e_deg_vec,'g','LineWidth', LS)
% hold off
% title('W vs Trim angles (wing-body)')
% xlabel('M (kg)')
% ylabel('\alpha_{trim} & \delta_{e_{trim}} (deg)')
% h_legend=legend('\alpha_{trim}','\delta_{e_{trim}}');
% set(h_legend, 'Location','Best','FontSize',TS)
% grid on
% if SAVE_FIGS==1
%     if dead_weight_small == 1
%         prefix = strcat('Cefiro3_Trim_longitudinal_PI_alpha_delta_e_real_deadweight_0_4');
%     else
%         prefix = strcat('Cefiro3_Trim_longitudinal_PI_alpha_delta_e_real_deadweight_3_2');
%     end
%     name   = strcat(prefix,sufix);
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end
% 
% Fig = Fig + 1;
% figure(Fig)
% plot(M_vec,SM_vec,'b','LineWidth', LS)
% title('W vs Static Margin')
% xlabel('M (kg)')
% ylabel('SM')
% grid on
% if SAVE_FIGS==1
%     if dead_weight_small == 1
%         prefix = strcat('Cefiro3_Trim_longitudinal_PI_M_vs_SM_deadweight_0_4');
%     else
%         prefix = strcat('Cefiro3_Trim_longitudinal_PI_M_vs_SM_deadweight_3_2');
%     end
%     name   = strcat(prefix,sufix);
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end
% 
% Fig = Fig + 1;
% figure(Fig)
% plot(M_vec,Xcg_vec,'b','LineWidth', LS)
% title('W vs X_{cg}')
% xlabel('M (kg)')
% ylabel('X_{cg}')
% grid on
% if SAVE_FIGS==1
%     if dead_weight_small == 1
%         prefix = strcat('Cefiro3_Trim_longitudinal_PI_M_vs_Xcg_deadweight_0_4');
%     else
%         prefix = strcat('Cefiro3_Trim_longitudinal_PI_M_vs_Xcg_deadweight_3_2');
%     end
%     name   = strcat(prefix,sufix);
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end

