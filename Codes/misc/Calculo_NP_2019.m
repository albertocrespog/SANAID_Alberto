function [NP,Propulsion] = Calculo_NP_2019(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,Body_Geo,...
    Design_criteria,Posicion_Palanca,SM_des,Data_ATM,XCG_data)

rho = Data_ATM.rho;
a = Data_ATM.a;

alpha_f = conditions.alpha_f;
beta_f = conditions.beta_f;
h = conditions.h;
V = conditions.V;

% Number of engines
n_eng = Prop_data.n_eng;

% XCG from design
x_XCG = XCG_data.x_XCG;
y_XCG = XCG_data.y_XCG;
z_XCG = XCG_data.z_XCG;

% Design Selected incidences
i_w1 = Design_criteria.i_w1;
i_w2 = Design_criteria.i_w2;

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

Height_x_fus = Body_Geo.Height;
Bisectriz_z_x_fus = Body_Geo.Bisectriz_z;
Centroid_z_x_fus = Body_Geo.Centroid_z;
Angle_fus_x_fus = Body_Geo.Angle_fus;
Angle_fus_interp = Body_Geo.Angle_fus_interp;
x_interp = Body_Geo.x_interp;

g = conv_UNITS.g;
in2m = conv_UNITS.in2m;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;

m_TOW = Weight_tier.m_TOW;
w_T0 = m_TOW*g;

NP.m_TO = m_TOW;

S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
cmac_w1 = Geo_tier.cmac_w1;

S_w2 = Geo_tier.S_w2; 
b_w2 = Geo_tier.b_w2;
cmac_w2 = Geo_tier.cmac_w2;

S_w1_e = Geo_tier.S_w1_e;
S_w2_e = Geo_tier.S_w2_e;
% b_w_e = Geo_tier.b_w1_e;
cmac_w1_e = Geo_tier.cmac_w1;
cmac_w2_e = Geo_tier.cmac_w2;
AR_w1_e = Geo_tier.AR_w1_e;
AR_w2_e = Geo_tier.AR_w2_e;
lambda_w1_e = Geo_tier.lambda_w1_e;
Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
lambda_w2_e = Geo_tier.lambda_w2_e;
Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;

S_ref = Geo_tier.S_ref;

% Position of the Xac for both CR and TOLD
x_xbar_w1_b = Geo_tier.x_xbar_w1;
x_xbar_w1_e = Geo_tier.x_xbar_w1;
xbar_w1 = Geo_tier.xbar_w1;
xbar_w1_e = Geo_tier.xbar_w1;
x_xbar_w2 = Geo_tier.x_xbar_w2;

Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
lambda_w1 = Geo_tier.lambda_w1;
l_xac_w1w2 = Geo_tier.l_xac_w1w2;

% x_h_LE = Geo_tier.x_w2_LE;
z_w2_LE = Geo_tier.z_w2_LE;
% z_v_LE = Geo_tier.z_w2_LE;

x_w1_LE = Geo_tier.x_w1_LE;
z_w1_LE = Geo_tier.z_w1_LE;

% Condiciones de vuelo
C_D0 = Aero_TH.CD0;
C_D1 = Aero_TH.CD1;
C_D2 = Aero_TH.CD2;

% Cruise
CLalpha_w1 = Aero.CL_alpha_w1_CR;
CLalpha_w1_e = CLalpha_w1;
CLalpha_w2 = Aero.CL_alpha_w2_CR;
CLalpha_w2_e = CLalpha_w2;

% Mach Number
M = V/a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CONFIGURACIÓN PUSHER. MODELO DE HELICÓPTERO%%%%%%%%%%%%%%%
D_prop = Prop_data.D_prop; 
R_heli = D_prop/2;
S_heli = (pi*R_heli^2);              %superficie de la hélice

% Velocidad inducida
% Tracción igual a la resistencia
q_inf = 0.5*rho*V^2;                      %presión dinámica inicial

%%%%% NUEVO MODELO PROPULSIVO %%%%
% [T,A,B,C] = propulsive_model_16(m_TO_1,q_inf,S_w,V,rho,D_heli,Aero_Polar,n_eng)

% Estimation of Desired Thrust
CL = w_T0/(q_inf*S_ref);
CD = C_D0 + C_D1*CL + C_D2*CL^2;
D = CD*q_inf*S_ref;
Fdes = D;
% [Propulsion] = get_EngineProperties_v2(V,h,Fdes,n_eng);
[Propulsion] = get_EngineProperties_v3(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,n_eng);
T_tot = Propulsion.Ti;
T_eng = Propulsion.Ti_eng;

% Estimation of data 
% v_i = -(1/2)*V + sqrt(1/4*V^2 + T_eng/(2*rho*S_heli));
% R_inf = R_heli*sqrt((V + v_i)/(V + 2*v_i));
% Propeller variables
v_i = Propulsion.v_i; % Induced velocity at prop disk
R_inf = Propulsion.R_inf; % Radius of prop wash at infinity
v_inf = Propulsion.v_inf; % induced velocity at propwash at infinity

z_zbar_w2 = Geo_tier.z_zbar_w2;
z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
x_xbar_w2 = Geo_tier.x_xbar_w2;
x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
cR_w1 = Geo_tier.cR_w1;
CD0_w1 = Aero_TH.CD0_w1;
alpha_CL_0_w1 = Aero.alpha_CL_0_w1_CR;

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

% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection

%la estela de los propulsores no afecta el HTP
S_w1_pw = Geo_tier.S_w1_pw;
S_w2_pw = Geo_tier.S_w2_pw;

S_w1_afe = S_w1_pw;
S_w2_afe = S_w2_pw;

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

%% Correction of Lift curve slope of w1 associated to pfuselage interaction
% Pendiente de sustentación del morro.
Fineness_Ratio = length_fus/w_Area_b_max;
%digitaliazacion figura PAMADI CAP3
x_f_k2_k1 = [4.,5.,6.,8.,10.,12.,14.,16.,18.,20.];
y_f_k2_k1 = [.77,.825,.865,.91,.94,.955,.965,.97,.973,.975];
f_k2_k1  = interp1(x_f_k2_k1,y_f_k2_k1,Fineness_Ratio,'spline');    
x1 = x_Area_b_max;

%% W1
aN_w1 = 2*f_k2_k1*(Area_b_max/S_w1);   %eq 3.26 PAMADI
KN_w1 = (aN_w1/CLalpha_w1_e_pw)*(S_w1/S_w1_e);  %eq 3.25 PAMADI    
CLalpha_fus = aN_w1;
KWB_w1 = 0.1714*(w_Area_b_max/b_w1)^2 + 0.8326*(w_Area_b_max/b_w1) + 0.9974; %eq 3.27 PAMADI
KBW_w1 = 0.7810*(w_Area_b_max/b_w1)^2 + 1.1976*(w_Area_b_max/b_w1) + 0.0088; %eq 3.28 PAMADI

a_bw_w1 = KBW_w1*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.34 PAMADI
a_wb_w1 = KWB_w1*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.33 PAMADI
a_WB_w1 = (KN_w1 + KWB_w1 + KBW_w1)*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.24 PAMADI
CLalpha_WB_w1 = a_WB_w1;
CLalpha_wb_w1 = CLalpha_WB_w1;

%% W2
aN_w2 = 2*f_k2_k1*(Area_b_max/S_w2);   %eq 3.26 PAMADI
KN_w2 = (aN_w2/CLalpha_w2_e_pw)*(S_w2/S_w2_e);  %eq 3.25 PAMADI    
CLalpha_fus = aN_w2;
KWB_w2 = 0.1714*(w_Area_b_max/b_w2)^2 + 0.8326*(w_Area_b_max/b_w2) + 0.9974; %eq 3.27 PAMADI
KBW_w2 = 0.7810*(w_Area_b_max/b_w2)^2 + 1.1976*(w_Area_b_max/b_w2) + 0.0088; %eq 3.28 PAMADI

a_bw_w2 = KBW_w1*CLalpha_w2_e_pw*(S_w2_e/S_w2); %eq 3.34 PAMADI
a_wb_w2 = KWB_w1*CLalpha_w2_e_pw*(S_w2_e/S_w2); %eq 3.33 PAMADI
a_WB_w2 = (KN_w2 + KWB_w2 + KBW_w1)*CLalpha_w2_e_pw*(S_w2_e/S_w2); %eq 3.24 PAMADI
CLalpha_WB_w2 = a_WB_w2;
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
% Variable independiente = bf_max / b
    
% Cálculo del centro aerodinamico del fuselaje
% Derivada del área:
dSdX = diff(Area_body)./diff(x_Area_body);
% l_N  = x_w_LE + mamparo_F_front;
l_N = x_Area_b_1_4; % estimation 
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
end
    
Xac_cre_WB = xbar_w1/cmac_w1_e; %eq 3.38 PAMADI
XAC_WB_cre = ((Xac_cre_N*aN_w1 + Xac_cre_WB*a_wb_w1 + Xac_cre_BW*a_bw_w1)/a_WB_w1)*cmac_w1_e; %eq 3.36 PAMADI
XAC_WB_LE = XAC_WB_cre*(cmac_w1_e/cmac_w1);                                        %eq 3.36 PAMADI
XAC_WB = XAC_WB_LE + x_w1_LE;                                                     %eq 3.36 PAMADI  

CL_alpha_w1w2b = CLalpha_wb_w1 + CLalpha_w2_e_pw*(downwash);
%eq PAMADI (sin canard)
X_NP = (CLalpha_wb_w1*XAC_WB + CLalpha_w2_e_pw*(downwash)*x_xbar_w2)/CL_alpha_w1w2b;

% Estimation of XCG for a desired SM
XCG_sm = X_NP - SM_des*cmac_w1;
% Estimation of SM for a desired XCG
SM_actual = (X_NP - x_XCG)/cmac_w1;

NP.CLalpha_wb_w1 = CLalpha_wb_w1;
NP.CLalpha_wb_w2 = CLalpha_wb_w2;
NP.CL_alpha_w1w2b = CL_alpha_w1w2b;
NP.XAC_WB = XAC_WB;
NP.x_xbar_w2 = x_xbar_w2;
NP.X_NP = X_NP;
NP.x_XCG = x_XCG;
NP.XCG_sm = XCG_sm;
NP.SM_des = SM_des;
NP.SM_actual = SM_actual;