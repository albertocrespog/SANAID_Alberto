function [PERFORMANCE_RESULTS,PERFORMANCE_ITER] = ...
    Calculo_Performance_Dic2019(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
    Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,V_Performance,Data_ATM,XCG_data)

% loop = 0; % Variable de control de ejecuciÃ³n
% while loop == 0;
% 
%     disp('""""""""""""""""""""""""""""""""""""""""""""""""""""""')
%     disp('"                                                    "')
%     disp('"  Performance Analysis                              "')
%     disp('"  Main Menu                                    "')
%     disp('"                                                    "')
%     disp('""""""""""""""""""""""""""""""""""""""""""""""""""""""')
%     disp(' ')
%     disp('[1]. Select Type of Aircrat')
%     disp('[2]. Introduce Number of Segments')
%     disp('[3]. Define type of segments')
%     disp('[4]. Define segments')
%     disp('[5]. Performance Analysis')
%     disp(' ')
%     disp('      [0]. Salir')
%     disp(' ')
%     output=input('	Select option: ');
% 
%     % Tolerancia de color en el seguimiento de lÃ­neas
%     colorTOL = 80;
% 
%     switch output
%         case 1
%             disp('[1]. EMERGENTIA')
%             disp('[2]. PEPIÑO XXL')
%             type_aircraft        = input('Select Type of aircraft: ','s');
%         case 2
%             num_segments        = input('Introduce number of segments: ','s');
%         case 3
%             disp('"""""""""""""""""""""""""""')
%             disp('      TYPE OF SEGMENTS     ')
%             disp('[1]. Conventional Take Off')
%             disp('[2]. VTOL Take Off')
%             disp('[3]. Climb')
%             disp('[4]. VTOL Climb')
%             disp('[5]. Cruise (Range)')
%             disp('[6]. Cruise (Endurance)')
%             disp('[7]. Descent')
%             disp('[8]. VTOL Descent')
%             disp('[9]. VTOL Landing')
%             disp('[10]. Landing')
%             loop = 0;
%             seg = 0;
%             while loop == 0
%                 segment = input('Introduce type of segment (type "q" to finish): ','s');
%                 seg = seg + 1;
%                 if strcmp(segment,'q')
%                     loop = 1;
%                 else
%                     CASE_iter{seg}.segment= segment;                    
%                 end
%             end
%         case 4
% 
%             switch output
%                 case 1 % [1]. Conventional Take Off
%                     segmento{1}.datos.altitude = 0;
%                     segmento{1}.datos.M = 0;
% 
%                 case 2
%                 case 3
%                 case 4
%                 case 5
% 
% 
% 
%         case 0
%             loop = 1;
%     end
% end



%                 for i=1:num_segments
%                 type_segments(i)        = input('Introduce type of segment: ','s');
            
            
% Calculo_Stability_Derivatives_June2019(Aero_Polar_PI, Aero_3D_PI, Aero_2D,...
%     Aero_3D_PI_CLalpha,Geo_calc_ITER_16,conv_UNITS,CG_PI_PD,V_CR,...
%     i_t,i_w,D_heli,NP_PI,CRUISE,PROPULSION,m_TO_1,...
%     Posicion_Palanca,DATA,Trim_ITER,TRIM_RESULTS,Stability_Pamadi,Fig,n_eng)

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
% % if XCG_FF == 1
% %     XCG from design
% %     x_XCG = NP.XCG_sm;
% %     y_XCG = XCG_data.z_XCG;
% %     z_XCG = XCG_data.z_XCG;
% % else
% %     XCG from design
% %     x_XCG = XCG_data.z_XCG;
% %     y_XCG = XCG_data.z_XCG;
% %     z_XCG = XCG_data.z_XCG;
% % end

x_XCG = XCG_data.x_XCG;
y_XCG = XCG_data.y_XCG;
z_XCG = XCG_data.z_XCG;

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
b_w1_e = Geo_tier.b_w1_e;
b_w2 = Geo_tier.b_w2;
b_w2_e = Geo_tier.b_w2_e;
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
y_d_T = Geo_tier.y_d_T;

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
Lambda_LE_w1_e = Geo_tier.Lambda_LE_w1_e;
Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
Lambda_LE_w2_e = Geo_tier.Lambda_LE_w2_e;
Lambda_c4_w2 = Geo_tier.Lambda_c4_w2;
Lambda_c2_w1 = Geo_tier.Lambda_c2_w1;
Lambda_c2_w2 = Geo_tier.Lambda_c2_w2;

dihedral_w1 = Geo_tier.dihedral_w1;
dihedral_w2 = Geo_tier.dihedral_w2;



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
Trim_ITER.x_XCG = x_XCG;
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

Delta_CM_Tprop = - T_avail*z_d_T/(q_inf*S_ref*cmac_w1);
CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;

CTx1 = T_avail/(q_inf*S_ref);
CMT1 = - (d_T/cmac_w1)*CTx1;

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
CM_alpha_w1w2b_fus = CM_alpha_w1w2b_CLalpha + CMTalpha + CMalpha_f; % Total CMalpha with fuselage contribution
SM_actual_fus = -CM_alpha_w1w2b_fus/CL_alpha_w1w2b; % Real SM with fuselage contribution

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
TRIM_RESULTS.trim_alpha = trim_alpha;
TRIM_RESULTS.trim_delta_e = trim_delta_e;
TRIM_RESULTS.trim_alpha_deg = trim_alpha_deg;
TRIM_RESULTS.trim_delta_e_deg = trim_delta_e_deg;
% Euler angle
theta1 = trim_alpha;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CDalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CD_alpha = (C_D1 + 2*C_D2*CL)*CL_alpha_w1w2b; 
CXalfa = CL-CD_alpha;
CD = C_D0 + C_D1*CL + C_D2*CL^2;  %polar aeronave
CZalfa = - CL_alpha_w1w2b - CD;
CMalfa = CM_alpha_w1w2b;

Stab_Der.CD_alpha = CD_alpha;
Stab_Der.CD = CD;
Stab_Der.CXalfa = CXalfa;
Stab_Der.CZalfa = CZalfa;
Stab_Der.CMalfa = CMalfa;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPEED DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cdu = 0; 
% Cxu = -2*CD - Cdu; 
% Clu = 0;                      %para vuelo subsónico  
% Czu = -2*CL - Clu - 2*m1*q0;
% Czu = -2*CL - Clu;
% CMu = 0;                      %0 para vuelo subsónico
% Cmu = 2*CM + CMu;                      %0 para vuelo subsónico
CDu = 0; 
CXu = -3*CD - CDu; 
% CXu = -2*CD - CDu; 
CLu = 0;                      %para vuelo subsónico  

Stab_Der.CDu = CDu;
Stab_Der.CXu = CXu;
Stab_Der.CLu = CLu;
    
% Czu = -2*CL - Clu - 2*m1*q0;
CZu = -2*CL - CLu;
% Cmu = 2*CM + CMu;                      %0 para vuelo subsónico
CMu = 0;                      %0 para vuelo subsónico

Stab_Der.CZu = CZu;
Stab_Der.CDu = CDu;
Stab_Der.CLu = CLu;
Stab_Der.CMu = CMu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Pitch Derivatives%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CXq = 0;
% Tail Volume Coefficient
V_h = ((x_xbar_w2 - x_XCG)/cmac_w1)*(S_w2_e/S_ref);
V_t = ((x_xbar_w1 - x_XCG)/cmac_w1)*(S_w2_e/S_ref);
V_vee = ((x_xbar_w1 - x_XCG)/cmac_w1)*(S_w2_e/S_ref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%CLq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_xbar_w1_e = Geo_tier.x_xbar_w1_e;
x_xbar_w2_e = Geo_tier.x_xbar_w2_e;

B_prandalt = sqrt(1 - (M^2)*(cos(Lambda_c4_w1))^2);     %eq 4.510 PAMADI
chi = (x_xbar_w1_e - x_XCG)/cmac_w1_e; %eq 4.490 PAMADI
% CLalpha_wb_w1 already includes (KWB_w1 + KBW_w1)
CLq_e = CLalpha_wb_w1*(0.5 + chi);   %eq 4.489 PAMADI
CLalpha_fus = 2*f_k2_k1*Area_b_max/(Vol_TOT^(2/3));  %eq 4.495 PAMADI
CLalpha_Bp = CLalpha_fus*(Vol_TOT^(2/3))/Area_b_max; %eq 4.494 PAMADI
CLq_B = 2*CLalpha_Bp*(1 + x_XCG/length_fus);           %eq 4.493 PAMADI

% W1 and w2 with no dynamic pressure correction
% CLq_wb_w1 = (KWB_w1 + KBW_w1)*((S_w1_e*cmac_w1_e)/(S_ref*cmac_w1))*CLq_e + ...
%     CLq_B*((length_fus*Area_b_max)/(S_w1*cmac_w1)); %eq 4.488 PAMADI
CLq_wb_w1 = ((S_w1_e*cmac_w1_e)/(S_ref*cmac_w1))*CLq_e + ...
    CLq_B*((length_fus*Area_b_max)/(S_w1*cmac_w1)); %eq 4.488 PAMADI
CLq_w2 = 2*CLalpha_w2_e*((x_xbar_w2_e - x_XCG)/cmac_w1); %eq 10.72 ROSKAM

% W1 and w2 with no dynamic pressure correction
CLq_wb_w1_pw = eta_w1_afe*CLq_wb_w1 + eta_w1_no_afe*CLq_wb_w1;
CLq_w2_pw = eta_w2_afe*CLq_w2 + eta_w2_no_afe*CLq_w2;

% Total Derivative CLq
CLq = CLq_wb_w1_pw + CLq_w2_pw;
CZq = -CLq;
% Storing DATA
Stab_Der.CLq_w2 = CLq_w2;
Stab_Der.CLq_wb_w1 = CLq_wb_w1;
Stab_Der.CLq_w2 = CLq_w2_pw;
Stab_Der.CLq_wb_w1 = CLq_wb_w1_pw;

Stab_Der.CXq = CXq;
Stab_Der.CLq = CLq;
Stab_Der.CZq = CZq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Cmq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Centroid fuselage
int_Sb1 = @(x) interp1(x_Area_body,Area_body,x).*(x);               %eq 4.514 PAMADI
x_c_fus = (1/Vol_TOT)*quad(@(x) int_Sb1(x),0,length_fus);           %eq 4.514 PAMADI
x_m = x_XCG;
xm_1 = x_m/length_fus;                      %eq 4.513 PAMADI
xc_1 = x_c_fus/length_fus;                  %eq 4.513 PAMADI
VB_1 = Vol_TOT/(Area_b_max*length_fus);     %eq 4.513 PAMADI

% Cálculo del centro aerodinamico del fuselaje
% Derivada del área:
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
% Ojo revisar por la contribución de la presión dinámica
CMq_wb_w1_part1 = (KWB_w1 + KBW_w1)*((cmac_w1_e/cmac_w1)^2)*CMq_e; %eq 4.502 PAMADI
CMq_wb_w1_part2 = CMq_B*(Area_b_max/S_w1)*((length_fus/cmac_w1)^2);                   %eq 4.502 PAMADI
CMq_wb_w1 = eta_w1_afe_S_w1_afe_S_ref*CMq_wb_w1_part1 + eta_w1_no_afe_S_w1_no_afe_S_ref*CMq_wb_w1_part1;
CMq_w2 = - CLq_w2*((x_xbar_w2 - x_XCG)/cmac_w1)^2;              %eq 10.78 ROSKAM

% W1 and w2 with no dynamic pressure correction
%% Revisar OJO - contribuición de presiópn dinámica en ala w1
CMq_wb_w1_pw = CMq_wb_w1 + CMq_wb_w1_part2;
CMq_w2_pw = eta_w2_afe_S_w2_afe_S_ref*CMq_w2 + eta_w2_no_afe_S_w2_no_afe_S_ref*CMq_w2;

% Total Derivative CLq
CMq = CMq_wb_w1_pw + CMq_w2_pw;

% Storing DATA
Stab_Der.CMq_w2 = CMq_w2;
Stab_Der.CMq_wb_w1 = CMq_wb_w1;
Stab_Der.CMq_w2 = CMq_w2_pw;
Stab_Der.CMq_wb_w1 = CMq_wb_w1_pw;
Stab_Der.CMq = CMq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Alpha punto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CXalfapunto = 0;              %vale 0 para vuelo subsónico
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% CL_alphapunto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = sqrt(1 - M^2);
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
CLalpha_B = 2*f_k2_k1*(Area_b_max/Vol_TOT^(2/3));      %eq 4.532 PAMADI
CLalpha_B_p = CLalpha_B*(Vol_TOT^(2/3)/(Area_b_max));    %eq 4.531 PAMADI
CL_alphapunto_B = 2*CLalpha_B_p*(Vol_TOT/(length_fus*Area_b_max));  %eq 4.530 PAMADI

% W1 and w2 with no dynamic pressure correction %eq 4.527 PAMADI
% Ojo revisar por la contribución de la presión dinámica
CL_alphapunto_wb_w1_part1 = (KWB_w1 + KBW_w1)*(cmac_w1_e/cmac_w1)*CL_alphapunto_e;    %eq 4.527 PAMADI
CL_alphapunto_wb_w1_part2 = CL_alphapunto_B*((length_fus*Area_b_max)/(S_ref*cmac_w1));      %eq 4.527 PAMADI
CL_alphapunto_wb_w1 = eta_w1_afe_S_w1_afe_S_ref*CL_alphapunto_wb_w1_part1 +...
    eta_w1_no_afe_S_w1_no_afe_S_ref*CL_alphapunto_wb_w1_part1;
% Modified version Pamadi to account for different dynamic pressure
l_w2 = ((x_xbar_w2 - x_XCG)/cmac_w1);
% W1 and w2 with no dynamic pressure correction
CL_alphapunto_w2 = 2*CLalpha_w2_e*l_w2*deps_dalpha; %eq 4.525 PAMADI

% W1 and w2 with no dynamic pressure correction
%% Revisar OJO - contribuición de presiópn dinámica en ala w1
CL_alphapunto_wb_w1_pw = CL_alphapunto_wb_w1 + CL_alphapunto_wb_w1_part2;
CL_alphapunto_w2_pw = eta_w2_afe_S_w2_afe_S_ref*CL_alphapunto_w2 + ...
    eta_w2_no_afe_S_w2_no_afe_S_ref*CL_alphapunto_w2;

% Total Derivative CLalphapunto
CLalphapunto = CL_alphapunto_wb_w1_pw + CL_alphapunto_w2_pw;
% Storing DATA
Stab_Der.CL_alphapunto_w2_pw = CL_alphapunto_w2_pw;
Stab_Der.CL_alphapunto_wb_w1_pw = CL_alphapunto_wb_w1_pw;
Stab_Der.CLalphapunto = CLalphapunto;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CM_alphapunto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CM_alphapunto_B = 2*CMalpha_p*((xc_1 - xm_1)/(1 - xm_1 - VB_1))*(Vol_TOT/(length_fus*Area_b_max)); %eq 4.543 PAMADI
% Only valid for tau < 4
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
% Ojo revisar por la contribución de la presión dinámica
CM_alphapunto_wb_w1_part1 = (KWB_w1 + KBW_w1)*(cmac_w1_e^2/cmac_w1^2)*CM_alphapunto_e;    %eq 4.539 PAMADI
CM_alphapunto_wb_w1_part2 = CM_alphapunto_B*((Area_b_max*length_fus^2)/(S_w1*cmac_w1^2));      %eq 4.539 PAMADI
CM_alphapunto_wb_w1 = eta_w1_afe_S_w1_afe_S_ref*CM_alphapunto_wb_w1_part1 +...
    eta_w1_no_afe_S_w1_no_afe_S_ref*CM_alphapunto_wb_w1_part1;
CM_alphapunto_w2 = - CL_alphapunto_w2*((x_xbar_w2 - x_XCG)/cmac_w1);      %eq 4.537 PAMADI

% W1 and w2 with no dynamic pressure correction
%% Revisar OJO - contribuición de presiópn dinámica en ala w1
CM_alphapunto_wb_w1_pw = CM_alphapunto_wb_w1 + CM_alphapunto_wb_w1_part2;
CM_alphapunto_w2_pw = eta_w2_afe_S_w2_afe_S_ref*CM_alphapunto_w2 + ...
    eta_w2_no_afe_S_w2_no_afe_S_ref*CM_alphapunto_w2;

% Total Derivative CMalphapunto
CMalphapunto = CM_alphapunto_wb_w1_pw + CM_alphapunto_w2_pw; %eq 4.538 PAMADI
% Storing DATA
Stab_Der.CM_alphapunto_w2_pw = CM_alphapunto_w2_pw;
Stab_Der.CM_alphapunto_wb_w1_pw = CM_alphapunto_wb_w1_pw;
Stab_Der.CMalphapunto = CMalphapunto;

CZalfapunto = -CLalphapunto; %OJO 
Stab_Der.CXalfapunto = CXalfapunto;
Stab_Der.CZalfapunto = CZalfapunto;

CZteta = -CL*sin(theta1);
CXteta = -CL*cos(theta1); 
CMteta = 0;

Stab_Der.CXteta=CXteta;
Stab_Der.CZteta=CZteta;
Stab_Der.CMteta=CMteta;

% Control derivatives 
CXdeltae = 0;
CDdeltae = -CXdeltae;
CZdeltae = -CL_delta_e;
CMdeltae = CM_delta_e;

Stab_Der.CL_delta_e = CL_delta_e;
Stab_Der.CXdeltae = CXdeltae;
Stab_Der.CDdeltae = CDdeltae;
Stab_Der.CZdeltae = CZdeltae;
Stab_Der.CMdeltae = CMdeltae;

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

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CTx1 %%%%%%%%%%%%%%%%%%%%%%%%%
CTx1 = CD;

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CTxu %%%%%%%%%%%%%%%%%%%%%%%%%
%% For future versions
alpha =0;
beta = 0;
% Associated to propulsive model
% [Propulsion] = get_EngineProperties_v3(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,n_eng)
D_prop = Prop_data.D_prop;
RPS = Propulsion.RPS;
d_CP_d_V = Propulsion.d_CP_d_V;
adimensional_P_Heli =rho*(RPS^3)*(D_prop^5); 
adimensional_P_Aero =q_inf*S_ref; 
adimensional_conv = adimensional_P_Aero/adimensional_P_Heli;
dcP_du = d_CP_d_V*adimensional_conv;
% ctxu = (1/(q_inf*qmet2qimp*S_w1*m22ft2))*(2*A_power*V + B_power)*W2pftsec;
CTxu = dcP_du - 2*CTx1;  %por helice a paso fijo
% CTxu = - 3*CTx1       %por helice a paso variable

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
Stab_Der.CM = CM;
Stab_Der.CL_alpha_w1w2b = CL_alpha_w1w2b;
Stab_Der.CM_alpha_w1w2b = CM_alpha_w1w2b;

Stab_Der.T = Propulsion.Ti;

% % Estimation through Fran's Formulation  if Stability_Pamadi = 0;
% % Estimation through Pamadi if Stability_Pamadi = 1;
% % Estimation through ROSKAM if Stability_Pamadi = 2;
% 
% Stability_Pamadi = 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stab_Dyn_Long = longitudinal_analysis(V,rho,S_w1,b_w1,cmac_w1,m_TOW,Stab_Der,...
%     q_inf,theta1,conv_UNITS,Stability_Pamadi,Weight_tier);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%COEFICIENTES Stab_Der ESTABILIDAD ESTÁTICA Stab_DerLATERAL.DIRECCIONAL %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% q_v = q_h_no_afe;
% eta_w2_no_afe = 1
% eta_w2_no_afe = q_v / q_inf;
% Arm from xcg to Xac_w2
l_xcg_w2 = x_xbar_w2 - x_XCG;
% Arm from Zcg to Zac_w2
l_zcg_w2 = z_zbar_w2 - z_XCG;
z_v = (l_zcg_w2*cos(trim_alpha) - l_xcg_w2*sin(trim_alpha));
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
modelo.general.mtow = m_TOW;
modelo.general.w_w0 = 1;     
modelo.general.Sref = S_w1;
modelo.general.Minf = M;
modelo.general.Vinf = V;
modelo.general.rhoinf = rho; 
modelo.general.Xcg = x_XCG;
modelo.general.Zcg = z_XCG;
modelo.general.CLa_WB = CLalpha_wb_w1;
modelo.general.downwash = deps_dalpha;

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
modelo.ala.Cla = CLalpha_wb_w1;
modelo.ala.CLa = CLalpha_wb_w1;
modelo.ala.CLa_we = CLalpha_wb_w1;
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

modelo.vertical.S = S_w2;
modelo.vertical.eta = eta_w2_no_afe;
modelo.vertical.Cla = CLalpha_wb_w2;
modelo.vertical.CLa = CLalpha_wb_w2;
modelo.vertical.b = b_w2;
modelo.vertical.Xca = x_xbar_w2;
modelo.vertical.Zca = l_zcg_w2;
modelo.vertical.AR = AR_w2;
modelo.vertical.TR = lambda_w2;
modelo.vertical.cm_c = 0.25*cmac_w2;
modelo.vertical.t_c = 0.12; %naca 0012

z_1R_y1_rudvtr = Geo_tier.z_1R_y1_rudvtr;
z_1R_y2_rudvtr = Geo_tier.z_1R_y2_rudvtr;

% modelo.vertical.y0_b2 = z_rudder_bottom;
modelo.vertical.y0_b2 = z_1R_y1_rudvtr;
% modelo.vertical.y1_b2 = z_rudder_top;
modelo.vertical.y1_b2 = z_1R_y2_rudvtr;

%%%%%%%%%%%%%%%%%%%%%derivadas en funcion de beta%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Cyb
% %Cybw
% Cyb_w_CL = (CL^2)*((6*tan(Lambda_c4_w1)*sin(Lambda_c4_w1))/(pi*AR_w*(AR_w + 4*cos(Lambda_c4_w1))));   
% Cyb_w_dihedro = -0.00573*abs(Sigma_w);
% Cyb_w = Cyb_w_CL + Cyb_w_dihedro;
% Stab_Der.Cyb_w_CL = Cyb_w_CL;
% Stab_Der.Cyb_w_dihedro = Cyb_w_dihedro;
% % Cyb fuselaje
% % Figure Smetana Fig 7
% Delta_y = 1.485 - 1;
% Delta_x = 1;
% m_figB6 = Delta_y/Delta_x; 
% Ki = 1 + m_figB6*(abs(z_ala_centerline/(h_Area_b_max/2)));
% CLalpha_fus = 2*f_k2_k1*Area_b_max/(Vol_TOT^(2/3));
% Cyb_fus = - Ki*CLalpha_fus*(Vol_TOT^(2/3))/S_w;       
% %Cybv
% % in_k_Cyb = b_v / 
%  k = 1.0; % Fig 3.75 Pamadi
% q_v = q_h_no_afe;
% eta_v = q_v / q_inf;
% % w_Area_b_max = Geo_calc_ITER_16.w_Area_b_max;
% % h_Area_b_max = Geo_calc_ITER_16.h_Area_b_max;
% sidewash = 0.724 + (3.06/(1+cos(Lambda_c4_w1)))*((S_v)/S_w) + ...
%     0.4*(-z_ala_centerline/h_Area_b_max) + 0.0009*AR_w;     
% Cyb_1v = -k*CLalpha_v*sidewash*((S_v)/S_w);                 
% Cyb_v = Cyb_1v;
% Cyb = Cyb_w + Cyb_fus + Cyb_v;                    

% ASPro Cyb
Stab_Der = getCybeta(modelo,Stab_Der);
Cyb_w = Stab_Der.Cyb_w;
Cyb_fus = Stab_Der.Cyb_fus;
Cyb_v = Stab_Der.Cyb_v;
Cyb = Stab_Der.Cyb;

% % Clb
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

% ASPro Clb
Stab_Der = getClbeta(modelo,trim_alpha,Stab_Der);
Clb_wb = Stab_Der.Clb_wb;
Clb_v = Stab_Der.Clb_v;
Clb = Stab_Der.Clb;

% % Cnb
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

% ASPro Cnb
Stab_Der = getCnbeta(modelo,Stab_Der);
% Stab_Der.Cnb_diedro = Cnb_diedro;
% Stab_Der.Cnb_flecha = Cnb_flecha;
Cnb_b = Stab_Der.Cnb_b;
Cnb_w = Stab_Der.Cnb_w;
Cnb_wb = Stab_Der.Cnb_wb; 
Cnb_v = Stab_Der.Cnb_v;
Cnb = Stab_Der.Cnb;

%%%%%%%%%%%%%%%derivadas propulsivas laterales-direccionales%%%%%%%%%%%%%%%
% % Data from NACA REport 640
% w_R_30 = 0.0525*2;
% w_R_60 = 0.073*2;
% w_R_90 = 0.045*2;
% Nprop = 1;
% D_prop = D_heli;
% CNalpha_p_KN = 0.1;
% KN = 262*w_R_30 + 262*w_R_60 + 135*w_R_90;
% dCN_dalpha = CNalpha_p_KN*(1 + 0.8*(KN/80.7) - 1);
% CyTb = (pi/4)*Nprop*(D_prop^2)*dCN_dalpha/S_w;
% % Thrust lines inclination angle 
% psi_T = 0;
% l_prop = (x_XCG-x_m_propeller)*cos(psi_T);
% CNTb = (pi/4)*Nprop*l_prop*(D_prop^2)*dCN_dalpha/(S_w*b_w);
% Stab_Der.CyTb = CyTb;
% Stab_Der.CNTb = CNTb;

%ASpro CyTb
Stab_Der = getCyTbeta(modelo,Stab_Der);
CyTb = Stab_Der.CyTb;

%Aspro CNTb
Stab_Der = getCnTbeta(modelo,Stab_Der);
CNTb = Stab_Der.CNTb;

%%%%%%%%%%%%%%%%%%%%%%%control alerones%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cydeltaa
% Cydeltaa=0;      

%ASpro Cydeltaa
Stab_Der = getCyda(modelo,Stab_Der);
Cydeltaa = Stab_Der.Cydeltaa;

% % Cldeltaa
% METODO_AIL = 1;
% if METODO_AIL == 0;
%     CLalpha_v_2D = 5.9656;
%     Cldelta_theory = 4.1666;
%     Cldelta_Cldelta_theory = 0.45;
%     K1 = 1.04;
%     K2 = 0.0714;
%     Cldeltaa = Cldelta_theory*Cldelta_Cldelta_theory*(CLalpha_v_2D/CLalpha_v)*K1*K2;
% elseif METODO_AIL == 1;
% %    Método II
%     tau_a = 0.45;
% %   ail_int = b_w/2 - b_ail;
%     ail_int = b_w/2 - 0.60;
%     ail_ext = b_w/2;
%     int_ca = @(x) cmac_w1.*x;
%     Cldeltaa = (2*CLalpha_we*tau_a/(S_w*b_w))*quad(@(x)int_ca(x),ail_int,ail_ext);
% end

%ASpro Cldelta
Stab_Der = getClda(modelo,Stab_Der);
Cldeltaa = Stab_Der.Cldeltaa;

% % Cndeltaa
% eta = ((b_w/2)- b_ail)/(b_w/2);
% K11=-0.15; 
% Cndeltaa = 2*K11*CL*Cldeltaa;

%ASpro Cndeltaa
Stab_Der = getCnda(modelo,Stab_Der);
Cndeltaa = Stab_Der.Cndeltaa;
%%%%%%%%%%%%%%%%%%%%%%%control rudder%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %derivadas en deltar
% % Tau is for surface relation
% SR_SV = (0.435*0.089)/(S_v/2);
% tau_r_S = 0.44;% 25% cuerda 
% % tau_r = 0.6; % 40% cuerda 
% % tau_r = 0.5; % 30% cuerda 
% % tau_r = 0.65;
% Cydeltar = tau_r_S*CLalpha_v*(S_v/S_w);
% Cldeltar = Cydeltar*(z_vt1/b_w);
% Cndeltar = -Cydeltar*((x_v_xbar - x_XCG)/b_w)*eta_h_no_afe;

%ASpro derivadas en deltar
Stab_Der = getCydr(modelo,Stab_Der);
Cydeltar = Stab_Der.Cydeltar;
Stab_Der = getCldr(modelo,trim_alpha,Stab_Der);
Cldeltar = Stab_Der.Cldeltar;
Stab_Der = getCndr(modelo,trim_alpha,Stab_Der);
Cndeltar = Stab_Der.Cndeltar;

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

save Stability_Derivatives_v3.mat Stab_Der
%save Stablity_Derivatites_16.mat Stab_Der
