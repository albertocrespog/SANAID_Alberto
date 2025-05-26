function [Trim_ITER_LAT] = Calculo_Trim_ITER_LAT_v4B(conv_UNITS,conditions_TRIM_lat,...
    Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,Storing_STABILITY_DATA_1,conditions,case_AC,OUTPUT_read_XLSX)

Geo_tier = Storing_GEO_DATA.Geo_tier;
Performance = Storing_AERO_DATA.Performance;
Aero = Storing_AERO_DATA.Aero;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;
Stab_Der = Storing_STABILITY_DATA_1.Stab_Der;

% Study conditions
direction_beta = conditions.conditions_TRIM_lat.direction_beta; 
beta = conditions_TRIM_lat.beta;
beta = direction_beta*beta; % corrects direction of beta
beta_vec = conditions_TRIM_lat.beta_vec;
beta_vec = direction_beta*conditions_TRIM_lat.beta_vec; % corrects direction of beta

direction_phi = conditions.conditions_TRIM_lat.direction_phi; 
phi = conditions_TRIM_lat.phi;
phi = direction_beta*phi; % corrects direction of beta
phi_vec = conditions_TRIM_lat.phi_vec;
phi_vec = direction_beta*conditions_TRIM_lat.phi_vec; % corrects direction of beta

% variable deltaaT
% deltaaT_i = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_TT_ail_min; % Minimum Trim Tab aileron deflection
% deltaaT_f = OUTPUT_read_XLSX.InputGeometry_Data_flags.delta_TT_ail_max; % Maximum Trim Tab aileron deflection
% N_Delta_deltaaT = OUTPUT_read_XLSX.Stability_flags.N_Delta_beta;
deltaaT_vec = conditions_TRIM_lat.deltaaT_vec;

K_y1_TT_ail_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_TT_ail_w1;
K_y2_TT_ail_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_TT_ail_w1;
cf_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_ail;

% Inertias
Ixx = Weight_tier.Ixx;
Iyy = Weight_tier.Iyy;
Izz = Weight_tier.Izz;
Ixz = Weight_tier.Ixz;

V = conditions_TRIM_lat.V_TO;
m_TOW = conditions_TRIM_lat.m_TOW;

% Performance
% V = Performance.V;
% V = V_TO;
rho = conditions_TRIM_lat.rho;

% Geometric data
S_ref = Geo_tier.S_ref;
S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
cmac_w1 = Geo_tier.cmac_w1;
% m_TOW = Weight_tier.m_TOW;

alpha_trim = Storing_STABILITY_DATA_1.TRIM_RESULTS.trim_alpha;

% Constants
g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
m2ft = conv_UNITS.m2ft;

% q_inf = Performance.q_inf;
q_inf = 0.5*rho*V^2;

T = Stab_Der.T;
CD = Stab_Der.CD;
D = q_inf*S_w1*CD;
gamma = asin((T-D)/(m_TOW*g));
gamma_deg = gamma*R2D;

% Conditions
% Bank angle zero
phi = conditions.conditions_TRIM_lat.phi_0; % defines the directión of bank angle (positive wing right down)
nhat = conditions.conditions_TRIM_lat.nhat; % nhat = 0; % For steady State & % nhat = 1; % For level turning flight

Trim_ITER_LAT_Viraje.gamma_deg = gamma_deg;

Cyb = Stab_Der.Cyb;
Clb = Stab_Der.Clb;
Cnb = Stab_Der.Cnb;

Cyr = Stab_Der.Cyr;
Clr = Stab_Der.Clr;
Cnr = Stab_Der.Cnr;

Cydeltaa = Stab_Der.Cydeltaa;
Cldeltaa = Stab_Der.Cldeltaa;
Cndeltaa = Stab_Der.Cndeltaa;

Cydeltar = Stab_Der.Cydeltar;
Cldeltar = Stab_Der.Cldeltar;
Cndeltar = Stab_Der.Cndeltar;

% Aileron Trim Tab
Cy_delta_a_Tab = Stab_Der.Cy_deltaa_Tab;
Cl_delta_a_Tab = Stab_Der.Cl_deltaa_Tab;
Cn_delta_a_Tab = Stab_Der.Cn_deltaa_Tab;
% Rudder Trim Tab
Cy_delta_r_Tab = Stab_Der.Cydeltar_Tab;
Cl_delta_r_Tab = Stab_Der.Cldeltar_Tab;
Cn_delta_r_Tab = Stab_Der.Cndeltar_Tab;

% Hinge linee derivatives
Ch_delta_a = Stab_Der.Ch_delta_a;
Ch_delta_r = Stab_Der.Ch_delta_r;
Ch_delta_at = Stab_Der.Ch_delta_at;
Ch_delta_rt = Stab_Der.Ch_delta_rt;

% Estimated from Literature
Chda = -0.4;
Chdr = -0.4;
ChdaT = -0.08;
% Estimated from SANAID
Chda = Ch_delta_a;
Chdr = Ch_delta_r;
ChdaT = Ch_delta_at;
% 
Chbeta_v = -0.05; % Estimation


Kra = 0; % Constant that depends on the mechanical instalation of a rudder-aileron interconection spring
% Medium Aircraft (e.g., Cessna Caravan, Beechcraft King Air): 3:1 to 4:1
Ga = 0.4*m2ft; % rad/ft -> rad/m, Stick to aileron gearing ratio
Gr = 1.5*m2ft; % rad/ft -> rad/m,

% Deflection of rudder trim TAB (no existent)
drT = 0;
daT = 20*D2R;

%% The SideSlip Equation
% Cybeta*beta+Cyda*da+CydaT*daT+Cydr*dr+CydrT*drT+Cond1+FTy
% a1*beta+a2*da+a3*dr+a4*daT+a5*drT+C1+FTy
a1 = Cyb;
a2 = Cydeltaa;
a3 = Cydeltar;
a4 = Cy_delta_a_Tab;
a5 = Cy_delta_r_Tab;
C1 = (((1-nhat)*(m_TOW*g)/(q_inf*S_w1)) + (nhat*Cyr*b_w1*g/(2*V^2)))*sin(phi);
a6 = get_side_force_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,conditions,conv_UNITS);
FTy = a6;

%% The Rolling Moment Equation
% Clbeta*beta+Clda*da+CldaT*daT+Cldr*dr+CldrT*drT+LTy+Cond2
% b1*beta+b2*da+b3*dr+b4*daT+b5*drT+C2+LTy
b1 = Clb;
b2 = Cldeltaa;
b3 = Cldeltar;
b4 = Cl_delta_a_Tab;
b5 = Cl_delta_r_Tab;
% Inertias wrt stability axis
Iyys = Iyy;
Ixxs = Ixx*((cos(alpha_trim))^2) + Izz*((sin(alpha_trim))^2) - Ixz*(sin(2*alpha_trim));
Izzs = Ixx*((sin(alpha_trim))^2) + Izz*((cos(alpha_trim))^2) + Ixz*(sin(2*alpha_trim));
Ixzs = 0.5*Ixx*sin(2*alpha_trim) - 0.5*Izz*cos(2*alpha_trim) + Ixz*((cos(alpha_trim))^2);
C2 = (nhat*(Iyys - Izzs)*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Clr*b_w1*g/(2*V^2));
%% Rolling Moments definition
% CAse 1: NO rolling movement
% case 2 % Asymetric loads in both wings
% case 3 % Asymetric load in Wing 1
% case 4 % Asymetric load in Wing 2
% Case 5 % Engine torque in center 
LTy = get_rolling_moment_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,conditions);

%% The Yawing Moment Equation
% Cnbeta*beta+Cnda*da+CndaT*daT+Cndr*dr+CndrT*drT+N*Ty+Cond3
% beta*c1+c2*da+c3*dr+c4*daT+c5*drT+C3+NTy
c1 = Cnb;
c2 = Cndeltaa;
c3 = Cndeltar;
c4 = Cn_delta_a_Tab;
c5 = Cn_delta_r_Tab;
C3 = ((-nhat*Ixzs*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Cnr*b_w1*g/(2*V^2)))*sin(phi);
%% Yawing Moments definition
% CAse 1: NO Yawing movement
% CAse 2: Yawing movement drag asymetry
% CAse 3: Yawing movement engine asymetry
% change yawing moment direction
% change_dir = -1;
% Yawing_maneuver = 1;
% conditions_TRIM_lat.Yawing_maneuver = Yawing_maneuver;
% N_moment = b_w1/2; % Worse case scenarios,
NTy =  get_yawing_moment_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,Storing_AERO_DATA,conditions,conv_UNITS,case_AC,OUTPUT_read_XLSX);

%% The Aileron Stick or Wheel Force Equation
% 2*Chda*Ga*da*k4+ChdaT*Ga*daT*k4-Ga*Kra*dr+Fa
% d2*da+d3*dr+d4*daT+Fa
% Aileron Surface
chord_ail = cf_ail*cmac_w1;
b_ail = ((K_y2_TT_ail_w1 - K_y1_TT_ail_w1)*b_w1);
S_ail = chord_ail*b_ail;
k4 = q_inf*S_ail*chord_ail*Ga;
d2 = Chda*Ga*k4;
d3 = -Ga*Kra;
d4 = ChdaT*Ga*k4;
%% The Rudder Pedal Force Equation

% Calculation of trim lateral conditions
beta = -(a2*b3*c4*daT+a2*b3*c5*drT-a2*b4*c3*daT-a2*b5*c3*drT-a3*b2*c4*daT-a3*b2*c5*drT+a3*b4*c2*daT+...
    a3*b5*c2*drT+a4*b2*c3*daT-a4*b3*c2*daT+a5*b2*c3*drT-a5*b3*c2*drT+C1*b2*c3-C1*b3*c2-C2*a2*c3+C2*a3*c2+...
    C3*a2*b3-C3*a3*b2+FTy*b2*c3-FTy*b3*c2-LTy*a2*c3+LTy*a3*c2+NTy*a2*b3-NTy*a3*b2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

deltaa = (a1*b3*c4*daT+a1*b3*c5*drT-a1*b4*c3*daT-a1*b5*c3*drT-a3*b1*c4*daT-a3*b1*c5*drT+a3*b4*c1*daT+a3*b5*c1*drT+...
    a4*b1*c3*daT-a4*b3*c1*daT+a5*b1*c3*drT-a5*b3*c1*drT+C1*b1*c3-C1*b3*c1-C2*a1*c3+C2*a3*c1+C3*a1*b3-C3*a3*b1+FTy*b1*c3-FTy*b3*c1-...
    LTy*a1*c3+LTy*a3*c1+NTy*a1*b3-NTy*a3*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1); 

Fa = (-a1*b2*c3*d4*daT+a1*b3*c2*d4*daT+a2*b1*c3*d4*daT-a2*b3*c1*d4*daT-a3*b1*c2*d4*daT+a3*b2*c1*d4*daT-...
    d2*a1*b3*c4*daT+d2*a1*b4*c3*daT+d2*a3*b1*c4*daT-d2*a3*b4*c1*daT-d2*a4*b1*c3*daT+d2*a4*b3*c1*daT+d3*a1*b2*c4*daT-...
    d3*a1*b4*c2*daT-d3*a2*b1*c4*daT+d3*a2*b4*c1*daT+d3*a4*b1*c2*daT-d3*a4*b2*c1*daT-a3*b5*c1*d2*drT+a5*b1*c2*d3*drT-...
    a5*b1*c3*d2*drT-a5*b2*c1*d3*drT+a5*b3*c1*d2*drT+a1*b2*c5*d3*drT-a1*b3*c5*d2*drT-a1*b5*c2*d3*drT+a1*b5*c3*d2*drT-...
    a2*b1*c5*d3*drT+a2*b5*c1*d3*drT+a3*b1*c5*d2*drT+C1*b1*c2*d3-C1*b1*c3*d2-C1*b2*c1*d3+C1*b3*c1*d2-C2*a1*c2*d3+...
    C2*a1*c3*d2+C2*a2*c1*d3-C2*a3*c1*d2+C3*a1*b2*d3-C3*a1*b3*d2-C3*a2*b1*d3+C3*a3*b1*d2+FTy*b1*c2*d3-FTy*b1*c3*d2-...
    FTy*b2*c1*d3+FTy*b3*c1*d2-LTy*a1*c2*d3+LTy*a1*c3*d2+LTy*a2*c1*d3-LTy*a3*c1*d2+NTy*a1*b2*d3-NTy*a1*b3*d2-NTy*a2*b1*d3+...
    NTy*a3*b1*d2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

deltar = -(a1*b2*c4*daT+a1*b2*c5*drT-a1*b4*c2*daT-a1*b5*c2*drT-a2*b1*c4*daT-a2*b1*c5*drT+a2*b4*c1*daT+...
    a2*b5*c1*drT+a4*b1*c2*daT-a4*b2*c1*daT+a5*b1*c2*drT-a5*b2*c1*drT+C1*b1*c2-C1*b2*c1-C2*a1*c2+C2*a2*c1+C3*a1*b2-...
    C3*a2*b1+FTy*b1*c2-FTy*b2*c1-LTy*a1*c2+LTy*a2*c1+NTy*a1*b2-NTy*a2*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

deltaa_deg = deltaa*R2D;
deltar_deg = deltar*R2D;
deltaaT_deg = daT*R2D;
deltarT_deg = drT*R2D;
beta_deg = beta*R2D;
phi_deg = phi*R2D;

Trim_ITER_LAT.deltaa_deg = deltaa_deg;
Trim_ITER_LAT.deltar_deg = deltar_deg;
Trim_ITER_LAT.deltaaT_deg = deltaaT_deg;
Trim_ITER_LAT.deltarT_deg = deltarT_deg;
Trim_ITER_LAT.beta_deg = beta_deg;
Trim_ITER_LAT.phi_deg = phi_deg;
Trim_ITER_LAT.Fa = Fa;
Trim_ITER_LAT.Fr = Fr;

% Variable phi study
for i=1:length(deltaaT_vec)
    q_inf = 0.5*rho*(V^2);
    daT = deltaaT_vec(i);
    phi = conditions.conditions_TRIM_lat.phi_0;
    
    C1 = (((1-nhat)*(m_TOW*g)/(q_inf*S_w1)) + (nhat*Cyr*b_w1*g/(2*V^2)))*sin(phi);
    C2 = (nhat*(Iyys - Izzs)*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Clr*b_w1*g/(2*V^2));
    C3 = ((-nhat*Ixzs*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Cnr*b_w1*g/(2*V^2)))*sin(phi);

    % Solutions

    % Calculation of trim lateral conditions
    beta_var(i) = -(a2*b3*c4*daT+a2*b3*c5*drT-a2*b4*c3*daT-a2*b5*c3*drT-a3*b2*c4*daT-a3*b2*c5*drT+a3*b4*c2*daT+...
        a3*b5*c2*drT+a4*b2*c3*daT-a4*b3*c2*daT+a5*b2*c3*drT-a5*b3*c2*drT+C1*b2*c3-C1*b3*c2-C2*a2*c3+C2*a3*c2+...
        C3*a2*b3-C3*a3*b2+FTy*b2*c3-FTy*b3*c2-LTy*a2*c3+LTy*a3*c2+NTy*a2*b3-NTy*a3*b2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

    deltaa_var(i) = (a1*b3*c4*daT+a1*b3*c5*drT-a1*b4*c3*daT-a1*b5*c3*drT-a3*b1*c4*daT-a3*b1*c5*drT+a3*b4*c1*daT+a3*b5*c1*drT+...
        a4*b1*c3*daT-a4*b3*c1*daT+a5*b1*c3*drT-a5*b3*c1*drT+C1*b1*c3-C1*b3*c1-C2*a1*c3+C2*a3*c1+C3*a1*b3-C3*a3*b1+FTy*b1*c3-FTy*b3*c1-...
        LTy*a1*c3+LTy*a3*c1+NTy*a1*b3-NTy*a3*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

    Fa_var(i) = (-a1*b2*c3*d4*daT+a1*b3*c2*d4*daT+a2*b1*c3*d4*daT-a2*b3*c1*d4*daT-a3*b1*c2*d4*daT+a3*b2*c1*d4*daT-...
        d2*a1*b3*c4*daT+d2*a1*b4*c3*daT+d2*a3*b1*c4*daT-d2*a3*b4*c1*daT-d2*a4*b1*c3*daT+d2*a4*b3*c1*daT+d3*a1*b2*c4*daT-...
        d3*a1*b4*c2*daT-d3*a2*b1*c4*daT+d3*a2*b4*c1*daT+d3*a4*b1*c2*daT-d3*a4*b2*c1*daT-a3*b5*c1*d2*drT+a5*b1*c2*d3*drT-...
        a5*b1*c3*d2*drT-a5*b2*c1*d3*drT+a5*b3*c1*d2*drT+a1*b2*c5*d3*drT-a1*b3*c5*d2*drT-a1*b5*c2*d3*drT+a1*b5*c3*d2*drT-...
        a2*b1*c5*d3*drT+a2*b5*c1*d3*drT+a3*b1*c5*d2*drT+C1*b1*c2*d3-C1*b1*c3*d2-C1*b2*c1*d3+C1*b3*c1*d2-C2*a1*c2*d3+...
        C2*a1*c3*d2+C2*a2*c1*d3-C2*a3*c1*d2+C3*a1*b2*d3-C3*a1*b3*d2-C3*a2*b1*d3+C3*a3*b1*d2+FTy*b1*c2*d3-FTy*b1*c3*d2-...
        FTy*b2*c1*d3+FTy*b3*c1*d2-LTy*a1*c2*d3+LTy*a1*c3*d2+LTy*a2*c1*d3-LTy*a3*c1*d2+NTy*a1*b2*d3-NTy*a1*b3*d2-NTy*a2*b1*d3+...
        NTy*a3*b1*d2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

    deltar_var(i) = -(a1*b2*c4*daT+a1*b2*c5*drT-a1*b4*c2*daT-a1*b5*c2*drT-a2*b1*c4*daT-a2*b1*c5*drT+a2*b4*c1*daT+...
        a2*b5*c1*drT+a4*b1*c2*daT-a4*b2*c1*daT+a5*b1*c2*drT-a5*b2*c1*drT+C1*b1*c2-C1*b2*c1-C2*a1*c2+C2*a2*c1+C3*a1*b2-...
        C3*a2*b1+FTy*b1*c2-FTy*b2*c1-LTy*a1*c2+LTy*a2*c1+NTy*a1*b2-NTy*a2*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

    % Solutions in degrees
    deltaa_deg_var(i) = deltaa_var(i)*R2D;
    deltar_deg_var(i) = deltar_var(i)*R2D;
    beta_deg_var(i) = beta_var(i)*R2D;
end

Trim_ITER_LAT.deltaaT_vec = deltaaT_vec;
Trim_ITER_LAT.phi_vec = phi_vec;
Trim_ITER_LAT.deltaa_var = deltaa_var;
Trim_ITER_LAT.deltar_var = deltar_var;
Trim_ITER_LAT.Fa_var = Fa_var;
Trim_ITER_LAT.beta_var = beta_var;

Trim_ITER_LAT.deltaa_deg_var = deltaa_deg_var;
Trim_ITER_LAT.deltar_deg_var = deltar_deg_var;
Trim_ITER_LAT.beta_deg_var = beta_deg_var;

% Variable phi study
for i=1:length(phi_vec)
    q_inf = 0.5*rho*(V^2);
    phi = phi_vec(i);
    nhat = 1;
    daT = 20*D2R;
    
    C1 = (((1-nhat)*(m_TOW*g)/(q_inf*S_w1)) + (nhat*Cyr*b_w1*g/(2*V^2)))*sin(phi);
    C2 = (nhat*(Iyys - Izzs)*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Clr*b_w1*g/(2*V^2));
    C3 = ((-nhat*Ixzs*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Cnr*b_w1*g/(2*V^2)))*sin(phi);

    % Solutions

    % Calculation of trim lateral conditions
    beta_var1(i) = -(a2*b3*c4*daT+a2*b3*c5*drT-a2*b4*c3*daT-a2*b5*c3*drT-a3*b2*c4*daT-a3*b2*c5*drT+a3*b4*c2*daT+...
        a3*b5*c2*drT+a4*b2*c3*daT-a4*b3*c2*daT+a5*b2*c3*drT-a5*b3*c2*drT+C1*b2*c3-C1*b3*c2-C2*a2*c3+C2*a3*c2+...
        C3*a2*b3-C3*a3*b2+FTy*b2*c3-FTy*b3*c2-LTy*a2*c3+LTy*a3*c2+NTy*a2*b3-NTy*a3*b2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

    deltaa_var1(i) = (a1*b3*c4*daT+a1*b3*c5*drT-a1*b4*c3*daT-a1*b5*c3*drT-a3*b1*c4*daT-a3*b1*c5*drT+a3*b4*c1*daT+a3*b5*c1*drT+...
        a4*b1*c3*daT-a4*b3*c1*daT+a5*b1*c3*drT-a5*b3*c1*drT+C1*b1*c3-C1*b3*c1-C2*a1*c3+C2*a3*c1+C3*a1*b3-C3*a3*b1+FTy*b1*c3-FTy*b3*c1-...
        LTy*a1*c3+LTy*a3*c1+NTy*a1*b3-NTy*a3*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

    Fa_var1(i) = (-a1*b2*c3*d4*daT+a1*b3*c2*d4*daT+a2*b1*c3*d4*daT-a2*b3*c1*d4*daT-a3*b1*c2*d4*daT+a3*b2*c1*d4*daT-...
        d2*a1*b3*c4*daT+d2*a1*b4*c3*daT+d2*a3*b1*c4*daT-d2*a3*b4*c1*daT-d2*a4*b1*c3*daT+d2*a4*b3*c1*daT+d3*a1*b2*c4*daT-...
        d3*a1*b4*c2*daT-d3*a2*b1*c4*daT+d3*a2*b4*c1*daT+d3*a4*b1*c2*daT-d3*a4*b2*c1*daT-a3*b5*c1*d2*drT+a5*b1*c2*d3*drT-...
        a5*b1*c3*d2*drT-a5*b2*c1*d3*drT+a5*b3*c1*d2*drT+a1*b2*c5*d3*drT-a1*b3*c5*d2*drT-a1*b5*c2*d3*drT+a1*b5*c3*d2*drT-...
        a2*b1*c5*d3*drT+a2*b5*c1*d3*drT+a3*b1*c5*d2*drT+C1*b1*c2*d3-C1*b1*c3*d2-C1*b2*c1*d3+C1*b3*c1*d2-C2*a1*c2*d3+...
        C2*a1*c3*d2+C2*a2*c1*d3-C2*a3*c1*d2+C3*a1*b2*d3-C3*a1*b3*d2-C3*a2*b1*d3+C3*a3*b1*d2+FTy*b1*c2*d3-FTy*b1*c3*d2-...
        FTy*b2*c1*d3+FTy*b3*c1*d2-LTy*a1*c2*d3+LTy*a1*c3*d2+LTy*a2*c1*d3-LTy*a3*c1*d2+NTy*a1*b2*d3-NTy*a1*b3*d2-NTy*a2*b1*d3+...
        NTy*a3*b1*d2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

    deltar_var1(i) = -(a1*b2*c4*daT+a1*b2*c5*drT-a1*b4*c2*daT-a1*b5*c2*drT-a2*b1*c4*daT-a2*b1*c5*drT+a2*b4*c1*daT+...
        a2*b5*c1*drT+a4*b1*c2*daT-a4*b2*c1*daT+a5*b1*c2*drT-a5*b2*c1*drT+C1*b1*c2-C1*b2*c1-C2*a1*c2+C2*a2*c1+C3*a1*b2-...
        C3*a2*b1+FTy*b1*c2-FTy*b2*c1-LTy*a1*c2+LTy*a2*c1+NTy*a1*b2-NTy*a2*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

    % Solutions in degrees
    deltaa_deg_var1(i) = deltaa_var1(i)*R2D;
    deltar_deg_var1(i) = deltar_var1(i)*R2D;
    beta_deg_var1(i) = beta_var1(i)*R2D;
end

Trim_ITER_LAT.phi_vec = phi_vec;
Trim_ITER_LAT.deltaa_var1 = deltaa_var1;
Trim_ITER_LAT.deltar_var1 = deltar_var1;
Trim_ITER_LAT.Fa_var1 = Fa_var1;
Trim_ITER_LAT.beta_var1 = beta_var1;

Trim_ITER_LAT.deltaa_deg_var1 = deltaa_deg_var1;
Trim_ITER_LAT.deltar_deg_var1 = deltar_deg_var1;
Trim_ITER_LAT.beta_deg_var1 = beta_deg_var1;

%% Variable Phi study and var V
V_VAR = conditions_TRIM_lat.V_VAR;
for i=1:length(phi_vec)
    for j=1:length(V_VAR)
        % ARE A FUNCTION OF SPEED
        V = V_VAR(j);
        q_inf = 0.5*rho*(V^2);
        phi = phi_vec(i);
        nhat = 1;
        daT = 20*D2R; 

        C1 = (((1-nhat)*(m_TOW*g)/(q_inf*S_w1)) + (nhat*Cyr*b_w1*g/(2*V^2)))*sin(phi);
        a6 = get_side_force_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,conditions,conv_UNITS);
        FTy = a6;
        C2 = (nhat*(Iyys - Izzs)*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Clr*b_w1*g/(2*V^2));
        LTy = get_rolling_moment_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,conditions);
        C3 = ((-nhat*Ixzs*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Cnr*b_w1*g/(2*V^2)))*sin(phi);
        NTy =  get_yawing_moment_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,Storing_AERO_DATA,conditions,conv_UNITS,case_AC,OUTPUT_read_XLSX);
        k4 = q_inf*S_ail*chord_ail*Ga;
        
        % SOLUTIONS
        % Calculation of trim lateral conditions
        beta_var2(i,j) = -(a2*b3*c4*daT+a2*b3*c5*drT-a2*b4*c3*daT-a2*b5*c3*drT-a3*b2*c4*daT-a3*b2*c5*drT+a3*b4*c2*daT+...
            a3*b5*c2*drT+a4*b2*c3*daT-a4*b3*c2*daT+a5*b2*c3*drT-a5*b3*c2*drT+C1*b2*c3-C1*b3*c2-C2*a2*c3+C2*a3*c2+...
            C3*a2*b3-C3*a3*b2+FTy*b2*c3-FTy*b3*c2-LTy*a2*c3+LTy*a3*c2+NTy*a2*b3-NTy*a3*b2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

        deltaa_var2(i,j) = (a1*b3*c4*daT+a1*b3*c5*drT-a1*b4*c3*daT-a1*b5*c3*drT-a3*b1*c4*daT-a3*b1*c5*drT+a3*b4*c1*daT+a3*b5*c1*drT+...
            a4*b1*c3*daT-a4*b3*c1*daT+a5*b1*c3*drT-a5*b3*c1*drT+C1*b1*c3-C1*b3*c1-C2*a1*c3+C2*a3*c1+C3*a1*b3-C3*a3*b1+FTy*b1*c3-FTy*b3*c1-...
            LTy*a1*c3+LTy*a3*c1+NTy*a1*b3-NTy*a3*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

        Fa_var2(i,j) = (-a1*b2*c3*d4*daT+a1*b3*c2*d4*daT+a2*b1*c3*d4*daT-a2*b3*c1*d4*daT-a3*b1*c2*d4*daT+a3*b2*c1*d4*daT-...
            d2*a1*b3*c4*daT+d2*a1*b4*c3*daT+d2*a3*b1*c4*daT-d2*a3*b4*c1*daT-d2*a4*b1*c3*daT+d2*a4*b3*c1*daT+d3*a1*b2*c4*daT-...
            d3*a1*b4*c2*daT-d3*a2*b1*c4*daT+d3*a2*b4*c1*daT+d3*a4*b1*c2*daT-d3*a4*b2*c1*daT-a3*b5*c1*d2*drT+a5*b1*c2*d3*drT-...
            a5*b1*c3*d2*drT-a5*b2*c1*d3*drT+a5*b3*c1*d2*drT+a1*b2*c5*d3*drT-a1*b3*c5*d2*drT-a1*b5*c2*d3*drT+a1*b5*c3*d2*drT-...
            a2*b1*c5*d3*drT+a2*b5*c1*d3*drT+a3*b1*c5*d2*drT+C1*b1*c2*d3-C1*b1*c3*d2-C1*b2*c1*d3+C1*b3*c1*d2-C2*a1*c2*d3+...
            C2*a1*c3*d2+C2*a2*c1*d3-C2*a3*c1*d2+C3*a1*b2*d3-C3*a1*b3*d2-C3*a2*b1*d3+C3*a3*b1*d2+FTy*b1*c2*d3-FTy*b1*c3*d2-...
            FTy*b2*c1*d3+FTy*b3*c1*d2-LTy*a1*c2*d3+LTy*a1*c3*d2+LTy*a2*c1*d3-LTy*a3*c1*d2+NTy*a1*b2*d3-NTy*a1*b3*d2-NTy*a2*b1*d3+...
            NTy*a3*b1*d2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

        deltar_var2(i,j) = -(a1*b2*c4*daT+a1*b2*c5*drT-a1*b4*c2*daT-a1*b5*c2*drT-a2*b1*c4*daT-a2*b1*c5*drT+a2*b4*c1*daT+...
            a2*b5*c1*drT+a4*b1*c2*daT-a4*b2*c1*daT+a5*b1*c2*drT-a5*b2*c1*drT+C1*b1*c2-C1*b2*c1-C2*a1*c2+C2*a2*c1+C3*a1*b2-...
            C3*a2*b1+FTy*b1*c2-FTy*b2*c1-LTy*a1*c2+LTy*a2*c1+NTy*a1*b2-NTy*a2*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

        % SOLUTIONS IN DEGREES
        deltaa_deg_var2(i,j) = deltaa_var2(i,j)*R2D;
        deltar_deg_var2(i,j) = deltar_var2(i,j)*R2D;
        beta_deg_var2(i,j) = beta_var2(i,j)*R2D;
    end
end
Trim_ITER_LAT.V_VAR = V_VAR;
Trim_ITER_LAT.deltaa_var2 = deltaa_var2;
Trim_ITER_LAT.deltar_var2 = deltar_var2;
Trim_ITER_LAT.Fa_var2 = Fa_var2;
Trim_ITER_LAT.beta_var2 = beta_var2;

Trim_ITER_LAT.deltaa_deg_var2 = deltaa_deg_var2;
Trim_ITER_LAT.deltar_deg_var2 = deltar_deg_var2;
Trim_ITER_LAT.beta_deg_var2 = beta_deg_var2;

%% Variable deltaaT study and var V
for i=1:length(deltaaT_vec)
    for j=1:length(V_VAR)
        % ARE A FUNCTION OF SPEED
        V = V_VAR(j);
        q_inf = 0.5*rho*(V^2);
        daT = deltaaT_vec(i);
        phi = conditions.conditions_TRIM_lat.phi_0;

        C1 = (((1-nhat)*(m_TOW*g)/(q_inf*S_w1)) + (nhat*Cyr*b_w1*g/(2*V^2)))*sin(phi);
        a6 = get_side_force_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,conditions,conv_UNITS);
        FTy = a6;
        C2 = (nhat*(Iyys - Izzs)*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Clr*b_w1*g/(2*V^2));
        LTy = get_rolling_moment_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,conditions);
        C3 = ((-nhat*Ixzs*(g^2)*tan(phi)*sin(phi))/(q_inf*S_w1*b_w1*V^2) + (nhat*Cnr*b_w1*g/(2*V^2)))*sin(phi);
        NTy =  get_yawing_moment_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,Storing_AERO_DATA,conditions,conv_UNITS,case_AC,OUTPUT_read_XLSX);
        k4 = q_inf*S_ail*chord_ail*Ga;
        
        % SOLUTIONS
        % Calculation of trim lateral conditions
        beta_var3(i,j) = -(a2*b3*c4*daT+a2*b3*c5*drT-a2*b4*c3*daT-a2*b5*c3*drT-a3*b2*c4*daT-a3*b2*c5*drT+a3*b4*c2*daT+...
            a3*b5*c2*drT+a4*b2*c3*daT-a4*b3*c2*daT+a5*b2*c3*drT-a5*b3*c2*drT+C1*b2*c3-C1*b3*c2-C2*a2*c3+C2*a3*c2+...
            C3*a2*b3-C3*a3*b2+FTy*b2*c3-FTy*b3*c2-LTy*a2*c3+LTy*a3*c2+NTy*a2*b3-NTy*a3*b2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

        deltaa_var3(i,j) = (a1*b3*c4*daT+a1*b3*c5*drT-a1*b4*c3*daT-a1*b5*c3*drT-a3*b1*c4*daT-a3*b1*c5*drT+a3*b4*c1*daT+a3*b5*c1*drT+...
            a4*b1*c3*daT-a4*b3*c1*daT+a5*b1*c3*drT-a5*b3*c1*drT+C1*b1*c3-C1*b3*c1-C2*a1*c3+C2*a3*c1+C3*a1*b3-C3*a3*b1+FTy*b1*c3-FTy*b3*c1-...
            LTy*a1*c3+LTy*a3*c1+NTy*a1*b3-NTy*a3*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

        Fa_var3(i,j) = (-a1*b2*c3*d4*daT+a1*b3*c2*d4*daT+a2*b1*c3*d4*daT-a2*b3*c1*d4*daT-a3*b1*c2*d4*daT+a3*b2*c1*d4*daT-...
            d2*a1*b3*c4*daT+d2*a1*b4*c3*daT+d2*a3*b1*c4*daT-d2*a3*b4*c1*daT-d2*a4*b1*c3*daT+d2*a4*b3*c1*daT+d3*a1*b2*c4*daT-...
            d3*a1*b4*c2*daT-d3*a2*b1*c4*daT+d3*a2*b4*c1*daT+d3*a4*b1*c2*daT-d3*a4*b2*c1*daT-a3*b5*c1*d2*drT+a5*b1*c2*d3*drT-...
            a5*b1*c3*d2*drT-a5*b2*c1*d3*drT+a5*b3*c1*d2*drT+a1*b2*c5*d3*drT-a1*b3*c5*d2*drT-a1*b5*c2*d3*drT+a1*b5*c3*d2*drT-...
            a2*b1*c5*d3*drT+a2*b5*c1*d3*drT+a3*b1*c5*d2*drT+C1*b1*c2*d3-C1*b1*c3*d2-C1*b2*c1*d3+C1*b3*c1*d2-C2*a1*c2*d3+...
            C2*a1*c3*d2+C2*a2*c1*d3-C2*a3*c1*d2+C3*a1*b2*d3-C3*a1*b3*d2-C3*a2*b1*d3+C3*a3*b1*d2+FTy*b1*c2*d3-FTy*b1*c3*d2-...
            FTy*b2*c1*d3+FTy*b3*c1*d2-LTy*a1*c2*d3+LTy*a1*c3*d2+LTy*a2*c1*d3-LTy*a3*c1*d2+NTy*a1*b2*d3-NTy*a1*b3*d2-NTy*a2*b1*d3+...
            NTy*a3*b1*d2)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

        deltar_var3(i,j) = -(a1*b2*c4*daT+a1*b2*c5*drT-a1*b4*c2*daT-a1*b5*c2*drT-a2*b1*c4*daT-a2*b1*c5*drT+a2*b4*c1*daT+...
            a2*b5*c1*drT+a4*b1*c2*daT-a4*b2*c1*daT+a5*b1*c2*drT-a5*b2*c1*drT+C1*b1*c2-C1*b2*c1-C2*a1*c2+C2*a2*c1+C3*a1*b2-...
            C3*a2*b1+FTy*b1*c2-FTy*b2*c1-LTy*a1*c2+LTy*a2*c1+NTy*a1*b2-NTy*a2*b1)/(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);

        % SOLUTIONS IN DEGREES
        deltaa_deg_var3(i,j) = deltaa_var3(i,j)*R2D;
        deltar_deg_var3(i,j) = deltar_var3(i,j)*R2D;
        beta_deg_var3(i,j) = beta_var3(i,j)*R2D;
    end
end
Trim_ITER_LAT.V_VAR = V_VAR;
Trim_ITER_LAT.deltaa_var3 = deltaa_var3;
Trim_ITER_LAT.deltar_var3 = deltar_var3;
Trim_ITER_LAT.Fa_var3 = Fa_var3;
Trim_ITER_LAT.beta_var3 = beta_var3;

Trim_ITER_LAT.deltaa_deg_var3 = deltaa_deg_var3;
Trim_ITER_LAT.deltar_deg_var3 = deltar_deg_var3;
Trim_ITER_LAT.beta_deg_var3 = beta_deg_var3;

