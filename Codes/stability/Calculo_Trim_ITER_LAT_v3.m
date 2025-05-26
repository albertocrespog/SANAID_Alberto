function [Trim_ITER_LAT] = Calculo_Trim_ITER_LAT_v3(conv_UNITS,conditions_TRIM_lat,...
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

% Constants
g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
% q_inf = Performance.q_inf;
q_inf = 0.5*rho*V^2;

T = Stab_Der.T;
CD = Stab_Der.CD;
D = q_inf*S_w1*CD;
gamma = asin((T-D)/(m_TOW*g));
gamma_deg = gamma*R2D;

Trim_ITER_LAT_Viraje.gamma_deg = gamma_deg;

Cyb = Stab_Der.Cyb;
Clb = Stab_Der.Clb;
Cnb = Stab_Der.Cnb;

Cydeltaa = Stab_Der.Cydeltaa;
Cldeltaa = Stab_Der.Cldeltaa;
Cndeltaa = Stab_Der.Cndeltaa;

Cydeltar = Stab_Der.Cydeltar;
Cldeltar = Stab_Der.Cldeltar;
Cndeltar = Stab_Der.Cndeltar;

a1 = Cyb;
a2 = Cydeltaa;
a3 = Cydeltar;
a4 = -m_TOW*g*cos(gamma)/(q_inf*S_w1);

b1 = Clb;
b2 = Cldeltaa;
b3 = Cldeltar;
%% Rolling Moments definition
% CAse 1: NO rolling movement
% case 2 % Asymetric loads in both wings
% case 3 % Asymetric load in Wing 1
% case 4 % Asymetric load in Wing 2
% Rolling_maneuver = 2;
% conditions_TRIM_lat.Rolling_maneuver = Rolling_maneuver;
% b4 = get_rolling_moment_trim_lateral(Geo_tier,rho,V_TO,conditions_TRIM_lat);
b4 = get_rolling_moment_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,conditions);

c1 = Cnb;
c2 = Cndeltaa;
c3 = Cndeltar;

%% Yawing Moments definition
% CAse 1: NO Yawing movement
% CAse 2: Yawing movement drag asymetry
% CAse 3: Yawing movement engine asymetry
% change yawing moment direction
% change_dir = -1;
% Yawing_maneuver = 1;
% conditions_TRIM_lat.Yawing_maneuver = Yawing_maneuver;
% N_moment = b_w1/2; % Worse case scenarios,
c4 = get_yawing_moment_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,Storing_AERO_DATA,conditions,conv_UNITS,case_AC,OUTPUT_read_XLSX);
% c4 = get_yawing_moment_trim_lateral(Geo_tier,rho,V_TO,conditions_TRIM_lat,Storing_AERO_DATA,N_moment,change_dir,conditions);

% Calculation of trim lateral conditions
deltaa = -(c3*b1*beta+c4*b3-b3*c1*beta-c3*b4)/...
    (-b3*c2+b2*c3);
deltar = (b2*c4-b2*c1*beta-b4*c2+c2*b1*beta)/...
    (-b3*c2+b2*c3);
phi = asin((b2*c3*a1*beta+b2*c4*a3-b2*a3*c1*beta-b4*c2*a3-c2*b3*a1*beta+b4*a2*c3+...
    a2*b3*c1*beta+a3*c2*b1*beta-c3*a2*b1*beta-b3*a2*c4)/(a4*(-b3*c2+b2*c3)));

deltaa_deg = deltaa*R2D;
deltar_deg = deltar*R2D;
phi_deg = phi*R2D;

Trim_ITER_LAT.deltaa_deg = deltaa_deg;
Trim_ITER_LAT.deltar_deg = deltar_deg;
Trim_ITER_LAT.phi_deg = phi_deg;

%% Variable Beta study
for i=1:length(beta_vec)
    q_inf = 0.5*rho*(V^2);
    a4 = -m_TOW*g*cos(gamma)/(q_inf*S_w1);
    beta = beta_vec(i);
    % Solutions
    deltaa_var(i) = -(c3*b1*beta+c4*b3-b3*c1*beta-c3*b4)/...
        (-b3*c2+b2*c3);
    deltar_var(i) = (b2*c4-b2*c1*beta-b4*c2+c2*b1*beta)/...
        (-b3*c2+b2*c3);
    phi_var(i) = asin((b2*c3*a1*beta+b2*c4*a3-b2*a3*c1*beta-b4*c2*a3-c2*b3*a1*beta+b4*a2*c3+...
        a2*b3*c1*beta+a3*c2*b1*beta-c3*a2*b1*beta-b3*a2*c4)/(a4*(-b3*c2+b2*c3)));
    % Solutions in degrees
    deltaa_deg_var(i) = deltaa_var(i)*R2D;
    deltar_deg_var(i) = deltar_var(i)*R2D;
    phi_deg_var(i) = phi_var(i)*R2D;
end
Trim_ITER_LAT.beta_vec = beta_vec;
Trim_ITER_LAT.deltaa_var = deltaa_var;
Trim_ITER_LAT.deltar_var = deltar_var;
Trim_ITER_LAT.phi_var = phi_var;
Trim_ITER_LAT.deltaa_deg_var = deltaa_deg_var;
Trim_ITER_LAT.deltar_deg_var = deltar_deg_var;
Trim_ITER_LAT.phi_deg_var = phi_deg_var;

%% Variable Beta study and var V
V_VAR = conditions_TRIM_lat.V_VAR;
for i=1:length(beta_vec)
    for j=1:length(V_VAR)
        q_inf = 0.5*rho*(V_VAR(j)^2);
        % aRE A FUNCTION OF sPEED
        b4 = get_rolling_moment_trim_lateral(Geo_tier,rho,V_VAR(j),conditions_TRIM_lat,conditions);
        c4 = get_yawing_moment_trim_lateral(Geo_tier,rho,V_VAR(j),conditions_TRIM_lat,Storing_AERO_DATA,conditions,conv_UNITS,case_AC,OUTPUT_read_XLSX);
        a4 = -m_TOW*g*cos(gamma)/(q_inf*S_w1);
        beta = beta_vec(i);
        % SOLUTIONS
        deltaa_var2(i) = -(c3*b1*beta+c4*b3-b3*c1*beta-c3*b4)/...
            (-b3*c2+b2*c3);
        deltar_var2(i) = (b2*c4-b2*c1*beta-b4*c2+c2*b1*beta)/...
            (-b3*c2+b2*c3);
        phi_var2(i) = asin((b2*c3*a1*beta+b2*c4*a3-b2*a3*c1*beta-b4*c2*a3-c2*b3*a1*beta+b4*a2*c3+...
            a2*b3*c1*beta+a3*c2*b1*beta-c3*a2*b1*beta-b3*a2*c4)/(a4*(-b3*c2+b2*c3)));
        % SOLUTIONS IN DEGREES
        deltaa_deg_var2(i) = deltaa_var2(i)*R2D;
        deltar_deg_var2(i) = deltar_var2(i)*R2D;
        phi_deg_var2(i) = phi_var2(i)*R2D;
    end
end
Trim_ITER_LAT.V_VAR = V_VAR;
Trim_ITER_LAT.deltaa_var2 = deltaa_var2;
Trim_ITER_LAT.deltar_var2 = deltar_var2;
Trim_ITER_LAT.phi_var2 = phi_var2;
Trim_ITER_LAT.deltaa_deg_var2 = deltaa_deg_var2;
Trim_ITER_LAT.deltar_deg_var2 = deltar_deg_var2;
Trim_ITER_LAT.phi_deg_var2 = phi_deg_var2;