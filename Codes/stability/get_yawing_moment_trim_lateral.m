function c4 = get_yawing_moment_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,Storing_AERO_DATA,conditions,conv_UNITS,case_AC,OUTPUT_read_XLSX)

% Yawing_maneuver = conditions_TRIM_lat.Yawing_maneuver;
Yawing_maneuver = conditions.conditions_TRIM_lat.Yawing_maneuver;

% V_TO = conditions_TRIM_lat.V_TO;
m_TOW = conditions_TRIM_lat.m_TOW;

% Performance
% V = V_TO;
rho = conditions_TRIM_lat.rho;

% Geometric data
S_ref = Geo_tier.S_ref;
S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
cmac_w1 = Geo_tier.cmac_w1;

% Constants
g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
q_inf = 0.5*rho*V^2;

% change_dir = conditions.conditions_TRIM_lat.change_dir;
change_dir_yaw = conditions.conditions_TRIM_lat.change_dir_yaw; 

% Number of missiles
% n_MSL = conditions.n_MSL; %0,1,2,3,4
N_moment = conditions.conditions_TRIM_lat.N_moment;

%% Yawing Moments definition
     % CAse 1: NO Yawing moment
        % CAse 2: Yawing movement drag asymetry different in each wing
        % CAse 3: Yawing movement drag asymetry  - only in wing 1
        % CAse 4: Yawing movement drag asymetry  - only in wing 2
        % CAse 5: Yawing movement engine asymetry
        % CAse 6: Yawing movement engine asymetry and drag asymmetry

switch Yawing_maneuver
    case 1 % CAse 1: NO Yawing moment
        % Yawing Moment due to engine asymetry
        NT1 = 0; 
        % Yawing Moment due to drag asymetry
        D_wing1 = 0;
        D_wing2 = 0;
        DeltaNDrag = -D_wing1*N_moment + D_wing2*N_moment; % Yawing Moment due to drag asymetry
    case 2 % CAse 2: Yawing movement drag asymetry different in each wing
        % Yawing Moment due to engine asymetry
        NT1 = 0; 
        % Yawing Moment due to drag asymetry
        DeltaDrag = get_DeltaDrag_Asymmetry(case_AC,conditions,OUTPUT_read_XLSX,Storing_AERO_DATA,Geo_tier,rho,V);       
        D_wing1 = DeltaDrag.D_wing1;
        D_wing2 = DeltaDrag.D_wing2;
        DeltaNDrag = -D_wing1*N_moment + D_wing2*N_moment; % Yawing Moment due to drag asymetry
    case 3 % CAse 3: Yawing movement drag asymetry  - only in wing 1 (Left)
        % Yawing Moment due to engine asymetry
        NT1 = 0; 
        % Yawing Moment due to drag asymetry
        DeltaDrag = get_DeltaDrag_Asymmetry(case_AC,conditions,OUTPUT_read_XLSX,Storing_AERO_DATA,Geo_tier,rho,V);
        D_wing1 = DeltaDrag.D_wing1;
        D_wing2 = 0;
        DeltaNDrag = -D_wing1*N_moment + D_wing2*N_moment; % Yawing Moment due to drag asymetry
    case 4 % CAse 4: Yawing movement drag asymetry  - only in wing 2
        % Yawing Moment due to engine asymetry
        NT1 = 0; 
        % Yawing Moment due to drag asymetry
        DeltaDrag = get_DeltaDrag_Asymmetry(case_AC,conditions,OUTPUT_read_XLSX,Storing_AERO_DATA,Geo_tier,rho,V);
        D_wing1 = 0;
        D_wing2 = DeltaDrag.D_wing2;
        DeltaNDrag = -D_wing1*N_moment + D_wing2*N_moment; % Yawing Moment due to drag asymetry
    case 5 % CAse 5: Yawing movement engine asymetry
        % Yawing Moment due to engine asymetry
        NT1 = 0; 
        % Yawing Moment due to drag asymetry
        DeltaDrag = get_DeltaDrag_Asymmetry(case_AC,conditions,OUTPUT_read_XLSX,Storing_AERO_DATA,Geo_tier,rho,V);
        D_wing1 = 0;
        D_wing2 = 0;
        DeltaNDrag = -D_wing1*N_moment + D_wing2*N_moment; % Yawing Moment due to drag asymetry
    case 6 % CAse 6: Yawing movement engine asymetry and drag asymmetry
        % Yawing Moment due to engine asymetry
        NT1 = 0; 
        % Yawing Moment due to drag asymetry
        DeltaDrag = get_DeltaDrag_Asymmetry(case_AC,conditions,OUTPUT_read_XLSX,Storing_AERO_DATA,Geo_tier,rho,V);
        D_wing1 = DeltaDrag.D_wing1;
        D_wing2 = DeltaDrag.D_wing2;
        DeltaNDrag = -D_wing1*N_moment + D_wing2*N_moment; % Yawing Moment due to drag asymetry
end

c4 = -change_dir_yaw*(NT1 + DeltaNDrag)/(q_inf*S_w1*b_w1);
