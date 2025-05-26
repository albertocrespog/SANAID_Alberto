function a6 = get_side_force_trim_lateral(Geo_tier,rho,V,conditions_TRIM_lat,conditions,conv_UNITS)

m_TOW = conditions_TRIM_lat.m_TOW;

% Performance
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

% % Rolling_maneuver = conditions_TRIM_lat.Rolling_maneuver;
% Rolling_maneuver = conditions.conditions_TRIM_lat.Rolling_maneuver;
% 
% change_dir_roll = conditions.conditions_TRIM_lat.change_dir_roll; 
% 
% sign1 = conditions.conditions_TRIM_lat.sign1;
% L_arm1 = conditions.conditions_TRIM_lat.L_arm1;
% L_weight1 = conditions.conditions_TRIM_lat.L_weight1;
% 
% sign2 = conditions.conditions_TRIM_lat.sign2;
% L_arm2 = conditions.conditions_TRIM_lat.L_arm2;
% L_weight2 = conditions.conditions_TRIM_lat.L_weight2;
% 
% 
% % Switching CASES
% switch Rolling_maneuver
%     case 1
%         LT = 0; % Rolling moment
%     case 2 % Asymetric loads in both wings
%         LT_1 = change_dir_roll*sign1*L_arm1*L_weight1; % Rolling moment BAT
%         LT_2 = change_dir_roll*sign2*L_arm2*L_weight2; % Rolling moment BAT
%         LT = LT_1 + LT_2; % Rolling moment
%     case 3 % Asymetric load in Wing 1
%         LT_1 = change_dir_roll*sign1*L_arm1*L_weight1; % Rolling moment BAT
%         LT = LT_1; % Rolling moment
%     case 4 % Asymetric load in Wing 2
%         LT_2 = change_dir_roll*sign2*L_arm2*L_weight2; % Rolling moment BAT
%         LT = LT_2; % Rolling moment  
%     case 5 % Rolling Moment Engine
%         LT_engine = conditions.conditions_TRIM_lat.LT_engine;
%         LT = LT_engine; % Rolling moment  
% end
FTy = 0;
% Geometric data
a6 = -FTy/(q_inf*S_w1);
