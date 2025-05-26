function conditions = initialization_STABILITY_Analysis(conditions,conv_UNITS)

m2ft = conv_UNITS.m2ft;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;

%% Initialize Stability Studies
% --------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALERT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All values here presented need to be defined at lest with this
% values which iomply that no studies are conducted, it is
% recommeded copy and past this at the beginning of the file and
% only modify in the sections if necessary to study variation of
% stability studies

%% Stability Studies
%% Asymmetric weight in each wing
% right wing
conditions.conditions_TRIM_lat.sign1 = 1; % Positive for right wing
conditions.conditions_TRIM_lat.L_arm1 = 0; % y location of mass in the wingspan
conditions.conditions_TRIM_lat.Z_arm1 = 0; % z location of mass in the wingspan
conditions.conditions_TRIM_lat.L_weight1 = 0; % Weight of mass in the wingspan
% left wing
conditions.conditions_TRIM_lat.sign2 = -1; % Negative for left wing
conditions.conditions_TRIM_lat.L_arm2 = 0; % y location of mass in the wingspan - opposite sing that right wing
conditions.conditions_TRIM_lat.Z_arm2 = 0; % z location of mass in the wingspan
conditions.conditions_TRIM_lat.L_weight2 = 0; % Weight of mass in the wingspan

%% Rolling Moments definition
% CAse 1: NO rolling movement external
% CAse 2: Rolling moment asymetric weight 1 different in each wing
% CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
% CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
conditions.conditions_TRIM_lat.Rolling_maneuver = 1;
conditions.conditions_TRIM_lat.change_dir_roll = 1;

%% Yawing Moments definition
% CAse 1: NO Yawing moment external
% CAse 2: Yawing movement drag asymetry different in each wing
% CAse 3: Yawing movement drag asymetry  - only in wing 1
% CAse 4: Yawing movement drag asymetry  - only in wing 2
% CAse 5: Yawing movement engine asymetry
% CAse 6: Yawing movement engine asymetry and drag asymmetry
% change yawing moment direction
conditions.conditions_TRIM_lat.Yawing_maneuver = 1;
conditions.conditions_TRIM_lat.change_dir_yaw = 1;
conditions.conditions_TRIM_lat.direction_beta = 1; % defines the directión of sideslip
% Additional Yawing moment
N_moment = 0; % Worse case scenarios % Actualizado
conditions.conditions_TRIM_lat.N_moment = N_moment;

%% Condition of initial Trim Lateral Directional
conditions.conditions_TRIM_lat.phi_0 = 0;
% defines the directión of bank angle (positive wing right down) For steady State & nhat = 1; For level turning flight
nhat = 0;
conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0;

%% From Lateral Trim conditions with Trim Tabs
% Assumptions for solving Trim-Lateral Directional equations with
% aileron trim-tab (daT) and rudder trim tab (daR) given by a
% constant value
conditions.conditions_TRIM_lat.daT = 0*D2R; % Fixed value
conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

% Assumptions for solving Trim-Lateral Directional equations with
% given Forces on aileron stick (Fa) and Rudder Peddals (in Newtons)
conditions.conditions_TRIM_lat.Fa = 0; % Fixed value
conditions.conditions_TRIM_lat.Fr = 0; % Fixed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALERT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------------------