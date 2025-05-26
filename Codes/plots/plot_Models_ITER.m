function PLOTTING_UAV = plot_Models_ITER(Geo_tier,Body_Geo,meshData,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX)

% Constants
f2m = 0.3048;
D2R = pi/180;
R2D = 180/pi;

n_eng = OUTPUT_read_XLSX.Propulsive_flags.propul(2);

% Stores Aircraft Configuration
% AC_CONFIGURATION.Control_surface = Control_surface;
AC_type = AC_CONFIGURATION.AC_type;
Engine_loc = AC_CONFIGURATION.Engine_loc;
Engine_conf = AC_CONFIGURATION.Engine_conf;

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Fuselage geometry
D_fus = Geo_tier.d_fus;
W_fus = Geo_tier.w_fus;
H_fus = Geo_tier.h_fus;
L_fus = Geo_tier.l_fus;

% Dummy so that it overwrites the variables
PLOTTING_UAV.dummy = 0;

% Aircraft type
switch AC_type
    case 1 % AC_type = 1 - flying wing
        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % VTP Geometry
        [PLOTTING_UAV, x_VTP_LE] = generate_VTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_w1_LE,x_VTP_LE];
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % w2 Geometry
        [PLOTTING_UAV, x_HTP_LE] = generate_HTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % VTP Geometry
        [PLOTTING_UAV, x_VTP_LE] = generate_VTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        
        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_w1_LE x_HTP_LE x_VTP_LE];

    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Canard Geometry
        [PLOTTING_UAV,x_can_LE] = generate_can_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % HTP Geometry
        [PLOTTING_UAV, x_HTP_LE] = generate_HTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % VTP Geometry
        [PLOTTING_UAV, x_VTP_LE] = generate_VTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_can_LE x_w1_LE x_HTP_LE x_VTP_LE];

    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % vee Geometry
        [PLOTTING_UAV, x_vee_LE] = generate_vee_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_w1_LE x_vee_LE];

    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Canard Geometry
        [PLOTTING_UAV,x_can_LE] = generate_can_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % vee Geometry
        [PLOTTING_UAV, x_vee_LE] = generate_vee_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_can_LE x_w1_LE x_vee_LE];
    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Canard Geometry
        [PLOTTING_UAV,x_can_LE] = generate_can_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % VTP Geometry
        [PLOTTING_UAV, x_VTP_LE] = generate_VTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_can_LE x_w1_LE x_VTP_LE];

    case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % vee Geometry
        [PLOTTING_UAV, x_vee_LE] = generate_vee_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % VTP Geometry
        [PLOTTING_UAV, x_VTP_LE] = generate_VTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_w1_LE x_vee_LE x_VTP_LE];

    case 8 % AC_type = 8 - 2 surface: wing + V-tail+V-tail
        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % vee Geometry
        [PLOTTING_UAV, x_vee_LE] = generate_vee_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % vee2 Geometry
        [PLOTTING_UAV, x_vee2_LE] = generate_vee2_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION);
        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_w1_LE x_vee_LE x_vee2_LE];
end

%--------------------------- FUSELAJE ---------------------------------
% Identifies the locations where to determine information of the fuselaje
% to ensure that that surfaces take into account the fuselaje
model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_fus);
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;

if OUTPUT_read_XLSX.Fuselage_flags.STL_reading_test == 1
    max_x_fus = max(max(model{1}))*SF_CAD;
    min_x_fus = min(min(model{1}))*SF_CAD;
    max_y_fus = max(max(model{2}))*SF_CAD;
    min_y_fus = min(min(model{2}))*SF_CAD;
    z_max_fus = max(max(model{3}))*SF_CAD;
    z_min_fus = min(min(model{3}))*SF_CAD;

    Delta_X_stl_fus = min_x_fus;
    Delta_Y_stl_fus = -(max_y_fus - (max_y_fus - min_y_fus)/2)/2;
    Delta_Z_stl_fus = 0;

else
    [points_fus,triangles,tri_norms] = import_stl_fast_v2(OUTPUT_read_XLSX.Fuselage_flags.STL_fus,1);
    max_x_fus = max(points_fus(:,1))*SF_CAD;
    min_x_fus = min(points_fus(:,1))*SF_CAD;
    max_y_fus = max(points_fus(:,2))*SF_CAD;
    min_y_fus = min(points_fus(:,2))*SF_CAD;
    z_max_fus = max(points_fus(:,3))*SF_CAD;
    z_min_fus = min(points_fus(:,3))*SF_CAD;

    Delta_X_stl_fus = points_fus(1,1);
    Delta_Y_stl_fus = points_fus(1,2);
    Delta_Y_stl_fus = 0;
    Delta_Z_stl_fus = points_fus(1,3);
end

Delta_Y_stl = Body_Geo.Delta_Z_stl;
Delta_Z_stl = Body_Geo.Delta_Y_stl;

x_max = Body_Geo.x_max;
y_max = Body_Geo.y_max;
z_max = Body_Geo.z_max;

% Location of aerodynamic Surfaces
for k = 1:length(x_loc)
    Y_vec(k) = interp1(x_max,y_max,x_loc(k),'spline');
    Z_vec(k) = interp1(x_max,z_max,x_loc(k),'spline');
end

PLOTTING_UAV.x_loc = x_loc;
PLOTTING_UAV.Y_vec = Y_vec;
PLOTTING_UAV.Z_vec = Z_vec;

% z_max_fus = max(max(meshData{3}));
% z_min_fus = min(min(meshData{3}));
PLOTTING_UAV.z_max_fus = z_max_fus;
PLOTTING_UAV.z_min_fus = z_min_fus;

% Propulsive_flags.l_eng = l_eng; %
% Engine Cylinder
eng_dia = OUTPUT_read_XLSX.Propulsive_flags.d_eng; %
eng_length = OUTPUT_read_XLSX.Propulsive_flags.l_eng; %

R_ENG = (eng_dia/2); % radii of disk
H_ENG = eng_length; %heith of disk
SD_ENG = 20; % Side count
[PLOT_ENG] = make_cylinder_special(R_ENG,H_ENG,SD_ENG);

% Engine Cylinder
nac_dia = OUTPUT_read_XLSX.Propulsive_flags.d_nc; %
nac_length = OUTPUT_read_XLSX.Propulsive_flags.l_nc; %

R_NAC = (nac_dia/2); % radii of disk
H_NAC = nac_length; %heith of disk
SD_NAC = 20; % Side count
[PLOT_NAC] = make_cylinder_special(R_NAC,H_NAC,SD_NAC);

%--------------------------- ENGINE ---------------------------------
% Engine Configuration
% Engine location
switch Engine_loc
    case 1 % Engine_loc = 1 - under wings
        % Generates de CYLINDER for the Left and right ENGINES
        % Reference point from the front point
        if n_eng == 2
            % Engine
            PLOTTING_UAV = generate_2engine1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_ENG);
            % Nacelle
            PLOTTING_UAV = generate_2nacelle1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_NAC);
            % Propeller disk
            PLOTTING_UAV = generate_2prop1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,AC_CONFIGURATION,Prop_data);
        elseif n_eng == 4
            % Engine and Nacelle first pair
            % Engine
            PLOTTING_UAV = generate_2engine1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_ENG);
            % Nacelle
            PLOTTING_UAV = generate_2nacelle1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_NAC);
            % Propeller disk
            PLOTTING_UAV = generate_2prop1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,AC_CONFIGURATION,Prop_data);
            % Engine and Nacelle second pair
            % Engine
            PLOTTING_UAV = generate_2engine2_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_ENG);
            % Nacelle
            PLOTTING_UAV = generate_2nacelle2_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_NAC);
            % Propeller disk
            PLOTTING_UAV = generate_2prop2_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,AC_CONFIGURATION,Prop_data);
        end
    case 2 % Engine_loc = 2 - fuselage front
        if n_eng > 1
            disp('For a front fuselage engine, only one engine can be defined. Please do correct Input Data');
        end
        % Location of engine and nacelle
%         X_LOC_ENG = - H_ENG/2;
%         X_LOC_NAC = - H_ENG/2;
%         X_LOC_PROP = - H_ENG/2;
        % Location of engine and nacelle
        X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1 + H_ENG/2;
        X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1 + H_NAC/2;
        X_LOC_PROP = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
        % Engine and Nacelle for single engine fuselage front
        % Engine
        PLOTTING_UAV = generate_1engine1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,X_LOC_ENG,PLOT_ENG);
        % Nacelle
        PLOTTING_UAV = generate_1nacelle1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,X_LOC_NAC,PLOT_NAC);
        % Propeller Disk
        PLOTTING_UAV = generate_1prop1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,AC_CONFIGURATION,X_LOC_PROP,Prop_data);
    case 3 % Engine_loc = 3 - fuselage rear
        if n_eng > 1
            disp('For a rear fuselage engine, only one engine can be defined. Please do correct Input Data');
        end

        % Location of engine and nacelle
        X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1 - H_ENG;
        X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1 - H_NAC;
        X_LOC_PROP = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
        % Engine and Nacelle for single engine fuselage front
        % Engine
        PLOTTING_UAV = generate_1engine1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,X_LOC_ENG,PLOT_ENG);
        % Nacelle
        PLOTTING_UAV = generate_1nacelle1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,X_LOC_NAC,PLOT_NAC);
        % Propeller Disk
        PLOTTING_UAV = generate_1prop1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,AC_CONFIGURATION,X_LOC_PROP,Prop_data);
    case 4 % Engine_loc = 4 - wingtips
        % Engine and Nacelle first pair
        % Engine
        PLOTTING_UAV = generate_2engine1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_ENG);
        % Nacelle
        PLOTTING_UAV = generate_2nacelle1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_NAC);
        % Propeller disk
        PLOTTING_UAV = generate_2prop1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,AC_CONFIGURATION,Prop_data);
    case 5 % Engine_loc = Engine_loc = 5 - wingtips for wing and canard configuration n_eng at each side
        % Generates de CYLINDER for the Left and right ENGINES
        % Engine and Nacelle first pair
        % Engine
        PLOTTING_UAV = generate_2engine1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_ENG);
        % Nacelle
        PLOTTING_UAV = generate_2nacelle1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_NAC);
        % Propeller disk
        PLOTTING_UAV = generate_2prop1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,AC_CONFIGURATION,Prop_data);
        % Engine and Nacelle second pair
        % Engine
        PLOTTING_UAV = generate_2engine2_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_ENG);
        % Nacelle
        PLOTTING_UAV = generate_2nacelle2_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_NAC);
        % Propeller disk
        PLOTTING_UAV = generate_2prop2_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,AC_CONFIGURATION,Prop_data);
end

PLOTTING_UAV.meshData = meshData; % FUSELAJE
PLOTTING_UAV.PLOT_ENG = PLOT_ENG; % DATA ENGINE
