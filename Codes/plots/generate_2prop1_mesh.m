function PLOTTING_UAV = generate_2prop1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,AC_CONFIGURATION,Prop_data)

Engine_conf = AC_CONFIGURATION.Engine_conf;

% Prop Disk
R = OUTPUT_read_XLSX.Propulsive_flags.D_prop/2; % radii of disk
H = 0.01; %heith of disk
SD = 20;
[PLOT_DISK] = make_cylinder_special(R,H,SD);
Delta_eng = 0.05;

eng_dia = OUTPUT_read_XLSX.Propulsive_flags.d_eng; %
eng_length = OUTPUT_read_XLSX.Propulsive_flags.l_eng; %

R_ENG = (eng_dia/2); % radii of disk
H_ENG = eng_length; %heith of disk

X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
Y_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
Z_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;

switch Engine_conf
    case 1 % Engine_conf = 1 - pusher prop
        % Generates de disk for the Left and right Rotors
        % Reference point from the front point
        X_LOC_PROP = X_LOC_ENG + H_ENG;
        Y_LOC_PROP = Y_LOC_ENG;
        Z_LOC_PROP = Z_LOC_ENG;
    case 2 % Engine_conf = 2 - puller prop
        % Generates de disk for the Left and right Rotors
        % Reference point from the front point
        X_LOC_PROP = X_LOC_ENG;
        Y_LOC_PROP = Y_LOC_ENG;
        Z_LOC_PROP = Z_LOC_ENG;
    case 3 % Engine_conf = 3 - turbine
        % Prop Disk
        R = OUTPUT_read_XLSX.Propulsive_flags.D_prop/2; % radii of disk
        H = 0.01; %heith of disk
        SD = 1;
        [PLOT_DISK] = make_cylinder_special(R,H,SD);
        Delta_eng = 0.05;        
        % Generates a really small disk for the Left and right Rotors
        % Reference point from the front point
        X_LOC_PROP = X_LOC_ENG;
        Y_LOC_PROP = Y_LOC_ENG;
        Z_LOC_PROP = Z_LOC_ENG;
end

PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

PLOTTING_UAV.Prop_R_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
PLOTTING_UAV.Prop_R_y_w1 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
PLOTTING_UAV.Prop_R_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

PLOTTING_UAV.Prop_R_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
PLOTTING_UAV.Prop_R_y_w1_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
PLOTTING_UAV.Prop_R_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

PLOTTING_UAV.PLOT_DISK = PLOT_DISK; % DATA Disk