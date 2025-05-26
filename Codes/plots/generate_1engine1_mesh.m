function PLOTTING_UAV = generate_1engine1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,X_LOC_ENG,PLOT_ENG)

% Generates de CYLINDER for the Left and right ENGINES
% Reference point from the front point
% X_LOC_ENG = - H_ENG/2;
Y_LOC_ENG = 0;
Z_LOC_ENG = 0;

PLOTTING_UAV.Eng_L_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
PLOTTING_UAV.Eng_L_y_w1 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;

PLOTTING_UAV.Eng_L_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
PLOTTING_UAV.Eng_L_y_w1_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;