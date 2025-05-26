function PLOTTING_UAV = generate_2engine1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_ENG)

% First Pair of engines
X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
Y_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
Z_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;

PLOTTING_UAV.Eng_L_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
PLOTTING_UAV.Eng_L_y_w1 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;

PLOTTING_UAV.Eng_L_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
PLOTTING_UAV.Eng_L_y_w1_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

PLOTTING_UAV.Eng_R_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
PLOTTING_UAV.Eng_R_y_w1 = PLOT_ENG.PatchData2_Y + Y_LOC_ENG;
PLOTTING_UAV.Eng_R_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;

PLOTTING_UAV.Eng_R_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
PLOTTING_UAV.Eng_R_y_w1_s = PLOT_ENG.PatchData1_Y + Y_LOC_ENG;
PLOTTING_UAV.Eng_R_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;
