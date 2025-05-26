function PLOTTING_UAV = generate_2nacelle1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_NAC)

% Nacelle
X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar1;
Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar1;
Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar1;

PLOTTING_UAV.Nac_L_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
PLOTTING_UAV.Nac_L_y_w1 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

PLOTTING_UAV.Nac_L_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
PLOTTING_UAV.Nac_L_y_w1_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;

PLOTTING_UAV.Nac_R_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
PLOTTING_UAV.Nac_R_y_w1 = PLOT_NAC.PatchData2_Y + Y_LOC_NAC;
PLOTTING_UAV.Nac_R_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

PLOTTING_UAV.Nac_R_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
PLOTTING_UAV.Nac_R_y_w1_s = PLOT_NAC.PatchData1_Y + Y_LOC_NAC;
PLOTTING_UAV.Nac_R_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
