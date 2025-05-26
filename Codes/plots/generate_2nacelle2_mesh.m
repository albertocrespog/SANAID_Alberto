function PLOTTING_UAV = generate_2nacelle2_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,PLOT_NAC)

% Nacelle
X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar2;
Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar2;
Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar2;

PLOTTING_UAV.Nac_L_x_w2 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
PLOTTING_UAV.Nac_L_y_w2 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
PLOTTING_UAV.Nac_L_z_w2 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

PLOTTING_UAV.Nac_L_x_w2_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
PLOTTING_UAV.Nac_L_y_w2_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
PLOTTING_UAV.Nac_L_z_w2_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;

PLOTTING_UAV.Nac_R_x_w2 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
PLOTTING_UAV.Nac_R_y_w2 = PLOT_NAC.PatchData2_Y + Y_LOC_NAC;
PLOTTING_UAV.Nac_R_z_w2 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

PLOTTING_UAV.Nac_R_x_w2_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
PLOTTING_UAV.Nac_R_y_w2_s = PLOT_NAC.PatchData1_Y + Y_LOC_NAC;
PLOTTING_UAV.Nac_R_z_w2_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
