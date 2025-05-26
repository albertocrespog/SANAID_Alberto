function PLOTTING_UAV = generate_1nacelle1_mesh(OUTPUT_read_XLSX,PLOTTING_UAV,X_LOC_NAC,PLOT_NAC)

% Generates de CYLINDER for the Left and right ENGINES
% Reference point from the front point
% Nacelle
% X_LOC_NAC = - H_ENG/2;
Y_LOC_NAC = 0;
Z_LOC_NAC = 0;

PLOTTING_UAV.Nac_L_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
PLOTTING_UAV.Nac_L_y_w1 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

PLOTTING_UAV.Nac_L_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
PLOTTING_UAV.Nac_L_y_w1_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
