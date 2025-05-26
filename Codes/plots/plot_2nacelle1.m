function plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac)

% Left Nacelle
Nac_L_x_w1 = PLOTTING_UAV.Nac_L_x_w1;
Nac_L_y_w1 = PLOTTING_UAV.Nac_L_y_w1;
Nac_L_z_w1 = PLOTTING_UAV.Nac_L_z_w1;
Nac_L_x_w1_s = PLOTTING_UAV.Nac_L_x_w1_s;
Nac_L_y_w1_s = PLOTTING_UAV.Nac_L_y_w1_s;
Nac_L_z_w1_s = PLOTTING_UAV.Nac_L_z_w1_s;

% Right Nacelle
Nac_R_x_w1 = PLOTTING_UAV.Nac_R_x_w1;
Nac_R_y_w1 = PLOTTING_UAV.Nac_R_y_w1;
Nac_R_z_w1 = PLOTTING_UAV.Nac_R_z_w1;
Nac_R_x_w1_s = PLOTTING_UAV.Nac_R_x_w1_s;
Nac_R_y_w1_s = PLOTTING_UAV.Nac_R_y_w1_s;
Nac_R_z_w1_s = PLOTTING_UAV.Nac_R_z_w1_s;

% Plots disck actuators Engines
pL_w1 = patch(Nac_L_x_w1,Nac_L_y_w1,Nac_L_z_w1,'y');
set(pL_w1,'FaceColor',color_nac,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_w1,'EraseMode','normal');
pL_w1_s = patch(Nac_L_x_w1_s,Nac_L_y_w1_s,Nac_L_z_w1_s,'y');
set(pL_w1_s,'FaceColor',color_nac,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_w1_s,'EraseMode','normal');

pR_w1 = patch(Nac_R_x_w1,Nac_R_y_w1,Nac_R_z_w1,'y');
set(pR_w1,'FaceColor',color_nac,'FaceLighting','phong','EdgeLighting','phong');
% set(pR_w1,'EraseMode','normal');
pR_w1_s = patch(Nac_R_x_w1_s,Nac_R_y_w1_s,Nac_R_z_w1_s,'y');
set(pR_w1_s,'FaceColor',color_nac,'FaceLighting','phong','EdgeLighting','phong');
% set(pR_w1_s,'EraseMode','normal');
