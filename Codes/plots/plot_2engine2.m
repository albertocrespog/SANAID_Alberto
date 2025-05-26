function plot_2engine2(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng)

% Second Pair of  Engines
% Left Engine
Eng_L_x_w2 = PLOTTING_UAV.Eng_L_x_w2;
Eng_L_y_w2 = PLOTTING_UAV.Eng_L_y_w2;
Eng_L_z_w2 = PLOTTING_UAV.Eng_L_z_w2;
Eng_L_x_w2_s = PLOTTING_UAV.Eng_L_x_w2_s;
Eng_L_y_w2_s = PLOTTING_UAV.Eng_L_y_w2_s;
Eng_L_z_w2_s = PLOTTING_UAV.Eng_L_z_w2_s;

% Right Engine
Eng_R_x_w2 = PLOTTING_UAV.Eng_R_x_w2;
Eng_R_y_w2 = PLOTTING_UAV.Eng_R_y_w2;
Eng_R_z_w2 = PLOTTING_UAV.Eng_R_z_w2;
Eng_R_x_w2_s = PLOTTING_UAV.Eng_R_x_w2_s;
Eng_R_y_w2_s = PLOTTING_UAV.Eng_R_y_w2_s;
Eng_R_z_w2_s = PLOTTING_UAV.Eng_R_z_w2_s;

% Plots disck actuators Engines
pL_w2 = patch(Eng_L_x_w2,Eng_L_y_w2,Eng_L_z_w2,'y');
set(pL_w2,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_can,'EraseMode','normal');
pL_w2_s = patch(Eng_L_x_w2_s,Eng_L_y_w2_s,Eng_L_z_w2_s,'y');
set(pL_w2_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_can_s,'EraseMode','normal');

pR_w2 = patch(Eng_R_x_w2,Eng_R_y_w2,Eng_R_z_w2,'y');
set(pR_w2,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
% set(pR_can,'EraseMode','normal');
pR_w2_s = patch(Eng_R_x_w2_s,Eng_R_y_w2_s,Eng_R_z_w2_s,'y');
set(pR_w2_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
% set(pR_can_s,'EraseMode','normal');