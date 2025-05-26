function plot_2prop2(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)

% Left Propeller Disk
Prop_L_x_w2 = PLOTTING_UAV.Prop_L_x_w2;
Prop_L_y_w2 = PLOTTING_UAV.Prop_L_y_w2;
Prop_L_z_w2 = PLOTTING_UAV.Prop_L_z_w2;
Prop_L_x_w2_s = PLOTTING_UAV.Prop_L_x_w2_s;
Prop_L_y_w2_s = PLOTTING_UAV.Prop_L_y_w2_s;
Prop_L_z_w2_s = PLOTTING_UAV.Prop_L_z_w2_s;

% Right Propeller Disk
Prop_R_x_w2 = PLOTTING_UAV.Prop_R_x_w2;
Prop_R_y_w2 = PLOTTING_UAV.Prop_R_y_w2;
Prop_R_z_w2 = PLOTTING_UAV.Prop_R_z_w2;
Prop_R_x_w2_s = PLOTTING_UAV.Prop_R_x_w2_s;
Prop_R_y_w2_s = PLOTTING_UAV.Prop_R_y_w2_s;
Prop_R_z_w2_s = PLOTTING_UAV.Prop_R_z_w2_s;

% Plots disck actuators wing 2
pL_w2 = patch(Prop_L_x_w2,Prop_L_y_w2,Prop_L_z_w2,'y');
set(pL_w2,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_w2,'EraseMode','normal');
pL_w2_s = patch(Prop_L_x_w2_s,Prop_L_y_w2_s,Prop_L_z_w2_s,'y');
set(pL_w2_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_w2_s,'EraseMode','normal');

pR_w2 = patch(Prop_R_x_w2,Prop_R_y_w2,Prop_R_z_w2,'y');
set(pR_w2,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
% set(pR_w2,'EraseMode','normal');
pR_w2_s = patch(Prop_R_x_w2_s,Prop_R_y_w2_s,Prop_R_z_w2_s,'y');
set(pR_w2_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
% set(pR_w2_s,'EraseMode','normal');
