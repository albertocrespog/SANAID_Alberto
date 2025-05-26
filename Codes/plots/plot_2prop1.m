function plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)

Prop_L_x_w1 = PLOTTING_UAV.Prop_L_x_w1;
Prop_L_y_w1 = PLOTTING_UAV.Prop_L_y_w1;
Prop_L_z_w1 = PLOTTING_UAV.Prop_L_z_w1;
Prop_L_x_w1_s = PLOTTING_UAV.Prop_L_x_w1_s;
Prop_L_y_w1_s = PLOTTING_UAV.Prop_L_y_w1_s;
Prop_L_z_w1_s = PLOTTING_UAV.Prop_L_z_w1_s;

Prop_R_x_w1 = PLOTTING_UAV.Prop_R_x_w1;
Prop_R_y_w1 = PLOTTING_UAV.Prop_R_y_w1;
Prop_R_z_w1 = PLOTTING_UAV.Prop_R_z_w1;
Prop_R_x_w1_s = PLOTTING_UAV.Prop_R_x_w1_s;
Prop_R_y_w1_s = PLOTTING_UAV.Prop_R_y_w1_s;
Prop_R_z_w1_s = PLOTTING_UAV.Prop_R_z_w1_s;

% Plots disck actuators wing 2
pL_w1 = patch(Prop_L_x_w1,Prop_L_y_w1,Prop_L_z_w1,'y');
set(pL_w1,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_w1,'EraseMode','normal');
pL_w1_s = patch(Prop_L_x_w1_s,Prop_L_y_w1_s,Prop_L_z_w1_s,'y');
set(pL_w1_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_w1_s,'EraseMode','normal');

pR_w1 = patch(Prop_R_x_w1,Prop_R_y_w1,Prop_R_z_w1,'y');
set(pR_w1,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
% set(pR_w1,'EraseMode','normal');
pR_w1_s = patch(Prop_R_x_w1_s,Prop_R_y_w1_s,Prop_R_z_w1_s,'y');
set(pR_w1_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
% set(pR_w1_s,'EraseMode','normal');
