function plot_1engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng)

% Single Engine
Eng_L_x_w1 = PLOTTING_UAV.Eng_L_x_w1;
Eng_L_y_w1 = PLOTTING_UAV.Eng_L_y_w1;
Eng_L_z_w1 = PLOTTING_UAV.Eng_L_z_w1;
Eng_L_x_w1_s = PLOTTING_UAV.Eng_L_x_w1_s;
Eng_L_y_w1_s = PLOTTING_UAV.Eng_L_y_w1_s;
Eng_L_z_w1_s = PLOTTING_UAV.Eng_L_z_w1_s;

% Plots disck actuators Engines
pL_w1 = patch(Eng_L_x_w1,Eng_L_y_w1,Eng_L_z_w1,'y');
set(pL_w1,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_w1,'EraseMode','normal');
pL_w1_s = patch(Eng_L_x_w1_s,Eng_L_y_w1_s,Eng_L_z_w1_s,'y');
set(pL_w1_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
% set(pL_w1_s,'EraseMode','normal');
