function MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC)

x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING
% Color
C_w1(:,:,1) = color_w1(1)*ones(size(x_mesh_w1_New));
C_w1(:,:,2) = color_w1(2)*ones(size(y_mesh_w1_New));
C_w1(:,:,3) = color_w1(3)*ones(size(z_mesh_w1_New));

mesh(x_mesh_w1_New,y_mesh_w1_New,z_mesh_w1_New,C_w1)

MESH_AC.x_mesh_w1_New = x_mesh_w1_New;
MESH_AC.y_mesh_w1_New = y_mesh_w1_New;
MESH_AC.z_mesh_w1_New = z_mesh_w1_New;
