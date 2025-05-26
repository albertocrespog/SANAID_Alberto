function MESH_AC = plot_w2(PLOTTING_UAV,color_w2,MESH_AC)

x_mesh_w2_New = PLOTTING_UAV.x_mesh_w2_New; % HTP
y_mesh_w2_New = PLOTTING_UAV.y_mesh_w2_New; % HTP
z_mesh_w2_New = PLOTTING_UAV.z_mesh_w2_New; % HTP
% Color
C_w2(:,:,3) = color_w2(1)*ones(size(x_mesh_w2_New));
C_w2(:,:,2) = color_w2(2)*ones(size(y_mesh_w2_New));
C_w2(:,:,3) = color_w2(3)*ones(size(z_mesh_w2_New));
mesh(x_mesh_w2_New,y_mesh_w2_New,z_mesh_w2_New,C_w2)

MESH_AC.x_mesh_w2_New = x_mesh_w2_New;
MESH_AC.y_mesh_w2_New = y_mesh_w2_New;
MESH_AC.z_mesh_w2_New = z_mesh_w2_New;