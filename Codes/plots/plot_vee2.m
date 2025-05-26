function MESH_AC = plot_vee2(PLOTTING_UAV,color_vee2,MESH_AC)

x_mesh_vee2_New = PLOTTING_UAV.x_mesh_vee2_New; % HTP
y_mesh_vee2_New = PLOTTING_UAV.y_mesh_vee2_New; % HTP
z_mesh_vee2_New = PLOTTING_UAV.z_mesh_vee2_New; % HTP
% Color
C_vee2(:,:,3) = color_vee2(1)*ones(size(x_mesh_vee2_New));
C_vee2(:,:,2) = color_vee2(2)*ones(size(y_mesh_vee2_New));
C_vee2(:,:,3) = color_vee2(3)*ones(size(z_mesh_vee2_New));
mesh(x_mesh_vee2_New,y_mesh_vee2_New,z_mesh_vee2_New,C_vee2)

MESH_AC.x_mesh_vee2_New = x_mesh_vee2_New;
MESH_AC.y_mesh_vee2_New = y_mesh_vee2_New;
MESH_AC.z_mesh_vee2_New = z_mesh_vee2_New;