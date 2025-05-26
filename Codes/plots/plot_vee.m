function MESH_AC = plot_vee(PLOTTING_UAV,color_vee,MESH_AC)

x_mesh_vee_New = PLOTTING_UAV.x_mesh_vee_New; % HTP
y_mesh_vee_New = PLOTTING_UAV.y_mesh_vee_New; % HTP
z_mesh_vee_New = PLOTTING_UAV.z_mesh_vee_New; % HTP
% Color
C_vee(:,:,3) = color_vee(1)*ones(size(x_mesh_vee_New));
C_vee(:,:,2) = color_vee(2)*ones(size(y_mesh_vee_New));
C_vee(:,:,3) = color_vee(3)*ones(size(z_mesh_vee_New));
mesh(x_mesh_vee_New,y_mesh_vee_New,z_mesh_vee_New,C_vee)

MESH_AC.x_mesh_vee_New = x_mesh_vee_New;
MESH_AC.y_mesh_vee_New = y_mesh_vee_New;
MESH_AC.z_mesh_vee_New = z_mesh_vee_New;