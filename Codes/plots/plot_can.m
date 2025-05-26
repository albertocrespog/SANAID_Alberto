function MESH_AC = plot_can(PLOTTING_UAV,color_can,MESH_AC)

x_mesh_can_New = PLOTTING_UAV.x_mesh_can_New; % HTP
y_mesh_can_New = PLOTTING_UAV.y_mesh_can_New; % HTP
z_mesh_can_New = PLOTTING_UAV.z_mesh_can_New; % HTP
% Color
C_can(:,:,3) = color_can(1)*ones(size(x_mesh_can_New));
C_can(:,:,2) = color_can(2)*ones(size(y_mesh_can_New));
C_can(:,:,3) = color_can(3)*ones(size(z_mesh_can_New));
mesh(x_mesh_can_New,y_mesh_can_New,z_mesh_can_New,C_can)

MESH_AC.x_mesh_can_New = x_mesh_can_New;
MESH_AC.y_mesh_can_New = y_mesh_can_New;
MESH_AC.z_mesh_can_New = z_mesh_can_New;