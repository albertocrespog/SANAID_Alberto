function MESH_AC = plot_HTP(PLOTTING_UAV,color_HTP,MESH_AC)

x_mesh_HTP_New = PLOTTING_UAV.x_mesh_HTP_New; % HTP
y_mesh_HTP_New = PLOTTING_UAV.y_mesh_HTP_New; % HTP
z_mesh_HTP_New = PLOTTING_UAV.z_mesh_HTP_New; % HTP
% Color
C_HTP(:,:,3) = color_HTP(1)*ones(size(x_mesh_HTP_New));
C_HTP(:,:,2) = color_HTP(2)*ones(size(y_mesh_HTP_New));
C_HTP(:,:,3) = color_HTP(3)*ones(size(z_mesh_HTP_New));
mesh(x_mesh_HTP_New,y_mesh_HTP_New,z_mesh_HTP_New,C_HTP)

MESH_AC.x_mesh_HTP_New = x_mesh_HTP_New;
MESH_AC.y_mesh_HTP_New = y_mesh_HTP_New;
MESH_AC.z_mesh_HTP_New = z_mesh_HTP_New;