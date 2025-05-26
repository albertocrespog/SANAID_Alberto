function [x_mesh_New, y_mesh_New, z_mesh_New] = get_DATA_New_Location(x_mesh_W,y_mesh_W,z_mesh_W,X,Y,Z,center_section)

if center_section==1
    y_mesh_New(1,:) = y_mesh_W(1,:) - Y;
    y_mesh_New(2,:) = y_mesh_W(2,:) - Y;
    y_mesh_New(3,:) = y_mesh_W(3,:) + Y;
    y_mesh_New(4,:) = y_mesh_W(4,:) + Y;
else
    y_mesh_New = y_mesh_W + Y;
end

x_mesh_New = x_mesh_W + X;
z_mesh_New = z_mesh_W + Z;

