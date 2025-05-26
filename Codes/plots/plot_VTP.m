function MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier)

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
twin_VTP = AC_CONFIGURATION.twin_VTP;

%--------------------------- VTP ---------------------------------
if VTP == 1
    if twin_VTP == 1
        % VTP 1
        y_offset_VTP = Geo_tier.y_offset_VTP;
        x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
        y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New;
        z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New - y_offset_VTP; % VTP 1
        z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New; % VTP 1

        % VTP 2
        x_mesh_VTP2_New = PLOTTING_UAV.x_mesh_VTP2_New; % VTP 2
        y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New; % VTP 2
        z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New + y_offset_VTP; % VTP 2
        z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New; % VTP 2
        % Assigns color
        C_vtp(:,:,3) = color_vtp(1)*ones(size(x_mesh_VTP1_New));
        C_vtp(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP1_New));
        C_vtp(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP1_New));

%         hm1 = mesh(x_mesh_VTP1_New,z_mesh_VTP1_New,y_mesh_VTP1_New,C_vtp);
%         %                 rotate(hm1, [1 0 0], 90)
%         hm2 = mesh(x_mesh_VTP2_New,z_mesh_VTP2_New,-y_mesh_VTP2_New,C_vtp);
        %                 rotate(hm2, [1 0 0], -90)
        mesh(x_mesh_VTP1_New,y_mesh_VTP1_New,z_mesh_VTP1_New,C_vtp)
        mesh(x_mesh_VTP2_New,y_mesh_VTP2_New,z_mesh_VTP2_New,C_vtp)


        MESH_AC.x_mesh_VTP1_New = x_mesh_VTP1_New;
        MESH_AC.y_mesh_VTP1_New = y_mesh_VTP1_New;
        MESH_AC.z_mesh_VTP1_New = z_mesh_VTP1_New;

        MESH_AC.x_mesh_VTP2_New = x_mesh_VTP2_New;
        MESH_AC.y_mesh_VTP2_New = y_mesh_VTP2_New;
        MESH_AC.z_mesh_VTP2_New = z_mesh_VTP2_New;
    else
        % Color
        x_mesh_VTP_New = PLOTTING_UAV.x_mesh_VTP_New; % VTP 1
        y_mesh_VTP_New = PLOTTING_UAV.y_mesh_VTP_New; % VTP 1
        z_mesh_VTP_New = PLOTTING_UAV.z_mesh_VTP_New; % VTP 1
        % Assigns color
        C_vtp(:,:,3) = color_vtp(1)*ones(size(x_mesh_VTP_New));
        C_vtp(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP_New));
        C_vtp(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP_New));
        % plots
        mesh(x_mesh_VTP_New,y_mesh_VTP_New,z_mesh_VTP_New,C_vtp)

        MESH_AC.x_mesh_VTP_New = x_mesh_VTP_New;
        MESH_AC.y_mesh_VTP_New = y_mesh_VTP_New;
        MESH_AC.z_mesh_VTP_New = z_mesh_VTP_New;
    end
end