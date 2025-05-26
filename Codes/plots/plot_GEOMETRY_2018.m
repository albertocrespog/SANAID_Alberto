% function plot_GEOMETRY_2018(Geo_tier,RESULTS,POST,N_dimension,Engine,SAVE_FIG,fig)
function Fig = plot_GEOMETRY_2022(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC,Performance,OUTPUT_read_XLSX)

% Aircraft type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
% Engine location
% Engine_loc = 1 - under wings
% Engine_loc = 2 - fuselage
% Engine_loc = 3 - wingtips
% Engine Configuration
% Engine_conf = 1 - pusher prop
% Engine_conf = 2 - puller prop
% Engine_conf = 3 - turbine

% Stores Aircraft Configuration
% AC_CONFIGURATION.Control_surface = Control_surface;
AC_type = AC_CONFIGURATION.AC_type;
Engine_loc = AC_CONFIGURATION.Engine_loc;
Engine_conf = AC_CONFIGURATION.Engine_conf;

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Plotting options
FS = Plot_Options.FS;
Video_3D = Plot_Options.Video_3D;
SAVE_FIGS = Plot_Options.SAVE_FIGS;

L_fus = Geo_tier.l_fus;
n_eng = Prop_data.n_eng;

Fig = Fig +1;
figure(Fig)
PLOTTING_UAV = plot_Models_ITER(Geo_tier,Body_Geo,meshData,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX);

%% plotting fuselage
meshData = PLOTTING_UAV.meshData; % FUSELAJE

% x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
% x_mesh_w2_New = PLOTTING_UAV.x_mesh_w2_New; % HTP
% % x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
% % x_mesh_VTP2_New = PLOTTING_UAV.x_mesh_VTP2_New; % VTP 2
%
% y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
% y_mesh_w2_New = PLOTTING_UAV.y_mesh_w2_New; % HTP
% % y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New; % VTP 1
% % y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New; % VTP 2
%
% z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING
% z_mesh_w2_New = PLOTTING_UAV.z_mesh_w2_New; % HTP
% % z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New; % VTP 1
% % z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New; % VTP 2

% patch(model{1}/1000,model{2}/1000,model{3}/1000,'b');
% axis equal; view(30,60);

%% Defines colors RGB
color_fus = COLOR_scheme.color_fus;
color_w1 = COLOR_scheme.color_w1;
color_w2 = COLOR_scheme.color_w2;
color_can = COLOR_scheme.color_can;
color_prop = COLOR_scheme.color_prop;
color_vtp = COLOR_scheme.color_vtp;
color_eng = COLOR_scheme.color_eng;
color_nac = COLOR_scheme.color_nac;
color_ac = COLOR_scheme.color_ac;

Delta_Z_stl = Body_Geo.Delta_Z_stl;
Delta_Y_stl = Body_Geo.Delta_Y_stl;

% Calculates the Aircraft Geometry
model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_ac);
ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;
[points_ac,triangles,tri norms] = import_stl_fast(OUTPUT_read_XLSX.Fuselage_flags.STL_ac,1);

max_x_ac = max(points_ac(:,1));
min_x_ac = min(points_ac(:,1));
max_y_ac = max(points_ac(:,2));
min_y_ac = min(points_ac(:,2));
max_z_ac = max(points_ac(:,3));
min_z_ac = min(points_ac(:,3));

Delta_X_stl = points_ac(1,1);
Delta_Y_stl = points_ac(1,2);
Delta_Z_stl = points_ac(1,3);

% Offsets
plotData{1} = (model{1}*ScaleFactor + Geo_tier.x_offset_CAD - min_x_ac).*SF_CAD;
plotData{2} = (model{2}*ScaleFactor - max_y_ac/2).*SF_CAD;
plotData{3} = (model{3}*ScaleFactor + Geo_tier.z_offset_CAD - Delta_Z_stl).*SF_CAD;
meshData_ac = plotData;
% Colors
C_ac(:,:,1) = color_ac(1)*ones(size(meshData_ac{1}));
C_ac(:,:,2) = color_ac(2)*ones(size(meshData_ac{1}));
C_ac(:,:,3) = color_ac(3)*ones(size(meshData_ac{1}));

if OUTPUT_read_XLSX.Fuselage_flags.AC_STL_Compare == 1;
    figure(Fig)
    mesh(meshData_ac{1},meshData_ac{2},meshData_ac{3},C_ac);
    case_AC = OUTPUT_read_XLSX.AC_Data_flags.case_AC;
    st = strcat(case_AC);
    title(st,'fontsize',FS)
    xlabel('y (m)')
    ylabel('z (m)')
    zlabel('z (m)')
    hold on
    Fig = Fig + 1;
end

% figure(Fig)
% if OUTPUT_read_XLSX.Fuselage_flags.AC_STL_Compare == 1
%     mesh(meshData_ac{1},meshData_ac{2},meshData_ac{3},C_ac);
% end
% st = strcat(case_AC);
% title(st,'fontsize',FS)
% xlabel('y (m)')
% ylabel('z (m)')
% zlabel('z (m)')
% hold on

%     hold on
%Plots Fuselage
% Asigns color to fuselage
model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_fus);
[points_fus,triangles,tri norms] = import_stl_fast(OUTPUT_read_XLSX.Fuselage_flags.STL_fus,1);

ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;

max_x_fus = max(points_fus(:,1));
min_x_fus = min(points_fus(:,1));
max_y_fus = max(points_fus(:,2));
min_y_fus = min(points_fus(:,2));
max_z_fus = max(points_fus(:,3));
min_z_fus = min(points_fus(:,3));

Delta_X_stl = points_fus(1,1);
Delta_Y_stl = points_fus(1,2);
Delta_Z_stl = points_fus(1,3);

plotData{1} = (model{1}*ScaleFactor + Geo_tier.x_offset_CAD - min_x_fus).*SF_CAD;
plotData{2} = (model{2}*ScaleFactor - max_y_fus/2).*SF_CAD;
plotData{3} = (model{3}*ScaleFactor + Geo_tier.z_offset_CAD - Delta_Z_stl).*SF_CAD;
meshData_fus = plotData;

C_fus(:,:,1) = color_fus(1)*ones(size(meshData_fus{1}));
C_fus(:,:,2) = color_fus(2)*ones(size(meshData_fus{1}));
C_fus(:,:,3) = color_fus(3)*ones(size(meshData_fus{1}));
% Plots fuselage

mesh(meshData_fus{1},meshData_fus{2},meshData_fus{3},C_fus);
case_AC = OUTPUT_read_XLSX.AC_Data_flags.case_AC;
st = strcat(case_AC);
title(st,'fontsize',FS)
xlabel('y (m)')
ylabel('z (m)')
zlabel('z (m)')
hold onr

% Aircraft type
switch AC_type
    case 1 % AC_type = 1 - flying wing
        plot_w1(PLOTTING_UAV,color_w1);
        x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
        y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
        z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING
        % Color
        C_w1(:,:,1) = color_w1(1)*ones(size(x_mesh_w1_New));
        C_w1(:,:,2) = color_w1(2)*ones(size(y_mesh_w1_New));
        C_w1(:,:,3) = color_w1(3)*ones(size(z_mesh_w1_New));
        mesh(x_mesh_w1_New,y_mesh_w1_New,z_mesh_w1_New,C_w1)

        %--------------------------- VTP ---------------------------------
        if VTP == 1
            if twin_VTP == 1

                % VTP 1
                y_offset_VTP = Geo_tier.y_offset_VTP
                x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
                y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New;
                z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New - y_offset_VTP; % VTP 1

                % VTP 2
                x_mesh_VTP2_New = PLOTTING_UAV.x_mesh_VTP2_New; % VTP 2
                y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New; % VTP 2
                %                 y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New - y_offset_VTP; % VTP 2
                z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New + y_offset_VTP; % VTP 2
                % Assigns color
                C_vtp(:,:,3) = color_vtp(1)*ones(size(x_mesh_VTP1_New));
                C_vtp(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP1_New));
                C_vtp(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP1_New));

                hm1 = mesh(x_mesh_VTP1_New,z_mesh_VTP1_New,y_mesh_VTP1_New,C_vtp)
                %                 rotate(hm1, [1 0 0], 90)
                hm2 = mesh(x_mesh_VTP2_New,z_mesh_VTP2_New,-y_mesh_VTP2_New,C_vtp)
                %                 rotate(hm2, [1 0 0], -90)
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
            end
        end

    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        %--------------------------- WING 1 ---------------------------------
        x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
        y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
        z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING
        % Color
        C_w1(:,:,1) = color_w1(1)*ones(size(x_mesh_w1_New));
        C_w1(:,:,2) = color_w1(2)*ones(size(y_mesh_w1_New));
        C_w1(:,:,3) = color_w1(3)*ones(size(z_mesh_w1_New));
        mesh(x_mesh_w1_New,y_mesh_w1_New,z_mesh_w1_New,C_w1)

        %--------------------------- WING 2 ---------------------------------
        x_mesh_w2_New = PLOTTING_UAV.x_mesh_w2_New; % HTP
        y_mesh_w2_New = PLOTTING_UAV.y_mesh_w2_New; % HTP
        z_mesh_w2_New = PLOTTING_UAV.z_mesh_w2_New; % HTP
        % Color
        C_w2(:,:,3) = color_w2(1)*ones(size(x_mesh_w2_New));
        C_w2(:,:,2) = color_w2(2)*ones(size(y_mesh_w2_New));
        C_w2(:,:,3) = color_w2(3)*ones(size(z_mesh_w2_New));
        mesh(x_mesh_w2_New,y_mesh_w2_New,z_mesh_w2_New,C_w2)

        %--------------------------- VTP ---------------------------------
        if VTP == 1
            if twin_VTP == 1

                % VTP 1
                y_offset_VTP = Geo_tier.y_offset_VTP
                x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
                y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New;
                z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New - y_offset_VTP; % VTP 1

                % VTP 2
                x_mesh_VTP2_New = PLOTTING_UAV.x_mesh_VTP2_New; % VTP 2
                y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New; % VTP 2
                %                 y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New - y_offset_VTP; % VTP 2
                z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New + y_offset_VTP; % VTP 2
                % Assigns color
                C_vtp(:,:,3) = color_vtp(1)*ones(size(x_mesh_VTP1_New));
                C_vtp(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP1_New));
                C_vtp(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP1_New));

                hm1 = mesh(x_mesh_VTP1_New,z_mesh_VTP1_New,y_mesh_VTP1_New,C_vtp)
                %                 rotate(hm1, [1 0 0], 90)
                hm2 = mesh(x_mesh_VTP2_New,z_mesh_VTP2_New,-y_mesh_VTP2_New,C_vtp)
                %                 rotate(hm2, [1 0 0], -90)
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
            end
        end

    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        %--------------------------- WING 1 ---------------------------------
        x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
        y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
        z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING
        % Color
        C_w1(:,:,1) = color_w1(1)*ones(size(x_mesh_w1_New));
        C_w1(:,:,2) = color_w1(2)*ones(size(y_mesh_w1_New));
        C_w1(:,:,3) = color_w1(3)*ones(size(z_mesh_w1_New));
        mesh(x_mesh_w1_New,y_mesh_w1_New,z_mesh_w1_New,C_w1)

        %--------------------------- WING 2 ---------------------------------
        x_mesh_w2_New = PLOTTING_UAV.x_mesh_w2_New; % HTP
        y_mesh_w2_New = PLOTTING_UAV.y_mesh_w2_New; % HTP
        z_mesh_w2_New = PLOTTING_UAV.z_mesh_w2_New; % HTP
        % Color
        C_w2(:,:,3) = color_w2(1)*ones(size(x_mesh_w2_New));
        C_w2(:,:,2) = color_w2(2)*ones(size(y_mesh_w2_New));
        C_w2(:,:,3) = color_w2(3)*ones(size(z_mesh_w2_New));
        mesh(x_mesh_w2_New,y_mesh_w2_New,z_mesh_w2_New,C_w2)

        %--------------------------- CANARD ---------------------------------
        x_mesh_can_New = PLOTTING_UAV.x_mesh_can_New; % HTP
        y_mesh_can_New = PLOTTING_UAV.y_mesh_can_New; % HTP
        z_mesh_can_New = PLOTTING_UAV.z_mesh_can_New; % HTP
        % Color
        C_can(:,:,3) = color_can(1)*ones(size(x_mesh_can_New));
        C_can(:,:,2) = color_can(2)*ones(size(y_mesh_can_New));
        C_can(:,:,3) = color_can(3)*ones(size(z_mesh_can_New));
        mesh(x_mesh_can_New,y_mesh_can_New,z_mesh_can_New,C_can)

        %--------------------------- VTP ---------------------------------
        if VTP == 1
            if twin_VTP == 1

                % VTP 1
                y_offset_VTP = Geo_tier.y_offset_VTP
                x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
                y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New;
                z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New - y_offset_VTP; % VTP 1

                % VTP 2
                x_mesh_VTP2_New = PLOTTING_UAV.x_mesh_VTP2_New; % VTP 2
                y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New; % VTP 2
                %                 y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New - y_offset_VTP; % VTP 2
                z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New + y_offset_VTP; % VTP 2
                % Assigns color
                C_vtp(:,:,3) = color_vtp(1)*ones(size(x_mesh_VTP1_New));
                C_vtp(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP1_New));
                C_vtp(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP1_New));

                hm1 = mesh(x_mesh_VTP1_New,z_mesh_VTP1_New,y_mesh_VTP1_New,C_vtp)
                %                 rotate(hm1, [1 0 0], 90)
                hm2 = mesh(x_mesh_VTP2_New,z_mesh_VTP2_New,-y_mesh_VTP2_New,C_vtp)
                %                 rotate(hm2, [1 0 0], -90)
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
            end
        end

    case 4 % AC_type = 4 - 2 surface: wing + V-tail

        %--------------------------- WING 1 ---------------------------------
        x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
        y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
        z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING

        x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
        y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
        z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING
        % Color
        C_w1(:,:,1) = color_w1(1)*ones(size(x_mesh_w1_New));
        C_w1(:,:,2) = color_w1(2)*ones(size(y_mesh_w1_New));
        C_w1(:,:,3) = color_w1(3)*ones(size(z_mesh_w1_New));
        mesh(x_mesh_w1_New,y_mesh_w1_New,z_mesh_w1_New,C_w1)

        %--------------------------- WING 2 ---------------------------------
        x_mesh_w2_New = PLOTTING_UAV.x_mesh_w2_New; % HTP
        y_mesh_w2_New = PLOTTING_UAV.y_mesh_w2_New; % HTP
        z_mesh_w2_New = PLOTTING_UAV.z_mesh_w2_New; % HTP
        % Color
        C_w2(:,:,3) = color_w2(1)*ones(size(x_mesh_w2_New));
        C_w2(:,:,2) = color_w2(2)*ones(size(y_mesh_w2_New));
        C_w2(:,:,3) = color_w2(3)*ones(size(z_mesh_w2_New));
        mesh(x_mesh_w2_New,y_mesh_w2_New,z_mesh_w2_New,C_w2)

    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        %--------------------------- WING 1 ---------------------------------
        x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
        y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
        z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING
        % Color
        C_w1(:,:,1) = color_w1(1)*ones(size(x_mesh_w1_New));
        C_w1(:,:,2) = color_w1(2)*ones(size(y_mesh_w1_New));
        C_w1(:,:,3) = color_w1(3)*ones(size(z_mesh_w1_New));
        mesh(x_mesh_w1_New,y_mesh_w1_New,z_mesh_w1_New,C_w1)

        %--------------------------- CANARD ---------------------------------
        x_mesh_can_New = PLOTTING_UAV.x_mesh_can_New; % HTP
        y_mesh_can_New = PLOTTING_UAV.y_mesh_can_New; % HTP
        z_mesh_can_New = PLOTTING_UAV.z_mesh_can_New; % HTP
        % Color
        C_can(:,:,3) = color_can(1)*ones(size(x_mesh_can_New));
        C_can(:,:,2) = color_can(2)*ones(size(y_mesh_can_New));
        C_can(:,:,3) = color_can(3)*ones(size(z_mesh_can_New));
        mesh(x_mesh_can_New,y_mesh_can_New,z_mesh_can_New,C_can)

        %--------------------------- WING 2 ---------------------------------
        x_mesh_w2_New = PLOTTING_UAV.x_mesh_w2_New; % HTP
        y_mesh_w2_New = PLOTTING_UAV.y_mesh_w2_New; % HTP
        z_mesh_w2_New = PLOTTING_UAV.z_mesh_w2_New; % HTP
        % Color
        C_w2(:,:,3) = color_w2(1)*ones(size(x_mesh_w2_New));
        C_w2(:,:,2) = color_w2(2)*ones(size(y_mesh_w2_New));
        C_w2(:,:,3) = color_w2(3)*ones(size(z_mesh_w2_New));
        mesh(x_mesh_w2_New,y_mesh_w2_New,z_mesh_w2_New,C_w2)

    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP

        %--------------------------- WING 1 ---------------------------------
        x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
        y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
        z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING
        % Color
        C_w1(:,:,1) = color_w1(1)*ones(size(x_mesh_w1_New));
        C_w1(:,:,2) = color_w1(2)*ones(size(y_mesh_w1_New));
        C_w1(:,:,3) = color_w1(3)*ones(size(z_mesh_w1_New));
        mesh(x_mesh_w1_New,y_mesh_w1_New,z_mesh_w1_New,C_w1)

        %--------------------------- CANARD ---------------------------------
        x_mesh_can_New = PLOTTING_UAV.x_mesh_can_New; % HTP
        y_mesh_can_New = PLOTTING_UAV.y_mesh_can_New; % HTP
        z_mesh_can_New = PLOTTING_UAV.z_mesh_can_New; % HTP
        % Color
        C_can(:,:,3) = color_can(1)*ones(size(x_mesh_can_New));
        C_can(:,:,2) = color_can(2)*ones(size(y_mesh_can_New));
        C_can(:,:,3) = color_can(3)*ones(size(z_mesh_can_New));
        mesh(x_mesh_can_New,y_mesh_can_New,z_mesh_can_New,C_can)

        %--------------------------- VTP ---------------------------------
        if VTP == 1
            if twin_VTP == 1

                % VTP 1
                y_offset_VTP = Geo_tier.y_offset_VTP
                x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
                y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New;
                z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New - y_offset_VTP; % VTP 1

                % VTP 2
                x_mesh_VTP2_New = PLOTTING_UAV.x_mesh_VTP2_New; % VTP 2
                y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New; % VTP 2
                %                 y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New - y_offset_VTP; % VTP 2
                z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New + y_offset_VTP; % VTP 2
                % Assigns color
                C_vtp(:,:,3) = color_vtp(1)*ones(size(x_mesh_VTP1_New));
                C_vtp(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP1_New));
                C_vtp(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP1_New));

                hm1 = mesh(x_mesh_VTP1_New,z_mesh_VTP1_New,y_mesh_VTP1_New,C_vtp)
                %                 rotate(hm1, [1 0 0], 90)
                hm2 = mesh(x_mesh_VTP2_New,z_mesh_VTP2_New,-y_mesh_VTP2_New,C_vtp)
                %                 rotate(hm2, [1 0 0], -90)
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
            end
        end
end

% Engine Configuration
% Engine location
switch Engine_loc
    case 1 % Engine_loc = 1 - under wings
        if n_eng < 2
            disp('For engines at the wing tips 2 engines need to be defined. Please do correct Input Data');
            pause
        end
        % Engine
        Eng_L_x_w1 = PLOTTING_UAV.Eng_L_x_w1;
        Eng_L_y_w1 = PLOTTING_UAV.Eng_L_y_w1;
        Eng_L_z_w1 = PLOTTING_UAV.Eng_L_z_w1;
        Eng_L_x_w1_s = PLOTTING_UAV.Eng_L_x_w1_s;
        Eng_L_y_w1_s = PLOTTING_UAV.Eng_L_y_w1_s;
        Eng_L_z_w1_s = PLOTTING_UAV.Eng_L_z_w1_s;

        Eng_R_x_w1 = PLOTTING_UAV.Eng_R_x_w1;
        Eng_R_y_w1 = PLOTTING_UAV.Eng_R_y_w1;
        Eng_R_z_w1 = PLOTTING_UAV.Eng_R_z_w1;
        Eng_R_x_w1_s = PLOTTING_UAV.Eng_R_x_w1_s;
        Eng_R_y_w1_s = PLOTTING_UAV.Eng_R_y_w1_s;
        Eng_R_z_w1_s = PLOTTING_UAV.Eng_R_z_w1_s;

        % Plots disck actuators Engines
        pL_w1 = patch(Eng_L_x_w1,Eng_L_y_w1,Eng_L_z_w1,'y');
        set(pL_w1,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1,'EraseMode','normal');
        pL_w1_s = patch(Eng_L_x_w1_s,Eng_L_y_w1_s,Eng_L_z_w1_s,'y');
        set(pL_w1_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1_s,'EraseMode','normal');

        pR_w1 = patch(Eng_R_x_w1,Eng_R_y_w1,Eng_R_z_w1,'y');
        set(pR_w1,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_w1,'EraseMode','normal');
        pR_w1_s = patch(Eng_R_x_w1_s,Eng_R_y_w1_s,Eng_R_z_w1_s,'y');
        set(pR_w1_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_w1_s,'EraseMode','normal');

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


        % Nacelle
        Nac_L_x_w1 = PLOTTING_UAV.Nac_L_x_w1;
        Nac_L_y_w1 = PLOTTING_UAV.Nac_L_y_w1;
        Nac_L_z_w1 = PLOTTING_UAV.Nac_L_z_w1;
        Nac_L_x_w1_s = PLOTTING_UAV.Nac_L_x_w1_s;
        Nac_L_y_w1_s = PLOTTING_UAV.Nac_L_y_w1_s;
        Nac_L_z_w1_s = PLOTTING_UAV.Nac_L_z_w1_s;

        Nac_R_x_w1 = PLOTTING_UAV.Nac_R_x_w1;
        Nac_R_y_w1 = PLOTTING_UAV.Nac_R_y_w1;
        Nac_R_z_w1 = PLOTTING_UAV.Nac_R_z_w1;
        Nac_R_x_w1_s = PLOTTING_UAV.Nac_R_x_w1_s;
        Nac_R_y_w1_s = PLOTTING_UAV.Nac_R_y_w1_s;
        Nac_R_z_w1_s = PLOTTING_UAV.Nac_R_z_w1_s;

        % Plots disck actuators Engines
        pL_w1 = patch(Nac_L_x_w1,Nac_L_y_w1,Nac_L_z_w1,'y');
        set(pL_w1,'FaceColor',color_nac,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1,'EraseMode','normal');
        pL_w1_s = patch(Nac_L_x_w1_s,Nac_L_y_w1_s,Nac_L_z_w1_s,'y');
        set(pL_w1_s,'FaceColor',color_nac,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1_s,'EraseMode','normal');

        pR_w1 = patch(Nac_R_x_w1,Nac_R_y_w1,Nac_R_z_w1,'y');
        set(pR_w1,'FaceColor',color_nac,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_w1,'EraseMode','normal');
        pR_w1_s = patch(Nac_R_x_w1_s,Nac_R_y_w1_s,Nac_R_z_w1_s,'y');
        set(pR_w1_s,'FaceColor',color_nac,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_w1_s,'EraseMode','normal');

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

    case 2 % Engine_loc = 2 - fuselage front
        % Engine
        Eng_L_x_w1 = PLOTTING_UAV.Eng_L_x_w1;
        Eng_L_y_w1 = PLOTTING_UAV.Eng_L_y_w1;
        Eng_L_z_w1 = PLOTTING_UAV.Eng_L_z_w1;
        Eng_L_x_w1_s = PLOTTING_UAV.Eng_L_x_w1_s;
        Eng_L_y_w1_s = PLOTTING_UAV.Eng_L_y_w1_s;
        Eng_L_z_w1_s = PLOTTING_UAV.Eng_L_z_w1_s;

        Prop_L_x_w1 = PLOTTING_UAV.Prop_L_x_w1;
        Prop_L_y_w1 = PLOTTING_UAV.Prop_L_y_w1;
        Prop_L_z_w1 = PLOTTING_UAV.Prop_L_z_w1;
        Prop_L_x_w1_s = PLOTTING_UAV.Prop_L_x_w1_s;
        Prop_L_y_w1_s = PLOTTING_UAV.Prop_L_y_w1_s;
        Prop_L_z_w1_s = PLOTTING_UAV.Prop_L_z_w1_s;

        % Plots disck actuators Engines
        pL_w1 = patch(Eng_L_x_w1,Eng_L_y_w1,Eng_L_z_w1,'y');
        set(pL_w1,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1,'EraseMode','normal');
        pL_w1_s = patch(Eng_L_x_w1_s,Eng_L_y_w1_s,Eng_L_z_w1_s,'y');
        set(pL_w1_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1_s,'EraseMode','normal');

    case 3 % Engine_loc = 3 - fuselage rear

        % Engine
        Eng_L_x_w1 = PLOTTING_UAV.Eng_L_x_w1;
        Eng_L_y_w1 = PLOTTING_UAV.Eng_L_y_w1;
        Eng_L_z_w1 = PLOTTING_UAV.Eng_L_z_w1;
        Eng_L_x_w1_s = PLOTTING_UAV.Eng_L_x_w1_s;
        Eng_L_y_w1_s = PLOTTING_UAV.Eng_L_y_w1_s;
        Eng_L_z_w1_s = PLOTTING_UAV.Eng_L_z_w1_s;

        Prop_L_x_w1 = PLOTTING_UAV.Prop_L_x_w1;
        Prop_L_y_w1 = PLOTTING_UAV.Prop_L_y_w1;
        Prop_L_z_w1 = PLOTTING_UAV.Prop_L_z_w1;
        Prop_L_x_w1_s = PLOTTING_UAV.Prop_L_x_w1_s;
        Prop_L_y_w1_s = PLOTTING_UAV.Prop_L_y_w1_s;
        Prop_L_z_w1_s = PLOTTING_UAV.Prop_L_z_w1_s;

        % Plots disck actuators Engines
        pL_w1 = patch(Eng_L_x_w1,Eng_L_y_w1,Eng_L_z_w1,'y');
        set(pL_w1,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1,'EraseMode','normal');
        pL_w1_s = patch(Eng_L_x_w1_s,Eng_L_y_w1_s,Eng_L_z_w1_s,'y');
        set(pL_w1_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1_s,'EraseMode','normal');

    case 4 % Engine_loc = 4 - wingtips
        if n_eng < 2
            disp('For engines at the wing tips 2 engines need to ve defined. Please do correct Input Data');
            pause
        end
        % Engine
        Eng_L_x_w1 = PLOTTING_UAV.Eng_L_x_w1;
        Eng_L_y_w1 = PLOTTING_UAV.Eng_L_y_w1;
        Eng_L_z_w1 = PLOTTING_UAV.Eng_L_z_w1;
        Eng_L_x_w1_s = PLOTTING_UAV.Eng_L_x_w1_s;
        Eng_L_y_w1_s = PLOTTING_UAV.Eng_L_y_w1_s;
        Eng_L_z_w1_s = PLOTTING_UAV.Eng_L_z_w1_s;

        Eng_R_x_w1 = PLOTTING_UAV.Eng_R_x_w1;
        Eng_R_y_w1 = PLOTTING_UAV.Eng_R_y_w1;
        Eng_R_z_w1 = PLOTTING_UAV.Eng_R_z_w1;
        Eng_R_x_w1_s = PLOTTING_UAV.Eng_R_x_w1_s;
        Eng_R_y_w1_s = PLOTTING_UAV.Eng_R_y_w1_s;
        Eng_R_z_w1_s = PLOTTING_UAV.Eng_R_z_w1_s;

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

        % Plots disck actuators Engines
        pL_w1 = patch(Eng_L_x_w1,Eng_L_y_w1,Eng_L_z_w1,'y');
        set(pL_w1,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1,'EraseMode','normal');
        pL_w1_s = patch(Eng_L_x_w1_s,Eng_L_y_w1_s,Eng_L_z_w1_s,'y');
        set(pL_w1_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1_s,'EraseMode','normal');

        pR_w1 = patch(Eng_R_x_w1,Eng_R_y_w1,Eng_R_z_w1,'y');
        set(pR_w1,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_w1,'EraseMode','normal');
        pR_w1_s = patch(Eng_R_x_w1_s,Eng_R_y_w1_s,Eng_R_z_w1_s,'y');
        set(pR_w1_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_w1_s,'EraseMode','normal');

    case 5 % Engine_loc = 5 - wingtips 2 pairs
        if n_eng < 4
            disp('For engines at the wing tips 2 engines need to ve defined. Please do correct Input Data');
            pause
        end

        % Wing Engine
        Eng_L_x_w1 = PLOTTING_UAV.Eng_L_x_w1;
        Eng_L_y_w1 = PLOTTING_UAV.Eng_L_y_w1;
        Eng_L_z_w1 = PLOTTING_UAV.Eng_L_z_w1;
        Eng_L_x_w1_s = PLOTTING_UAV.Eng_L_x_w1_s;
        Eng_L_y_w1_s = PLOTTING_UAV.Eng_L_y_w1_s;
        Eng_L_z_w1_s = PLOTTING_UAV.Eng_L_z_w1_s;

        Eng_R_x_w1 = PLOTTING_UAV.Eng_R_x_w1;
        Eng_R_y_w1 = PLOTTING_UAV.Eng_R_y_w1;
        Eng_R_z_w1 = PLOTTING_UAV.Eng_R_z_w1;
        Eng_R_x_w1_s = PLOTTING_UAV.Eng_R_x_w1_s;
        Eng_R_y_w1_s = PLOTTING_UAV.Eng_R_y_w1_s;
        Eng_R_z_w1_s = PLOTTING_UAV.Eng_R_z_w1_s;

        % Plots disck actuators Engines
        pL_w1 = patch(Eng_L_x_w1,Eng_L_y_w1,Eng_L_z_w1,'y');
        set(pL_w1,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1,'EraseMode','normal');
        pL_w1_s = patch(Eng_L_x_w1_s,Eng_L_y_w1_s,Eng_L_z_w1_s,'y');
        set(pL_w1_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_w1_s,'EraseMode','normal');

        pR_w1 = patch(Eng_R_x_w1,Eng_R_y_w1,Eng_R_z_w1,'y');
        set(pR_w1,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_w1,'EraseMode','normal');
        pR_w1_s = patch(Eng_R_x_w1_s,Eng_R_y_w1_s,Eng_R_z_w1_s,'y');
        set(pR_w1_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_w1_s,'EraseMode','normal');

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

        % Second Pair of  Engines
        Eng_L_x_w2 = PLOTTING_UAV.Eng_L_x_w2;
        Eng_L_y_w2 = PLOTTING_UAV.Eng_L_y_w2;
        Eng_L_z_w2 = PLOTTING_UAV.Eng_L_z_w2;
        Eng_L_x_w2_s = PLOTTING_UAV.Eng_L_x_w2_s;
        Eng_L_y_w2_s = PLOTTING_UAV.Eng_L_y_w2_s;
        Eng_L_z_w2_s = PLOTTING_UAV.Eng_L_z_w2_s;

        Eng_R_x_w2 = PLOTTING_UAV.Eng_R_x_w2;
        Eng_R_y_w2 = PLOTTING_UAV.Eng_R_y_w2;
        Eng_R_z_w2 = PLOTTING_UAV.Eng_R_z_w2;
        Eng_R_x_w2_s = PLOTTING_UAV.Eng_R_x_w2_s;
        Eng_R_y_w2_s = PLOTTING_UAV.Eng_R_y_w2_s;
        Eng_R_z_w2_s = PLOTTING_UAV.Eng_R_z_w2_s;

        % Plots disck actuators Engines
        pL_w2 = patch(Eng_L_x_w2,Eng_L_y_w2,Eng_L_z_w2,'y');
        set(pL_w2,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_can,'EraseMode','normal');
        pL_w2_s = patch(Eng_L_x_w2_s,Eng_L_y_w2_s,Eng_L_z_w2_s,'y');
        set(pL_w2_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pL_can_s,'EraseMode','normal');

        pR_w2 = patch(Eng_R_x_w2,Eng_R_y_w2,Eng_R_z_w2,'y');
        set(pR_w2,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_can,'EraseMode','normal');
        pR_w2_s = patch(Eng_R_x_w2_s,Eng_R_y_w2_s,Eng_R_z_w2_s,'y');
        set(pR_w2_s,'FaceColor',color_eng,'FaceLighting','phong','EdgeLighting','phong');
        % set(pR_can_s,'EraseMode','normal');

        Prop_L_x_w2 = PLOTTING_UAV.Prop_L_x_w2;
        Prop_L_y_w2 = PLOTTING_UAV.Prop_L_y_w2;
        Prop_L_z_w2 = PLOTTING_UAV.Prop_L_z_w2;
        Prop_L_x_w2_s = PLOTTING_UAV.Prop_L_x_w2_s;
        Prop_L_y_w2_s = PLOTTING_UAV.Prop_L_y_w2_s;
        Prop_L_z_w2_s = PLOTTING_UAV.Prop_L_z_w2_s;

        Prop_R_x_w2 = PLOTTING_UAV.Prop_R_x_w2;
        Prop_R_y_w2 = PLOTTING_UAV.Prop_R_y_w2;
        Prop_R_z_w2 = PLOTTING_UAV.Prop_R_z_w2;
        Prop_R_x_w2_s = PLOTTING_UAV.Prop_R_x_w2_s;
        Prop_R_y_w2_s = PLOTTING_UAV.Prop_R_y_w2_s;
        Prop_R_z_w2_s = PLOTTING_UAV.Prop_R_z_w2_s;

end

% Engine Configuration
switch Engine_conf
    case 1 % Engine_conf = 1 - pusher prop
        % Depending the number of engines
        if n_eng == 1
            % Plots disck actuators wing 2
            pL_w1 = patch(Prop_L_x_w1,Prop_L_y_w1,Prop_L_z_w1,'y');
            set(pL_w1,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1,'EraseMode','normal');
            pL_w1_s = patch(Prop_L_x_w1_s,Prop_L_y_w1_s,Prop_L_z_w1_s,'y');
            set(pL_w1_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1_s,'EraseMode','normal');
        elseif n_eng == 2
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
        elseif n_eng == 4
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

            % Plots disck actuators wing 2
            pL_w2 = patch(Prop_L_x_w2,Prop_L_y_w2,Prop_L_z_w2,'y');
            set(pL_w2,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1,'EraseMode','normal');
            pL_w2_s = patch(Prop_L_x_w2_s,Prop_L_y_w2_s,Prop_L_z_w2_s,'y');
            set(pL_w2_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1_s,'EraseMode','normal');

            pR_w2 = patch(Prop_R_x_w2,Prop_R_y_w2,Prop_R_z_w2,'y');
            set(pR_w2,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pR_w1,'EraseMode','normal');
            pR_w2_s = patch(Prop_R_x_w2_s,Prop_R_y_w2_s,Prop_R_z_w2_s,'y');
            set(pR_w2_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pR_w1_s,'EraseMode','normal');
end
    case 2 % Engine_conf = 2 - puller prop

        % Depending the number of engines
        if n_eng == 1
            % Plots disck actuators wing 2
            pL_w1 = patch(Prop_L_x_w1,Prop_L_y_w1,Prop_L_z_w1,'y');
            set(pL_w1,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1,'EraseMode','normal');
            pL_w1_s = patch(Prop_L_x_w1_s,Prop_L_y_w1_s,Prop_L_z_w1_s,'y');
            set(pL_w1_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1_s,'EraseMode','normal');

        elseif n_eng == 2
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
        elseif n_eng == 4
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

            % Plots disck actuators wing 2
            pL_w2 = patch(Prop_L_x_w2,Prop_L_y_w2,Prop_L_z_w2,'y');
            set(pL_w2,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1,'EraseMode','normal');
            pL_w2_s = patch(Prop_L_x_w2_s,Prop_L_y_w2_s,Prop_L_z_w2_s,'y');
            set(pL_w2_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1_s,'EraseMode','normal');

            pR_w2 = patch(Prop_R_x_w2,Prop_R_y_w2,Prop_R_z_w2,'y');
            set(pR_w2,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pR_w1,'EraseMode','normal');
            pR_w2_s = patch(Prop_R_x_w2_s,Prop_R_y_w2_s,Prop_R_z_w2_s,'y');
            set(pR_w2_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pR_w1_s,'EraseMode','normal');
        end
    case 3 % Engine_conf = 3 - turbine ASSUMES A REALLY SMALL DISK
        % Depending the number of engines
        if n_eng == 1
            % Plots disck actuators wing 2
            pL_w1 = patch(Prop_L_x_w1,Prop_L_y_w1,Prop_L_z_w1,'y');
            set(pL_w1,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1,'EraseMode','normal');
            pL_w1_s = patch(Prop_L_x_w1_s,Prop_L_y_w1_s,Prop_L_z_w1_s,'y');
            set(pL_w1_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1_s,'EraseMode','normal');

        elseif n_eng == 2
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
        elseif n_eng == 4
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

            % Plots disck actuators wing 2
            pL_w2 = patch(Prop_L_x_w2,Prop_L_y_w2,Prop_L_z_w2,'y');
            set(pL_w2,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1,'EraseMode','normal');
            pL_w2_s = patch(Prop_L_x_w2_s,Prop_L_y_w2_s,Prop_L_z_w2_s,'y');
            set(pL_w2_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pL_w1_s,'EraseMode','normal');

            pR_w2 = patch(Prop_R_x_w2,Prop_R_y_w2,Prop_R_z_w2,'y');
            set(pR_w2,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pR_w1,'EraseMode','normal');
            pR_w2_s = patch(Prop_R_x_w2_s,Prop_R_y_w2_s,Prop_R_z_w2_s,'y');
            set(pR_w2_s,'FaceColor',color_prop,'FaceLighting','phong','EdgeLighting','phong');
            % set(pR_w1_s,'EraseMode','normal');
        end
end

grid on
hold off
% axis equal
% grid on

%% Defines limits of plots axis
z_max_fus = PLOTTING_UAV.z_max_fus;
z_min_fus = PLOTTING_UAV.z_min_fus;
% Creates Margin to all sizes of plots around the airplane
Delta_Plot = 0.05;
b_w1 = Geo_tier.b_w1;
MAX_X = L_fus + Delta_Plot;
MAX_Y = b_w1 + Delta_Plot;
MAX_Z = z_max_fus + Delta_Plot;
x_min_PLOT = -Delta_Plot;
x_max_PLOT = MAX_X;
y_min_PLOT = -MAX_Y/2 ;
y_max_PLOT = MAX_Y/2 ;
z_min_PLOT = (z_min_fus);
z_max_PLOT = MAX_Z;
axis([x_min_PLOT x_max_PLOT y_min_PLOT y_max_PLOT z_min_PLOT z_max_PLOT]);
axis equal
PLOTTING_DATA = PLOTTING_UAV;
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat(num2str(Performance.V),'_m_s');
    name   = strcat(prefix,st);
    saveas(gcf,name,'fig');
    saveas(gcf,name,'pdf');
    saveas(gcf,name,'bmp');
    saveas(gcf,name,'png');
end

if Video_3D == 1
    switch case_AC
        case 1 % case_AC = 1 - EMERGENTIA 1:1
            VID_TXT = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        case 2 % case_AC = 2 - EMERGENTIA 1:2
            VID_TXT = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        case 3 % case_AC = 3 - PEPIÑO XXL
            VID_TXT = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    end
    OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
    CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],VID_TXT,OptionZ)
end

% model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_ac);
% ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
% plotData{1} = model{1}*ScaleFactor + Geo_tier.x_offset_CAD;
% plotData{2} = model{2}*ScaleFactor;
% plotData{3} = model{3}*ScaleFactor + Geo_tier.z_offset_CAD;
% meshData = plotData;
Fig = Fig + 1;