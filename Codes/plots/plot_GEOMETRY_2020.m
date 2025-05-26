% function that plots geometry
function Fig = plot_GEOMETRY_2020(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC,Performance)
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
PLOTTING_UAV = plot_Models_VTOL_ITER(Geo_tier,Body_Geo,meshData,Prop_data,AC_CONFIGURATION);

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
color_prop = COLOR_scheme.color_prop;
color_vtp = COLOR_scheme.color_vtp;
color_eng = COLOR_scheme.color_eng;

figure(Fig)
% meshData = stl2matlab('ProVant.stl');
% Asigns color to fuselage
C_fus(:,:,1) = color_fus(1)*ones(size(meshData{1}));
C_fus(:,:,2) = color_fus(2)*ones(size(meshData{1}));
C_fus(:,:,3) = color_fus(3)*ones(size(meshData{1}));

% Plots fuselage
mesh(meshData{1},meshData{2},meshData{3},C_fus);

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        st = strcat('EMERGENTIA UAV 100%');
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        st = strcat('EMERGENTIA UAV 50%');
    case 3 % case_AC = 3 - PEPIÑO XXL
        st = strcat('PEPIÑO XXL UAV');
    case 4 % case_AC = 4 - COMMERCIAL
        st = strcat('COMMERCIAL AIRPLANE');
    case 5 % case_AC = 5 - WIG
        st = strcat('WIG AIRPLANE');
    case 6 % case_AC = 6 - CERVERA
        st = strcat('CERVERA');
end

title(st,'fontsize',FS)
xlabel('y (m)')
ylabel('z (m)')
zlabel('z (m)')
hold on

% Aircraft type
switch AC_type
    case 1 % AC_type = 1 - flying wing
        x_mesh_w1_New = PLOTTING_UAV.x_mesh_w1_New; % WING
        y_mesh_w1_New = PLOTTING_UAV.y_mesh_w1_New; % WING
        z_mesh_w1_New = PLOTTING_UAV.z_mesh_w1_New; % WING
        % Color
        C_w1(:,:,1) = color_w1(1)*ones(size(x_mesh_w1_New));
        C_w1(:,:,2) = color_w1(2)*ones(size(y_mesh_w1_New));
        C_w1(:,:,3) = color_w1(3)*ones(size(z_mesh_w1_New));
        mesh(x_mesh_w1_New,y_mesh_w1_New,z_mesh_w1_New,C_w1)
        
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
                y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New ;
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
                x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
                y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New; % VTP 1
                z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New; % VTP 1
                % Assigns color
                C_vtp(:,:,3) = color_vtp(1)*ones(size(x_mesh_VTP1_New));
                C_vtp(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP1_New));
                C_vtp(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP1_New));
                % plots
                mesh(x_mesh_VTP1_New,y_mesh_VTP1_New,z_mesh_VTP1_New,C_vtp)
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
        %--------------------------- VTP ---------------------------------
        if n_VTP == 2
            x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
            y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New; % VTP 1
            z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New; % VTP 1
            x_mesh_VTP2_New = PLOTTING_UAV.x_mesh_VTP2_New; % VTP 2
            y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New; % VTP 2
            z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New; % VTP 2
            % Assigns color
            C_vtp(:,:,3) = color_vtp(1)*ones(size(x_mesh_VTP1_New));
            C_vtp(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP1_New));
            C_vtp(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP1_New));
            mesh(x_mesh_VTP1_New,y_mesh_VTP1_New,z_mesh_VTP1_New,C_vtp)
            mesh(x_mesh_VTP2_New,y_mesh_VTP2_New,z_mesh_VTP2_New,C_vtp)
        elseif n_VTP == 1
            % Color
            x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
            y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New; % VTP 1
            z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New; % VTP 1
            % Assigns color
            C_vtp(:,:,3) = color_vtp(1)*ones(size(x_mesh_VTP1_New));
            C_vtp(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP1_New));
            C_vtp(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP1_New));
            % plots
            mesh(x_mesh_VTP1_New,y_mesh_VTP1_New,z_mesh_VTP1_New,C_vtp)
        end
     
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        
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
        %--------------------------- WING 2 ---------------------------------
        x_mesh_w2_New = PLOTTING_UAV.x_mesh_w2_New; % HTP
        y_mesh_w2_New = PLOTTING_UAV.y_mesh_w2_New; % HTP
        z_mesh_w2_New = PLOTTING_UAV.z_mesh_w2_New; % HTP
        % Color
        C_w2(:,:,3) = color_w2(1)*ones(size(x_mesh_w2_New));
        C_w2(:,:,2) = color_w2(2)*ones(size(y_mesh_w2_New));
        C_w2(:,:,3) = color_w2(3)*ones(size(z_mesh_w2_New));
        mesh(x_mesh_w2_New,y_mesh_w2_New,z_mesh_w2_New,C_w2)
end

% Engine Configuration
% Engine location
switch Engine_loc
    case 1 % Engine_loc = 1 - under wings
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
        end
    case 3 % Engine_conf = 3 - turbine        
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
    switch case_AC
        case 1 % case_AC = 1 - EMERGENTIA 1:1
            prefix = strcat('ProVant4_');
        case 2 % case_AC = 2 - EMERGENTIA 1:2
            prefix = strcat('ProVant4_');
        case 3 % case_AC = 3 - PEPIÑO XXL
            prefix = strcat('PEPINOXXL_');
    end
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
            VID_TXT = strcat('ProVant4');
        case 2 % case_AC = 2 - EMERGENTIA 1:2
            VID_TXT = strcat('ProVant4');
        case 3 % case_AC = 3 - PEPIÑO XXL
            VID_TXT = strcat('PEPINOXXL');
    end
    OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
    CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],VID_TXT,OptionZ)
end

Fig = Fig + 1;