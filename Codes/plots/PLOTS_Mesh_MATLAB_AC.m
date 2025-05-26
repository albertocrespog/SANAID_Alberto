function PLOTS_Mesh_MATLAB_AC(OUTPUT_read_XLSX,PLOTTING_UAV,COLOR_scheme,AC_type,Engine_loc,AC_CONFIGURATION,Geo_tier)

%% Defines colors RGB
color_fus = COLOR_scheme.color_fus;
color_w1 = COLOR_scheme.color_w1;
color_HTP = COLOR_scheme.color_HTP;
color_vee = COLOR_scheme.color_vee;
color_vee2 = COLOR_scheme.color_vee2;
color_can = COLOR_scheme.color_can;
color_prop = COLOR_scheme.color_prop;
color_vtp = COLOR_scheme.color_vtp;
color_eng = COLOR_scheme.color_eng;
color_nac = COLOR_scheme.color_nac;
color_ac = COLOR_scheme.color_ac;

% Plots Fuselage
% Asigns color to fuselage
model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_fus);

if OUTPUT_read_XLSX.Fuselage_flags.STL_reading_test == 1
    max_x_fus = max(max(model{1}));
    min_x_fus = min(min(model{1}));
    max_y_fus = max(max(model{2}));
    min_y_fus = min(min(model{2}));
    z_max_fus = max(max(model{3}));
    z_min_fus = min(min(model{3}));
    
    Delta_X_stl_fus = min_x_fus;
    Delta_Y_stl_fus = -(max_y_fus - (max_y_fus - min_y_fus)/2)/2;
    Delta_Z_stl_fus = 0;
    dummy=0;
else
    [points_fus,triangles,tri_norms] = import_stl_fast_v2(OUTPUT_read_XLSX.Fuselage_flags.STL_fus,1);
    max_x_fus = max(points_fus(:,1));
    min_x_fus = min(points_fus(:,1));
    max_y_fus = max(points_fus(:,2));
    min_y_fus = min(points_fus(:,2));
    z_max_fus = max(points_fus(:,3));
    z_min_fus = min(points_fus(:,3));

    Delta_X_stl_fus = points_fus(1,1);
    Delta_Y_stl_fus = points_fus(1,2);
    Delta_Y_stl_fus = 0;
    Delta_Z_stl_fus = points_fus(1,3);

    % Calculate stats (min, max, and deltas for each axis)
    points = points_fus;
    min_vals = min(points);
    max_vals = max(points);
    deltas = points - min(points);

    dummy=0;
end

ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;

% plotData{1} = (model{1}*ScaleFactor - Delta_X_stl_fus).*SF_CAD;
% plotData{2} = (model{2}*ScaleFactor - Delta_Y_stl_fus).*SF_CAD;
% % plotData{3} = (model{3}*ScaleFactor + Geo_tier.z_offset_CAD).*SF_CAD;
% plotData{3} = (model{3}*ScaleFactor - Delta_Z_stl_fus).*SF_CAD;

plotData{1} = (model{1}*ScaleFactor).*SF_CAD;
plotData{2} = (model{2}*ScaleFactor).*SF_CAD;
% plotData{3} = (model{3}*ScaleFactor + Geo_tier.z_offset_CAD).*SF_CAD;
plotData{3} = (model{3}*ScaleFactor).*SF_CAD;

meshData_fus = plotData;

C_fus(:,:,1) = color_fus(1)*ones(size(meshData_fus{1}));
C_fus(:,:,2) = color_fus(2)*ones(size(meshData_fus{1}));
C_fus(:,:,3) = color_fus(3)*ones(size(meshData_fus{1}));
% Plots fuselage

% First Mesh Fuselage
mesh(meshData_fus{1},meshData_fus{2},meshData_fus{3},C_fus);
case_AC = OUTPUT_read_XLSX.AC_Data_flags.case_AC;
% st = strcat(case_AC);
% title(st,'fontsize',FS)
% title('Aircrat Original and SANAID Re-Construction')
% xlabel('y (m)')
% ylabel('z (m)')
% zlabel('z (m)')
% hold on

% Initializes
MESH_AC.dummy = 0;
% Aircraft type
switch AC_type
    case 1 % AC_type = 1 - flying wing
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);        

    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- HTP ------------------------------
        MESH_AC = plot_HTP(PLOTTING_UAV,color_HTP,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);        

    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- CANARD ---------------------------------
        MESH_AC = plot_can(PLOTTING_UAV,color_HTP,MESH_AC);
        %--------------------------- HTP ------------------------------
        MESH_AC = plot_HTP(PLOTTING_UAV,color_w2,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);        

    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- Vee ------------------------------
        MESH_AC = plot_vee(PLOTTING_UAV,color_vee,MESH_AC);

    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- CANARD ---------------------------------
        MESH_AC = plot_can(PLOTTING_UAV,color_can,MESH_AC);
        %--------------------------- Vee ------------------------------
        MESH_AC = plot_vee(PLOTTING_UAV,color_vee,MESH_AC);
        
    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- CANARD ---------------------------------
        MESH_AC = plot_can(PLOTTING_UAV,color_can,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);  


    case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- Vee ------------------------------
        MESH_AC = plot_vee(PLOTTING_UAV,color_vee,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);  

    case 8 % AC_type = 8 - 2 surface: wing + V-tail+V-tail

        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- Vee ------------------------------
        MESH_AC = plot_vee(PLOTTING_UAV,color_vee,MESH_AC);
        %--------------------------- Vee2 ---------------------------------
        MESH_AC = plot_vee2(PLOTTING_UAV,color_vee2,MESH_AC);

end

% % Engine Configuration
% % Engine location
% switch Engine_loc    
%     case 1 % Engine_loc = 1 - under wings
%         if n_eng < 2
%             disp('For engines at the wing tips 2 engines need to be defined. Please do correct Input Data');
%             pause
%         end
%        if n_eng == 2
%            % Engine
%            plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
%            % Nacelle
%            plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
%            % Propeller disk
%            plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
%         elseif n_eng == 4
%             % Engine and Nacelle first pair
%            % Engine
%            plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
%            % Nacelle
%            plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
%            % Propeller disk
%            plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
%            % Engine and Nacelle second pair
%            % Engine
%            plot_2engine2(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
%            % Nacelle
%            plot_2nacelle2(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
%            % Propeller disk
%            plot_2prop2(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop);
%         end
%     case 2 % Engine_loc = 2 - fuselage front
% 
%            % Engine
%            plot_1engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
%            % Nacelle
%            plot_1nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
%            % Propeller disk
%            plot_1prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
% 
%     case 3 % Engine_loc = 3 - fuselage rear
% 
%            % Engine
%            plot_1engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
%            % Nacelle
%            plot_1nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
%            % Propeller disk
%            plot_1prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
% 
%     case 4 % Engine_loc = 4 - wingtips
%         if n_eng < 2
%             disp('For engines at the wing tips 2 engines need to ve defined. Please do correct Input Data');
%             pause
%         end
%         % Engine
%         plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
%         % Nacelle
%         plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
%         % Propeller disk
%         plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
% 
%     case 5 % Engine_loc = 5 - wingtips 2 pairs
%         if n_eng < 4
%             disp('For engines at the wing tips 2 engines need to be defined. Please do correct Input Data');
%             pause
%         end
%         if n_eng == 2
%             % Engine
%             plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
%             % Nacelle
%             plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
%             % Propeller disk
%             plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
%         elseif n_eng == 4
%             % Engine and Nacelle first pair
%             % Engine
%             plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
%             % Nacelle
%             plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
%             % Propeller disk
%             plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
%             % Engine and Nacelle second pair
%             % Engine
%             plot_2engine2(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
%             % Nacelle
%             plot_2nacelle2(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
%             % Propeller disk
%             plot_2prop2(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop);
%         end
% end