function PLOTS_Real_AC(OUTPUT_read_XLSX,PLOTTING_UAV,COLOR_scheme,Body_Geo)


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

% Calculates the Aircraft Geometry
model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_ac);

ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;
SF_CAD_AC = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD_AC;

if OUTPUT_read_XLSX.Fuselage_flags.STL_reading_test == 1
    max_x_ac = max(max(model{1}));
    min_x_ac = min(min(model{1}));
    max_y_ac = max(max(model{2}));
    min_y_ac = min(min(model{2}));
    max_z_ac = max(max(model{3}));
    min_z_ac = min(min(model{3}));

    Delta_X_stl = min_x_ac;
    Delta_Y_stl = -(max_y_ac - (max_y_ac - min_y_ac)/2)/2;
    Delta_Z_stl = 0;
    dummy=0;

else
    [points_ac,triangles,tri_norms] = import_stl_fast_v2(OUTPUT_read_XLSX.Fuselage_flags.STL_ac,1);
    max_x_ac = max(points_ac(:,1));
    min_x_ac = min(points_ac(:,1));
    max_y_ac = max(points_ac(:,2));
    min_y_ac = min(points_ac(:,2));
    max_z_ac = max(points_ac(:,3));
    min_z_ac = min(points_ac(:,3));

    Delta_X_stl = points_ac(1,1);
    Delta_Y_stl = points_ac(1,2);
    Delta_Z_stl = points_ac(1,3);

    % Calculate stats (min, max, and deltas for each axis)
    points = points_ac;
    min_vals = min(points);
    max_vals = max(points);
    deltas = points - min(points);

%     [points, triangles, tri_norms, stats] = import_stl_fast_v3(OUTPUT_read_XLSX.Fuselage_flags.STL_ac,1);
% 
% % Acceder a los resultados
% disp('Valores mínimos:');
% disp(stats.min);
% 
% disp('Valores máximos:');
% disp(stats.max);
% 
% disp('Desviaciones (Delta):');
% disp(stats.deltas);
    dummy=0;
end
% 
% plotData{1} = (model{1}*ScaleFactor - Delta_X_stl).*SF_CAD;
% plotData{2} = (model{2}*ScaleFactor - Delta_Y_stl).*SF_CAD;
% plotData{3} = (model{3}*ScaleFactor - Delta_Z_stl).*SF_CAD;

% 
% Delta_Y_stl = Body_Geo.Delta_Y_stl;
% Delta_Z_stl = Body_Geo.Delta_Z_stl;
% Delta_X_stl = 0;

% Generates the Plots 
plotData{1} = (model{1}*ScaleFactor).*SF_CAD_AC;
plotData{2} = (model{2}*ScaleFactor).*SF_CAD_AC;
plotData{3} = (model{3}*ScaleFactor).*SF_CAD_AC;
meshData_ac = plotData;

% Colors
C_ac(:,:,1) = color_ac(1)*ones(size(meshData_ac{1}));
C_ac(:,:,2) = color_ac(2)*ones(size(meshData_ac{1}));
C_ac(:,:,3) = color_ac(3)*ones(size(meshData_ac{1}));
% mesh(meshData_ac{1},meshData_ac{2},meshData_ac{3},C_ac);
% 
% 
% plotData{1} = (model{1}*ScaleFactor).*SF_CAD_AC - Delta_X_stl;
% plotData{2} = (model{2}*ScaleFactor).*SF_CAD_AC - Delta_Y_stl;
% plotData{3} = (model{3}*ScaleFactor).*SF_CAD_AC - Delta_Z_stl;
% meshData_ac = plotData;
% 
% % Colors
% C_ac(:,:,1) = 0.6*color_ac(1)*ones(size(meshData_ac{1}));
% C_ac(:,:,2) = 0.8*color_ac(2)*ones(size(meshData_ac{1}));
% C_ac(:,:,3) = 0.9*color_ac(3)*ones(size(meshData_ac{1}));
% mesh(meshData_ac{1},meshData_ac{2},meshData_ac{3},C_ac);
% pause

% figure(Fig)
% Fig = Fig + 1;
%% PLOTS Sequence
mesh(meshData_ac{1},meshData_ac{2},meshData_ac{3},C_ac);

% case_AC = OUTPUT_read_XLSX.AC_Data_flags.case_AC;
% st = strcat(case_AC);
% title(st,'fontsize',FS)
% title('Aircrat Original CAD')
% xlabel('y (m)')
% ylabel('z (m)')
% zlabel('z (m)')