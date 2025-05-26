function [Fig] = Generates_Plots_Lateral_Trim4D(Trim_ITER_LAT4D,...
    Geo_tier,Plot_Options,conv_UNITS,conditions_TRIM_lat,OUTPUT_read_XLSX,Fig,filenameS)

g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;

% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;
MATLAB_in = Plot_Options.MATLAB_in;

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;
 
% Defines the Target for the plots
fname = filenameS.filename_Plots;

% deltaa_deg = Trim_ITER_LAT4D.deltaa_deg;
% deltar_deg = Trim_ITER_LAT4D.deltar_deg;
% daT_deg = Trim_ITER_LAT4D.daT_deg;
% drT_deg = Trim_ITER_LAT4D.drT_deg;
% Fa = Trim_ITER_LAT4D.Fa; 
% Fr = Trim_ITER_LAT4D.Fr; 
% beta_deg = Trim_ITER_LAT4D.beta_deg;
% phi_vec = conditions_TRIM_lat.phi_vec;
% phi = conditions_TRIM_lat.phi;

% 3D Figures for variation phi and V
V_vec = conditions_TRIM_lat.V_VAR;
phi_vec = conditions_TRIM_lat.phi_vec;
% deltaaT_vec = conditions_TRIM_lat.deltaaT_vec;

% deltaa_deg_var1 = Trim_ITER_LAT4D.deltaa_deg_var1;
% deltar_deg_var1 = Trim_ITER_LAT4D.deltar_deg_var1;
% daT_deg_var1 = Trim_ITER_LAT4D.daT_deg_var1;
% drT_deg_var1 = Trim_ITER_LAT4D.drT_deg_var1;
% beta_deg_var1 = Trim_ITER_LAT4D.beta_deg_var1;

% Variable phi and V
deltaa_deg2 = Trim_ITER_LAT4D.deltaa_deg_var2;
deltar_deg2 = Trim_ITER_LAT4D.deltar_deg_var2;
daT_deg_var2 = Trim_ITER_LAT4D.daT_deg_var2;
drT_deg_var2 = Trim_ITER_LAT4D.drT_deg_var2;
beta_deg2 = Trim_ITER_LAT4D.beta_deg_var2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Variable phi and V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plots for result of variable trim conditions
deltaa_deg2_vec_gs = deltaa_deg2;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,phi_vec*R2D,deltaa_deg2_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\delta_{a} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{a}  (deg)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltaa_vs_V&phi_5eqs_sol_zeroFAFR');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_deltaa_plot = linspace(min(min(deltaa_deg2_vec_gs)),max(max(deltaa_deg2_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,phi_vec*R2D,deltaa_deg2_vec_gs,vect_cc_deltaa_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('\delta_{a} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{a}  (deg)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltaa_vs_V&phi_contour_5eqs_sol_zeroFAFR');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end


%% Plots for result of variable trim conditions
deltar_deg2_vec_gs = deltar_deg2;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,phi_vec*R2D,deltar_deg2_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\delta_{r} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{r}  (deg)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltar_vs_V&phi_5eqs_sol_zeroFAFR');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_deltar_plot = linspace(min(min(deltar_deg2_vec_gs)),max(max(deltar_deg2_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,phi_vec*R2D,deltar_deg2_vec_gs,vect_cc_deltar_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('\delta_{r} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{r}  (deg)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltar_vs_V&phi_contour_5eqs_sol_zeroFAFR');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end


%% Plots for result of variable trim conditions
daT_deg_var2_deg_vec_gs = daT_deg_var2;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,phi_vec*R2D,daT_deg_var2_deg_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\delta_{aT} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{aT}  (N)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('daT_vs_V&phi_5eqs_sol_zeroFAFR');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_daT2_plot = linspace(min(min(daT_deg_var2_deg_vec_gs)),max(max(daT_deg_var2_deg_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,phi_vec*R2D,daT_deg_var2_deg_vec_gs,vect_cc_daT2_plot');
clabel(pdot_c,h_pdot_c)
hold on

% Add specific contour line for constant value (e.g., const_value)
const_value = -8; % Define the constant value you want to highlight
[C, h_const] = contour(V_vec, phi_vec*R2D, daT_deg_var2_deg_vec_gs, [const_value const_value], 'LineColor', 'r', 'LineWidth', 2);
% Enhance the plot appearance
shading interp;
grid on;
% Optionally, add a label to the constant contour line
clabel(C, h_const, 'FontSize', 12, 'Color', 'r');

title('\delta_{aT} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{aT}','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('daT_vs_V&phi_contour_5eqs_sol_zeroFAFR');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end




% Plots for result of variable trim conditions
drT_deg_var2_deg_vec_gs = drT_deg_var2;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,phi_vec*R2D,drT_deg_var2_deg_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\delta_{rT} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{rT}  (N)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('drT_vs_V&phi_5eqs_sol_zeroFAFR');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_drT2_plot = linspace(min(min(drT_deg_var2_deg_vec_gs)),max(max(drT_deg_var2_deg_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,phi_vec*R2D,drT_deg_var2_deg_vec_gs,vect_cc_drT2_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('\delta_{rT} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{rT}','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('drT_vs_V&phi_contour_5eqs_sol_zeroFAFR');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end




