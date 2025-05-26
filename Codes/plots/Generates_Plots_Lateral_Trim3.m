function [Fig] = Generates_Plots_Lateral_Trim3(Trim_ITER_LAT3,...
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
 
% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filenamec = 'Figs';
% filename_DATA = strcat(filename,filenamed);
% filename_Plots = strcat(filename,filenamec);
fname = filenameS.plots;

% Datos para aceleraciones
V_vec = Trim_ITER_LAT3.V_vec;
deltaa_vec = Trim_ITER_LAT3.deltaa_vec;
deltar_vec = Trim_ITER_LAT3.deltar_vec;
pdot_vec = Trim_ITER_LAT3.pdot_vec;
rdot_vec = Trim_ITER_LAT3.rdot_vec;

%% Plots for result of variable trim conditions
pdot_vec_gs = pdot_vec;
Fig = Fig + 1;
figure(Fig)
surf(deltaa_vec*R2D,V_vec,pdot_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('dp/dt (rads/s^2) vs V and aileron deflection','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\delta_{a} (deg)','FontSize',FS)
zlabel('dp/dt  (rads/s^2)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('pdot_vs_V&deltaa');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_pdot_plot = linspace(min(min(pdot_vec_gs)),max(max(pdot_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(deltaa_vec*R2D,V_vec,pdot_vec_gs,vect_cc_pdot_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('dp/dt (rads/s^2) vs V and aileron deflection','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\delta_{a} (deg)','FontSize',FS)
zlabel('dp/dt (rads/s^2)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('pdot_vs_V&deltaa_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

rdot_vec_gs = rdot_vec;

%% Plots for result of variable trim conditions
Fig = Fig + 1;
figure(Fig)
surf(deltar_vec*R2D,V_vec,rdot_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
title('dr/dt (rads/s^2) vs V and rudder deflection','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\delta_{r} (deg)','FontSize',FS)
zlabel('dr/dt (rads/s^2)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('rdot_vs_V&deltar');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_rdot_plot = linspace(min(min(rdot_vec_gs)),max(max(rdot_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(deltar_vec*R2D,V_vec,rdot_vec_gs,vect_cc_rdot_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('dr/dt (rads/s^2) vs V and rudder deflection','FontSize',FS)
ylabel('Velocity (m/s)','FontSize',FS)
xlabel('\delta_{r} (deg)','FontSize',FS)
zlabel('dr/dt (rads/s^2)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('rdot_vs_V&deltar_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            SAVE_types(fname,name,gca,gcf); 
end
