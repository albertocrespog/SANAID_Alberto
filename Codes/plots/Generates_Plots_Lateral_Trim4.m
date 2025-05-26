function [Fig] = Generates_Plots_Lateral_Trim4(Trim_ITER_LAT4,...
    Geo_tier,Plot_Options,conv_UNITS,conditions_TRIM_lat,OUTPUT_read_XLSX,Fig)

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
 
deltaa_deg = Trim_ITER_LAT4.deltaa_deg;
deltar_deg = Trim_ITER_LAT4.deltar_deg;
deltaaT_deg = Trim_ITER_LAT4.deltaaT_deg; 
beta_deg = Trim_ITER_LAT4.beta_deg;

phi_vec = conditions_TRIM_lat.phi_vec;
phi = conditions_TRIM_lat.phi;

% beta_vec = conditions_TRIM_lat.beta_vec;
% beta = conditions_TRIM_lat.beta;

deltaa_deg_var = Trim_ITER_LAT4.deltaa_deg_var;
deltar_deg_var = Trim_ITER_LAT4.deltar_deg_var;
deltaaT_deg_var = Trim_ITER_LAT4.deltaaT_deg_var;
beta_deg_var = Trim_ITER_LAT4.beta_deg_var;

Fig = Fig + 1;
figure(Fig)
plot(phi_vec*R2D,deltaa_deg_var,'b','LineWidth', LS)
hold on
plot(phi_vec*R2D,deltar_deg_var,'r','LineWidth', LS)
plot(phi_vec*R2D,deltaaT_deg_var,'m','LineWidth', LS)
plot(phi_vec*R2D,beta_deg_var,'g','LineWidth', LS)
plot(phi*R2D,deltaa_deg,'bo','LineWidth', LS*1.5)
plot(phi*R2D,deltar_deg,'ro','LineWidth', LS*1.5)
plot(phi*R2D,deltaaT_deg,'mo','LineWidth', LS*1.5)
plot(phi*R2D,beta_deg,'go','LineWidth', LS*1.5)
hold off
title('Lateral Trim Conditions - Sidewind and Drag Asymmetry')
xlabel('\phi (deg)')
ylabel('trims angles (deg)')

if MATLAB_in == 1
    legend('\delta_a','\delta_r','\delta_{rT}','\beta','\delta_{a-TO}','\delta_{r-TO}','\delta_{rT-TO}','\beta_{TO}')
else
    h_legend=legend('\delta_a','\delta_r','\delta_{rT}','\beta','\delta_{a-TO}','\delta_{r-TO}','\delta_{rT-TO}','\beta_{TO}');
    set(h_legend, 'Location','Best','FontSize',LFS)
end
    
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('da_dr_drT_beta&YawingAsymmetry');
    name   = strcat(prefix,st);
    fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

% 3D Figures for variation phi and V
V_vec = conditions_TRIM_lat.V_VAR;
phi_vec = conditions_TRIM_lat.phi_vec;
deltaa_deg = Trim_ITER_LAT4.deltaa_deg_var2;
deltar_deg = Trim_ITER_LAT4.deltar_deg_var2;
deltaaT_deg = Trim_ITER_LAT4.deltaaT_deg_var2;
beta_deg = Trim_ITER_LAT4.beta_deg_var2;

%% Plots for result of variable trim conditions
deltaa_deg_vec_gs = deltaa_deg;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,phi_vec*R2D,deltaa_deg_vec_gs)
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
    st = strcat('deltaa_vs_V&phi');
    name   = strcat(prefix,st);
    fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_deltaa_plot = linspace(min(min(deltaa_deg_vec_gs)),max(max(deltaa_deg_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,phi_vec*R2D,deltaa_deg_vec_gs,vect_cc_deltaa_plot');
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
    st = strcat('deltaa_vs_V&beta_contour');
    name   = strcat(prefix,st);
    fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end


%% Plots for result of variable trim conditions
deltar_deg_vec_gs = deltar_deg;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,phi_vec*R2D,deltar_deg_vec_gs)
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
    st = strcat('deltar_vs_V&phi');
    name   = strcat(prefix,st);
    fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_deltar_plot = linspace(min(min(deltar_deg_vec_gs)),max(max(deltar_deg_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,phi_vec*R2D,deltar_deg_vec_gs,vect_cc_deltar_plot');
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
    st = strcat('deltar_vs_V&phi_contour');
    name   = strcat(prefix,st);
    fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end




%% Plots for result of variable trim conditions
deltaaT_deg_vec_gs = deltaaT_deg;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,phi_vec*R2D,deltaaT_deg_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\delta_{aT} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{aT}  (deg)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltaaT_vs_V&phi');
    name   = strcat(prefix,st);
    fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_deltaaT_plot = linspace(min(min(deltaaT_deg_vec_gs)),max(max(deltaaT_deg_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,phi_vec*R2D,deltaaT_deg_vec_gs,vect_cc_deltaaT_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('\delta_{aT} vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\delta_{aT}  (deg)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltaaT_vs_V&phi_contour');
    name   = strcat(prefix,st);
    fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end



%% Plots for result of variable trim conditions
beta_deg_vec_gs = beta_deg;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,phi_vec*R2D,beta_deg_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\beta vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\beta  (deg)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('beta_vs_V&phi');
    name   = strcat(prefix,st);
    fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_beta_plot = linspace(min(min(beta_deg_vec_gs)),max(max(beta_deg_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,phi_vec*R2D,beta_deg_vec_gs,vect_cc_beta_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('\beta vs V and \phi ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\phi (deg)','FontSize',FS)
zlabel('\beta  (deg)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('beta_vs_V&phi_contour');
    name   = strcat(prefix,st);
    fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end
