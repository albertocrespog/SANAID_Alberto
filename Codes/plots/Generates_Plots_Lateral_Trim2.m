function [Fig] = Generates_Plots_Lateral_Trim2(Trim_ITER_LAT2,...
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

deltaa_deg = Trim_ITER_LAT2.deltaa_deg;
deltar_deg = Trim_ITER_LAT2.deltar_deg;
phi_deg = Trim_ITER_LAT2.phi_deg;
beta_vec = conditions_TRIM_lat.beta_vec;
beta = conditions_TRIM_lat.beta;

deltaa_deg_var = Trim_ITER_LAT2.deltaa_deg_var;
deltar_deg_var = Trim_ITER_LAT2.deltar_deg_var;
phi_deg_var = Trim_ITER_LAT2.phi_deg_var;

Fig = Fig + 1;
figure(Fig)
plot(beta_vec*R2D,deltaa_deg_var,'b','LineWidth', LS)
hold on
plot(beta_vec*R2D,deltar_deg_var,'r','LineWidth', LS)
plot(beta_vec*R2D,phi_deg_var,'g','LineWidth', LS)
plot(beta*R2D,deltaa_deg,'bo','LineWidth', LS*1.5)
plot(beta*R2D,deltar_deg,'ro','LineWidth', LS*1.5)
plot(beta*R2D,phi_deg,'go','LineWidth', LS*1.5)
hold off
title('Lateral Trim Conditions - Sidewind and Drag Asymmetry')
xlabel('\beta (deg)')
ylabel('trims angles (deg)')

if MATLAB_in == 1
    legend('\delta_a','\delta_r','\phi','\delta_{a-TO}','\delta_{r-TO}','\phi_{TO}')
else
    h_legend=legend('\delta_a','\delta_r','\phi','\delta_{a-TO}','\delta_{r-TO}','\phi_{TO}');
    set(h_legend, 'Location','Best','FontSize',LFS)
end
    
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('da_dr_phi&YawingAsymmetry');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

% 3D Figures for variation Beta and V
V_vec = conditions_TRIM_lat.V_VAR;
beta_vec = conditions_TRIM_lat.beta_vec;
deltaa_deg = Trim_ITER_LAT2.deltaa_deg_var3;
deltar_deg = Trim_ITER_LAT2.deltar_deg_var3;
phi_deg = Trim_ITER_LAT2.phi_deg_var3;

%% Plots for result of variable trim conditions
deltaa_deg_vec_gs = deltaa_deg;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,beta_vec*R2D,deltaa_deg_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\delta_{a} vs V and \beta ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\beta$ (deg)','FontSize',FS)
zlabel('\delta_{a}  (deg)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltaa_vs_V&beta');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_deltaa_plot = linspace(min(min(deltaa_deg_vec_gs)),max(max(deltaa_deg_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,beta_vec*R2D,deltaa_deg_vec_gs,vect_cc_deltaa_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('\delta_{a} vs V and \beta ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\beta$ (deg)','FontSize',FS)
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
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end


%% Plots for result of variable trim conditions
deltar_deg_vec_gs = deltar_deg;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,beta_vec*R2D,deltar_deg_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\delta_{r} vs V and \beta ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\beta (deg)','FontSize',FS)
zlabel('\delta_{a}  (deg)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltar_vs_V&beta');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_deltar_plot = linspace(min(min(deltar_deg_vec_gs)),max(max(deltar_deg_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,beta_vec*R2D,deltar_deg_vec_gs,vect_cc_deltar_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('\delta_{r} vs V and \beta ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\beta (deg)','FontSize',FS)
zlabel('\delta_{r}  (deg)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltar_vs_V&beta_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end


%% Plots for result of variable trim conditions
phi_deg_vec_gs = phi_deg;
Fig = Fig + 1;
figure(Fig)
surf(V_vec,beta_vec*R2D,phi_deg_vec_gs)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\phi vs V and \beta ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\beta (deg)','FontSize',FS)
zlabel('phi  (deg)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('phi_vs_V&beta');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end

Fig = Fig + 1;
figure(Fig)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_phi_plot = linspace(min(min(phi_deg_vec_gs)),max(max(phi_deg_vec_gs)),N_contour_lines);
[pdot_c,h_pdot_c] = contourf(V_vec,beta_vec*R2D,phi_deg_vec_gs,vect_cc_phi_plot');
clabel(pdot_c,h_pdot_c)
hold on
shading interp
grid on
title('\phi vs V and \beta ','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('\beta (deg)','FontSize',FS)
zlabel('\phi  (deg)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('phi_vs_V&beta_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
end


