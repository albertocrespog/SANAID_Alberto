%%%%%% PRUEBA PLOTS %%%%%%
function [Fig] = Generates_Plots_Performance_v5(Storing_PERFORMANCE_DATA,Plot_Options,OUTPUT_read_XLSX,Restrictions,Fig,filenameS)
%close all

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;
% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;
MATLAB_in = Plot_Options.MATLAB_in;

% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filenamec = 'Figs';
% filename_DATA = strcat(filename,filenamed);
% filename_Plots = strcat(filename,filenamec);
fname = filenameS.plots;

% Plot options
Plot_OptionsB = Storing_PERFORMANCE_DATA.Plot_Options;
h_VAR = Plot_OptionsB.h_VAR;
V_VAR = Plot_OptionsB.V_VAR;
N_h_VAR_perfo = Plot_OptionsB.N_h_VAR_perfo;
N_V_VAR_perfo = Plot_OptionsB.N_V_VAR_perfo;

Restrictions = Storing_PERFORMANCE_DATA.Restrictions_var_h_V;
% V_stall = Restrictions.V_stall;
% V_TO = Restrictions.V_TO;
% V_min_ope =  Restrictions.V_min_ope;
for i=1:length(V_VAR)
    for j=1:length(h_VAR)
        
        Time_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_24{i,j}.Total_Datos.tiempo_total/60;
        Distance_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_24{i,j}.Total_Datos.distancia_total/1000;
        Wempty_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_24{i,j}.Weights_AP.m_TOW;
        W_final_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_24{i,j}.Weights_AP.m_F;

        % energia_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_4{i,j}.datos(1).segmento.energia;
        % Power_PLOTS(i,j) = energia_PLOTS(i,j)./Time_PLOTS(i,j);

        % palanca_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_4{i,j}.datos(1).segmento.palanca;
        % RPM_max(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_4{i,j}.Storing_PROPULSION_DATA.Prop_data.RPM_max;
        % RPM_PLOTS(i,j) = palanca_PLOTS(i,j).*RPM_max(i,j);

        % Thrust_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_4{i,j}.datos(1).segmento.empuje;
        % L_D_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_4{i,j}.datos(1).segmento.L_D;
        % % etamp_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_4{i,j}.datos(1).segmento.etamp;
        % CL_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_4{i,j}.datos(1).segmento.CL;
        % CD_PLOTS(i,j) = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_4{i,j}.datos(1).segmento.CD;

        % Power_PLOTS(i,j) = Thrust_PLOTS(i,j).*V_VAR(j)./etamp_PLOTS(i,j);

        h_VAR_MAT(i,j) = h_VAR(j);
        V_VAR_MAT(i,j) = V_VAR(i);
        % V_stall(j,i) = Restrictions_var_V_m{i,j}.V_stall;
        % V_TO(j,i) = Restrictions_var_V_m{i,j}.V_TO;
        % V_min_ope(j,i) = Restrictions_var_V_m{i,j}.V_min_ope;

        V_stall(i,j) = Restrictions{i,j}.V_stall;
        V_TO(i,j) = Restrictions{i,j}.V_TO;
        V_min_ope(i,j) = Restrictions{i,j}.V_min_ope;
        % D_prop = Restrictions{i,j}.D_prop*100/2.54;
    end
end
       
%% Plots for Mission Time
Fig = Fig + 1;
figure(Fig)
% subplot (1,2,1)
[x1 y1] = contour(V_VAR_MAT,h_VAR_MAT,Time_PLOTS,[0 0]); % x1 will be your "two-column matrix"
V_alpha_cero = x1(2,2:end);
M_alpha_cero = x1(1,2:end);
surf(h_VAR,V_VAR,Time_PLOTS)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Mission Time (min) vs V and Altitude','FontSize',FS)
xlabel('Altitude (m)','FontSize',FS)
ylabel('Velocity (m/s)','FontSize',FS)
zlabel('Time Mission (min)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('time_mission_vs_h&V');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end

%% Plot contours
Fig = Fig + 1;
figure(Fig)
N_contour_lines = 20; % Number of contour lines
vect_cc_time_mission_plot = linspace(min(min(Time_PLOTS')),max(max(Time_PLOTS')),N_contour_lines);
% [alpha_c,h_alpha_c] = contourf(V_VAR,m_bat_VAR,Time_PLOTS,vect_cc_time_mission_plot');
[alpha_c,h_alpha_c] = contourf(h_VAR,V_VAR,Time_PLOTS,vect_cc_time_mission_plot');
clabel(alpha_c,h_alpha_c)
hold on
plot(h_VAR,V_stall(1,:),'m--','LineWidth',LS)
plot(h_VAR,V_TO(1,:),'o-.','LineWidth',LS)
plot(h_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
legend('t_{mission}(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest')
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
grid on
Title = strcat('Mission Time (min) vs V and Altitude');
title(Title,'FontSize',FS)
xlabel('Altitude (m)','FontSize',FS)
ylabel('Velocity (m/s)','FontSize',FS)
zlabel('Time Mission (min)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('time_mission_vs_h&V_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
end

%% Plots for thresult of variable trim conditions
Fig = Fig + 1;
figure(Fig)
% subplot (1,2,1)
[x1 y1] = contour(V_VAR_MAT,h_VAR_MAT,Distance_PLOTS,[0 0]); % x1 will be your "two-column matrix"
V_alpha_cero = x1(2,2:end);
M_alpha_cero = x1(1,2:end);
surf(h_VAR,V_VAR,Distance_PLOTS)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Distance Mission (km) vs V and Altitude','FontSize',FS)
xlabel('Altitude (m)','FontSize',FS)
ylabel('Velocity (m/s)','FontSize',FS)
zlabel('Distance Mission (km)','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('x_mission_vs_h&V');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end

%% Plot contours
Fig = Fig + 1;
figure(Fig)
N_contour_lines = 20; % Number of contour lines
vect_cc_time_mission_plot = linspace(min(min(Distance_PLOTS')),max(max(Distance_PLOTS')),N_contour_lines);
[alpha_c,h_alpha_c] = contourf(h_VAR,V_VAR,Distance_PLOTS,vect_cc_time_mission_plot');
clabel(alpha_c,h_alpha_c)
hold on
plot(h_VAR,V_stall(1,:),'m--','LineWidth',LS)
plot(h_VAR,V_TO(1,:),'o-.','LineWidth',LS)
plot(h_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
legend('x_{mission}(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest')
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
grid on
title('Distance Mission (km) vs V and Altitude','FontSize',FS)
xlabel('Altitude (m)','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
zlabel('Distance Mission (km)','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('x_mission_vs_h&V_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
end

% %% Plots for Electric Power
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,1)
% [x1 y1] = contour(V_VAR_MAT,m_VAR_MAT,Power_PLOTS/1000,[0 0]); % x1 will be your "two-column matrix"
% V_alpha_cero = x1(2,2:end);
% M_alpha_cero = x1(1,2:end);
% surf(V_VAR,m_bat_VAR,Power_PLOTS/1000)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
% %    rotate3d on % allows you to rotate the graph with the mouse
% colorbar % adds the color bar to the right of the graph
% title('Electric Power vs V and battery mass','FontSize',FS)
% ylabel('batery mass (kg)','FontSize',FS)
% xlabel('Velocity (m/s)','FontSize',FS)
% zlabel('P_{electric} (KW)','FontSize',FS)
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('Pe_vs_V&mbat');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     % Unified plot Options
%     SAVE_types(fname,name,gca,gcf);
% end
% 
% %% Plot contours
% Fig = Fig + 1;
% figure(Fig)
% N_contour_lines = 20; % Number of contour lines
% vect_cc_time_mission_plot = linspace(min(min(Power_PLOTS'/1000)),max(max(Power_PLOTS'/1000)),N_contour_lines);
% [alpha_c,h_alpha_c] = contourf(m_bat_VAR,V_VAR,Power_PLOTS'/1000,vect_cc_time_mission_plot');
% clabel(alpha_c,h_alpha_c)
% hold on
% plot(m_bat_VAR,V_stall(:,1),'m--','LineWidth',LS)
% plot(m_bat_VAR,V_TO(:,1),'o-.','LineWidth',LS)
% plot(m_bat_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% legend('P_{e-mission}(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest')
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% grid on
% title('Electric Power vs V and battery mass','FontSize',FS)
% xlabel('batery mass (kg)','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% zlabel('P_{electric} (KW)','FontSize',FS)
% hold off
% sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% % caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('Pe_vs_V&mbat');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%             % Unified plot Options             
%             SAVE_types(fname,name,gca,gcf); 
% %             saveas(gcf,fullfile(fname, name),'fig');
% %             saveas(gcf,fullfile(fname, name),'pdf');
% %             saveas(gcf,fullfile(fname, name),'bmp');
% %             saveas(gcf,fullfile(fname, name),'png');
% end
% 
% %% Plots for Energy Mission
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,1)
% [x1 y1] = contour(V_VAR_MAT,m_VAR_MAT,energia_PLOTS/1000,[0 0]); % x1 will be your "two-column matrix"
% V_alpha_cero = x1(2,2:end);
% M_alpha_cero = x1(1,2:end);
% surf(V_VAR,m_bat_VAR,energia_PLOTS/1000)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
% %    rotate3d on % allows you to rotate the graph with the mouse
% colorbar % adds the color bar to the right of the graph
% title('Energy Consummed vs V and battery mass','FontSize',FS)
% ylabel('batery mass (kg)','FontSize',FS)
% xlabel('Velocity (m/s)','FontSize',FS)
% zlabel('Energy Consummed (KW)','FontSize',FS)
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('Energy_mission_vs_V&mbat');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     % Unified plot Options
%     SAVE_types(fname,name,gca,gcf);
% end
% 
% %% Plot contours
% Fig = Fig + 1;
% figure(Fig)
% N_contour_lines = 20; % Number of contour lines
% vect_cc_time_mission_plot = linspace(min(min(energia_PLOTS'/1000)),max(max(energia_PLOTS'/1000)),N_contour_lines);
% [alpha_c,h_alpha_c] = contourf(m_bat_VAR,V_VAR,energia_PLOTS'/1000,vect_cc_time_mission_plot');
% clabel(alpha_c,h_alpha_c)
% hold on
% plot(m_bat_VAR,V_stall(:,1),'m--','LineWidth',LS)
% plot(m_bat_VAR,V_TO(:,1),'o-.','LineWidth',LS)
% plot(m_bat_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% legend('E_{mission}(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest')
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% grid on
% title('Energy Consummed vs V and battery mass','FontSize',FS)
% xlabel('batery mass (kg)','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% zlabel('Energy Consummed (KW)','FontSize',FS)
% hold off
% sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% % caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('Energy_mission_vs_V&mbat_contour');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%             % Unified plot Options             
%             SAVE_types(fname,name,gca,gcf); 
% %             saveas(gcf,fullfile(fname, name),'fig');
% %             saveas(gcf,fullfile(fname, name),'pdf');
% %             saveas(gcf,fullfile(fname, name),'bmp');
% %             saveas(gcf,fullfile(fname, name),'png');
% end
% 
% %% Plots for RPMs
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,1)
% [x1 y1] = contour(V_VAR_MAT,m_VAR_MAT,RPM_PLOTS,[0 0]); % x1 will be your "two-column matrix"
% V_alpha_cero = x1(2,2:end);
% M_alpha_cero = x1(1,2:end);
% surf(V_VAR,m_bat_VAR,RPM_PLOTS)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
% %    rotate3d on % allows you to rotate the graph with the mouse
% colorbar % adds the color bar to the right of the graph
% title('RPM vs V and battery mass','FontSize',FS)
% ylabel('batery mass (kg)','FontSize',FS)
% xlabel('Velocity (m/s)','FontSize',FS)
% zlabel('RPM (-)','FontSize',FS)
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('RPM_mission_vs_V&mbat');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     % Unified plot Options
%     SAVE_types(fname,name,gca,gcf);
% end
% 
% %% Plot contours
% Fig = Fig + 1;
% figure(Fig)
% N_contour_lines = 20; % Number of contour lines
% vect_cc_time_mission_plot = linspace(min(min(RPM_PLOTS')),max(max(RPM_PLOTS')),N_contour_lines);
% [alpha_c,h_alpha_c] = contourf(m_bat_VAR,V_VAR,RPM_PLOTS',vect_cc_time_mission_plot');
% clabel(alpha_c,h_alpha_c)
% hold on
% plot(m_bat_VAR,V_stall(:,1),'m--','LineWidth',LS)
% plot(m_bat_VAR,V_TO(:,1),'o-.','LineWidth',LS)
% plot(m_bat_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% legend('RPM_{mission}(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest')
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% grid on
% title('RPM vs V and battery mass','FontSize',FS)
% xlabel('batery mass (kg)','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% zlabel('RPM (-)','FontSize',FS)
% hold off
% sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% % caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('RP_mission_vs_V&mbat_contour');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%             % Unified plot Options             
%             SAVE_types(fname,name,gca,gcf); 
% %             saveas(gcf,fullfile(fname, name),'fig');
% %             saveas(gcf,fullfile(fname, name),'pdf');
% %             saveas(gcf,fullfile(fname, name),'bmp');
% %             saveas(gcf,fullfile(fname, name),'png');
% end
% 
% %% Plots for Throttle position
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,1)
% [x1 y1] = contour(V_VAR_MAT,m_VAR_MAT,palanca_PLOTS,[0 0]); % x1 will be your "two-column matrix"
% V_alpha_cero = x1(2,2:end);
% M_alpha_cero = x1(1,2:end);
% surf(V_VAR,m_bat_VAR,palanca_PLOTS)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
% %    rotate3d on % allows you to rotate the graph with the mouse
% colorbar % adds the color bar to the right of the graph
% title('Throttle vs V and battery mass','FontSize',FS)
% ylabel('batery mass (kg)','FontSize',FS)
% xlabel('Velocity (m/s)','FontSize',FS)
% zlabel('Throttle (-)','FontSize',FS)
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('deltaT_mission_vs_V&mbat');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     % Unified plot Options
%     SAVE_types(fname,name,gca,gcf);
% end
% 
% %% Plot contours
% Fig = Fig + 1;
% figure(Fig)
% N_contour_lines = 20; % Number of contour lines
% vect_cc_time_mission_plot = linspace(min(min(palanca_PLOTS')),max(max(palanca_PLOTS')),N_contour_lines);
% [alpha_c,h_alpha_c] = contourf(m_bat_VAR,V_VAR,palanca_PLOTS',vect_cc_time_mission_plot');
% clabel(alpha_c,h_alpha_c)
% hold on
% plot(m_bat_VAR,V_stall(:,1),'m--','LineWidth',LS)
% plot(m_bat_VAR,V_TO(:,1),'o-.','LineWidth',LS)
% plot(m_bat_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% legend('\delta_{T_{mission}}(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest')
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% grid on
% title('Throttle vs V and battery mass','FontSize',FS)
% xlabel('batery mass (kg)','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% zlabel('Throttle (-)','FontSize',FS)
% hold off
% sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% % caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('deltaT_mission_vs_V&mbat_contour');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%             % Unified plot Options             
%             SAVE_types(fname,name,gca,gcf); 
% %             saveas(gcf,fullfile(fname, name),'fig');
% %             saveas(gcf,fullfile(fname, name),'pdf');
% %             saveas(gcf,fullfile(fname, name),'bmp');
% %             saveas(gcf,fullfile(fname, name),'png');
% end
% 
% 
% 
% %% Plots for moto-propulsife efficiency
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,1)
% [x1 y1] = contour(V_VAR_MAT,m_VAR_MAT,etamp_PLOTS,[0 0]); % x1 will be your "two-column matrix"
% V_alpha_cero = x1(2,2:end);
% M_alpha_cero = x1(1,2:end);
% surf(V_VAR,m_bat_VAR,etamp_PLOTS)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
% %    rotate3d on % allows you to rotate the graph with the mouse
% colorbar % adds the color bar to the right of the graph
% title('Propulsive Efficiency vs V and battery mass','FontSize',FS)
% ylabel('batery mass (kg)','FontSize',FS)
% xlabel('Velocity (m/s)','FontSize',FS)
% zlabel('Propulsive Efficiency (-)','FontSize',FS)
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('etaMP_mission_vs_V&mbat');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     % Unified plot Options
%     SAVE_types(fname,name,gca,gcf);
% end
% 
% %% Plot contours
% Fig = Fig + 1;
% figure(Fig)
% N_contour_lines = 20; % Number of contour lines
% vect_cc_time_mission_plot = linspace(min(min(etamp_PLOTS')),max(max(etamp_PLOTS')),N_contour_lines);
% [alpha_c,h_alpha_c] = contourf(m_bat_VAR,V_VAR,etamp_PLOTS',vect_cc_time_mission_plot');
% clabel(alpha_c,h_alpha_c)
% hold on
% plot(m_bat_VAR,V_stall(:,1),'m--','LineWidth',LS)
% plot(m_bat_VAR,V_TO(:,1),'o-.','LineWidth',LS)
% plot(m_bat_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% legend('x_{mission}(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest')
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% grid on
% title('Propulsive Efficiency vs V and battery mass','FontSize',FS)
% xlabel('batery mass (kg)','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% zlabel('Propulsive Efficiency (-)','FontSize',FS)
% hold off
% sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% % caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('etaMP_mission_vs_V&mbat_contour');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%             % Unified plot Options             
%             SAVE_types(fname,name,gca,gcf); 
% %             saveas(gcf,fullfile(fname, name),'fig');
% %             saveas(gcf,fullfile(fname, name),'pdf');
% %             saveas(gcf,fullfile(fname, name),'bmp');
% %             saveas(gcf,fullfile(fname, name),'png');
% end
% 
% %% Plots for Mission L/D
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,1)
% [x1 y1] = contour(V_VAR_MAT,m_VAR_MAT,L_D_PLOTS,[0 0]); % x1 will be your "two-column matrix"
% V_alpha_cero = x1(2,2:end);
% M_alpha_cero = x1(1,2:end);
% surf(V_VAR,m_bat_VAR,L_D_PLOTS)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
% %    rotate3d on % allows you to rotate the graph with the mouse
% colorbar % adds the color bar to the right of the graph
% title('L/D vs V and battery mass','FontSize',FS)
% ylabel('batery mass (kg)','FontSize',FS)
% xlabel('Velocity (m/s)','FontSize',FS)
% zlabel('L/D (-)','FontSize',FS)
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('L_D_mission_vs_V&mbat');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     % Unified plot Options
%     SAVE_types(fname,name,gca,gcf);
% end
% 
% %% Plot contours
% Fig = Fig + 1;
% figure(Fig)
% N_contour_lines = 20; % Number of contour lines
% vect_cc_time_mission_plot = linspace(min(min(L_D_PLOTS')),max(max(L_D_PLOTS')),N_contour_lines);
% [alpha_c,h_alpha_c] = contourf(m_bat_VAR,V_VAR,L_D_PLOTS',vect_cc_time_mission_plot');
% clabel(alpha_c,h_alpha_c)
% hold on
% plot(m_bat_VAR,V_stall(:,1),'m--','LineWidth',LS)
% plot(m_bat_VAR,V_TO(:,1),'o-.','LineWidth',LS)
% plot(m_bat_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% legend('x_{mission}(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest')
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% grid on
% title('L/D vs V and battery mass','FontSize',FS)
% xlabel('batery mass (kg)','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% zlabel('L/D (-)','FontSize',FS)
% hold off
% sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% % caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('L_D_mission_vs_V&mbat_contour');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%             % Unified plot Options             
%             SAVE_types(fname,name,gca,gcf); 
% %             saveas(gcf,fullfile(fname, name),'fig');
% %             saveas(gcf,fullfile(fname, name),'pdf');
% %             saveas(gcf,fullfile(fname, name),'bmp');
% %             saveas(gcf,fullfile(fname, name),'png');
% end


