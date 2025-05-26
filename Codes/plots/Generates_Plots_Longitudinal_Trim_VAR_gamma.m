function [Fig,M_alpha_cero,V_alpha_cero] = Generates_Plots_Longitudinal_Trim_VAR_gamma(TRIM_RESULTS_var,Trim_ITER,Restrictions_var_V_XCG,...
    Stab_Der_var_V_XCG,Geo_tier,Plot_Options,Fig,OUTPUT_read_XLSX,x_XCG_VAR,filenameS,conv_UNITS)

% g = conv_UNITS.g;
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

% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filenamec = 'Figs';
% filename_DATA = strcat(filename,filenamed);
% filename_Plots = strcat(filename,filenamec);
fname = filenameS.plots;

V_VAR = Plot_Options.V_VAR;
gamma_VAR = Plot_Options.gamma_VAR*R2D;
       
for i=1:length(gamma_VAR)
    for j=1:length(V_VAR)
        trim_alpha_deg_plot(i,j) = TRIM_RESULTS_var{i,j}.trim_alpha_deg;
        trim_delta_e_deg_plot(i,j) = TRIM_RESULTS_var{i,j}.trim_delta_e_deg;
        trim_delta_T(i,j) = TRIM_RESULTS_var{i,j}.delta_T;
        CL_w1(i,j) = Trim_ITER{i,j}.CL_w1;
%        CL_HTP(i,j) = Trim_ITER{i,j}.CL_HTP;
        CL_needed(i,j) = Trim_ITER{i,j}.CL_needed;
        SM(i,j) = TRIM_RESULTS_var{i,j}.SM;

        trim_CM_alpha_ac(i,j) = TRIM_RESULTS_var{i,j}.CM_alpha_ac;
        trim_CM0_ac(i,j) = TRIM_RESULTS_var{i,j}.CM0_ac;
        
        V_stall(i,j) = Restrictions_var_V_XCG{i,j}.V_stall;
        V_TO(i,j) = Restrictions_var_V_XCG{i,j}.V_TO;
        V_min_ope(i,j) = Restrictions_var_V_XCG{i,j}.V_min_ope;
        gamma_VAR_MAT(i,j) = gamma_VAR(i);
        V_VAR_MAT(i,j) = V_VAR(j);
        x_XCG_MAT(i,j) = x_XCG_VAR(i,j);
    end
end


[x1 y1] = contour(V_VAR_MAT,gamma_VAR_MAT,trim_alpha_deg_plot,[0 0]); % x1 will be your "two-column matrix"
M_alpha_cero = x1(2,2:end);
V_alpha_cero = x1(1,2:end);


%% Plots for result of variable trim conditions
Fig = Fig + 1;
figure(Fig)
% subplot (1,2,1)
surf(gamma_VAR,V_VAR,trim_alpha_deg_plot')
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
%like hot, o r graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\alpha_{trim} vs V and \gamma','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\gamma (deg)','FontSize',FS)
zlabel('\alpha trim','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('alpha_vs_V&gamma');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    SAVE_types(fname,name,gca,gcf);
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
end


%% Plot contours
Fig = Fig + 1;
figure(Fig)
N_contour_lines = 20; % Number of contour lines
vect_cc_trim_alpha_deg_plot = linspace(min(min(trim_alpha_deg_plot)),max(max(trim_alpha_deg_plot)),N_contour_lines);
[alpha_c,h_alpha_c] = contourf(gamma_VAR,V_VAR,trim_alpha_deg_plot',vect_cc_trim_alpha_deg_plot);
clabel(alpha_c,h_alpha_c)
hold on
plot(gamma_VAR,V_stall(:,1),'m--','LineWidth',LS)
plot(gamma_VAR,V_TO(:,1),'o-.','LineWidth',LS)
plot(gamma_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest')
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
grid on
title('\alpha_{trim} vs V and \gamma','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\gamma (deg)','FontSize',FS)
zlabel('\alpha trim','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('alpha_vs_V&gamma_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
% subplot (1,2,1)
surf(gamma_VAR,V_VAR,trim_delta_e_deg_plot')
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
title('\delta_{trim} vs V and \gamma','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\gamma (deg)','FontSize',FS)
zlabel('\delta_{trim}','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('delta_e_vs_V&gamma');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
% subplot (1,2,2)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_trim_delta_deg_plot = linspace(min(min(trim_delta_e_deg_plot)),max(max(trim_delta_e_deg_plot)),N_contour_lines);
[delta_c,h_delta_c] = contourf(gamma_VAR,V_VAR,trim_delta_e_deg_plot',vect_cc_trim_delta_deg_plot);
clabel(delta_c,h_delta_c)
hold on
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
plot(gamma_VAR,V_stall(:,1),'m--','LineWidth',LS)
plot(gamma_VAR,V_TO(:,1),'o-.','LineWidth',LS)
plot(gamma_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
shading interp
grid on
title('\delta_{trim} vs V and \gamma','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\gamma (deg)','FontSize',FS)
zlabel('\delta_{trim}','FontSize',FS)
hold off

if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('delta_e_vs_V&gamma_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
% subplot (1,2,1)
surf(gamma_VAR,V_VAR,trim_delta_T')
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
title('\delta_{T} vs V and \gamma','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\gamma (deg)','FontSize',FS)
zlabel('\delta_{T}','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltaT_vs_V&gamma');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
% subplot (1,2,2)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_trim_deltaT_plot = linspace(min(min(trim_delta_T)),max(max(trim_delta_T)),N_contour_lines);
[deltaT_c,h_deltaT_c] = contourf(gamma_VAR,V_VAR,trim_delta_T',vect_cc_trim_deltaT_plot');
clabel(deltaT_c,h_deltaT_c)
% colormap(flipud(colormap('gray')))
shading interp
grid on
hold on
plot(gamma_VAR,V_stall(:,1),'m--','LineWidth',LS)
plot(gamma_VAR,V_TO(:,1),'o-.','LineWidth',LS)
plot(gamma_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
title('\delta_{T} vs V and \gamma','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\gamma (deg)','FontSize',FS)
zlabel('\delta_{T}','FontSize',FS)
hold off
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('deltaT_vs_V&gamma_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end



SM_values = SM * 100;
isConstantZData = all(SM_values(:) == SM_values(1)); % Chequeo de si todos los valores son iguales
if isConstantZData
    Message_plot1 = 'The ZData (SM*100) is constant. Contour plot will not be generated.';
    disp(Message_plot1)
else
    Fig = Fig + 1;
    figure(Fig)
    % subplot (1,2,1)
    surf(gamma_VAR,V_VAR,SM'*100)
    shading interp
    colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
    colorbar % adds the color bar to the right of the graph
    title('SM (%) vs V and \gamma','FontSize',FS)
    ylabel('Velocity  (m/s)','FontSize',FS)
    xlabel('\gamma (deg)','FontSize',FS)
    zlabel('SM (%)','FontSize',FS)
    if SAVE_FIGS==1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('SM_vs_V&gamma');
        name   = strcat(prefix,st);
        % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
        % Unified plot Options
        SAVE_types(fname,name,gca,gcf);
        %             saveas(gcf,fullfile(fname, name),'fig');
        %             saveas(gcf,fullfile(fname, name),'pdf');
        %             saveas(gcf,fullfile(fname, name),'bmp');
        %             saveas(gcf,fullfile(fname, name),'png');

    end

    Fig = Fig + 1;
    figure(Fig)
    % subplot (1,2,2)
    %% Plot contours
    N_contour_lines = 20; % Number of contour lines
    vect_cc_trim_SM_plot = linspace(min(min(SM*100)),max(max(SM*100)),N_contour_lines);
    [delta_SM_c,h_XCG_c] = contourf(gamma_VAR,V_VAR,SM*100',vect_cc_trim_SM_plot');
    clabel(delta_SM_c,h_XCG_c)
    % colormap(flipud(colormap('gray')))
    shading interp
    grid on
    hold on
    plot(gamma_VAR,V_stall(:,1),'m--','LineWidth',LS)
    plot(gamma_VAR,V_TO(:,1),'o-.','LineWidth',LS)
    plot(gamma_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
    plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
    colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
    colorbar % adds the color bar to the right of the graph
    legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
    title('SM vs V and \gamma','FontSize',FS)
    ylabel('Velocity  (m/s)','FontSize',FS)
    xlabel('\gamma (deg)','FontSize',FS)
    zlabel('Static Margin (SM)','FontSize',FS)
    hold off
    if SAVE_FIGS==1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('SM_vs_V&gamma');
        name   = strcat(prefix,st);
        % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
        % Unified plot Options
        SAVE_types(fname,name,gca,gcf);
        %             saveas(gcf,fullfile(fname, name),'fig');
        %             saveas(gcf,fullfile(fname, name),'pdf');
        %             saveas(gcf,fullfile(fname, name),'bmp');
        %             saveas(gcf,fullfile(fname, name),'png');

    end
end


isConstantZData = all(trim_CM_alpha_ac(:) == trim_CM_alpha_ac(1)); % Chequeo de si todos los valores son iguales
if isConstantZData
    Message_plot1 = 'The ZData (trim_CM_alpha_ac) is constant. Contour plot will not be generated.';
    disp(Message_plot1)
    % warning('The ZData (trim_CM_alpha_ac) is constant. Contour plot will not be generated.');
else
    Fig = Fig + 1;
    figure(Fig)
    % subplot (1,2,1)
    surf(gamma_VAR,V_VAR,trim_CM_alpha_ac')
    shading interp
    colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
    colorbar % adds the color bar to the right of the graph
    title('C_{M_{\alpha}} vs V and \gamma vs V and \gamma','FontSize',FS)
    ylabel('Velocity  (m/s)','FontSize',FS)
    xlabel('\gamma (deg)','FontSize',FS)
    zlabel('C_{M_{\alpha}}','FontSize',FS)
    if SAVE_FIGS==1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('CMalpha_vs_V&gamma');
        name   = strcat(prefix,st);
        % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
        % Unified plot Options
        SAVE_types(fname,name,gca,gcf);
        %             saveas(gcf,fullfile(fname, name),'fig');
        %             saveas(gcf,fullfile(fname, name),'pdf');
        %             saveas(gcf,fullfile(fname, name),'bmp');
        %             saveas(gcf,fullfile(fname, name),'png');

    end

    Fig = Fig + 1;
    figure(Fig)
    % subplot (1,2,2)
    %% Plot contours
    N_contour_lines = 20; % Number of contour lines
    vect_cc_trim_CMalpha_plot = linspace(min(min(trim_CM_alpha_ac)),max(max(trim_CM_alpha_ac)),N_contour_lines);
    [delta_CMalpha_c,h_CMalpha_c] = contourf(gamma_VAR,V_VAR,trim_CM_alpha_ac',vect_cc_trim_CMalpha_plot');
    clabel(delta_CMalpha_c,h_CMalpha_c)
    % colormap(flipud(colormap('gray')))
    shading interp
    grid on
    hold on
    plot(gamma_VAR,V_stall(:,1),'m--','LineWidth',LS)
    plot(gamma_VAR,V_TO(:,1),'o-.','LineWidth',LS)
    plot(gamma_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
    plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
    colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
    colorbar % adds the color bar to the right of the graph
    legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
    title('C_{M_{\alpha}} vs V and \gamma','FontSize',FS)
    ylabel('Velocity  (m/s)','FontSize',FS)
    xlabel('\gamma (deg)','FontSize',FS)
    zlabel('C_{M_{\alpha}}','FontSize',FS)
    hold off
    if SAVE_FIGS==1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('CMalpha_vs_V&gamma_contour');
        name   = strcat(prefix,st);
        % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
        % Unified plot Options
        SAVE_types(fname,name,gca,gcf);
        %             saveas(gcf,fullfile(fname, name),'fig');
        %             saveas(gcf,fullfile(fname, name),'pdf');
        %             saveas(gcf,fullfile(fname, name),'bmp');
        %             saveas(gcf,fullfile(fname, name),'png');

    end
end



isConstantZData = all(trim_CM0_ac(:) == trim_CM0_ac(1)); % Chequeo de si todos los valores son iguales
if isConstantZData
    Message_plot1 = 'The ZData (trim_CM0_ac) is constant. Contour plot will not be generated.';
    disp(Message_plot1)
    % warning('The ZData (trim_CM0_ac) is constant. Contour plot will not be generated.');
else
    Fig = Fig + 1;
    figure(Fig)
    % subplot (1,2,1)
    surf(gamma_VAR,V_VAR,trim_CM0_ac')
    shading interp
    colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
    colorbar % adds the color bar to the right of the graph
    title('C_{M_{0}} vs V and \gamma vs V and \gamma','FontSize',FS)
    ylabel('Velocity  (m/s)','FontSize',FS)
    xlabel('\gamma (deg)','FontSize',FS)
    zlabel('C_{M_{0}}','FontSize',FS)
    if SAVE_FIGS==1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('CM0_vs_V&gamma');
        name   = strcat(prefix,st);
        % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
        % Unified plot Options
        SAVE_types(fname,name,gca,gcf);
        %             saveas(gcf,fullfile(fname, name),'fig');
        %             saveas(gcf,fullfile(fname, name),'pdf');
        %             saveas(gcf,fullfile(fname, name),'bmp');
        %             saveas(gcf,fullfile(fname, name),'png');

    end

Fig = Fig + 1;
    figure(Fig)
    % subplot (1,2,2)
    %% Plot contours
    N_contour_lines = 20; % Number of contour lines
    vect_cc_trim_CM0_plot = linspace(min(min(trim_CM0_ac)),max(max(trim_CM0_ac)),N_contour_lines);
    [delta_CM0_c,h_CM0_c] = contourf(gamma_VAR,V_VAR,trim_CM0_ac',vect_cc_trim_CM0_plot');
    clabel(delta_CM0_c,h_CM0_c)
    % colormap(flipud(colormap('gray')))
    shading interp
    grid on
    hold on
    plot(gamma_VAR,V_stall(:,1),'m--','LineWidth',LS)
    plot(gamma_VAR,V_TO(:,1),'o-.','LineWidth',LS)
    plot(gamma_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
    plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
    colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
    colorbar % adds the color bar to the right of the graph
    legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
    title('C_{M_{0}} vs V and \gamma','FontSize',FS)
    ylabel('Velocity  (m/s)','FontSize',FS)
    xlabel('\gamma (deg)','FontSize',FS)
    zlabel('C_{M_{0}}','FontSize',FS)
    hold off
    if SAVE_FIGS==1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('CM0_vs_V&gamma_contour');
        name   = strcat(prefix,st);
        % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
        % Unified plot Options
        SAVE_types(fname,name,gca,gcf);
        %             saveas(gcf,fullfile(fname, name),'fig');
        %             saveas(gcf,fullfile(fname, name),'pdf');
        %             saveas(gcf,fullfile(fname, name),'bmp');
        %             saveas(gcf,fullfile(fname, name),'png');

    end
end

Fig = Fig + 1;
figure(Fig)
% subplot (1,2,1)
surf(gamma_VAR,V_VAR,CL_w1')
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
title('C_{L_{w1}} vs V and \gamma vs V and \gamma','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\gamma (deg)','FontSize',FS)
zlabel('C_{L_{w1}}','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('SM_vs_V&gamma');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
% subplot (1,2,2)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_trim_CL_w1_plot = linspace(min(min(CL_w1)),max(max(CL_w1)),N_contour_lines);
[delta_CL_w1_c,h_CL_w1_c] = contourf(gamma_VAR,V_VAR,CL_w1',vect_cc_trim_CL_w1_plot');
clabel(delta_CL_w1_c,h_CL_w1_c)
% colormap(flipud(colormap('gray')))
shading interp
grid on
hold on
plot(gamma_VAR,V_stall(:,1),'m--','LineWidth',LS)
plot(gamma_VAR,V_TO(:,1),'o-.','LineWidth',LS)
plot(gamma_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
title('C_{L_{w1}} vs V and \gamma','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\gamma (deg)','FontSize',FS)
zlabel('C_{L_{w1}}','FontSize',FS)
hold off
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('CLw1_vs_V&gamma');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,2)
% %% Plot contours
% N_contour_lines = 20; % Number of contour lines
% vect_cc_trim_CL_HTP_plot = linspace(min(min(CL_HTP)),max(max(CL_HTP)),N_contour_lines);
% [delta_CL_HTP_c,h_CL_HTP_c] = contourf(gamma_VAR,V_VAR,CL_HTP,vect_cc_trim_CL_HTP_plot');
% clabel(delta_CL_HTP_c,h_CL_HTP_c)
% % colormap(flipud(colormap('gray')))
% shading interp
% grid on
% hold on
% plot(gamma_VAR,V_stall(1,:),'m--','LineWidth',LS)
% plot(gamma_VAR,V_TO(1,:),'o-.','LineWidth',LS)
% plot(gamma_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
% title('C_{L_{HTP}} vs V and \gamma','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('\gamma (deg)','FontSize',FS)
% zlabel('C_{L_{HTP}}','FontSize',FS)
% hold off
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('CLHTP_vs_V&gamma');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end

Fig = Fig + 1;
figure(Fig)
% subplot (1,2,2)
%% Plot contours
N_contour_lines = 20; % Number of contour lines
vect_cc_trim_CL_needed_plot = linspace(min(min(CL_needed)),max(max(CL_needed)),N_contour_lines);
[delta_CL_needed_c,h_CL_needed_c] = contourf(gamma_VAR,V_VAR,CL_needed',vect_cc_trim_CL_needed_plot');
clabel(delta_CL_needed_c,h_CL_needed_c)
% colormap(flipud(colormap('gray')))
shading interp
grid on
hold on
plot(gamma_VAR,V_stall(:,1),'m--','LineWidth',LS)
plot(gamma_VAR,V_TO(:,1),'o-.','LineWidth',LS)
plot(gamma_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
title('C_{L_{needed}} vs V and \gamma','FontSize',FS)
ylabel('Velocity  (m/s)','FontSize',FS)
xlabel('\gamma (deg)','FontSize',FS)
zlabel('C_{L_{needed}}','FontSize',FS)
hold off
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('CLneeded_vs_V&gamma');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
    %             saveas(gcf,fullfile(fname, name),'fig');
    %             saveas(gcf,fullfile(fname, name),'pdf');
    %             saveas(gcf,fullfile(fname, name),'bmp');
    %             saveas(gcf,fullfile(fname, name),'png');

end

isConstantZData = all(x_XCG_MAT(:) == x_XCG_MAT(1)); % Chequeo de si todos los valores son iguales
if isConstantZData
    Message_plot1 = 'The ZData (x_XCG_MAT) is constant. Contour plot will not be generated.';
    disp(Message_plot1)
    % warning('The ZData (x_XCG_MAT) is constant. Contour plot will not be generated.');
else
    Fig = Fig + 1;
    figure(Fig)
    % subplot (1,2,1)
    surf(gamma_VAR,V_VAR,x_XCG_MAT')
    shading interp
    colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
    colorbar % adds the color bar to the right of the graph
    title('$x_{CG}$ (m) vs V and \gamma','FontSize',FS)
    ylabel('Velocity  (m/s)','FontSize',FS)
    xlabel('\gamma (deg)','FontSize',FS)
    zlabel('x_{CG} (m)','FontSize',FS)
    if SAVE_FIGS==1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('XCG_vs_V&gamma');
        name   = strcat(prefix,st);
        % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
        % Unified plot Options
        SAVE_types(fname,name,gca,gcf);
        %             saveas(gcf,fullfile(fname, name),'fig');
        %             saveas(gcf,fullfile(fname, name),'pdf');
        %             saveas(gcf,fullfile(fname, name),'bmp');
        %             saveas(gcf,fullfile(fname, name),'png');

    end


Fig = Fig + 1;
    figure(Fig)
    % subplot (1,2,2)
    %% Plot contours
    N_contour_lines = 20; % Number of contour lines
    vect_cc_trim_XCG_plot = linspace(min(min(x_XCG_MAT)),max(max(x_XCG_MAT)),N_contour_lines);
    [delta_XCG_c,h_XCG_c] = contourf(gamma_VAR,V_VAR,x_XCG_MAT',vect_cc_trim_XCG_plot');
    clabel(delta_XCG_c,h_XCG_c)
    % colormap(flipud(colormap('gray')))
    shading interp
    grid on
    hold on
    plot(gamma_VAR,V_stall(:,1),'m--','LineWidth',LS)
    plot(gamma_VAR,V_TO(:,1),'o-.','LineWidth',LS)
    plot(gamma_VAR,V_min_ope(:,1),'^:','LineWidth',LS)
    plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
    colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
    colorbar % adds the color bar to the right of the graph
    legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
    title('$x_{CG}$ vs V and \gamma','FontSize',FS)
    ylabel('Velocity  (m/s)','FontSize',FS)
    xlabel('\gamma (deg)','FontSize',FS)
    zlabel('Center of Gravity (m)','FontSize',FS)
    hold off
    if SAVE_FIGS==1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('XCG_vs_V&gamma');
        name   = strcat(prefix,st);
        % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
        % Unified plot Options
        SAVE_types(fname,name,gca,gcf);
        %             saveas(gcf,fullfile(fname, name),'fig');
        %             saveas(gcf,fullfile(fname, name),'pdf');
        %             saveas(gcf,fullfile(fname, name),'bmp');
        %             saveas(gcf,fullfile(fname, name),'png');

    end
end
