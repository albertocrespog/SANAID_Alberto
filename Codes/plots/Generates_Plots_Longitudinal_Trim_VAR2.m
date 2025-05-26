function [Fig,M_alpha_cero,V_alpha_cero] = Generates_Plots_Longitudinal_Trim_VAR2(TRIM_RESULTS_var,Trim_ITER,Restrictions_var_V_XCG,...
    Stab_Der_var_V_XCG,Geo_tier,Plot_Options,Plot_options2,Fig,OUTPUT_read_XLSX,Aero,filenameS,conv_UNITS)

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;
% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;
MATLAB_in = Plot_Options.MATLAB_in;


m2ft = conv_UNITS.m2ft;
knot2mps = conv_UNITS.knot2mps;
kg2lb = conv_UNITS.kg2lb;
g = conv_UNITS.g;

% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filenamec = 'Figs';
% filename_DATA = strcat(filename,filenamed);
% filename_Plots = strcat(filename,filenamec);
fname = filenameS.plots;

V_VAR = Plot_options2.V_VAR;
n_VAR = Plot_options2.n_VAR;
       
for i=1:length(n_VAR)
    for j=1:length(V_VAR)
        trim_alpha_deg_plot(i,j) = TRIM_RESULTS_var{i,j}.trim_alpha_deg;
        trim_delta_e_deg_plot(i,j) = TRIM_RESULTS_var{i,j}.trim_delta_e_deg;
        trim_delta_T(i,j) = TRIM_RESULTS_var{i,j}.delta_T;
        CL_w1(i,j) = Trim_ITER{i,j}.CL_w1;
%         CL_HTP(i,j) = Trim_ITER{i,j}.CL_HTP;
        CL_needed(i,j) = Trim_ITER{i,j}.CL_needed;
        SM(i,j) = TRIM_RESULTS_var{i,j}.SM;

        trim_CM_alpha_ac(i,j) = TRIM_RESULTS_var{i,j}.CM_alpha_ac;
        trim_CM0_ac(i,j) = TRIM_RESULTS_var{i,j}.CM0_ac;
        
        V_stall(i,j) = Restrictions_var_V_XCG{i,j}.V_stall;
        V_TO(i,j) = Restrictions_var_V_XCG{i,j}.V_TO;
        V_min_ope(i,j) = Restrictions_var_V_XCG{i,j}.V_min_ope;
        n_VAR_MAT(i,j) = n_VAR(i);
        V_VAR_MAT(i,j) = V_VAR(j);
    end
end

CL_alpha_ac(1,1) = TRIM_RESULTS_var{1,1}.CL_alpha_ac;
CL0_ac(1,1) = TRIM_RESULTS_var{1,1}.CL0_ac;

%% Generates the envelope of V-n diagram
W = OUTPUT_read_XLSX.Weights_flags.MTOW_true*9.81;
S = Geo_tier.S_w1;
h = 0;
[Temp,rho,p,a]=atmos_inter(h);

% C_L_max_w1_CR = Aero.CL_max_w1_CR;
% alpha_max_w1_CR = Aero.alpha_max_w1_CR;

% C_L_max_CR = Aero.CL_max_w1_CR;
alpha_max_w1_CR = Aero.alpha_max_w1_CR;
C_L_max_CR = CL0_ac + CL_alpha_ac*alpha_max_w1_CR*pi/180;

C_L_max_neg_CR = -C_L_max_CR*0.75; 
C_D0_CR    = Aero.Polar.C_D0; 
C_D1_CR    = Aero.Polar.C_D1; 
C_D2_CR    = Aero.Polar.C_D2;  

n_max = OUTPUT_read_XLSX.Stability_flags.n_max;
n_min = OUTPUT_read_XLSX.Stability_flags.n_min;

% Vs = sqrt(W/(0.5*rho*C_L_max_CR*S));
Vs = min(min(V_stall(:,:)));
% Vs_neg = sqrt(W/(0.5*rho*abs(C_L_max_neg_CR)*S));
Vs_neg = Vs*(sqrt(1/0.75));

W_2lbs = W*kg2lb/g;
S_2ft = S*m2ft^2;
W_S = W_2lbs/S_2ft;

Kc = compute_kc(W_S);
Vc_knots = Kc*sqrt(W_S);
Vc = Vc_knots*knot2mps;
Vd = 1.25*Vc;
Va = Vs*sqrt(n_max);
Vg = Vs*sqrt(abs(n_min));

n_vec_pos = linspace(0,n_max,100);
n_vec_neg = linspace(0,n_min,100);
V_stall_pos = sqrt(n_vec_pos*2*W/(S*rho*C_L_max_CR));
V_stall_neg = sqrt(n_vec_neg*2*W/(S*rho*C_L_max_neg_CR));

origen = [0 0];
Stall_pos = [V_stall_pos' n_vec_pos'];
Stall_neg = [V_stall_neg' n_vec_neg'];
A = [Va n_max];
C = [Vc n_max];
D = [Vd n_max];
G = [Vg n_min];
F = [Vc n_min];
E = [Vd n_min];

% env = [origen; Stall_pos; A; C; D; E; F; G; origen];
env = [A; C; D; E; F; G];
env_stall_pos = [origen; Stall_pos];
env_stall_neg = [origen; Stall_neg];

perdida = [Vs n_max; Vs n_min];


%% Plots for result of variable trim conditions
Fig = Fig + 1;
figure(Fig)
% subplot (1,2,1)
[x1 y1] = contour(n_VAR_MAT,V_VAR_MAT,trim_alpha_deg_plot,[0 0]); % x1 will be your "two-column matrix"
V_alpha_cero = x1(2,2:end);
M_alpha_cero = x1(1,2:end);
surf(V_VAR,n_VAR,trim_alpha_deg_plot)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('\alpha_{trim} vs V and mass','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('n (gs)','FontSize',FS)
zlabel('\alpha trim','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('alpha_vs_V&n');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end



%% Plot contours
Fig = Fig + 1;
figure(Fig)
N_contour_lines = 30; % Number of contour lines
vect_cc_trim_alpha_deg_plot = linspace(min(min(trim_alpha_deg_plot)),max(max(trim_alpha_deg_plot)),N_contour_lines);
[alpha_c,h_alpha_c] = contourf(V_VAR,n_VAR,trim_alpha_deg_plot,vect_cc_trim_alpha_deg_plot');
clabel(alpha_c,h_alpha_c)
hold on
plot(V_stall(:,1),n_VAR,'m--','LineWidth',LS)
plot(V_TO(:,1),n_VAR,'o-.','LineWidth',LS)
plot(V_min_ope(:,1),n_VAR,'^:','LineWidth',LS)
% plot(V_alpha_cero,M_alpha_cero,'k*-.','LineWidth',LS)
% Generates Envelope
plot(env(:,1),env(:,2), 'k-o','LineWidth',1.5*LS)
plot(env_stall_pos(:,1),env_stall_pos(:,2), 'k-','LineWidth',1.5*LS)
plot(env_stall_neg(:,1),env_stall_neg(:,2), 'k-','LineWidth',1.5*LS)
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','Location','northwest','FontSize',LFS)
shading interp

colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
grid on
title('\alpha_{trim} vs V and mass','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('n (gs)','FontSize',FS)
zlabel('\alpha trim','FontSize',FS)
hold off
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(trim_alpha_deg_plot)) max(max(trim_alpha_deg_plot))]); %%%%%%%%%%%%%%%
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('alpha_vs_V&n_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end

Fig = Fig + 1;
figure(Fig)
% subplot (1,2,1)
surf(V_VAR,n_VAR,trim_delta_e_deg_plot)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
title('\delta_{trim} vs V and mass','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('n (gs)','FontSize',FS)
zlabel('\delta_{trim}','FontSize',FS)
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('delta_e_vs_V&n');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end

Fig = Fig + 1;
figure(Fig)
% subplot (1,2,2)
%% Plot contours
N_contour_lines = 30; % Number of contour lines
vect_cc_trim_delta_deg_plot = linspace(min(min(trim_delta_e_deg_plot)),max(max(trim_delta_e_deg_plot)),N_contour_lines);
[delta_c,h_delta_c] = contourf(V_VAR,n_VAR,trim_delta_e_deg_plot,vect_cc_trim_delta_deg_plot');
clabel(delta_c,h_delta_c)
hold on
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
colorbar % adds the color bar to the right of the graph
plot(V_stall(:,1),n_VAR,'m--','LineWidth',LS)
plot(V_TO(:,1),n_VAR,'o-.','LineWidth',LS)
plot(V_min_ope(:,1),n_VAR,'^:','LineWidth',LS)
% plot(V_alpha_cero,M_alpha_cero,'k*-.','LineWidth',LS)
% Generates Envelope
plot(env(:,1),env(:,2), 'k-o','LineWidth',1.5*LS)
plot(env_stall_pos(:,1),env_stall_pos(:,2), 'k-','LineWidth',1.5*LS)
plot(env_stall_neg(:,1),env_stall_neg(:,2), 'k-','LineWidth',1.5*LS)
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','Location','northwest','FontSize',LFS)
shading interp
grid on
title('\delta_{trim} vs V and mass','FontSize',FS)
xlabel('Velocity  (m/s)','FontSize',FS)
ylabel('n (gs)','FontSize',FS)
zlabel('\delta_{trim}','FontSize',FS)
hold off

if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('delta_e_vs_V&n_contour');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end

% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,1)
% surf(n_VAR,V_VAR,trim_delta_T)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% title('\delta_{T} vs V and mass','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('n (gs)','FontSize',FS)
% zlabel('\delta_{T}','FontSize',FS)
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('deltaT_vs_V&m');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end
% 
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,2)
% %% Plot contours
% N_contour_lines = 20; % Number of contour lines
% vect_cc_trim_deltaT_plot = linspace(min(min(trim_delta_T)),max(max(trim_delta_T)),N_contour_lines);
% [deltaT_c,h_deltaT_c] = contourf(n_VAR,V_VAR,trim_delta_T,vect_cc_trim_deltaT_plot');
% clabel(deltaT_c,h_deltaT_c)
% % colormap(flipud(colormap('gray')))
% shading interp
% grid on
% hold on
% plot(n_VAR,V_stall(1,:),'m--','LineWidth',LS)
% plot(n_VAR,V_TO(1,:),'o-.','LineWidth',LS)
% plot(n_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
% title('\delta_{T} vs V and mass','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('n (gs)','FontSize',FS)
% zlabel('\delta_{T}','FontSize',FS)
% hold off
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('deltaT_vs_V&m_contour');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end

% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,1)
% surf(n_VAR,V_VAR,SM*100)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% title('SM (%) vs V and mass','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('n (gs)','FontSize',FS)
% zlabel('SM (%)','FontSize',FS)
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('SM_vs_V&m');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end

% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,2)
% %% Plot contours
% N_contour_lines = 20; % Number of contour lines
% vect_cc_trim_SM_plot = linspace(min(min(SM)),max(max(SM)),N_contour_lines);
% [delta_SM_c,h_SM_c] = contourf(n_VAR,V_VAR,SM,vect_cc_trim_SM_plot');
% clabel(delta_SM_c,h_SM_c)
% % colormap(flipud(colormap('gray')))
% shading interp
% grid on
% hold on
% plot(n_VAR,V_stall(1,:),'m--','LineWidth',LS)
% plot(n_VAR,V_TO(1,:),'o-.','LineWidth',LS)
% plot(n_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
% title('SM vs V and mass','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('n (gs)','FontSize',FS)
% zlabel('Static Margin (SM)','FontSize',FS)
% hold off
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('SM_vs_V&m');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end

% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,2)
% %% Plot contours
% N_contour_lines = 20; % Number of contour lines
% vect_cc_trim_CMalpha_plot = linspace(min(min(trim_CM_alpha_ac)),max(max(trim_CM_alpha_ac)),N_contour_lines);
% [delta_CMalpha_c,h_CMalpha_c] = contourf(n_VAR,V_VAR,trim_CM_alpha_ac,vect_cc_trim_CMalpha_plot');
% clabel(delta_CMalpha_c,h_CMalpha_c)
% % colormap(flipud(colormap('gray')))
% shading interp
% grid on
% hold on
% plot(n_VAR,V_stall(1,:),'m--','LineWidth',LS)
% plot(n_VAR,V_TO(1,:),'o-.','LineWidth',LS)
% plot(n_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
% title('C_{M_{\alpha}} vs V and mass','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('n (gs)','FontSize',FS)
% zlabel('C_{M_{\alpha}}','FontSize',FS)
% hold off
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('CMalpha_vs_V&m');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end
% 
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,2)
% %% Plot contours
% N_contour_lines = 20; % Number of contour lines
% vect_cc_trim_CM0_plot = linspace(min(min(trim_CM0_ac)),max(max(trim_CM0_ac)),N_contour_lines);
% [delta_CM0_c,h_CM0_c] = contourf(n_VAR,V_VAR,trim_CM0_ac,vect_cc_trim_CM0_plot');
% clabel(delta_CM0_c,h_CM0_c)
% % colormap(flipud(colormap('gray')))
% shading interp
% grid on
% hold on
% plot(n_VAR,V_stall(1,:),'m--','LineWidth',LS)
% plot(n_VAR,V_TO(1,:),'o-.','LineWidth',LS)
% plot(n_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
% title('C_{M_{0}} vs V and mass','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('n (gs)','FontSize',FS)
% zlabel('C_{M_{0}}','FontSize',FS)
% hold off
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('CM0_vs_V&m');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end
% 
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,2)
% %% Plot contours
% N_contour_lines = 20; % Number of contour lines
% vect_cc_trim_CL_w1_plot = linspace(min(min(CL_w1)),max(max(CL_w1)),N_contour_lines);
% [delta_CL_w1_c,h_CL_w1_c] = contourf(n_VAR,V_VAR,CL_w1,vect_cc_trim_CL_w1_plot');
% clabel(delta_CL_w1_c,h_CL_w1_c)
% % colormap(flipud(colormap('gray')))
% shading interp
% grid on
% hold on
% plot(n_VAR,V_stall(1,:),'m--','LineWidth',LS)
% plot(n_VAR,V_TO(1,:),'o-.','LineWidth',LS)
% plot(n_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
% title('C_{L_{w1}} vs V and mass','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('n (gs)','FontSize',FS)
% zlabel('C_{L_{w1}}','FontSize',FS)
% hold off
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('CLw1_vs_V&m');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end
% 
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,2)
% %% Plot contours
% N_contour_lines = 20; % Number of contour lines
% vect_cc_trim_CL_HTP_plot = linspace(min(min(CL_HTP)),max(max(CL_HTP)),N_contour_lines);
% [delta_CL_HTP_c,h_CL_HTP_c] = contourf(n_VAR,V_VAR,CL_HTP,vect_cc_trim_CL_HTP_plot');
% clabel(delta_CL_HTP_c,h_CL_HTP_c)
% % colormap(flipud(colormap('gray')))
% shading interp
% grid on
% hold on
% plot(n_VAR,V_stall(1,:),'m--','LineWidth',LS)
% plot(n_VAR,V_TO(1,:),'o-.','LineWidth',LS)
% plot(n_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
% title('C_{L_{HTP}} vs V and mass','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('n (gs)','FontSize',FS)
% zlabel('C_{L_{HTP}}','FontSize',FS)
% hold off
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('CLHTP_vs_V&m');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end
% 
% Fig = Fig + 1;
% figure(Fig)
% % subplot (1,2,2)
% %% Plot contours
% N_contour_lines = 20; % Number of contour lines
% vect_cc_trim_CL_needed_plot = linspace(min(min(CL_needed)),max(max(CL_needed)),N_contour_lines);
% [delta_CL_needed_c,h_CL_needed_c] = contourf(n_VAR,V_VAR,CL_needed,vect_cc_trim_CL_needed_plot');
% clabel(delta_CL_needed_c,h_CL_needed_c)
% % colormap(flipud(colormap('gray')))
% shading interp
% grid on
% hold on
% plot(n_VAR,V_stall(1,:),'m--','LineWidth',LS)
% plot(n_VAR,V_TO(1,:),'o-.','LineWidth',LS)
% plot(n_VAR,V_min_ope(1,:),'^:','LineWidth',LS)
% plot(M_alpha_cero,V_alpha_cero,'k*-.','LineWidth',LS)
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% colorbar % adds the color bar to the right of the graph
% legend('\alpha(V,m)','V_{stall}','1.2 V_{stall}','1.3 V_{stall}','V_{\alpha=0^\circ}','Location','northwest','FontSize',LFS)
% title('C_{L_{needed}} vs V and mass','FontSize',FS)
% ylabel('Velocity  (m/s)','FontSize',FS)
% xlabel('n (gs)','FontSize',FS)
% zlabel('C_{L_{needed}}','FontSize',FS)
% hold off
% if SAVE_FIGS==1
%     prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
%     st = strcat('CLneeded_vs_V&m');
%     name   = strcat(prefix,st);
%     % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
%     saveas(gca, fullfile(fname, name), 'jpeg');
%     saveas(gcf,fullfile(fname, name),'fig');
%     saveas(gcf,fullfile(fname, name),'pdf');
%     saveas(gcf,fullfile(fname, name),'bmp');
%     saveas(gcf,fullfile(fname, name),'png');
% end
