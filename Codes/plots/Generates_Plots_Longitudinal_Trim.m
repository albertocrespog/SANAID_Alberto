function [Fig] = Generates_Plots_Longitudinal_Trim(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var,Trim_ITER_var,Geo_tier,Plot_Options,OUTPUT_read_XLSX,Fig,filenameS)

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

% x_XCG_fwd = TRIM_RESULTS.x_XCG_fwd;
% x_XCG_rwd = TRIM_RESULTS.x_XCG_rwd;
% N_x_XCG_VAR = 100;
% x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);
x_XCG_VAR = Plot_Options.x_XCG_VAR;

% SM_des = Trim_ITER.SM_des;
% % Generation of variable to plot
for i=1:length(x_XCG_VAR)
    CM_alpha_w1w2b_PLOT(i) = TRIM_RESULTS_var{i}.CM_alpha_ac;
    CM0_w1w2b_PLOT(i) = TRIM_RESULTS_var{i}.CM;
    X_NP_PLOT(i) = TRIM_RESULTS_var{i}.X_NP;
    SM_real_PLOT(i) = TRIM_RESULTS_var{i}.SM;
    %     SM_actual_PLOT(i) = TRIM_RESULTS_var{i}.SM_actual;
    trim_alpha_deg_PLOT(i) = TRIM_RESULTS_var{i}.trim_alpha_deg;
    trim_delta_e_deg_PLOT(i) = TRIM_RESULTS_var{i}.trim_delta_e_deg;
    X_NP_PLOT(i) = TRIM_RESULTS_var{i}.X_NP;
end

SM_des = Trim_ITER.SM_des;
% Estimation of Real Static Marging satisfying CM_alpha = 0
X_AC_w1Bw2_real = interp1(CM_alpha_w1w2b_PLOT,x_XCG_VAR,0,'spline');
% Fot
for i=1:length(x_XCG_VAR)
    SM_real_CMalpha_PLOT(i) = (X_AC_w1Bw2_real - x_XCG_VAR(i))/Geo_tier.cmac_w1;
    SM_des_PLOT(i) = SM_des;
end

Fig = Fig + 1;
figure(Fig)
plot(x_XCG_VAR,CM_alpha_w1w2b_PLOT)
title('C_{M_\alpha} vs X_{CG}','FontSize',FS)
xlabel('X_{CG}','FontSize',FS)
ylabel('C_{M_\alpha}','FontSize',FS)
grid on
hold off
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('CMalpha_vs_XCG');
    name   = strcat(prefix,st);
    % % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
    %             saveas(gcf,fullfile(fname, name),'fig');
    %             saveas(gcf,fullfile(fname, name),'pdf');
    %             saveas(gcf,fullfile(fname, name),'bmp');
    %             saveas(gcf,fullfile(fname, name),'png');
    
end

Fig = Fig + 1;
figure(Fig)
%     plot(x_XCG_VAR,SM_actual_PLOT,'k')
plot(x_XCG_VAR,SM_real_PLOT,'b-.')
hold on
plot(x_XCG_VAR,SM_real_CMalpha_PLOT,'r:')
plot(x_XCG_VAR,SM_des_PLOT,'r','LineWidth',LS)
title('Static Margin vs X_{CG}','FontSize',FS)
xlabel('X_{CG}','FontSize',FS)
ylabel('SM','FontSize',FS)
legend('SM = -C_{M_{\alpha_{w1&B&w2&}}}/C_{L_{\alpha_{w1&B&w2&}}}',...
    'SM = (X_{NP}-x_{CG})/cmac_{w_1}',...
    'SM = (X_{ac_{w1Bw2}}-X_{CG})/c_{ac_{w_1}} with X_{ac_{w1Bw2}} from C_{M_{\alpha_{w1&B&w2&}}} = 0',...
    'SM_{des}')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('SM_vs_XCG');
    name   = strcat(prefix,st);
    % % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
    %             saveas(gcf,fullfile(fname, name),'fig');
    %             saveas(gcf,fullfile(fname, name),'pdf');
    %             saveas(gcf,fullfile(fname, name),'bmp');
    %             saveas(gcf,fullfile(fname, name),'png');
    
end

Fig = Fig + 1;
figure(Fig)
plot(x_XCG_VAR,CM0_w1w2b_PLOT,'LineWidth',LS)
title('C_{M_0} vs X_{CG}','FontSize',FS)
xlabel('X_{CG}','FontSize',FS)
ylabel('C_{M_0}','FontSize',FS)
grid on
hold off
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('CM0_vs_XCG');
    name   = strcat(prefix,st);
    % % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
    %             saveas(gcf,fullfile(fname, name),'fig');
    %             saveas(gcf,fullfile(fname, name),'pdf');
    %             saveas(gcf,fullfile(fname, name),'bmp');
    %             saveas(gcf,fullfile(fname, name),'png');
    
end

Fig = Fig + 1;
figure(Fig)
plot(x_XCG_VAR,trim_alpha_deg_PLOT)
title('\alpha_{trim} vs X_{CG}','FontSize',FS)
xlabel('X_{CG}','FontSize',FS)
ylabel('\alpha_{trim}','FontSize',FS)
grid on
hold off
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('alpha_vs_XCG');
    name   = strcat(prefix,st);
    % % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
    %             saveas(gcf,fullfile(fname, name),'fig');
    %             saveas(gcf,fullfile(fname, name),'pdf');
    %             saveas(gcf,fullfile(fname, name),'bmp');
    %             saveas(gcf,fullfile(fname, name),'png');
    
end

Fig = Fig + 1;
figure(Fig)
plot(x_XCG_VAR,trim_delta_e_deg_PLOT)
title('\delta_{trim} vs X_{CG}','FontSize',FS)
xlabel('X_{CG}','FontSize',FS)
ylabel('\delta_{trim}','FontSize',FS)
grid on
hold off
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('delta_vs_XCG');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
    %             saveas(gcf,fullfile(fname, name),'fig');
    %             saveas(gcf,fullfile(fname, name),'pdf');
    %             saveas(gcf,fullfile(fname, name),'bmp');
    %             saveas(gcf,fullfile(fname, name),'png');
    
end
