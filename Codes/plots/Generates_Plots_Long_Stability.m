function [Fig] = Generates_Plots_Long_Stability(TRIM_RESULTS,Trim_ITER,TRIM_RESULTS_var,Geo_tier,Plot_Options,Fig)

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;
% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;
MATLAB_in = Plot_Options.MATLAB_in;

x_XCG_fwd = TRIM_RESULTS.x_XCG_fwd;
x_XCG_rwd = TRIM_RESULTS.x_XCG_rwd;
N_x_XCG_VAR = 100;
x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);

SM_des = Trim_ITER.SM_des;
% Generation of variable to plot
for i=1:N_x_XCG_VAR
    CM_alpha_w1w2b_PLOT(i) = TRIM_RESULTS_var{i}.CM_alpha_w1w2b;
    CM0_w1w2b_PLOT(i) = TRIM_RESULTS_var{i}.CM;
    X_NP_PLOT(i) = TRIM_RESULTS_var{i}.X_NP;
    SM_real_PLOT(i) = TRIM_RESULTS_var{i}.SM_real;
    SM_actual_PLOT(i) = TRIM_RESULTS_var{i}.SM_actual;
end

% Estimation of Real Static Marging satisfying CM_aalpha = 0
X_AC_w1Bw2_real = interp1(CM_alpha_w1w2b_PLOT,x_XCG_VAR,0,'spline');
% Fot 
for i=1:N_x_XCG_VAR
    SM_real_CMalpha_PLOT(i) = (X_AC_w1Bw2_real - x_XCG_VAR(i))/Geo_tier.cmac_w1;
    SM_des_PLOT(i) = SM_des;
end

PRINT_PLOTS_STABILITY_SM = Plot_Options.PRINT_PLOTS_STABILITY_SM;

%  PRINT_PLOTS_STABILITY_SM
if PRINT_PLOTS_STABILITY_SM ==1
    
    Fig = Fig + 1;
    figure(Fig)
    plot(x_XCG_VAR,CM_alpha_w1w2b_PLOT)
    title('C_{M_\alpha} vs X_{CG}','FontSize',FS)
    xlabel('X_{CG}','FontSize',FS)
    ylabel('C_{M_\alpha}','FontSize',FS)
    grid on
    hold off
    if SAVE_FIGS==1
        name   = strcat(prefix,'CMalpha_vs_XCG');    
        saveas(gcf,name,'fig');
        saveas(gca,name,'epsc');
        saveas(gcf,name,'pdf');
        saveas(gcf,name,'bmp');
        saveas(gcf,name,'png');
    end
    
    Fig = Fig + 1;
    figure(Fig)
    plot(x_XCG_VAR,SM_actual_PLOT,'k')
    hold on
    plot(x_XCG_VAR,SM_real_PLOT,'b-.')
    plot(x_XCG_VAR,SM_real_CMalpha_PLOT,'r:')
    plot(x_XCG_VAR,SM_des_PLOT,'r','LineWidth',LS*2)
    title('Static Margin vs X_{CG}','FontSize',FS)
    xlabel('X_{CG}','FontSize',FS)
    ylabel('SM','FontSize',FS)
    legend('SM = -C_{M_{\alpha_{w1&B&w2&}}}/C_{L_{\alpha_{w1&B&w2&}}}',...
        'SM = (X_{NP}-x_{CG})/cmac_{w_1}',...
        'SM = (X_{ac_{w1Bw2}}-X_{CG})/c_{ac_{w_1}} with X_{ac_{w1Bw2}} from C_{M_{\alpha_{w1&B&w2&}}} = 0',...
        'SM_{des}')
    grid on
    if SAVE_FIGS==1
        name   = strcat(prefix,'SM_vs_XCG');    
        saveas(gcf,name,'fig');
        saveas(gca,name,'epsc');
        saveas(gcf,name,'pdf');
        saveas(gcf,name,'bmp');
        saveas(gcf,name,'png');
    end
end