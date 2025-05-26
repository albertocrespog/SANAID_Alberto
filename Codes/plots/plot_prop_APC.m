function [Fig] = plot_prop_APC(Data_P,Plot_Options,Fig,prefix,VECTOR_Prop)

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;

MATLAB_in = Plot_Options.MATLAB_in;

% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
% PRINT_PLOTS_XFLR5 = Plot_Options.PRINT_PLOTS_XFLR5;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;

Fig = Fig + 1;
figure(Fig)
for i=1:length(VECTOR_Prop)
    plot(Data_P(VECTOR_prop(i)).J,Data_P(VECTOR_prop(i)).Ct,mark_Type{i},'LineWidth', LS); hold on
    leg1{i} = strcat(mark_legend{VECTOR_prop(i)});
end
title('C_T vs J','FontSize',FS)
xlabel('J','FontSize',FS)
xlabel('C_T','FontSize',FS)
if MATLAB_in == 1
    legend(leg1)
else
    h_legend = legend(leg1)
    set(h_legend, 'Location','Best','FontSize',LFS)
end
hold off
grid on
if SAVE_FIGS==1
    name   = strcat(prefix,'CT_vs_J');
    saveas(gcf,name,'fig');
    saveas(gca,name,'epsc');
    saveas(gcf,name,'pdf');
    saveas(gcf,name,'bmp');
    saveas(gcf,name,'png');
end
Fig = Fig + 1;
   
figure(Fig)
for i=1:length(VECTOR_Prop)
    plot(Data_P(VECTOR_prop(i)).J,Data_P(VECTOR_prop(i)).Cp,mark_Type{i},'LineWidth', LS); hold on
    leg2{i} = strcat(mark_legend{VECTOR_prop(i)});
end
title('C_P vs J','FontSize',FS)
xlabel('J','FontSize',FS)
xlabel('C_P','FontSize',FS)
if MATLAB_in == 1
    legend(leg2)
else
    h_legend = legend(leg1)
    set(h_legend, 'Location','Best','FontSize',LFS)
end
hold off
grid on
if SAVE_FIGS==1
    name   = strcat(prefix,'CP_vs_J');
    saveas(gcf,name,'fig');
    saveas(gca,name,'epsc');
    saveas(gcf,name,'pdf');
    saveas(gcf,name,'bmp');
    saveas(gcf,name,'png');
end
Fig = Fig + 1;

figure(Fig)
for i=1:length(VECTOR_Prop)
    plot(Data_P(VECTOR_prop(i)).J,Data_P(VECTOR_prop(i)).Pe,mark_Type{i},'LineWidth', LS); hold on
    leg3{i} = strcat(mark_legend{VECTOR_prop(i)});
end
title('\eta_{m} vs J','FontSize',FS)
xlabel('J','FontSize',FS)
xlabel('\eta_{m}','FontSize',FS)
if MATLAB_in == 1
    legend(leg1)
else
    h_legend = legend(leg3)
    set(h_legend, 'Location','Best','FontSize',LFS)
end
hold off
grid on
if SAVE_FIGS==1
    name   = strcat(prefix,'eta_m_vs_J');
    saveas(gcf,name,'fig');
    saveas(gca,name,'epsc');
    saveas(gcf,name,'pdf');
    saveas(gcf,name,'bmp');
    saveas(gcf,name,'png');
end
Fig = Fig + 1;
