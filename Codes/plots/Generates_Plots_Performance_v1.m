function [Fig] = Generates_Plots_Performance_v1(Geo_tier,Plot_Options,conv_UNITS,Fig,...
    Segments, output_caracteristics,Weights_AP_var,fuel_total_var,tiempo_total_var,distancia_total_var,W_var,datos_var,case_AC,Plots_performance,Post_processing_PERFORMANCE)

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        prefix = strcat('ProVant4_');
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        prefix = strcat('ProVant4_');
    case 3 % case_AC = 3 - PEPIÑO XXL
        prefix = strcat('PEPINOXXL_');
    case 4 % case_AC = 3 - PEPIÑO XXL
        prefix = strcat('Comercial_');
    case 5 % case_AC = 3 - PEPIÑO XXL
        prefix = strcat('WIG_');
    case 6 % case_AC = 3 - PEPIÑO XXL
        prefix = strcat('CERVERA_');     
end

CI = Plots_performance.CI;

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS*0.75;
% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;
MATLAB_in = Plot_Options.MATLAB_in;

N_V_VAR_perf = Plot_Options.N_V_VAR_perf;
V_VAR_perf = Plot_Options.V_VAR_perf;

Wp_VAR_perf = Plot_Options.Wp_VAR_perf;
N_Wp_VAR_perf = Plot_Options.N_Wp_VAR_perf;
N_Wp_plot = Plot_Options.N_Wp_plot;

kg2Tm = conv_UNITS.kg2Tm;
m2km = conv_UNITS.m2km;
seg2hrs = conv_UNITS.seg2hrs;
km2nm = conv_UNITS.km2nm;
l2gal = conv_UNITS.l2gal;

m_f_PLOT =     Plots_performance.m_f_PLOT;
tiempo_total_PLOT = Plots_performance.tiempo_total_PLOT;
distancia_total_PLOT = Plots_performance.distancia_total_PLOT;
Cost_fuel = Plots_performance.Cost_fuel;
DOC_PLOT = Plots_performance.DOC_PLOT;
CAPM_PLOT = Plots_performance.CAPM_PLOT;
distancia_PLOT = Plots_performance.distancia_PLOT;
L_D_PLOT = Plots_performance.L_D_PLOT;
tiempo_PLOT = Plots_performance.tiempo_PLOT;
CL_PLOT = Plots_performance.CL_PLOT;
palanca_PLOT = Plots_performance.palanca_PLOT;

% CI = 0; % kg/s
% density_fuel =  0.809; % kg/l
% cost_fuel = 161.69;%  cts/gal - 7 Feb 2020
% % Generation of variable to plot
% for i=1:N_V_VAR_perf
%     for j=1:N_Wp_VAR_perf
%         m_f_PLOT(i,j) = Weights_AP_var{i,j}.m_f;
%         tiempo_total_PLOT(i,j) = tiempo_total_var{i,j};
%         distancia_total_PLOT(i,j) = distancia_total_var{i,j};
%         Cost_fuel(i,j) = (m_f_PLOT(i,j)/density_fuel)*l2gal*cost_fuel;
%         DOC_PLOT(i,j) = (tiempo_total_PLOT(i,j)*CI + Cost_fuel(i,j));
%         CAPM_PLOT(i,j) = DOC_PLOT(i,j)/((distancia_total_PLOT(i,j)*m2km)*(Wp_VAR_perf(j)*kg2Tm));
%     end
% end
%
% for j=1:length(datos_var{i}.crucero.tiempo)
%     distancia_PLOT(j) = datos_var{i}.crucero.distancia(j);
% end
%
% % Generation of variable to plot for Cruise
% for i=1:N_V_VAR_perf
%     for j=1:length(datos_var{i,N_Wp_plot}.crucero.L_D)
%         L_D_PLOT(i,j) = datos_var{i,N_Wp_plot}.crucero.L_D(j);
%         tiempo_PLOT(i,j) = datos_var{i,N_Wp_plot}.crucero.tiempo(j);
%         CL_PLOT(i,j)  = datos_var{i,N_Wp_plot}.crucero.CL(j);
%         palanca_PLOT(i,j) = datos_var{i,N_Wp_plot}.crucero.palanca(j);
%     end
% end

if Post_processing_PERFORMANCE == 1
    
    Fig = Fig + 1;
    figure(Fig)
    for i=1:N_Wp_VAR_perf
        plot(V_VAR_perf,m_f_PLOT(:,i),mark_Type{i},'LineWidth', LS);
        leg1{i} = strcat('W_{payload}= ',num2str(Wp_VAR_perf(i)),'kg');
        hold on
    end
    grid on
    title('m_{fuel} vs. V and W_{payload}','FontSize',FS)
    xlabel('V (m/s)','FontSize',FS)
    ylabel('m_{fuel} - (kg)','FontSize',FS)
    if MATLAB_in == 1
        legend(leg1)
    else
        h_legend=legend(leg1);
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    hold off
    grid on
    if SAVE_FIGS==1
        name   = strcat(prefix,'mf_vs_W_pl_&_V_trends');
        saveas(gcf,name,'fig');
        saveas(gca,name,'epsc');
        saveas(gcf,name,'pdf');
        saveas(gcf,name,'bmp');
        saveas(gcf,name,'png');
    end
    
    Fig = Fig + 1;
    figure(Fig)
    for i=1:N_Wp_VAR_perf
        plot(V_VAR_perf,tiempo_total_PLOT(:,i)*seg2hrs,mark_Type{i},'LineWidth', LS);
        leg1{i} = strcat('W_{payload}= ',num2str(Wp_VAR_perf(i)),'kg');
        hold on
    end
    grid on
    title('time vs. V and W_{payload}','FontSize',FS)
    xlabel('V (m/s)','FontSize',FS)
    ylabel('Time (hrs)','FontSize',FS)
    if MATLAB_in == 1
        legend(leg1)
    else
        h_legend=legend(leg1);
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    hold off
    grid on
    if SAVE_FIGS==1
        name   = strcat(prefix,'mf_vs_W_pl_&_V');
        saveas(gcf,name,'fig');
        saveas(gca,name,'epsc');
        saveas(gcf,name,'pdf');
        saveas(gcf,name,'bmp');
        saveas(gcf,name,'png');
    end
    
    cts2Kdollar = 1/(100*1000);
    Fig = Fig + 1;
    figure(Fig)
    for i=1:N_Wp_VAR_perf
        plot(V_VAR_perf,Cost_fuel(:,i)*cts2Kdollar,mark_Type{i},'LineWidth', LS);
        leg1{i} = strcat('W_{payload}= ',num2str(Wp_VAR_perf(i)),'kg');
        hold on
    end
    grid on
    title('Fuel Cost vs. V and W_{payload}','FontSize',FS)
    xlabel('V (m/s)','FontSize',FS)
    ylabel('Cost Fuel (K$)','FontSize',FS)
    if MATLAB_in == 1
        legend(leg1)
    else
        h_legend=legend(leg1);
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    hold off
    grid on
    if SAVE_FIGS==1
        name   = strcat(prefix,'Cost_Fuel_vs_W_pl_&_V');
        saveas(gcf,name,'fig');
        saveas(gca,name,'epsc');
        saveas(gcf,name,'pdf');
        saveas(gcf,name,'bmp');
        saveas(gcf,name,'png');
    end
    
    Fig = Fig + 1;
    figure(Fig)
    for i=1:N_Wp_VAR_perf
        plot(V_VAR_perf,DOC_PLOT(:,i)*cts2Kdollar,mark_Type{i},'LineWidth', LS);
        leg1{i} = strcat('W_{payload}= ',num2str(Wp_VAR_perf(i)),'kg, and CI = ',num2str(CI),' Kg/seg');
        hold on
    end
    grid on
    title('DOC vs. V and W_{payload}','FontSize',FS)
    xlabel('V (m/s)','FontSize',FS)
    ylabel('DOC Fuel (K$)','FontSize',FS)
    if MATLAB_in == 1
        legend(leg1)
    else
        h_legend=legend(leg1);
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    hold off
    grid on
    if SAVE_FIGS==1
        name   = strcat(prefix,'DOC_vs_W_pl_&_V');
        saveas(gcf,name,'fig');
        saveas(gca,name,'epsc');
        saveas(gcf,name,'pdf');
        saveas(gcf,name,'bmp');
        saveas(gcf,name,'png');
    end
    
    Fig = Fig + 1;
    figure(Fig)
    for i=1:N_Wp_VAR_perf
        plot(V_VAR_perf,CAPM_PLOT(:,i),mark_Type{i},'LineWidth', LS);
        leg1{i} = strcat('W_{payload}= ',num2str(Wp_VAR_perf(i)),'kg, and CI = ',num2str(CI),' Kg/seg');
        hold on
    end
    grid on
    title('CAPM vs. V and W_{payload}','FontSize',FS)
    xlabel('V (m/s)','FontSize',FS)
    ylabel('Cost per Available Payload Mile (CAPM) (cts$/km ton)','FontSize',FS)
    if MATLAB_in == 1
        legend(leg1)
    else
        h_legend=legend(leg1);
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    hold off
    grid on
    if SAVE_FIGS==1
        name   = strcat(prefix,'CAPM_vs_W_pl_&_V');
        saveas(gcf,name,'fig');
        saveas(gca,name,'epsc');
        saveas(gcf,name,'pdf');
        saveas(gcf,name,'bmp');
        saveas(gcf,name,'png');
    end
    
    if N_Wp_plot>1
        Fig = Fig + 1;
        figure(Fig)
        surf(Wp_VAR_perf,V_VAR_perf,m_f_PLOT)
        shading interp
        colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
        %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
        %    rotate3d on % allows you to rotate the graph with the mouse
        colorbar % adds the color bar to the right of the graph
        title('m_f vs W_{payload} and V')
        ylabel('V (m/s)')
        xlabel('W_{payload} - (kg)')
        zlabel('m_f - (kg)')
        grid on
        hold off
        if SAVE_FIGS==1
            name   = strcat(prefix,'mf_vs_W_pl_&_V');
            saveas(gcf,name,'fig');
            saveas(gca,name,'epsc');
            saveas(gcf,name,'pdf');
            saveas(gcf,name,'bmp');
            saveas(gcf,name,'png');
        end
        
        % Fig = Fig + 1;
        % figure(Fig)
        % plot(V_VAR_perf,tiempo_total_PLOT)
        % title('time vs V','FontSize',FS)
        % xlabel('V (m/s)','FontSize',FS)
        % ylabel('time (sec)','FontSize',FS)
        % grid on
        % hold off
        % if SAVE_FIGS==1
        %     name   = strcat(prefix,'time_vs_V');
        %     saveas(gcf,name,'fig');
        %     saveas(gca,name,'epsc');
        %     saveas(gcf,name,'pdf');
        %     saveas(gcf,name,'bmp');
        %     saveas(gcf,name,'png');
        % end
        
        
        Fig = Fig + 1;
        figure(Fig)
        surf(distancia_PLOT,V_VAR_perf,L_D_PLOT)
        shading interp
        colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
        %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
        %    rotate3d on % allows you to rotate the graph with the mouse
        colorbar % adds the color bar to the right of the graph
        Title2 = strcat('L/D vs distance and V, for W_{payload}= ',num2str(Wp_VAR_perf(N_Wp_plot)),'kg');
        title(Title2)
        ylabel('V (m/s)')
        xlabel('x (Km)')
        zlabel('L/D')
        grid on
        hold off
        if SAVE_FIGS==1
            name   = strcat(prefix,'time_vs_L_D');
            saveas(gcf,name,'fig');
            saveas(gca,name,'epsc');
            saveas(gcf,name,'pdf');
            saveas(gcf,name,'bmp');
            saveas(gcf,name,'png');
        end
        
        Fig = Fig + 1;
        figure(Fig)
        surf(distancia_PLOT*m2km,V_VAR_perf,CL_PLOT)
        shading interp
        colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
        %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
        %    rotate3d on % allows you to rotate the graph with the mouse
        colorbar % adds the color bar to the right of the graph
        Title3 = strcat('C_L vs distance and V, for W_{payload}= ',num2str(Wp_VAR_perf(N_Wp_plot)),'kg');
        title(Title3)
        ylabel('V (m/s)')
        xlabel('x (Km)')
        zlabel('C_L')
        grid on
        hold off
        if SAVE_FIGS==1
            name   = strcat(prefix,'time_vs_CL');
            saveas(gcf,name,'fig');
            saveas(gca,name,'epsc');
            saveas(gcf,name,'pdf');
            saveas(gcf,name,'bmp');
            saveas(gcf,name,'png');
        end
        
        Fig = Fig + 1;
        figure(Fig)
        surf(distancia_PLOT,V_VAR_perf,palanca_PLOT)
        shading interp
        colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
        %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
        %    rotate3d on % allows you to rotate the graph with the mouse
        colorbar % adds the color bar to the right of the graph
        Title2 = strcat('\delta_T vs distance and V, for W_{payload}= ',num2str(Wp_VAR_perf(N_Wp_plot)),'kg');
        title(Title2)
        ylabel('V (m/s)')
        xlabel('x (Km)')
        zlabel('\delta_T')
        grid on
        hold off
        if SAVE_FIGS==1
            name   = strcat(prefix,'time_vs_deltaT');
            saveas(gcf,name,'fig');
            saveas(gca,name,'epsc');
            saveas(gcf,name,'pdf');
            saveas(gcf,name,'bmp');
            saveas(gcf,name,'png');
        end
        
        Fig = Fig + 1;
        figure(Fig)
        surf(distancia_PLOT,V_VAR_perf,tiempo_PLOT*seg2hrs)
        shading interp
        colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
        %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
        %    rotate3d on % allows you to rotate the graph with the mouse
        colorbar % adds the color bar to the right of the graph
        Title1 = strcat('Time vs distance and V, for W_{payload}= ',num2str(Wp_VAR_perf(N_Wp_plot)),'kg');
        title(Title1)
        % title('Time vs distance and V, and ')
        ylabel('V (m/s)')
        xlabel('x (Km)')
        zlabel('Time (hrs)')
        grid on
        hold off
        if SAVE_FIGS==1
            name   = strcat(prefix,'time_vs_deltaT');
            saveas(gcf,name,'fig');
            saveas(gca,name,'epsc');
            saveas(gcf,name,'pdf');
            saveas(gcf,name,'bmp');
            saveas(gcf,name,'png');
        end
    end
end
end


