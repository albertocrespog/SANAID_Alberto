function [Fig] = plot_XFLR5_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix,OUTPUT_read_XLSX,AC_CONFIGURATION,filenameS)

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filenamec = 'Figs';
% filename_DATA = strcat(filename,filenamed);
% filename_Plots = strcat(filename,filenamec);
% fname = filenameS.filename_Plots;
fname = filenameS.plots;

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;
% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
% PRINT_PLOTS_XFLR5 = Plot_Options.PRINT_PLOTS_XFLR5;
% PRINT_PLOTS_AERO = Plot_Options.PRINT_PLOTS_AERO;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;
% VECTOR = Plot_Options.VECTOR;
VECTOR_XFLR5 = Plot_Options.VECTOR_XFLR5;
MATLAB_in = Plot_Options.MATLAB_in;

index_w1 =  OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1;
index_HTP =  OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_HTP;
index_vee =  OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee;
index_vee2 =  OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee2;
index_VTP =  OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP;

compare = VECTOR_XFLR5.compare;
if compare == 1
    VEC{1} = VECTOR_XFLR5.v1;
elseif compare == 2
    VEC{1} = VECTOR_XFLR5.v1;
    VEC{2} = VECTOR_XFLR5.v2;
elseif compare == 3
    VEC{1} = VECTOR_XFLR5.v1;
    VEC{2} = VECTOR_XFLR5.v2;
    VEC{3} = VECTOR_XFLR5.v3;
end

for n=1:length(VEC)
    if n == 1 || n == 2
        VECTOR = VEC{n};
        if OUTPUT_read_XLSX.PLOT_flags.plot_individuals ==1; %
            Fig = Fig + 1;
            figure(Fig)
            for i=1:length(VECTOR)
                plot(DATA_Ae(VECTOR(i)).alpha,DATA_Ae(VECTOR(i)).CL,mark_Type{i},'LineWidth', LS); hold on
                leg1{i} = strcat(mark_legend{VECTOR(i)});
            end
            title('C_L vs \alpha','FontSize',FS)
            xlabel('alpha (deg)','FontSize',FS)
            ylabel('C_L','FontSize',FS)
            if MATLAB_in == 1
                legend(leg1)
            else
                h_legend = legend(leg1);
                set(h_legend, 'Location','Best','FontSize',LFS)
            end
            hold off
            grid on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('CL_vs_alpha');
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
            for i=1:length(VECTOR)
                plot(DATA_Ae(VECTOR(i)).CD,DATA_Ae(VECTOR(i)).CL,mark_Type{i},'LineWidth', LS); hold on
                leg1{i} = strcat(mark_legend{VECTOR(i)});
            end
            %     plot(CD1,DATA_Ae(VECTOR(1)).CL,mark_Type{3},'LineWidth', LS); hold on
            %         leg1{3} = strcat('CD = CDv + CDi - XFLR5');
            %     plot(CD2,DATA_Ae(VECTOR(1)).CL,mark_Type{4},'LineWidth', LS); hold on
            %         leg1{4} = strcat('MacCrowing Correction to CDi OGE');
            
            title('C_L vs C_D','FontSize',FS)
            xlabel('C_D','FontSize',FS)
            ylabel('C_L','FontSize',FS)
            if MATLAB_in == 1
                legend(leg1)
            else
                h_legend = legend(leg1);
                set(h_legend, 'Location','Best','FontSize',LFS)
            end
            hold off
            grid on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('CL_vs_CD');
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
            for i=1:length(VECTOR)
                plot(DATA_Ae(VECTOR(i)).alpha,DATA_Ae(VECTOR(i)).Cm,mark_Type{i},'LineWidth', LS); hold on
                leg1{i} = strcat(mark_legend{VECTOR(i)});
            end
            title('C_M vs \alpha','FontSize',FS)
            xlabel('alpha (deg)','FontSize',FS)
            ylabel('C_M','FontSize',FS)
            if MATLAB_in == 1
                legend(leg1)
            else
                h_legend = legend(leg1);
                set(h_legend, 'Location','Best','FontSize',LFS)
            end
            hold off
            grid on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('CM_vs_alpha');
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
        
        %% Colect All graphs
        Fig = Fig + 1;
        figure(Fig)
        subplot(1,3,1)
        for i=1:length(VECTOR)
            plot(DATA_Ae(VECTOR(i)).alpha,DATA_Ae(VECTOR(i)).CL,mark_Type{i},'LineWidth', LS); hold on
            leg1{i} = strcat(mark_legend{VECTOR(i)});
        end
        title('C_L vs \alpha','FontSize',FS)
        xlabel('alpha (deg)','FontSize',FS)
        ylabel('C_L','FontSize',FS)
        if MATLAB_in == 1
            legend(leg1)
        else
            h_legend = legend(leg1);
            set(h_legend, 'Location','Best','FontSize',LFS)
        end
        hold off
        grid on
        
        subplot(1,3,2)
        for i=1:length(VECTOR)
            plot(DATA_Ae(VECTOR(i)).CD,DATA_Ae(VECTOR(i)).CL,mark_Type{i},'LineWidth', LS); hold on
            leg1{i} = strcat(mark_legend{VECTOR(i)});
        end
        title('C_L vs C_D','FontSize',FS)
        xlabel('C_D','FontSize',FS)
        ylabel('C_L','FontSize',FS)
        if MATLAB_in == 1
            legend(leg1)
        else
            h_legend = legend(leg1);
            set(h_legend, 'Location','Best','FontSize',LFS)
        end
        hold off
        grid on
        
        subplot(1,3,3)
        for i=1:length(VECTOR)
            plot(DATA_Ae(VECTOR(i)).alpha,DATA_Ae(VECTOR(i)).Cm,mark_Type{i},'LineWidth', LS); hold on
            leg1{i} = strcat(mark_legend{VECTOR(i)});
        end
        title('C_M vs \alpha','FontSize',FS)
        xlabel('alpha (deg)','FontSize',FS)
        ylabel('C_M','FontSize',FS)
        if MATLAB_in == 1
            legend(leg1)
        else
            h_legend = legend(leg1);
            set(h_legend, 'Location','Best','FontSize',LFS)
        end
        hold off
        grid on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CLCDCM_vs_alpha');
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
        subplot(1,3,1)
        for i=1:length(VECTOR)
            CLCD = DATA_Ae(VECTOR(i)).CL./DATA_Ae(VECTOR(i)).CD;
            plot(DATA_Ae(VECTOR(i)).alpha,CLCD,mark_Type{i},'LineWidth', LS); hold on
            leg1{i} = strcat(mark_legend{VECTOR(i)});
        end
        title('C_L/C_D vs \alpha','FontSize',FS)
        xlabel('alpha (deg)','FontSize',FS)
        ylabel('C_L/C_D','FontSize',FS)
        if MATLAB_in == 1
            legend(leg1)
        else
            h_legend = legend(leg1);
            set(h_legend, 'Location','Best','FontSize',LFS)
        end
        hold off
        grid on
        
        figure(Fig)
        subplot(1,3,2)
        for i=1:length(VECTOR)
            CL1_2_CD = (DATA_Ae(VECTOR(i)).CL).^(1/2)./DATA_Ae(VECTOR(i)).CD;
            plot(DATA_Ae(VECTOR(i)).alpha,CL1_2_CD,mark_Type{i},'LineWidth', LS); hold on
            leg1{i} = strcat(mark_legend{VECTOR(i)});
        end
        title('C_L^{1/2}/C_D vs \alpha','FontSize',FS)
        xlabel('alpha (deg)','FontSize',FS)
        ylabel('C_L^{1/2}/C_D','FontSize',FS)
        if MATLAB_in == 1
            legend(leg1)
        else
            h_legend = legend(leg1);
            set(h_legend, 'Location','Best','FontSize',LFS)
        end
        hold off
        grid on

        figure(Fig)
        subplot(1,3,3)
        for i=1:length(VECTOR)
            CL3_2_CD = (DATA_Ae(VECTOR(i)).CL).^(3/2)./DATA_Ae(VECTOR(i)).CD;
            plot(DATA_Ae(VECTOR(i)).alpha,CL3_2_CD,mark_Type{i},'LineWidth', LS); hold on
            leg1{i} = strcat(mark_legend{VECTOR(i)});
        end
        title('C_L^{3/2}/C_D vs \alpha','FontSize',FS)
        xlabel('alpha (deg)','FontSize',FS)
        ylabel('C_L^{3/2}/C_D','FontSize',FS)
        if MATLAB_in == 1
            legend(leg1)
        else
            h_legend = legend(leg1);
            set(h_legend, 'Location','Best','FontSize',LFS)
        end
        hold off
        grid on

    end
    
    %% plots different values for VTP
    if n == 3
        VECTOR = VEC{n};
        if VTP == 1
            if OUTPUT_read_XLSX.PLOT_flags.plot_individuals ==1; %
                Fig = Fig + 1;
                figure(Fig)
                for i=1:length(VECTOR)
                    plot(DATA_Ae(VECTOR(i)).beta,DATA_Ae(VECTOR(i)).CY,mark_Type{i},'LineWidth', LS); hold on
                    leg1{i} = strcat(mark_legend{VECTOR(i)});
                end
                title('C_Y vs \beta','FontSize',FS)
                xlabel('beta (deg)','FontSize',FS)
                ylabel('C_Y','FontSize',FS)
                if MATLAB_in == 1
                    legend(leg1)
                else
                    h_legend = legend(leg1);
                    set(h_legend, 'Location','Best','FontSize',LFS)
                end
                hold off
                grid on
                if SAVE_FIGS==1
                    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                    st = strcat('CY_vs_beta');
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
                for i=1:length(VECTOR)
                    plot(DATA_Ae(VECTOR(i)).beta,DATA_Ae(VECTOR(i)).Cl,mark_Type{i},'LineWidth', LS); hold on
                    leg1{i} = strcat(mark_legend{VECTOR(i)});
                end
                title('C_l vs \beta','FontSize',FS)
                xlabel('beta (deg)','FontSize',FS)
                ylabel('C_l','FontSize',FS)
                if MATLAB_in == 1
                    legend(leg1)
                else
                    h_legend = legend(leg1);
                    set(h_legend, 'Location','Best','FontSize',LFS)
                end
                hold off
                grid on
                if SAVE_FIGS==1
                    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                    st = strcat('Cl_vs_beta');
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
                for i=1:length(VECTOR)
                    plot(DATA_Ae(VECTOR(i)).beta,DATA_Ae(VECTOR(i)).Cn,mark_Type{i},'LineWidth', LS); hold on
                    leg1{i} = strcat(mark_legend{VECTOR(i)});
                end
                title('C_n vs \beta','FontSize',FS)
                xlabel('beta (deg)','FontSize',FS)
                ylabel('C_n','FontSize',FS)
                if MATLAB_in == 1
                    legend(leg1)
                else
                    h_legend = legend(leg1);
                    set(h_legend, 'Location','Best','FontSize',LFS)
                end
                hold off
                grid on
                if SAVE_FIGS==1
                    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                    st = strcat('Cn_vs_beta');
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
                for i=1:length(VECTOR)
                    plot(DATA_Ae(VECTOR(i)).beta,DATA_Ae(VECTOR(i)).CD,mark_Type{i},'LineWidth', LS); hold on
                    leg1{i} = strcat(mark_legend{VECTOR(i)});
                end
                title('C_D vs \beta','FontSize',FS)
                xlabel('beta (deg)','FontSize',FS)
                ylabel('C_D','FontSize',FS)
                if MATLAB_in == 1
                    legend(leg1)
                else
                    h_legend = legend(leg1);
                    set(h_legend, 'Location','Best','FontSize',LFS)
                end
                hold off
                grid on
                if SAVE_FIGS==1
                    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                    st = strcat('CD_vs_beta');
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
            
            %% Colect All graphs
            Fig = Fig + 1;
            figure(Fig)
            subplot(2,2,1)
            for i=1:length(VECTOR)
                plot(DATA_Ae(VECTOR(i)).beta,DATA_Ae(VECTOR(i)).CY,mark_Type{i},'LineWidth', LS); hold on
                leg1{i} = strcat(mark_legend{VECTOR(i)});
            end
            title('C_Y vs \beta','FontSize',FS)
            xlabel('beta (deg)','FontSize',FS)
            ylabel('C_Y','FontSize',FS)
            if MATLAB_in == 1
                legend(leg1)
            else
                h_legend = legend(leg1);
                set(h_legend, 'Location','Best','FontSize',LFS)
            end
            hold off
            grid on
            
            subplot(2,2,2)
            for i=1:length(VECTOR)
                plot(DATA_Ae(VECTOR(i)).beta,DATA_Ae(VECTOR(i)).Cl,mark_Type{i},'LineWidth', LS); hold on
                leg1{i} = strcat(mark_legend{VECTOR(i)});
            end
            title('C_l vs \beta','FontSize',FS)
            xlabel('beta (deg)','FontSize',FS)
            ylabel('C_l','FontSize',FS)
            if MATLAB_in == 1
                legend(leg1)
            else
                h_legend = legend(leg1);
                set(h_legend, 'Location','Best','FontSize',LFS)
            end
            hold off
            grid on
            
            subplot(2,2,3)
            for i=1:length(VECTOR)
                plot(DATA_Ae(VECTOR(i)).beta,DATA_Ae(VECTOR(i)).Cn,mark_Type{i},'LineWidth', LS); hold on
                leg1{i} = strcat(mark_legend{VECTOR(i)});
            end
            title('C_n vs \beta','FontSize',FS)
            xlabel('beta (deg)','FontSize',FS)
            ylabel('C_n','FontSize',FS)
            if MATLAB_in == 1
                legend(leg1)
            else
                h_legend = legend(leg1);
                set(h_legend, 'Location','Best','FontSize',LFS)
            end
            hold off
            grid on
            
            subplot(2,2,4)
            for i=1:length(VECTOR)
                plot(DATA_Ae(VECTOR(i)).beta,DATA_Ae(VECTOR(i)).CD,mark_Type{i},'LineWidth', LS); hold on
                leg1{i} = strcat(mark_legend{VECTOR(i)});
            end
            title('C_D vs \beta','FontSize',FS)
            xlabel('beta (deg)','FontSize',FS)
            ylabel('C_D','FontSize',FS)
            if MATLAB_in == 1
                legend(leg1)
            else
                h_legend = legend(leg1);
                set(h_legend, 'Location','Best','FontSize',LFS)
            end
            hold off
            grid on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('C_LC_nCD_vs_beta');
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
    end
    
    
end