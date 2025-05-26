function [Fig] = plot_XFLR5_Polar_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix,OUTPUT_read_XLSX,AC_CONFIGURATION,filenameS)

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
%         if OUTPUT_read_XLSX.PLOT_flags.plot_individuals ==1; %
            Fig = Fig + 1;
            figure(Fig)
            for i=1:length(VECTOR)
                plot(DATA_Ae(VECTOR(i)).CD,DATA_Ae(VECTOR(i)).CL,mark_Type{i},'LineWidth', LS); hold on
                leg1{i} = strcat(mark_legend{VECTOR(i)});
            end
            leg1{i+1} = strcat('poly-fit Polar');
            plot(DATA_PL.CD_poly,DATA_Ae(VECTOR(i)).CL,mark_Type{i+1},'LineWidth', LS)
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
                name   = strcat(prefix,st)
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