function [Fig] = Generates_Plots_Lateral_Trim(Trim_ITER_LAT,...
    Geo_tier,Plot_Options,conv_UNITS,conditions_TRIM_lat,OUTPUT_read_XLSX,Fig,filenameS)

g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;

% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS
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
% 
deltaa_deg = Trim_ITER_LAT.deltaa_deg;
deltar_deg = Trim_ITER_LAT.deltar_deg;
phi_deg = Trim_ITER_LAT.phi_deg;
beta_vec = conditions_TRIM_lat.beta_vec;
beta = conditions_TRIM_lat.beta;

deltaa_deg_var = Trim_ITER_LAT.deltaa_deg_var;
deltar_deg_var = Trim_ITER_LAT.deltar_deg_var;
phi_deg_var = Trim_ITER_LAT.phi_deg_var;

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
title('Lateral Trim Conditions - Sidewind')
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
    st = strcat('da_dr_phi&m');
    name   = strcat(prefix,st);
    % % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end