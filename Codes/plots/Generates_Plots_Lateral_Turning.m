function [Fig] = Generates_Plots_Lateral_Turning(Trim_ITER_LAT_Viraje,Geo_tier,...
    Plot_Options,conv_UNITS,conditions_TRIM_turning,OUTPUT_read_XLSX,Fig,filenameS)

g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

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

% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;
MATLAB_in = Plot_Options.MATLAB_in;

phi = conditions_TRIM_turning.phi;
phi_vec = conditions_TRIM_turning.phi_vec;
n_viraje = conditions_TRIM_turning.n_viraje;

% Ixx = Weight_tier.Ixx;
% Iyy = Weight_tier.Iyy;
% Izz = Weight_tier.Izz;
% Ixz = Weight_tier.Ixz;
% 
% % Conversion units and constans
% g = conv_UNITS.g;
% R2D = conv_UNITS.R2D;
% D2R = conv_UNITS.D2R;
% ft2m = conv_UNITS.ft2m;
% m2ft = conv_UNITS.m2ft;
% W2hp = conv_UNITS.W2hp;
% mps2ftps = conv_UNITS.mps2ftps;
% kg2lb = conv_UNITS.kg2lb;
% lb2kg = conv_UNITS.lb2kg;
% ftpm2mps = conv_UNITS.ftpm2mps;
% in2m = conv_UNITS.in2m;
% W2pftsec = conv_UNITS.W2pftsec;
% m22ft2 = conv_UNITS.m22ft2;
% rho_SI2rho_IMP = conv_UNITS.rho_SI2rho_IMP;
% qmet2qimp = conv_UNITS.qmet2qimp;
% N2lbf = conv_UNITS.N2lbf;

beta_viraje_var = Trim_ITER_LAT_Viraje.beta_viraje_var;
deltaa_viraje_var = Trim_ITER_LAT_Viraje.deltaa_viraje_var;
deltar_viraje_var = Trim_ITER_LAT_Viraje.deltar_viraje_var;

beta_deg_viraje_var = Trim_ITER_LAT_Viraje.beta_deg_viraje_var;
deltaa_deg_viraje_var = Trim_ITER_LAT_Viraje.deltaa_deg_viraje_var;
deltar_deg_viraje_var = Trim_ITER_LAT_Viraje.deltar_deg_viraje_var;

n_viraje_var = Trim_ITER_LAT_Viraje.n_viraje_var;
psi_dot_viraje_var = Trim_ITER_LAT_Viraje.psi_dot_viraje_var;
R_t_viraje_var = Trim_ITER_LAT_Viraje.R_t_viraje_var;

Fig = Fig + 1;
figure(Fig)
plot(phi_vec*R2D,deltaa_deg_viraje_var,'b','LineWidth', LS)
hold on
plot(phi_vec*R2D,deltar_deg_viraje_var,'r','LineWidth', LS)
plot(phi_vec*R2D,beta_deg_viraje_var,'g','LineWidth', LS)
hold off
title('Lateral Trim Conditions - Level Turn')
xlabel('\phi (deg)')
ylabel('\delta_a, \delta_r, \beta - (deg)')
if MATLAB_in == 1
    legend(leg1)
else
    h_legend=legend('\delta_a','\delta_r','\beta');
    set(h_legend, 'Location','Best','FontSize',LFS)
end
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('da_dr_phi');
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
plot(phi_vec*R2D,n_viraje_var,'b','LineWidth', LS)
title('Load Factor vs \phi - Level Turn')
xlabel('\phi (deg)')
ylabel('n')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('n');
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
plot(phi_vec*R2D,psi_dot_viraje_var*R2D,'b','LineWidth', LS)
title('Turn Speed vs \phi - Level Turn')
xlabel('\phi (deg)')
ylabel('\psi-dot (deg/seg)')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('psi');
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
plot(phi_vec*R2D,R_t_viraje_var,'b','LineWidth', LS)
title('Turn Radius vs \phi - Level Turn')
xlabel('\phi (deg)')
ylabel('R_t (m)')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('TurnRadius_phi_lateral_LevelTurn');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end
