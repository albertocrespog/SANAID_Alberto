function [mark_Type COLOR_scheme Plot_Options] = Generate_Plot_Options(prefix,mark_legend,VECTOR_XFLR5,Plot_Options,OUTPUT_read_XLSX)
%% Defiones options for the plots
% size of plot letter
% LS = 2; % Line size
% FS = 12; % Text Font size
% LFS = 8; % Legend Font Size
% Fig = 0;
% Video_3D = 0; % Saves video 3D
% SAVE_FIGS = 0; % saves the plots: fig, jpg and pdf

LS = OUTPUT_read_XLSX.PLOT_flags.LS; % Line size	LS
FS = OUTPUT_read_XLSX.PLOT_flags.FS; % Text Font size
LFS = OUTPUT_read_XLSX.PLOT_flags.LFS; % Line size	LS
Video_3D = OUTPUT_read_XLSX.PLOT_flags.Video_3D; % Legend Font Size
SAVE_FIGS = OUTPUT_read_XLSX.PLOT_flags.SAVE_FIGS; % Line size	LS

% Store DATA
Plot_Options.LS = LS;
Plot_Options.FS = FS;
Plot_Options.LFS = LFS;
% Plot_Options.Fig = Fig;
Plot_Options.Video_3D = Video_3D;
Plot_Options.SAVE_FIGS = SAVE_FIGS;
Plot_Options.mark_legend = mark_legend;
Plot_Options.VECTOR_XFLR5 = VECTOR_XFLR5;

PL = 1;

%% Defines color scheme for 3D plots
color_fus = [0.2 0.2 0.8]; % fuselage
color_w1 = [0.5 0.5 0.8]; % w1
color_w2 = [0.2 0.6 0.8]; % w2
color_HTP = [0.15 0.6 0.8]; % HTP
color_vee = [0.25 0.65 0.75]; % vee
color_vee2 = [0.25 0.65 0.85]; % vee2
color_can = [0.6 0.2 0.8]; % can
color_vtp = [0.4 0.4 0.8]; % vtp
color_prop = [0.9 0.9 0.1]; % Propeller 
color_eng = [0.7 0.4 0.3]; % engine
color_nac = [0.5 0.3 0.2]; % nacelle
color_ac = [0.2000 0.2000 0.2000]; % aicraft

% Storing DATA
COLOR_scheme.color_fus = color_fus;
COLOR_scheme.color_w1 = color_w1;
COLOR_scheme.color_w2 = color_w2;
COLOR_scheme.color_HTP = color_HTP;
COLOR_scheme.color_vee = color_vee;
COLOR_scheme.color_vee2 = color_vee2;
COLOR_scheme.color_can = color_can;
COLOR_scheme.color_vtp = color_vtp;
COLOR_scheme.color_prop = color_prop;
COLOR_scheme.color_eng = color_eng;
COLOR_scheme.color_nac = color_nac;
COLOR_scheme.color_ac = color_ac;

%% Defines the Plotting options, different style lines
mrk_Law1  = 'r-';
mrk_Law2  = 'b-.';
mrk_Law3  = 'g--';
mrk_Law4  = 'k:';
mrk_Law5  = 'm-.';
mrk_Law6  = 'r:';
mrk_Law7  = 'b-+';
mrk_Law8  = 'g-s';
mrk_Law9  = 'k-d';
mrk_Law10 = 'r-^';
mrk_Law11 = 'b-v';
mrk_Law12 = 'g->';
mrk_Law13 = 'k-<';
mrk_Law14 = 'm:';
mrk_Law15 = 'b:';
mrk_Law16 = 'g:';
mrk_Law17 = 'k:';
mrk_Law18 = 'r:';
mrk_Law19  = 'b-+';
mrk_Law20  = 'g-s';
mrk_Law21  = 'k-d';
mrk_Law22  = 'r-^';
mrk_Law23  = 'b-v';
mrk_Law24  = 'g->';

mark_Type   = {mrk_Law1,mrk_Law2,mrk_Law3,mrk_Law4,mrk_Law5,...
                mrk_Law6,mrk_Law7,mrk_Law8,mrk_Law9,mrk_Law10,...
                mrk_Law11,mrk_Law12,mrk_Law13,mrk_Law14,mrk_Law15,...
                mrk_Law16,mrk_Law17,mrk_Law18,mrk_Law19,mrk_Law20,...
                mrk_Law21,mrk_Law22,mrk_Law23,mrk_Law24};
            
Plot_Options.mark_Type = mark_Type;
