function OUTPUT_read_XLSX = initialization_PLOTS_in_MISSIONS(OUTPUT_read_XLSX);


%% PLOTS
% plot_perfo1 = Prints PLOTS PERFORMANCE STUDY
% plot_perfo2 = Prints plots of Performance for Variable V and mass - Electric
% plot_perfo3 = Prints plots of Performance for Variable h and V
% plot_perfo4 = Prints plots of Performance Glide for Variable h and V
% plot_perfo5 = Prints plots of Performance Glide for Variable h Max
% Initializes
OUTPUT_read_XLSX.PLOT_flags.plot_perfo1 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_perfo2 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_perfo3 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_perfo4 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_perfo5 = 0;
% Initializes Stability Plots
OUTPUT_read_XLSX.PLOT_flags.plot_stability1 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability2 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability3 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability4 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability5 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability6 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability7 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability8 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability9 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability10 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability11 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability12 = 0;
OUTPUT_read_XLSX.PLOT_flags.plot_stability13 = 0;

% initializes Optimization
OUTPUT_read_XLSX.STUDY_flags.Perfo1 = 0; % Variable Study
OUTPUT_read_XLSX.STUDY_flags.Perfo2 = 0; % Variable Study
OUTPUT_read_XLSX.STUDY_flags.Perfo3 = 0; % Variable Study
OUTPUT_read_XLSX.STUDY_flags.Perfo4 = 0; % Variable Study
OUTPUT_read_XLSX.STUDY_flags.Perfo5 = 0; % Variable Study
OUTPUT_read_XLSX.STUDY_flags.Perfo6 = 0; % Variable Study
OUTPUT_read_XLSX.STUDY_flags.Perfo7 = 0; % Variable Study