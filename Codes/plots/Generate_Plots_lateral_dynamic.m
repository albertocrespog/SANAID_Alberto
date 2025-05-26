function [Fig] = Generate_Plots_lateral_dynamic(Plot_Options,OUTPUT_read_XLSX,conv_UNITS,Fig,Stab_Dyn_LatDir,filenameS)

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;
% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
MATLAB_in = Plot_Options.MATLAB_in;

% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filenamec = 'Figs';
% filename_DATA = strcat(filename,filenamed);
% filename_Plots = strcat(filename,filenamec);
fname = filenameS.plots;

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

%% Loads data
% saved_file = 'Lateral.mat';
% vars = whos('-file', saved_file);
% load(saved_file, vars(1:2).name)
Matrix_lat     = Stab_Dyn_LatDir.Matrix_lat;
Deflection_lat = Stab_Dyn_LatDir.Deflection_lat;
A_latdir = Matrix_lat;
B_latdir = Deflection_lat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%RESPUESTAS LONGITUDINALES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=0;
C1=[1 0 0 0 0];
C2=[0 1 0 0 0];
C3=[0 0 1 0 0];
C4=[0 0 0 1 0];
C5=[0 0 0 0 1];

Matrix=A_latdir;
Deflection=B_latdir;

sys_beta=ss(Matrix,Deflection,C1,D);
sys_p=ss(Matrix,Deflection,C2,D);
sys_r=ss(Matrix,Deflection,C3,D);
sys_phi=ss(Matrix,Deflection,C4,D);
sys_psi=ss(Matrix,Deflection,C5,D);

% Time period ploting Dutch Roll Mode
tf1 = OUTPUT_read_XLSX.Stability_flags.tf1_lat;

% Perturbation in side slip angle
Dbeta = OUTPUT_read_XLSX.Stability_flags.Dbeta;
% Perturbation in roll rate
Dp = OUTPUT_read_XLSX.Stability_flags.Dp;
% Perturbation in yaw rate
Dr = OUTPUT_read_XLSX.Stability_flags.Dr;
% Perturbation in bank angle
Dphi = OUTPUT_read_XLSX.Stability_flags.Dphi;

%% Figures
Fig = Fig + 1;
figure(Fig)
x0=[Dbeta 0 0 0 0];

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*R2D,'LineWidth', LS)
title(['Lateral-Directional Response for the airplane to a change in \beta of ',num2str(Dbeta*R2D) 'deg'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*R2D,'LineWidth', LS)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*R2D,'LineWidth', LS)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*R2D,'LineWidth', LS)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*R2D,'LineWidth', LS)
ylabel('\psi, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_latdir_pert_Dbeta');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end

Fig = Fig + 1;
figure(Fig)
x0=[0 Dp 0 0 0];

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*R2D,'LineWidth', LS)
title(['Lateral-Directional Response for the airplane to a change in p of ',num2str(Dp*R2D) 'deg/s'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*R2D,'LineWidth', LS)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*R2D,'LineWidth', LS)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*R2D,'LineWidth', LS)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*R2D,'LineWidth', LS)
ylabel('\psi, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_latdir_pert_Dp');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end


Fig = Fig + 1;
figure(Fig)
x0=[0 0 Dr 0 0];
% tf1=20;

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*R2D,'LineWidth', LS)
title(['Lateral-Directional Response for the airplane to a change in r of ',num2str(Dr*R2D) 'deg/s'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*R2D,'LineWidth', LS)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*R2D,'LineWidth', LS)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*R2D,'LineWidth', LS)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*R2D,'LineWidth', LS)
ylabel('\psi, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_latdir_pert_Dr');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end

Fig = Fig + 1;
figure(Fig)
x0=[0 0 0 Dphi 0];
% tf1=20;

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*R2D,'LineWidth', LS)
title(['Lateral-Directional Response for the airplane to a change in \phi of ',num2str(Dphi*R2D) 'deg'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*R2D,'LineWidth', LS)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*R2D,'LineWidth', LS)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*R2D,'LineWidth', LS)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*R2D,'LineWidth', LS)
ylabel('\psi, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_latdir_pert_Dphi');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end

Fig = Fig + 1;
figure(Fig)
Dbeta=10*pi/180; 
Dp=20*pi/180; 
x0=[Dbeta Dp 0 0 0];
% tf1=20;

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*R2D,'LineWidth', LS)
title(['Lateral-Directional Response for the airplane to a change in \beta of ',num2str(Dbeta*R2D) 'deg and p of', num2str(Dp*R2D),'deg/s'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*R2D,'LineWidth', LS)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*R2D,'LineWidth', LS)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*R2D,'LineWidth', LS)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*R2D,'LineWidth', LS)
ylabel('\psi, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_latdir_pert_DbetaDp');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    % Unified plot Options
    SAVE_types(fname,name,gca,gcf);
end

