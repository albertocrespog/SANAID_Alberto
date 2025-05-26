function [Fig] = Generate_Plots_longitudinal_dynamic(Plot_Options,OUTPUT_read_XLSX,conv_UNITS,Fig, Stab_Dyn_Long,filenameS)

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
fname = filenameS.plots;ts;

% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
MATLAB_in = Plot_Options.MATLAB_in;

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

%% Loads data
% saved_file = 'Longitudinal.mat';
% vars = whos('-file', saved_file);
% load(saved_file, vars(1:4).name)
% A_lon = Matrix_lon;
% B_lon = Deflection_lon;
% alpha = theta1;
% U1 = u1;
Matrix_lon     = Stab_Dyn_Long.Matrix_lon;
Deflection_lon = Stab_Dyn_Long.Deflection_lon;
u1             = Stab_Dyn_Long.u1;

A_lon = Matrix_lon;
B_lon = Deflection_lon;
U1 = u1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%RESPUESTAS LONGITUDINALES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 0;
C1=[1 0 0 0];
C2=[0 1 0 0];
C3=[0 0 1 0];
C4=[0 0 0 1];

Matrix=A_lon;
Deflection=B_lon;

sys_speed=ss(Matrix,Deflection,C1,D);
sys_alpha=ss(Matrix,Deflection,C2,D);
sys_q=ss(Matrix,Deflection,C3,D);
sys_theta=ss(Matrix,Deflection,C4,D);

% Time period ploting Phugoid Mode
tf1 = OUTPUT_read_XLSX.Stability_flags.tf1_long;
% Time period ploting Short Period Mode
tf2 =OUTPUT_read_XLSX.Stability_flags.tf2_long;

% Determine the perturbations
% Perturbation in forward speed velocity (percentage of trim velocity)
Du = OUTPUT_read_XLSX.Stability_flags.Du;
Du = Du*U1;
% Perturbation in angle of attack
Dalpha = OUTPUT_read_XLSX.Stability_flags.Dalpha;
% Perturbation in pitch rate
Dq = OUTPUT_read_XLSX.Stability_flags.Dq;
% Perturbation in pitch angle
Dtheta = OUTPUT_read_XLSX.Stability_flags.Dtheta;

%% Figures
Fig = Fig + 1;
figure(Fig)
x0=[Du 0 0 0];
subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed,'LineWidth', LS)
title(['Longitudinal Response for the airplane to a change in Speed of ',num2str(Du) 'm/s'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*R2D,'LineWidth', LS)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*R2D,'LineWidth', LS)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*R2D,'LineWidth', LS)
ylabel('q, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_long_pert_Du');
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
x0=[0 Dalpha 0 0];
subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed,'LineWidth', LS)
title(['Longitudinal Response for the airplane to a change in Angle of Attack of ',num2str(Dalpha*180/pi) '\circ'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*R2D,'LineWidth', LS)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*R2D,'LineWidth', LS)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*R2D,'LineWidth', LS)
ylabel('q, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_long_pert_Dalpha');
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
x0=[-Du Dalpha 0 0];
subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed,'LineWidth', LS)
title(['Longitudinal Response for the airplane to a \Deltau=',num2str(-Du),'m/s and a \Delta\alpha=',num2str(Dalpha*180/pi),'\circ'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*R2D,'LineWidth', LS)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*R2D,'LineWidth', LS)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*R2D,'LineWidth', LS)
ylabel('q, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_long_pert_DuDalpha');
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
Dq=10*pi/180; % change in speed m/s
x0=[0 0 Dq 0];
subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed,'LineWidth', LS)
title(['Longitudinal Response for the airplane to a change in Pitch Rate of ',num2str(Dq*180/pi) '\circ/sec'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*R2D,'LineWidth', LS)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*R2D,'LineWidth', LS)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*R2D,'LineWidth', LS)
ylabel('q, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_long_pert_Dq');
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
Dtheta=3*pi/180;
x0=[0 0 0 Dtheta];

subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed,'LineWidth', LS)
title(['Longitudinal Response for the airplane to a change in Pitch Angle of ',num2str(Dtheta*180/pi) '\circ'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*R2D,'LineWidth', LS)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*R2D,'LineWidth', LS)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*R2D,'LineWidth', LS)
ylabel('q, deg/s')
grid on
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('plot_long_pert_Dtheta');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end
