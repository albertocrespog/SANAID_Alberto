function Plots_longitudinal

saved_file = 'Longitudinal.mat';
vars = whos('-file', saved_file);
load(saved_file, vars(1:4).name)
A_lon = Matrix_lon;
B_lon = Deflection_lon;
alpha = theta1;
U1 = u1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%RESPUESTAS LONGITUDINALES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=0;
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

tf1=50;
tf2=0.1;

figure(1)
Du=u1*0.20; % change in speed m/s
x0=[Du 0 0 0];
% tf1=100;
% tf2=0.1;

subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed)
title(['Longitudinal Response for the airplane to a change in Speed of ',num2str(Du) 'm/s'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*180/pi)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*180/pi)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*180/pi)
ylabel('q, deg/s')
grid on

figure(2)
Dalpha=10*pi/180; 
x0=[0 Dalpha 0 0];
% tf1=100;
% tf2=0.1;

subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed)
title(['Longitudinal Response for the airplane to a change in Angle of Attack of ',num2str(Dalpha*180/pi) '\circ'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*180/pi)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*180/pi)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*180/pi)
ylabel('q, deg/s')
grid on


figure(3)
x0=[-Du Dalpha 0 0];
% tf1=100;
% tf2=0.1;

subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed)
title(['Longitudinal Response for the airplane to a \Deltau=',num2str(-Du),'m/s and a \Delta\alpha=',num2str(Dalpha*180/pi),'\circ'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*180/pi)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*180/pi)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*180/pi)
ylabel('q, deg/s')
grid on


figure(4)
Dq=10*pi/180; % change in speed m/s
x0=[0 0 Dq 0];
% tf1=90;
% tf2=0.1;

subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed)
title(['Longitudinal Response for the airplane to a change in Pitch Rate of ',num2str(Dq*180/pi) '\circ/sec'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*180/pi)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*180/pi)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*180/pi)
ylabel('q, deg/s')
grid on

figure(5)
Dtheta=3*pi/180;
x0=[0 0 0 Dtheta];
% tf1=100;
% tf2=0.1;

subplot(4,1,1)
[y_speed,t_speed,x_speed]=initial(sys_speed,x0,tf1);
plot(t_speed,y_speed)
title(['Longitudinal Response for the airplane to a change in Pitch Angle of ',num2str(Dq*180/pi) '\circ'])
ylabel('\Deltau, m/s')
grid on

subplot(4,1,2)
[y_alpha,t_alpha,x_alpha]=initial(sys_alpha,x0,tf1);
plot(t_alpha,y_alpha*180/pi)
ylabel('\Delta\alpha, deg')
grid on

subplot(4,1,3)
[y_theta,t_theta,x_theta]=initial(sys_theta,x0,tf1);
plot(t_theta,y_theta*180/pi)
ylabel('\Delta\theta, deg')
grid on

subplot(4,1,4)
[y_q,t_q,x_q]=initial(sys_q,x0,tf1);
plot(t_q,y_q*180/pi)
ylabel('q, deg/s')
grid on

