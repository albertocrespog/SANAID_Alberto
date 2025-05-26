function Generate_Plots_laterall_dynamic

saved_file = 'Lateral.mat';
vars = whos('-file', saved_file);
load(saved_file, vars(1:2).name)
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

figure(1)
Dbeta=10*pi/180; 
x0=[Dbeta 0 0 0 0];
tf1=20;

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*180/pi)
title(['Lateral-Directional Response for the airplane to a change in \beta of ',num2str(Dbeta*180/pi) 'deg'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*180/pi)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*180/pi)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*180/pi)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*180/pi)
ylabel('\psi, deg/s')
grid on


figure(2)
Dp=20*pi/180; 
x0=[0 Dp 0 0 0];
tf1=20;

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*180/pi)
title(['Lateral-Directional Response for the airplane to a change in p of ',num2str(Dp*180/pi) 'deg/s'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*180/pi)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*180/pi)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*180/pi)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*180/pi)
ylabel('\psi, deg/s')
grid on

figure(3)
Dr=20*pi/180; 
x0=[0 0 Dr 0 0];
tf1=20;

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*180/pi)
title(['Lateral-Directional Response for the airplane to a change in r of ',num2str(Dr*180/pi) 'deg/s'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*180/pi)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*180/pi)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*180/pi)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*180/pi)
ylabel('\psi, deg/s')
grid on

figure(4)
Dphi=20*pi/180; 
x0=[0 0 0 Dphi 0];
tf1=20;

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*180/pi)
title(['Lateral-Directional Response for the airplane to a change in \phi of ',num2str(Dphi*180/pi) 'deg'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*180/pi)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*180/pi)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*180/pi)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*180/pi)
ylabel('\psi, deg/s')
grid on


figure(5)
Dbeta=10*pi/180; 
Dp=20*pi/180; 
x0=[Dbeta Dp 0 0 0];
tf1=20;

subplot(5,1,1)
[y_beta,t_beta,x_beta]=initial(sys_beta,x0,tf1);
plot(t_beta,y_beta*180/pi)
title(['Lateral-Directional Response for the airplane to a change in \beta of ',num2str(Dbeta*180/pi) 'deg and p of', num2str(Dp*180/pi),'deg/s'])
ylabel('\Delta \beta, deg')
grid on

subplot(5,1,2)
[y_p,t_p,x_p]=initial(sys_p,x0,tf1);
plot(t_p,y_p*180/pi)
ylabel('\Delta p, deg/s')
grid on

subplot(5,1,3)
[y_r,t_r,x_r]=initial(sys_r,x0,tf1);
plot(t_r,y_r*180/pi)
ylabel('\Delta r deg/s')
grid on

subplot(5,1,4)
[y_phi,t_phi,x_phi]=initial(sys_phi,x0,tf1);
plot(t_phi,y_phi*180/pi)
ylabel('\phi, deg/s')
grid on

subplot(5,1,5)
[y_psi,t_psi,x_psi]=initial(sys_psi,x0,tf1);
plot(t_psi,y_psi*180/pi)
ylabel('\psi, deg/s')
grid on
