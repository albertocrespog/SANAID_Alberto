% Function that generate the states and controls signals
function [x,delta] = Generate_States(Data_Trim,Data_Der,h,V)

x(1) = V;
x(2) = h;

alpha = Data_Trim.alpha_1;
beta  = 0;
gamma = 0;

phi      = 0; % approximation
theta    = gamma - alpha;
psi      = 0;

x(6) = phi;
x(7) = theta;
x(8) = psi;

% Angular velocitis in the x,y and z axis
p        = 0;
q        = 0;
r        = 0;

x(9) = p;
x(10) = q;
x(11) = r;

% Position of the x,y and z axis
x_p      = 0;
y_p      = 0;
z_p      = 2;

x(12) = x_p;
x(13) = y_p;
x(14) = z_p;

% Aerodynamic component of the velocity in the x,y and z axis
u = V*cos(alpha)*cos(beta);
v = V*sin(beta);
w = V*sin(alpha)*cos(beta);

x(15) = u;
x(16) = v;
x(17) = w;

delta_e = Data_Trim.delta_e_1;
delta_T = 1;
delta_a = 0; 
delta_r = 0;

delta.delta_e = delta_e;
delta.delta_T = delta_T;
delta.delta_a = delta_a;
delta.delta_r = delta_r;