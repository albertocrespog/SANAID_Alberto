function Data_const = get_const_1214(DATA_int,Data_Der,Data_ATM,h_trim,V_trim,n_trim)

R2D=180/pi;
D2R=pi/180;

% Data around the equilibrium point
% Data_ATM = get_Atmospheric_Cefiro(h);
rho = Data_ATM.rho;
a_speed = Data_ATM.a_speed;
M = V_trim/a_speed;
Q = 0.5*rho*V_trim^2;

g  = Data_Der.g;
b = Data_Der.b;
S  = Data_Der.S;
Iy = Data_Der.Iy;
Ix = Data_Der.Ix;
Iz = Data_Der.Iz;
Ixz = Data_Der.Ixz;
m  = Data_Der.m;

% Recalls the values of the derivatives
Cyb        = Data_Der.Cyb;
Clb        = Data_Der.Clb;
Cnb        = Data_Der.Cnb;
Cyr        = Data_Der.Cyr;
Clr        = Data_Der.Clr;
Cnr        = Data_Der.Cnr;
Cyp        = Data_Der.Cyp;
Clp        = Data_Der.Clp;
Cnp        = Data_Der.Cnp;
Cydeltaa   = Data_Der.Cydeltaa;
Cldeltaa   = Data_Der.Cldeltaa;
Cndeltaa   = Data_Der.Cndeltaa;
Cydeltar   = Data_Der.Cydeltar;
Cldeltar   = Data_Der.Cldeltar;
Cndeltar   = Data_Der.Cndeltar;

alpha1 = Data_Der.alpha1;
theta1 = alpha1;

% Load Factor
phi_ref   = acos(1/n_trim);
RT = (V_trim^2)*tan(phi_ref)/g;

Q1 = (g/V_trim)*(n_trim - 1/n_trim);
Q_ref = Q1;

R1 = (g*(sin(phi_ref)))/(V_trim);
R_ref = R1;

psi_dot1 = g*tan(phi_ref)/V_trim;

C_der_ad = b/(2*V_trim);
Cyr1 = C_der_ad*Cyr;
Clr1 = C_der_ad*Clr;
Cnr1 = C_der_ad*Cnr;

Cyp1 = C_der_ad*Cyp;
Clp1 = C_der_ad*Clp*1.5;
Cnp1 = C_der_ad*Cnp;

Data_const.Cyr1 = Cyr1;
Data_const.Clr1 = Clr1;
Data_const.Cnr1 = Cnr1;

Data_const.Cyp1 = Cyp1;
Data_const.Clp1 = Clp1;
Data_const.Cnp1 = Cnp1;

% Adimensionalize moments of inertia 
I1 = (-Ixz^2 + Ix*Iy - Ix^2)/(Ixz^2 - Ix*Iz);
I2 = (Ixz*Iz - Ixz*Iy + Ixz*Ix)/(Ixz^2 - Ix*Iz);
I5 = -(Ixz)/(Ixz^2 - Ix*Iz);
I7 = -(Ix)/(Ixz^2 - Ix*Iz);

I3 = (-Ixz^2 - Iz^2 - Iz*Iy)/(Ix*Iz-Ixz^2);
I4 = (-Ixz*Iy + Ixz*Ix  + Ixz*Iz)/(Ix*Iz-Ixz^2);
I6 = (Iz)/(Ix*Iz-Ixz^2);
I8 = (Ixz)/(Ix*Iz-Ixz^2);

I9 = (Iz - Ix)/Iy;
I10 = Ixz/Iy;

l1 = Q*S*b*Clb;
l2 = Q*S*b*Clp1; 
l3 = Q*S*b*Clr1; 
l4 = Q*S*b*Cldeltaa; 
l5 = Q*S*b*Cldeltar; 

Data_const.l1 = l1;
Data_const.l2 = l2;
Data_const.l3 = l3;
Data_const.l4 = l4;
Data_const.l5 = l5;

n1 = Q*S*b*Cnb;
n2 = Q*S*b*Cnp1; 
n3 = Q*S*b*Cnr1; 
n4 = Q*S*b*Cndeltaa; 
n5 = Q*S*b*Cndeltar; 

Data_const.n1 = n1;
Data_const.n2 = n2;
Data_const.n3 = n3;
Data_const.n4 = n4;
Data_const.n5 = n5;

y1 = Q*S*Cyb/m;
y2 = Q*S*Cyp1/m;
y3 = Q*S*Cyr1/m;
y4 = Q*S*Cydeltar/m;
y5 = Q*S*Cydeltaa/m;

Data_const.y1 = y1;
Data_const.y2 = y2;
Data_const.y3 = y3;
Data_const.y4 = y4;
Data_const.y5 = y5;

Data_const.I1 = I1;
Data_const.I2 = I2;
Data_const.I3 = I3;
Data_const.I4 = I4;
Data_const.I5 = I5;
Data_const.I6 = I6;
Data_const.I7 = I7;
Data_const.I8 = I8;

Data_const.I9 = I9;
Data_const.I10 = I10;

b1 = I5*cos(alpha1) + I7*sin(alpha1);
b2 = -I5*sin(alpha1) + I7*cos(alpha1);
b3 = I6*cos(alpha1) + I8*sin(alpha1);
b4 = -I6*sin(alpha1) + I8*cos(alpha1);

Data_const.b1 = b1;
Data_const.b2 = b2;
Data_const.b3 = b3;
Data_const.b4 = b4;

v1 = V_trim*cos(alpha1);
v2 = V_trim*sin(alpha1);
v3 = Q*S*Cyb/m;
v4 = Q*S*Cyp1/m;
v5 = Q*S*Cyr1/m;
v6 = Q*S*Cydeltar/m;
v7 = Q*S*Cydeltaa/m;
v8 = g*cos(theta1);

Data_const.v1 = v1;
Data_const.v2 = v2;
Data_const.v3 = v3;
Data_const.v4 = v4;
Data_const.v5 = v5;
Data_const.v6 = v6;
Data_const.v7 = v7;
Data_const.v8 = v8;

beta1 = v1/V_trim;
beta2 = v2/V_trim;
beta3 = v3/V_trim;
beta4 = v4/V_trim;
beta5 = v5/V_trim;
beta6 = v6/V_trim;
beta7 = v7/V_trim;
beta8 = v8/V_trim;

Data_const.beta1 = beta1;
Data_const.beta2 = beta2;
Data_const.beta3 = beta3;
Data_const.beta4 = beta4;
Data_const.beta5 = beta5;
Data_const.beta6 = beta6;
Data_const.beta7 = beta7;
Data_const.beta8 = beta8;

r1 = Q*S*b*(b1*Clb + b2*Cnb);
r2 = Q*S*b*(b1*Clp1 + b2*Cnp1) + I1*Q_ref;
r3 = Q*S*b*(b1*Clr1 + b2*Cnr1) + I2*Q_ref;
r4 = Q*S*b*(b1*Cldeltaa + b2*Cndeltaa);
r5 = Q*S*b*(b1*Cldeltar + b2*Cndeltar);

Data_const.r1 = r1;
Data_const.r2 = r2;
Data_const.r3 = r3;
Data_const.r4 = r4;
Data_const.r5 = r5;

p1 = Q*S*b*(b3*Clb + b4*Cnb);
p2 = Q*S*b*(b3*Clp1 + b4*Cnp1) + I3*Q_ref;
p3 = Q*S*b*(b3*Clr1 + b4*Cnr1) + I4*Q_ref;
p4 = Q*S*b*(b3*Cldeltaa + b4*Cndeltaa);
p5 = Q*S*b*(b3*Cldeltar + b4*Cndeltar);

Data_const.p1 = p1;
Data_const.p2 = p2;
Data_const.p3 = p3;
Data_const.p4 = p4;
Data_const.p5 = p5;

d1 = -(r5*p1 - r1*p5)/(p4*r5 - p5*r4);
d2 = -(r5*p3 - r3*p5)/(p4*r5 - p5*r4);
d3 = -(r4*p1 - r1*p4)/(p4*r5 - p5*r4);
d4 = -(r4*p3 - r3*p4)/(p4*r5 - p5*r4);

Data_const.d1 = d1;
Data_const.d2 = d2;
Data_const.d3 = d3;
Data_const.d4 = d4;