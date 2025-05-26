function [F Mm DynVar] = ForcesMoments(x,h,V,Data_ATM,Wind,Data_Der,delta,Propulsion,Geo_tier)

%--------------------------States & Control--------------------------------
V       = x(1);
h       = x(2);

alpha    = x(3);
beta     = x(4);
gamma    = x(5);

alpha_dot = 0;
beta_dot = 0;

phi      = x(6);
theta    = x(7);
psi      = x(8);

p        = x(9);
q        = x(10);
r        = x(11);

x_p      = x(12);
y_p      = x(13);
z_p      = x(14);

u_p      = x(15);
v_p      = x(16);
w_p      = x(17);

delta_e = delta.delta_e;
delta_T = delta.delta_T;
delta_a = delta.delta_a;
delta_r = delta.delta_r;

%--------------------------Dynamic model--------------------------------
% Velocities considering the wind component
u = V*cos(alpha)*cos(beta);
v = V*sin(beta);
w = V*sin(alpha)*cos(beta);

%--------------------------Atmospheric variable--------------------------------
% Data_ATM = get_Atmospheric_Cefiro(h);
rho = Data_ATM.rho;
Q = 0.5*rho*V^2;
Wx = Wind.Wx;
Wy = Wind.Wy;
Wz = Wind.Wz;

%--------------------------Derivatives--------------------------------
CL_0        = Data_Der.CL0;
CM_0        = Data_Der.CM0;
CD_0        = Data_Der.Cd0;
CL_alpha    = Data_Der.CLalpha;
CM_alpha    = Data_Der.CMalpha;
CL_delta_e  = Data_Der.CLdelta_e;
CM_delta_e  = Data_Der.CMdelta_e;

CX_alpha = Data_Der.Cxalfa;
CZ_alpha = Data_Der.Czalfa;
% CM_alpha = Data_Der.Cmalfa;

CX_u = Data_Der.Cxu;
CZ_u = Data_Der.Czu;
CM_u = Data_Der.Cmu;

CX_q = Data_Der.Cxq;
CZ_q = Data_Der.Czq;
CM_q        = Data_Der.Cmq;

CX_alpha_dot = Data_Der.Cxalfapunto;
CZ_alpha_dot = Data_Der.Czalfapunto;
CM_alpha_dot = Data_Der.Cmalfapunto;

CX_delta_e = Data_Der.Cxdeltae;
CZ_delta_e = Data_Der.Czdeltae;
CM_delta_e = Data_Der.Cmdeltae;

CX_theta = Data_Der.Cxteta;
CZ_theta = Data_Der.Czteta;
CM_theta = Data_Der.Cmteta;

CY_beta = Data_Der.Cyb;
CN_beta = Data_Der.Cnb;
CL_beta = Data_Der.Clb;

CY_delta_a = Data_Der.Cydeltaa;
CN_delta_a = Data_Der.Cndeltaa;
CL_delta_a = Data_Der.Cldeltaa;

CY_delta_r = Data_Der.Cydeltar;
CN_delta_r = Data_Der.Cndeltar;
CL_delta_r = Data_Der.Cldeltar;

CY_p = Data_Der.Cyp;
CL_p = Data_Der.Clp;
CN_p = Data_Der.Cnp;
CY_r = Data_Der.Cyr;
CL_r = Data_Der.Clr;
CN_r = Data_Der.Cnr;

CY_beta_dot = Data_Der.Cybpunto;
CL_beta_dot = Data_Der.Clbpunto;
CN_beta_dot = Data_Der.Cnbpunto;

k1  = Data_Der.k1;
k2  = Data_Der.k2;
b  = Data_Der.b;
S  = Data_Der.S;
c  = Data_Der.c;
m  = Data_Der.m;
g  = Data_Der.g;

Ix =  Data_Der.Ix;
Iy =  Data_Der.Iy;
Iz =  Data_Der.Iz;
Ixz = Data_Der.Ixz;

%--------------------------Constants--------------------------------

V_ref = V;
X_adi = (c/(2*V_ref)); % Adimensinalize longitudinal derivatives
Y_adi = (b/(2*V_ref)); % Adimensinalize lateral-directional derivatives

CL     = CL_0 + CL_alpha*alpha + CL_delta_e*delta_e;
CD     = CD_0 + k1*CL + k2*CL^2;

% sigma = sigmoide(alpha);
% 
% C_L = (1-sigma)*(C_L0 + C_Lalpha*alpha) + sigma*(2*sign(alpha)*sin(alpha)^2*...
% cos(alpha)) + 0.5*C_Lq*c*Q/(Uinf+1e-16) + C_Ldelta_e*delta_e;
% C_D = C_D0 + (1-sigma)*((C_L0 + C_Lalpha*alpha)^2/(pi*e*AR)) + sigma*(2*...
% sign(alpha)*sin(alpha)^3) + 0.5*C_Dq*c*Q/(Uinf+1e-16) + C_Dbeta1*beta + ...
% C_Dbeta2*beta^2 + C_Ddelta_e*delta_e;

%% WARNING!!!
% Assume initial, need to be changed with actual derivatives!!!
CX_alpha = CD;
CZ_alpha = CL_0 + CL_alpha*alpha;

CX = CX_u*u + CX_alpha*alpha + CX_theta*theta + CX_q*X_adi*q + CX_delta_e*delta_e;
CY     = CY_beta*beta + CY_beta_dot*Y_adi*beta_dot + CY_p*Y_adi*p + CY_r*Y_adi*r + CY_delta_a*delta_a + CY_delta_r*delta_r;
CZ = CZ_u*u + CZ_alpha*alpha + CZ_alpha_dot*X_adi*alpha_dot + CZ_theta*theta + CZ_q*X_adi*q + CZ_delta_e*delta_e;

CL_bar = CL_beta*beta + CL_beta_dot*Y_adi*beta_dot + CL_p*Y_adi*p + CL_r*Y_adi*r + CL_delta_a*delta_a + CL_delta_r*delta_r;
CM     = CM_0 + CM_alpha*X_adi*alpha + CM_theta*theta + CM_q*q + CM_delta_e*delta_e;
CN     = CN_beta*beta + CN_beta_dot*Y_adi*beta_dot + CN_p*Y_adi*p + CN_r*Y_adi*r + CN_delta_a*delta_a + CN_delta_r*delta_r;

D = Q*S*CD;
L = Q*S*CL;
Y = Q*S*CY;

Xs = Q*S*CX;
Ys = Y;
Zs = Q*S*CZ;

T_prop = Propulsion.Ti;
Q_prop = Propulsion.Qi;
CT = T_prop/(Q*S);

T = T_prop;

CT_x_1 = CT;
CT_x_u = -3*CT_x_1;
CT_X_alpha = 0;

CT_x_u = CT_x_u + 2*CT_x_1;
CT_x = CT_x_u*(u/V) + CT_X_alpha*alpha;
CT_y = 0;
CT_z = 0;

dT = 0;
CM_t_u = - CT_x_u*(dT/c);
CM_T1 = 0;
CM_T_alpha = 0;
CM_T_u = CM_t_u + 2*CM_T1;
CM_T_x = CM_T_u*u/V + CM_T_alpha*alpha;

CL_T_y = 0;
CN_T_y = 0;

CN_T_beta = 0;
CN_T_z = CN_T_beta*beta;

T_x_s = Q*S*CT_x;
T_y_s = Q*S*CT_y;
T_z_s = Q*S*CT_z;

L_T_s = Q*S*b*CL_T_y;
M_T_s = Q*S*c*CM_T_x;
N_T_s = Q*S*b*CN_T_z;

X_a = cos(alpha)*Xs - sin(alpha)*Zs;
Y_a = Ys;
Z_a = sin(alpha)*Xs + cos(alpha)*Zs;

T_x_a = cos(alpha)*T_x_s - sin(alpha)*T_z_s;
T_y_a = T_y_s;
T_z_a = sin(alpha)*T_x_s + cos(alpha)*T_z_s;

Cc = -CX*sin(beta) + CY*cos(beta);
C = Q*S*Cc;

L_bar = Q*S*b*CL_bar;
M     = Q*S*c*CM;
N     = Q*S*b*CN;

% Conversion de stabilida da cuerpo
L_a = cos(alpha)*L_bar - sin(alpha)*N;
M_a = M;
N_a = sin(alpha)*L_bar + cos(alpha)*N;

L_T_a = cos(alpha)*L_T_s - sin(alpha)*N_T_s;
M_T_a = M_T_s;
N_T_a = sin(alpha)*L_T_s + cos(alpha)*N_T_s;

%--------------------------Forces--------------------------------
F_x = X_a + T_x_a;
F_y = Y_a + T_y_a;
F_z = Z_a + T_z_a;

F.F_x = F_x;
F.F_y = F_y;
F.F_z = F_z;

%--------------------------Moments--------------------------------
% WARNING!!!!
% Need to include the gyroscopic effects and the aerodyamic torque contributions for rotating engines
% Gyroscopic Moments
% Matrix of angular velocity of the system in body axis
% Mw = [0 - r q; r 0 -p; -q p 0];
% % Angular momentum
% h_rotor_x = I_r*omega_x;
% h_rotor_y = I_r*omega_x;
% h_rotor_z = I_r*omega_x;

L_prop_gyro = 0;
M_prop_gyro = 0;
N_prop_gyro = 0;

% Aerodynamic Torque
L_prop_tau = 0;
M_prop_tau = 0;
N_prop_tau = 0;

L_prop = L_prop_tau + L_prop_tau;
M_prop = M_prop_tau + M_prop_tau;
N_prop = N_prop_tau + N_prop_tau;

% Momoments in Body Axis
L_bar = L_a + L_T_a + L_prop;
M = M_a + M_T_a + M_prop;
N = N_a + N_T_a + N_prop;

Mm_x = L_bar;
Mm_y = M;
Mm_z = N;

Mm.Mm_x = Mm_x;
Mm.Mm_y = Mm_y;
Mm.Mm_z = Mm_z;

DynVar.F = F;
DynVar.Mm = Mm;

%--------------------------Saving in structure--------------------------------
DynVar.u = u;
DynVar.v = v;
DynVar.w = w;

DynVar.CL     = CL;
DynVar.CD     = CD;
DynVar.CY     = CY;

DynVar.CX = CX;
DynVar.CZ = CZ;

DynVar.CL_bar = CL_bar;
DynVar.CM     = CM;
DynVar.CN     = CN;

DynVar.D = D;
DynVar.L = L;
DynVar.Y = Y;

DynVar.Xs = Xs;
DynVar.Ys = Ys;
DynVar.Zs = Zs;

DynVar.CT = CT;

DynVar.T = T;

DynVar.CTx = CT_x;
DynVar.CTy = CT_y;
DynVar.CTz = CT_z;

DynVar.dT = dT;
DynVar.CMTx = CM_T_x;
DynVar.CLTy = CL_T_y;
DynVar.CNTz = CN_T_z;

DynVar.Txs = T_x_s;
DynVar.Tys = T_y_s;
DynVar.Tzs = T_z_s;

DynVar.MTs = M_T_s;
DynVar.LTs = L_T_s;
DynVar.NTs = N_T_s;

DynVar.MTa = M_T_a;
DynVar.LTa = L_T_a;
DynVar.NTa = N_T_a;

DynVar.Xa = X_a;
DynVar.Ya = Y_a;
DynVar.Za = Z_a;

DynVar.Txa = T_x_a;
DynVar.Tya = T_y_a;
DynVar.Tza = T_z_a;

DynVar.Ma = M_a;
DynVar.La = L_a;
DynVar.Na = N_a;

DynVar.Cc = Cc;
DynVar.C = C;

DynVar.L_bar = L_bar;
DynVar.M     = M;
DynVar.N     = N;

% DynVar.udot = udot;
% DynVar.vdot = vdot;
% DynVar.wdot = wdot;