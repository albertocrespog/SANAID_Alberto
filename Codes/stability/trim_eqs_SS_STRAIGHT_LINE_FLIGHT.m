function F = trim_eqs_SS_STRAIGHT_LINE_FLIGHT(x,CD_delta,CL_delta,CM_delta,AC_CONFIGURATION,Aero_TH,Aero,Prop_data,Geo_tier, conditions, Performance, conv_UNITS, OUTPUT_read_XLSX,Stab_Der_parts,Stab_Der,Weight_tier)

phiT = 0;
psiT = 0;
gamma = conditions.gamma; 
phi = conditions.phi;

z_d_T = Geo_tier.z_d_T;
x_d_T = Geo_tier.x_d_T; % Positive for engine behind Xcg : x_d_T = x_eng_xbar - x_XCG; 
% y_d_T = Geo_tier.y_d_T;

S_ref = Geo_tier.S_ref;
V = conditions.V;
rho = Performance.rho;
q_inf = 0.5*rho*V^2;
cmac_w1 = Geo_tier.cmac_w1;

m_TOW = conditions.m_TOW;
g = conv_UNITS.g;

C_D0 = Aero.Polar.C_D0;
CL0_ac = Stab_Der_parts.CL0_ac;
CM0_ac = Stab_Der.CM0_ac;

CD_alpha_ac = Stab_Der.CD_alpha;
CL_alpha_ac = Stab_Der_parts.CL_alpha_ac;
CM_alpha_ac = Stab_Der.CM_alpha_ac;

% Variables independientes
alpha = x(1); % Ángulo de ataque
delta = x(2); % Deflexión de control
T = x(3);     % Empuje
delta_T = x(4); % Posición de la palanca de empuje

% Obtener propiedades del motor (T es función de deltaT)
[Propulsion] = get_propulsion_TRIM(AC_CONFIGURATION,Aero_TH,Aero,Prop_data,Geo_tier, conditions, Performance, conv_UNITS, OUTPUT_read_XLSX,T);

deltaT_computed = Propulsion.delta_T;

% Ecuaciones del sistema
% F(1) = CL_needed - CL0_ac - CL_alpha_ac*alpha - CL_delta*delta - T*sin(alpha);
F(1) = -m_TOW*g*cos(gamma)*cos(phi) + (CL0_ac + CL_alpha_ac*alpha + CL_delta*delta)*q_inf*S_ref - T*(-cos(phiT)*cos(psiT)*sin(alpha) + sin(phiT)*cos(alpha));
% F(2) = CM0_needed + CM_alpha_ac*alpha + CM_delta*delta + T*z_d_T;
F(2) = (CM0_ac + CM_alpha_ac*alpha + CM_delta*delta)*q_inf*S_ref*cmac_w1 + T*(cos(phiT)*cos(psiT)*z_d_T - sin(phiT)*x_d_T);
% F(3) = T - (CD0_ac + CD_alpha_ac*alpha + CD_delta*delta)/cos(alpha);
F(3) = m_TOW*g*sin(gamma) + (C_D0 + CD_alpha_ac*alpha + CD_delta*delta)*q_inf*S_ref - T*(cos(phiT)*cos(psiT)*cos(alpha) + sin(phiT)*sin(alpha));
F(4) = delta_T - deltaT_computed; % Relación entre deltaT y T

end