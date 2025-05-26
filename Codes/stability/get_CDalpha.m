function Stab_Der = get_CDalpha(Aero, conditions, Stab_Der, conv_UNITS, Performance, Geo_tier)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CDalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selects the source of the polar model
C_D0 = Aero.Polar.C_D0;
C_D1 = Aero.Polar.C_D1;
C_D2 = Aero.Polar.C_D2;

n = conditions.n;
m_TOW = conditions.m_TOW;
g = conv_UNITS.g;
V = conditions.V;
S_ref = Geo_tier.S_ref;
% n = 3.8;
% n = 1;
phi_n = n-1;
q = (g/V)*phi_n;
rho = Performance.rho;
q_inf = 0.5*rho*V^2;
q_ref = (2*V/0.469);

% Pitch Derivatives for Maneu ver static stability
CLq = Stab_Der.CLq;

w_T0 = m_TOW*g;
CL_needed = n*w_T0/(q_inf*S_ref) - CLq*q/q_ref;
CL_alpha_ac = Stab_Der.CL_alpha_ac;

CD_alpha = (C_D1 + 2*C_D2*CL_needed)*CL_alpha_ac;

Stab_Der.CD_alpha = CD_alpha;
end
