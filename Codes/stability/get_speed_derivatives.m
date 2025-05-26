function [Stab_Der] = get_speed_derivatives(Stab_Der,conditions,Performance,Aero, conv_UNITS,Geo_tier)
% CL = Stab_Der.CL;
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
CL = CL_needed;
V = conditions.V;
a = Performance.a;
Mach = V/a;
CD = Stab_Der.CD;
C_D1 = Aero.Polar.C_D1;
C_D2 = Aero.Polar.C_D2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPEED DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CLu = (Mach^2*CL)/(1 - Mach^2);                 %para vuelo subsónico 
    CZu = -2*CL - CLu ;
    % Czu = -2*CL - Clu - 2*m1*q0;
%     CDu = 0;
    CDu = (C_D1 + 2*C_D2*CL)*CLu; 
    CXu = -2*CD - CDu;  % Pamadi 4.432
    CMu = 0;                      %0 para vuelo subsónico
    % CMu = 2*CM + CMu;                      %0 para vuelo subsónico
    
    Stab_Der.CDu = CDu;
    Stab_Der.CZu = CZu;
    Stab_Der.CXu = CXu;
    Stab_Der.CLu = CLu;
    Stab_Der.CMu = CMu;
end
