function Stab_Der = getClbdot(modelo,alpha,Stab_Der)
b_w       = modelo.ala.b;

Zca_v     = modelo.vertical.Zca;
Xca_v     = modelo.vertical.Xca;

Xcg       = modelo.general.Xcg;

Cy_betaDot   = Stab_Der.Cybpunto;

%% Cl_betaDot
Cl_betaDot          = Cy_betaDot*(Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w;

if isnan(Cl_betaDot)
    Cl_betaDot = 0;
end

Stab_Der.Clbpunto = Cl_betaDot;
end