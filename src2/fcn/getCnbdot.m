function Stab_Der = getCnbdot(modelo,alpha,Stab_Der)
b_w       = modelo.ala.b;

Zca_v     = modelo.vertical.Zca;
Xca_v     = modelo.vertical.Xca;

Xcg       = modelo.general.Xcg;

Cy_betaDot   = Stab_Der.Cybpunto;

%% Cn_betaDot

Cn_betaDot          = -Cy_betaDot*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))/b_w;

if isnan(Cn_betaDot)
    Cn_betaDot = 0;
end
Stab_Der.Cnbpunto = Cn_betaDot;

end