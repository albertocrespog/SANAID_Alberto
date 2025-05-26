
% Cálculo preliminar de la derivada de estabilidad C_n_dr.
function Stab_Der = getCndr(modelo, alpha, Stab_Der)
%% DATOS REQUERIDOS DEL MODELO
y0_b2   = modelo.vertical.y0_b2;
y1_b2   = modelo.vertical.y1_b2;
AR      = modelo.vertical.AR;
TR      = modelo.vertical.TR;
cf_c    = modelo.vertical.cm_c;
t_c     = modelo.vertical.t_c;
CLa     = modelo.vertical.CLa;
Cla     = modelo.vertical.Cla;
S_v     = modelo.vertical.S;
z_v     = modelo.vertical.Zca;
Zca     = modelo.vertical.Zca;
Xca     = modelo.vertical.Xca;
x_XCG   = modelo.general.Xcg;
z_XCG   = modelo.general.Zcg;
l_v     = modelo.vertical.Xca - modelo.general.Xcg;
rend    = modelo.vertical.eta;

Sref    = modelo.general.Sref;
b       = modelo.ala.b;


%% CÁLCULO DE Cy_dr (pag 461, 2151 PDF)

k_prima                 = 1;
alphaCL_alphaCl         = alphaCL_alphaCl_calc(cf_c,AR);
Cldelta_Cldeltatheory   = Cldelta_Cldeltatheory_calc(cf_c,Cla);
Cldelta_theory          = Cldeltatheory_calc(cf_c,t_c);
K_b                     = Kb_calc(y0_b2, y1_b2, TR);

Cy_dr   = K_b*CLa*S_v/Sref*rend*k_prima/Cla*alphaCL_alphaCl*Cldelta_Cldeltatheory*Cldelta_theory;
if isnan(Cy_dr)
    Cy_dr = 0;
end

%% CÁLCULO DE LA DERIVADA DE CONTROL
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 462, 2152 PDF)

Cn_dr   = -((Zca - z_XCG)*sin(alpha)+ (Xca - x_XCG)*cos(alpha))/b*Cy_dr;

Stab_Der.Cndeltar = Cn_dr;
end