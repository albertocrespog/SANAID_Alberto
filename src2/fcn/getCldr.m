
% Cálculo preliminar de la derivada de estabilidad C_l_dr.
function Stab_Der = getCldr(modelo, alpha, Stab_Der)
%% DATOS REQUERIDOS DEL MODELO
b       = modelo.ala.b;
y0_b2   = modelo.vertical.y0_b2;
y1_b2   = modelo.vertical.y1_b2;
AR      = modelo.vertical.AR;
TR      = modelo.vertical.TR;
cf_c    = modelo.vertical.cm_c;
t_c     = modelo.vertical.t_c;
CLa     = modelo.vertical.CLa;
Cla     = modelo.vertical.Cla;
S_v     = modelo.vertical.S;
Zca     = modelo.vertical.Zca;
Xca     = modelo.vertical.Xca;
x_XCG = modelo.general.Xcg;
z_XCG = modelo.general.Zcg;
% l_v     = modelo.vertical.Xca - modelo.general.Xcg;


Sref    = modelo.general.Sref;
rend    = 0.95;
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
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)

Cl_dr   = ((Zca - z_XCG)*cos(alpha)-(Xca - x_XCG)*sin(alpha))/b*Cy_dr;

Stab_Der.Cldeltar = Cl_dr;
end