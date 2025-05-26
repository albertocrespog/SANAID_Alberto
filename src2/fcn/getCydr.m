
% Cálculo preliminar de la derivada de estabilidad C_y_dr.
function Stab_Der = getCydr(modelo,Stab_Der)
%% DATOS REQUERIDOS DEL MODELO
y0_b2   = modelo.vertical.y0_b2;
y1_b2   = modelo.vertical.y1_b2;
AR      = modelo.vertical.AR;
TR      = modelo.vertical.TR;
cf_c    = modelo.vertical.cm_c;
t_c     = modelo.vertical.t_c;
Cla     = modelo.vertical.Cla;
CLa     = modelo.vertical.CLa;
S_v     = modelo.vertical.S;

Sref    = modelo.general.Sref;
rend    = 0.95;

%% CÁLCULO DE LA DERIVADA DE ESTABILIDAD
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)
%% CÁLCULO DE LA DERIVADA DE ESTABILIDAD
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)
k_prima                 = 1;
% The three dimensional ruddervator effectiveness parameter is determined from Figure 8.53 in Airplane Design Part VI 
% and is a function of the V-tail aspect ratio, and ruddervator chord to V-tail chord ratio:
alphaCL_alphaCl         = alphaCL_alphaCl_calc(cf_c,AR);
% Calculo de Cldelta_Cldeltatheory: Correction Factor for Plain Flap Lift. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.15 (pag 230, 1920 PDF))
Cldelta_Cldeltatheory   = Cldelta_Cldeltatheory_calc(cf_c,Cla);
% Calculo de Cldeltatheory: Lift Effectiveness. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.14 (pag 228, 1918 PDF)
Cldelta_theory          = Cldeltatheory_calc(cf_c,t_c);
% Calculo de Kb: Effect of the taper ratio and flap span. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.4.2, Fig. 8.52 (pag 260, 1950 PDF))
K_b                     = Kb_calc(y0_b2, y1_b2, TR);

% The change in sideslip due to ruddervator directional deflection 
beta_delta_r = -K_b*Cldelta_Cldeltatheory*Cldelta_theory*(k_prima/Cla)*alphaCL_alphaCl;

Cy_dr   = -rend*CLa*(S_v/Sref)*beta_delta_r;
% Cy_dr   = K_b*CLa*S_v/Sref*rend*k_prima/Cla*alphaCL_alphaCl*Cldelta_Cldeltatheory*Cldelta_theory

if isnan(Cy_dr)
    Cy_dr = 0;
end
Stab_Der.Cydeltar = Cy_dr;
end