
% Cálculo preliminar de la derivada de estabilidad C_n_da.

function Stab_Der = getCnda(modelo,Stab_Der)
%% DATOS REQUERIDOS DEL MODELO
Minf    = modelo.general.Minf;
W       = modelo.general.mtow*modelo.general.w_w0*9.8065;    
qinf    = modelo.general.qinf;
Sref    = modelo.general.Sref;

yf_b    = modelo.ala.y1_b2;
y0_b    = modelo.ala.y0_b2;
b       = modelo.ala.b;
AR      = modelo.ala.AR;
TR      = modelo.ala.TR;
LAMc4   = modelo.ala.LAMc4;
ca_c    = modelo.ala.cm_c;
t_c     = modelo.ala.t_c;
Cla     = modelo.ala.Cla;
CLa     = modelo.ala.CLa;

CL      = W/qinf/Sref; 




%% CALCULO DE Cn_da
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 448, 2138 PDF)
Ka      = Ka_calc(yf_b, AR, TR);

Cn_da   = 2*Ka*CL*Cl_da;

Stab_Der.Cndeltaa = Cn_da;
end
