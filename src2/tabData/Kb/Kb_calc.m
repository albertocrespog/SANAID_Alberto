%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   14 de Agosto de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de Kb: Effect of the taper ratio and flap span. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.4.2, Fig. 8.52 (pag 260, 1950 PDF))
function Kb = Kb_calc(y1_b2, y2_b2, TR)

data = load('Kb.dat');

TR_pos = [0,0.5,1];
[~,kk] = min(abs(TR-TR_pos));
nPuntos = 5;
kk_i = (kk-1)*nPuntos+1;
kk_f = kk_i + (nPuntos - 1);
data = data(kk_i:kk_f,:);
Kb = interp1(data(:,1),data(:,2),y2_b2,'pchip') - interp1(data(:,1),data(:,2),y1_b2,'pchip');

end