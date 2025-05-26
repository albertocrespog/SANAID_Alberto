
% Calculo de Cldeltatheory: Lift Effectiveness. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.14 (pag 228, 1918 PDF)

function Cldeltatheory = Cldeltatheory_calc(cf_c,t_c)

t_c_pos = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.15];
[~,kk] = min(abs(t_c_pos - t_c));

nPuntos = 9;
kk_i = (kk-1)*nPuntos+1;
kk_f = kk_i + (nPuntos - 1);
data = load('Cl_delta_theory.dat');
data = data(kk_i:kk_f,:);

Cldeltatheory = interp1(data(:,1),data(:,2),cf_c,'pchip');
end