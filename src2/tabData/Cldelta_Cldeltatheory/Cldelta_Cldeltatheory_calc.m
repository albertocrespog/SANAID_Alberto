
% Calculo de Cldelta_Cldeltatheory: Correction Factor for Plain Flap Lift. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.15 (pag 230, 1920 PDF))

function Cldelta_Cldeltatheory = Cldelta_Cldeltatheory_calc(cf_c,Cla)

Cla_Clatheory_pos = 0.7:0.02:1;
% Cla_Clatheory = Cla/2/pi;  
Cla_Clatheory = Cla;
[~,kk] = min(abs(Cla_Clatheory - Cla_Clatheory_pos));

nPuntos = 4;
kk_i = (kk-1)*nPuntos+1;
kk_f = kk_i + (nPuntos - 1);
data = load('Cldelta_Cldeltatheory.dat');
data = data(kk_i:kk_f,:);

Cldelta_Cldeltatheory = interp1(data(:,1),data(:,2),cf_c,'pchip');

end