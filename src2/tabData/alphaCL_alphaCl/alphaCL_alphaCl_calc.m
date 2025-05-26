%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   14 de Agosto de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de alphaCL_alphaCl: Effect of the aspect ratio and flap-chord ratio on the 3-d flap effectiveness. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.4.2, Fig. 8.53 (pag 261, 1951 PDF))
function alphaCL_alphaCl = alphaCL_alphaCl_calc(cf_c,AR)

data1 = load('alphaCL_alphaCl_1.dat');
alpha = interp1(data1(:,1),data1(:,2),cf_c,'pchip');
alpha_pos = 0.1:0.1:1;
[~,kk] = min(abs(alpha-alpha_pos));
nPuntos = 7;
kk_i = (kk-1)*nPuntos+1;
kk_f = kk_i + (nPuntos - 1);
data2 = load('alphaCL_alphaCl_2.dat');
data2 = data2(kk_i:kk_f,:);
alphaCL_alphaCl = interp1(data2(:,1),data2(:,2),AR,'pchip');

end