%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   28 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de K_MLAM: Compressibility correction to sweep. AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.21 (pag 394, 2084 PDF))

function K_MLAM = K_MLAM_calc(M, AR, LAMc2)
nPoint = 27; % Numero de puntos por curva
par = AR/cos(LAMc2);
par_pos = [10 8 6 5 4 3 2];
[~, kk] = min(abs(par - par_pos));
kk_i = (kk-1)*nPoint+1;
kk_f = kk_i + (nPoint-1);
data = load('K_MLAM.dat');
data = data(kk_i:kk_f,:);
K_MLAM = interp1(data(:,1),data(:,2),M*cos(LAMc2),'pchip');

end