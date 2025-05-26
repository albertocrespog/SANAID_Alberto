%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   28 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de Clbeta_CL_AR: Wing aspect ratio contribution AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.23 (pag 394, 2084 PDF))
function Clbeta_CL_AR = Clbeta_CL_AR_calc(AR, lambda)
nPoint = 11; % Numero de puntos por curva
par_pos = [0 0.5 1];
[~, kk] = min(abs(lambda - par_pos));
kk_i = (kk-1)*nPoint+1;
kk_f = kk*nPoint;
data = load('Clbeta_CL_A.dat');
data = data(kk_i:kk_f,:);
if AR > data(end,1)
    AR = data(end,1);
elseif AR < data(1,1)
    AR = data(1,1);
end
Clbeta_CL_AR = interp1(data(:,1),data(:,2),AR,'pchip');
end