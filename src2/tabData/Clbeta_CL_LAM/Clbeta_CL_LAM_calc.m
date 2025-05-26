%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   28 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de Clbeta_CL_LAM: Wing sweep contribution. AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.2 (pag 393, 2083 PDF))
function Clbeta_CL_LAM = Clbeta_CL_LAM_calc(lambda, AR, LAMc2)
nPoint = 23; % Numero de puntos por curva
par_pos = [0 0.5 1];
[~, k] = min(abs(lambda - par_pos));

switch k
    case 1
        data = load('Clbeta_CL_0.dat');
        par2_pos = [1 2 4 6 8];
        [~, kk] = min(abs(AR - par2_pos));
        kk_i = (kk-1)*nPoint+1;
        kk_f = kk*nPoint;
        data = data(kk_i:kk_f,:);
    case 2
        data = load('Clbeta_CL_05.dat');
        par2_pos = [1 2 4 6 8];
        [~, kk] = min(abs(AR - par2_pos));
        kk_i = (kk-1)*nPoint+1;
        kk_f = kk*nPoint;
        data = data(kk_i:kk_f,:);
    case 3
        data = load('Clbeta_CL_1.dat');
        par2_pos = [1 1.5 2 4 6];
        [~, kk] = min(abs(AR - par2_pos));
        kk_i = (kk-1)*nPoint+1;
        kk_f = kk*nPoint;
        data = data(kk_i:kk_f,:);
end
Clbeta_CL_LAM = interp1(data(:,1),data(:,2),LAMc2,'pchip');
end