%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   28 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de Clbeta_tau: Wing dihedral angle effect. AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 395, 2085 PDF)

function Clbeta_die = Clbeta_die_calc(lambda, AR, LAMc2)
nPoint = 10; % Numero de puntos por curva
par_pos = [0 0.5 1];
[~, k] = min(abs(lambda - par_pos));
par_pos2 = [0 40 60]*pi/180;
[~,kk] = min(abs(abs(LAMc2) - par_pos2));
kk_i = (kk-1)*nPoint+1;
kk_f = kk*nPoint;

switch k
    case 1
        data = load('Clbeta_tau_0.dat');
        data = data(kk_i:kk_f,:);
    case 2
        data = load('Clbeta_tau_05.dat');
        data = data(kk_i:kk_f,:);
    case 3
        data = load('Clbeta_tau_1.dat');
        data = data(kk_i:kk_f,:);
end
if AR>8
    AR = 8;
end
Clbeta_die = interp1(data(:,1),data(:,2),AR,'pchip');