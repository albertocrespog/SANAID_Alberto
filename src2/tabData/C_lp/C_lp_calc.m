%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   12 de Agosto de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de C_lp: Gráfica tema 8 diapositiva 102 Cálculo de aeronaves

function Clp = C_lp_calc(lam,AR,LAMc4,Minf)
beta = sqrt(1-Minf^2);
lambda_pos = [0, 0.25, 0.5, 1];
[~,ii] = min(abs(lam-lambda_pos));
ARcorr = beta*AR;
ARcorr_pos = [1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10];
[~,ii2] = min(abs(ARcorr - ARcorr_pos));

switch ii
    case 1
        nPuntos = 9;
        data = load('Clp_0.dat');
    case 2
        nPuntos = 10;
        data = load('Clp_025.dat');
    case 3
        nPuntos = 10;
        data = load('Clp_050.dat');
    case 4
        nPuntos = 9;
        data = load('Clp_1.dat');
end

kk_i = (ii2-1)*nPuntos+1;
kk_f = kk_i + (nPuntos-1);
data = data(kk_i:kk_f,:);

LAMe = atan(1/beta*tan(LAMc4));
Clp = -interp1(data(:,1),data(:,2),LAMe,'pchip');
end


