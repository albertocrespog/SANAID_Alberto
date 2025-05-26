%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   13 de Octubre de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de Clr_CL_M0: Pamadi, Bandu. Fig. 4.28, pag. 417

function Clr_CL_M0 = Clr_CL_M0_calc(lambda, AR, LAMc4)
nPoint1 = 10; % Numero de puntos por curva
nPoint2 = 2;
par_pos = [0 0.25 0.5 1];
[~, k] = min(abs(lambda - par_pos));
k_i = (k-1)*nPoint1+1;
k_f = k*nPoint1;
par_pos2 = [0 15 30 45 60]*pi/180;
[~,kk] = min(abs(abs(LAMc4) - par_pos2));
kk_i = (kk-1)*nPoint2+1;
kk_f = kk*nPoint2;
if AR > 10
    AR = 10;
end

data1 = load('Clr_CL_M0_1.dat');
data1 = data1(k_i:k_f,:);
val1 = interp1(data1(:,1),data1(:,2),AR,'pchip');

data2 = load('Clr_CL_M0_2.dat');
data2 = data2(kk_i:kk_f,:);
Clr_CL_M0 = interp1(data2(:,1),data2(:,2),val1,'linear');

end

