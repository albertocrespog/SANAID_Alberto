%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   27 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de AvHB_AvB: Ratio of the vertical tail aspect ratio in the 
% presence of the fuselage and the horizontal tail to that in the presence 
% of the fuselage alone.
%
function AvHB_AvB = AvHB_AvB_calc(xac_h, bv, zac_h, cv)
par = xac_h/cv;
par_pos = [0.5 0.6 0.7 0.8];
[~, k] = min(abs(par - par_pos));
switch k
    case 1
        data = load('AvHB_AvB_05.dat');
    case 2
        data = load('AvHB_AvB_06.dat');
    case 3
        data = load('AvHB_AvB_07.dat');
    case 4
        data = load('AvHB_AvB_08.dat');
end

eval = zac_h/bv;
if eval > max(data(:,1))
    eval = max(data(:,1));
elseif eval < min(data(:,1))
    eval = min(data(:,1));
end

AvHB_AvB = interp1(data(:,1),data(:,2),eval,'pchip');
end