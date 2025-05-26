%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   27 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de AvB_Av: Ratio of the vertical tail aspect ratio in the 
% presence of the fuselage to that of the isolated vertical tail.

function AvB_Av = AvB_Av_calc(b,d,lambda)
% Se tienen dos curvas, una para lambda=1 y otra para la lambda <= 0.6.
% Se comprueba la curva mas adecuada.
if abs(lambda - 1) < abs(lambda - 0.6)
    data = load('AvB_Av_1.dat');
else
    data = load('AvB_Av_06.dat');
end

eval = b/d;

if eval > max(data(:,1))
    eval = max(data(:,1));
elseif eval < min(data(:,1))
    eval = min(data(:,1));
end

AvB_Av = interp1(data(:,1),data(:,2),eval,'pchip');

end