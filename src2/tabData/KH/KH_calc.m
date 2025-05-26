%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   27 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de KH: Factor that accounts for the relative size of the 
% horizontal and the vertical tail.

function KH = KH_calc(St, Sv)

data = load('KH.dat');
KH = interp1(data(:,1),data(:,2),St/Sv,'pchip');

end