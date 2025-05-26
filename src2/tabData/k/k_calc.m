%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   27 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de k: Factor empleado en la estimacion de las derivadas de
% estabilidad direccional frente a viento cruzado. AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.12 (pag 385, 2075 PDF))

function k = k_calc(bv, d1)
data = load('k.dat');
if bv/d1 > 6
    k = data(end,2);
elseif bv/d1 < 0
    k = data(1,2);
else
    k = interp1(data(:,1),data(:,2),bv/d1,'linear');
end
