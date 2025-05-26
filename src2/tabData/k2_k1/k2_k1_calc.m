%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   17 de Agosto de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de k2-k1: Reduced mass factor.
% Performance, Stability, Dynamics and Control; Pamadi, Bandu. Pag. 170, Figure 3.6

function k2_k1 = k2_k1_calc(l_fus, D_max)
    data = load('k2_k1.dat');
    k2_k1 = interp1(data(:,1), data(:,2), l_fus/D_max, 'pchip');
end