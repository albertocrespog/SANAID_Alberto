%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   4 de Septiembre de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de xac_cre_BW
% Performance, Stability, Dynamics and Control; Pamadi, Bandu. Pag. 193,
% Figure 3.20

function xac_cre_BW0 = xac_cre_BW0_calc(ARexp_w, TRexp_w, LAM_w)
    par     = 0.25*ARexp_w*(1+TRexp_w)*tan(LAM_w);
    data = load('xac_cre_BW0.dat');
    xac_cre_BW0 = interp1(data(:,1), data(:,2), par, 'linear');
end