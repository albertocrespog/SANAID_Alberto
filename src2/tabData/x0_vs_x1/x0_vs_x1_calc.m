%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   17 de Agosto de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de x0/l_fus: Body Station where flow becomes viscous. AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.10 (pag 385, 2075 PDF))

function x0_lfus = x0_vs_x1_calc(x1, l_fus)
    x0_lfus = 0.378 + 0.527*x1/l_fus;

%     data = load('../digitalized_graph/x0_vs_x1.dat');
%     x0_lfus = interp1(data(:,1), data(:,2), x1/l_fus,'pchip');
end