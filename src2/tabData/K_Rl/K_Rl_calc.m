%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   28 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de K_Rl: factor dependent on Reynold's Number. AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.29 (pag 400, 2090 PDF))

function K_Rl = K_Rl_calc(Re)
data = load('K_Rl.dat');
K_Rl =interp1(data(:,1),data(:,2),Re,'linear');
end