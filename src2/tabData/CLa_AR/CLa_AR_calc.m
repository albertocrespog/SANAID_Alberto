%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   27 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de CLa_AR: Pendiente de la curva de sustentacion por unidad de
% relación de aspecto. El valor de la pendiente esta dado en (1/rad).

function CLa_AR = CLa_AR_calc(AR,beta,LAMc2,k)
param = AR/k*(beta^2 + (tan(LAMc2))^2);
data = load('CLa_AR.dat');

if param > max(data(:,1))
    param = max(data(:,1));
elseif param < min(data(:,1))
    param = min(data(:,1));
end

CLa_AR = interp1(data(:,1),data(:,2),param,'pchip');

end