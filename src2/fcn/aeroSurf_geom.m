function [c_r, c_t, x_ac, y_ac, z_ac, LAMc4, LAMc2, LAM, MISC] = aeroSurf_geom(S, b, y, c_y, le_y,diedro)
l_ad   = y(2:end) - y(1:(end-1));          % Envergadura de los paneles
off_ad = le_y(2:end) - le_y(1:(end-1));    % Desplazamiento en x de cada panel
%b   = 2*y(end);

%% Cálculo de las areas de cada panel y cuerdas en raiz y punta
Si_ad       = 0.5*l_ad.*(c_y(1:(end-1)) + c_y(2:(end)));

c_r         = S/b/sum(Si_ad);
c_t         = c_r*c_y(end);

c_y     = c_y*c_r;
le_y    = le_y*c_r;
y       = y*b/2;
l       = l_ad*b/2;
off     = off_ad*c_r;

Si      = b/2*c_r*Si_ad;


%% Cálculo del centro de cada panel. Fuente: http://www.efunda.com/math/areas/trapezoid.cfm
% X_c     = (2*c_y(2:end).*off + c_y(2:end).^2 + off.*c_y(1:(end-1))... 
%           + c_y(1:(end-1)).*c_y(2:end) + c_y(1:(end-1)).^2)./(c_y(2:end)...
%           + c_y(1:(end-1)))/3;
cm_i    = Si./l;
Y_c     = l.*(2*c_y(2:end) + c_y(1:(end-1)))./(c_y(2:end)+c_y(1:(end-1)))/3;

%% Cálculo de las coordenadas del centro aerodinámico del ala. Fuente: http://naca.central.cranfield.ac.uk/reports/1942/naca-report-751.pdf
%X_c     = le_y(1:(end-1)) + X_c;
X_c     = le_y(1:(end-1)) + off/2 + cm_i*0.25;
Y_c     = y(1:(end-1)) + Y_c;

x_ac    = 2/S*sum(Si.*X_c);
y_ac    = 2/S*sum(Si.*Y_c);
z_ac    = 2/S*sum(Si.*Y_c*tan(diedro));

%% Cálculo de la flecha
LAM     = atan(le_y(end)/y(end));
LAMc4   = atan((le_y(end) + 0.25*(c_y(end) - c_y(1)))/y(end));
LAMc2   = atan((le_y(end) + 0.5*(c_y(end) - c_y(1)))/y(end));

LAM_deg = LAM*180/pi;
LAMc4_deg = LAMc4*180/pi;
LAMc2_deg = LAMc2*180/pi;

MISC.LAM_deg = LAM_deg;
MISC.LAMc4_deg = LAMc4_deg;
MISC.LAMc2_deg = LAMc2_deg;

end