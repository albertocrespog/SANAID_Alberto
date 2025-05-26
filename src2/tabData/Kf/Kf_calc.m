%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   28 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de Kf: Fuselage correction factor. AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 394, 2084 PDF))

function Kf = Kf_calc(lt, b, AR, LAMc2)
% if AR/cos(LAMc2)> 8
%     par_pos = 8;
% elseif AR/cos(LAMc2)< 4
%     par_pos = 4;
% end
%     
% nPoint = 16; % Numero de puntos por curva
% par = AR/cos(LAMc2);
% par_pos = [8 7 6 5.5 5 4.5 4];
% [~, kk] = min(abs(par - par_pos));
% kk_i = (kk-1)*nPoint+1;
% kk_f = kk_i + nPoint;
% data = load('Kf.dat');
% data = data(kk_i:kk_f,:);
% Kf = interp1(data(:,1),data(:,2),lt/b,'pchip');
% end

%%  -----------------------------------------------------------------------
%   MODIFIED EMERGENTIA
%   06 de Febrero de 2020
%   Creado por Sergio Esteban
%%   ----------------------------------------------------------------------
% Calculo de Kf: The fuselage correction factor is obtained from Figure 10.22 in Airplane Design Part VI 
% and is a function of the lifting surface aspect ratio, lifting surface mid-chord sweep angle, the length of the fuselage and the lifting surface span:
% DATCOM FIGURE 5.2.2.1-26

Y = [4.0,4.5,5.0,5.5,6.0,7.0,8.0];
X = [0.0,.2,.4,.6,.8,1.0,1.2,1.4,1.6];
Z = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,.99,.97;
     1.0,1.0,1.0,1.0,1.0,1.0,.98,.948,.911;
     1.0,1.0,1.0,1.0,.997,.971,.933,.883,.827;
     1.0,1.0,1.0,.991,.963,.922,.870,.811,.746;
     1.0,1.0,.995,.970,.932,.884,.829,.764,.695;
     1.0,1.0,.977,.944,.899,.845,.780,.715,.641;
     1.0,.985,.960,.921,.870,.812,.745,.670,.592];

Yq = AR/cos(LAMc2);
if Yq > 8;
    Yq = 8;
elseif Yq < 4
    Yq = 4;
end
Xq = lt/b;

if Xq > max(X)
    % Extrapolate using the two last points in X
    slope = (Z(:,end) - Z(:,end-1)) / (X(end) - X(end-1));
    Z_extrapolated = Z(:,end) + slope * (Xq - X(end));
    Kf = interp1(Y, Z_extrapolated, Yq, 'linear', 'extrap');
else
    % Interpolate within the bounds
    Kf = interp2(X, Y, Z, Xq, Yq);
end

