%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   28 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de K_Mtau: Compressibility correction to dihedral. AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 396, 2086 PDF)

function K_Mdie = K_Mdie_calc(M, AR, LAMc2)
% nPoint = 13; % Numero de puntos por curva
% par = AR/cos(LAMc2);
% par_pos = [10 8 6 4 2];
% [~, kk] = min(abs(par - par_pos));
% kk_i = (kk-1)*nPoint+1;
% kk_f = kk_i + nPoint;
% data = load('K_Mtau.dat');
% data = data(kk_i:kk_f,:);
% K_Mdie = interp1(data(:,1),data(:,2),M*cos(LAMc2),'pchip');
% end

%%  -----------------------------------------------------------------------
%   MODIFIED EMERGENTIA
%   06 de Febrero de 2020
%   Creado por Sergio Esteban
%%   ----------------------------------------------------------------------
% Calculo de K_Mdie: The compressibility correction to the lifting surface dihedral effect can be found from 
% Figure 10.25 in Airplane Design Part VI and is a function of the steady state flight Mach number, the lifting 
% surface aspect ratio, and the lifting surface half chord sweep angle:

% DATCOM FIGURE 5.1.2.1-30-A
Y = [2., 4., 6., 8., 10.];
X = [0.0,.2,.4,.5,.6,.7,.8,.9,.95];
Z =       [1.0, 1.01, 1.018, 1.02, 1.023, 1.03, 1.04, 1.05, 1.057,;
       1.0, 1.012, 1.03, 1.045, 1.06, 1.085, 1.118, 1.16, 1.19,;  
     1.0, 1.015, 1.045, 1.07, 1.1, 1.14, 1.197, 1.27, 1.33,; 
     1.0, 1.018, 1.05, 1.085, 1.125, 1.19, 1.26, 1.39, 1.485,;
     1.0, 1.02, 1.058, 1.097, 1.148, 1.215, 1.325, 1.495, 1.635];

Yq = AR/cos(LAMc2);
if Yq > 10;
    Yq = 10;
elseif Yq < 2
    Yq = 2;
end
Xq = M*cos(LAMc2);
K_Mdie = interp2(X,Y,Z,Xq,Yq);
