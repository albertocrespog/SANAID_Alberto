%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   26 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Proporciona el coeficiente Kn empleado en el calculo de la contribución
% del fuselaje a la derivada de estabilidad C_n_beta. Este factor introduce
% el efecto de la interferencia ala-fuselaje y es funcion de la geometria y
% de la posicion del centro de gravedad.
% 
% Las graficas digitalizadas se pueden encontrar en la seccion 3, pagina 
% 272 del libro Perfomance, Stability, Dynamics and Control; Autor: Pamadi, Bandu N.

%%

function K_N = K_N_calc(l_fus, S_bSide, Xcg, h, h1, h2, bf_max)
%% ARGUMENTOS DE ENTRADA
% l_fus: longitud del fuselaje
% S_bSide: area lateral del fuselaje
% Xcg: posicion del centro de graveda medida desde el morro
% h: altura maxima del fuselaje
% h1: altura del fuselaje a 1/4 de la longitud del fuselaje
% h2: altura del fuselaje a 3/4 de la longitud del fuselaje
% bf_max: anchura maxima del fuselaje

%% CARGA DE LAS TRES GRAFICAS
% Todas las gráficas están formadas por rectas y = a + bx. 
% coefX_pos =   [a1     b1; 
%                a2     b2;
%                .       .
%                an     bn];

kn1 = load('kn1.dat');
[M,~] = size(kn1);

% Rectas paralelas => b igual para todas
coef1_pos(1,:) = ([1 kn1(1,1); 1 kn1(2,1)]\[kn1(1,2); kn1(2,2)])';
coef1_pos(2:(M-1),1) = kn1(3:M,2) - coef1_pos(1,2)*kn1(3:M,1);
coef1_pos(:,2) = coef1_pos(1,2);


kn2 = load('kn2.dat');
[M,~] = size(kn2);

% Rectas que pasan por el origen => a igual para todas
coef2_pos(1,:) = ([1 kn2(1,1); 1 kn2(2,1)]\[kn2(1,2); kn2(2,2)])';
coef2_pos(2:(M-1),2) = (kn2(3:M,2) - coef2_pos(1,1))./kn2(3:M,1);
coef2_pos(:,1) = coef2_pos(1,1);


kn3 = load('kn3.dat');
[M,~] = size(kn3);
kn3_offset = 4.70339e-04;

% Rectas que pasan por el origen => a igual para todas
coef3_pos(1,:) = ([1 kn3(1,1); 1 kn3(2,1)]\[kn3(1,2); kn3(2,2)])';
coef3_pos(2:(M-1),2) = (kn3(3:M,2) - coef3_pos(1,1))./kn3(3:M,1);
coef3_pos(:,1) = coef3_pos(1,1);

%% CALCULO DEL COEFICIENTE
% Primera grafica
par1 = l_fus^2/S_bSide;
par1_pos = [20 14 10 8 7 6 5 4 3 2.5];
[~,k_par1] = min(abs(par1-par1_pos));
coef1 = coef1_pos(k_par1,:);
x1 = Xcg/l_fus;
y1 = coef1(1) + coef1(2)*x1;

% Segunda grafica
par2 = sqrt(h1/h2);
par2_pos = 1.6:-0.2:0.8;
[~,k_par2] = min(abs(par2-par2_pos));
coef2 = coef2_pos(k_par2,:);
y2 = coef2(1) + coef2(2)*y1;

% Tercera grafica
par3 = h/bf_max;
par3_pos = [0.5 0.6 0.8 1 2];
[~,k_par3] = min(abs(par3-par3_pos));
coef3 = coef3_pos(k_par3,:);
K_N = coef3(1) + coef3(2)*y2-kn3_offset;

end