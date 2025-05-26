close all
clear all
% Cargar los datos de las curvas de CP
load('curvas_cp.mat', 'curvas_cp');

% Ajuste de modelos polinómicos
n_curvas_cp = length(curvas_cp);
grados_pol = 4; % Grado del polinomio
modelos_cp = cell(1, n_curvas_cp);
grados_pol_red = 3; % Grado del polinomio
modelos_cp = cell(1, n_curvas_cp);

figure(1);
hold on;
title('Curvas CP y modelos polinómicos');
xlabel('J');
ylabel('CP');

for i = 1:n_curvas_cp
    x = curvas_cp{i}(:, 1);
    y = curvas_cp{i}(:, 2);
    p = polyfit(x, y, grados_pol);
    p_red = polyfit(x, y, grados_pol_red);
    modelos_cp{i} = p;
    modelos_cp_red{i} = p_red;
    
    % Muestra la curva ajustada
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);
    plot(x, y, 'o'); % Puntos digitalizados
    plot(x_fit, y_fit, 'LineWidth', 2); % Curva ajustada
end
grid on

legend('Data \beta=15', 'Poly \beta=25','Data \beta=20', 'Poly \beta=20','Data \beta=25', 'Poly \beta=25','Data \beta=30', 'Poly \beta=30','Data \beta=35', 'Poly \beta=35','Data \beta=40', 'Poly \beta=40','Data \beta=45', 'Poly \beta=45');
hold off;

% Guardar los modelos polinómicos
save('modelos_cp.mat', 'modelos_cp');
save('modelos_cp_red.mat', 'modelos_cp_red');


% Cargar los datos de las curvas de CP
load('curvas_ct.mat', 'curvas_ct');

% Ajuste de modelos polinómicos
n_curvas_ct = length(curvas_ct);
grados_pol = 4; % Grado del polinomio
grados_pol_red = 2; % Grado del polinomio
modelos_ct = cell(1, n_curvas_ct);


figure(2);
hold on;
title('Curvas CT y modelos polinómicos');
xlabel('J');
ylabel('CT');

for i = 1:n_curvas_ct
    x = curvas_ct{i}(:, 1);
    y = curvas_ct{i}(:, 2);
    p = polyfit(x, y, grados_pol);
    p_red = polyfit(x, y, grados_pol_red);
    modelos_ct{i} = p;
    modelos_ct_red{i} = p_red;
    
    % Muestra la curva ajustada
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);
    plot(x, y, 'o'); % Puntos digitalizados
    plot(x_fit, y_fit, 'LineWidth', 2); % Curva ajustada
end

legend('Data \beta=15', 'Poly \beta=25','Data \beta=20', 'Poly \beta=20','Data \beta=25', 'Poly \beta=25','Data \beta=30', 'Poly \beta=30','Data \beta=35', 'Poly \beta=35','Data \beta=40', 'Poly \beta=40','Data \beta=45', 'Poly \beta=45');
hold off;
grid on

% Guardar los modelos polinómicos
save('modelos_ct.mat', 'modelos_ct');
save('modelos_ct_red.mat', 'modelos_ct_red');