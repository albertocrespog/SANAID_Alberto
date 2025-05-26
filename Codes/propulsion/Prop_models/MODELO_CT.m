close all
clear all
% Cargar los datos de las curvas de CP
load('curvas_cp.mat', 'curvas_cp');

% Ajuste de modelos polinómicos
n_curvas_cp = length(curvas_cp);
grados_pol = 4; % Grado del polinomio
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
    modelos_cp{i} = p;
    
    % Muestra la curva ajustada
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);
    plot(x, y, 'o'); % Puntos digitalizados
    plot(x_fit, y_fit, 'LineWidth', 2); % Curva ajustada
end
grid on

legend('Datos', 'Modelo polinómico');
hold off;

% Guardar los modelos polinómicos
save('modelos_cp.mat', 'modelos_cp');


% Cargar los datos de las curvas de CP
load('curvas_ct.mat', 'curvas_ct');

% Ajuste de modelos polinómicos
n_curvas_ct = length(curvas_ct);
grados_pol = 4; % Grado del polinomio
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
    modelos_ct{i} = p;
    
    % Muestra la curva ajustada
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);
    plot(x, y, 'o'); % Puntos digitalizados
    plot(x_fit, y_fit, 'LineWidth', 2); % Curva ajustada
end

legend('Datos', 'Modelo polinómico','1','1');
hold off;
grid on

% Guardar los modelos polinómicos
save('modelos_ct.mat', 'modelos_ct');