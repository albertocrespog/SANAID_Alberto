% Digitalización de la gráfica

% Carga de la imagen de la gráfica
img = imread('Imagen1.jpg');
imshow(img);
hold on;

% Define las curvas de CP para cada ángulo de hélice
n_curvas = 7; % Número de curvas (una para cada ángulo de hélice)
curvas_cp = cell(1, n_curvas);

% Digitaliza los puntos de cada curva
for i = 1:n_curvas
    disp(['Seleccione los puntos de la curva ' num2str(i)]);
    [x, y] = ginput; % Digitalización manual de puntos
    curvas_cp{i} = [x, y];
    plot(x, y, 'o'); % Muestra los puntos seleccionados
end

% Ajuste de modelos polinómicos
grados_pol = 3; % Grado del polinomio
modelos_cp = cell(1, n_curvas);

for i = 1:n_curvas
    x = curvas_cp{i}(:, 1);
    y = curvas_cp{i}(:, 2);
    p = polyfit(x, y, grados_pol);
    modelos_cp{i} = p;
    
    % Muestra la curva ajustada
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'LineWidth', 2);
end

% Interpolación de CT (Ejemplo para una curva)
% Aquí deberías digitalizar las curvas de CT de manera similar y crear modelos de interpolación

% Asumiendo que ya tienes las curvas de CT digitalizadas en curvas_ct
curvas_ct = cell(1, n_curvas);
modelos_ct = cell(1, n_curvas);

for i = 1:n_curvas
    disp(['Seleccione los puntos de la curva CT ' num2str(i)]);
    [x, y] = ginput; % Digitalización manual de puntos
    curvas_ct{i} = [x, y];
    plot(x, y, '--'); % Muestra los puntos seleccionados
    
    % Crear interpolación
    x = curvas_ct{i}(:, 1);
    y = curvas_ct{i}(:, 2);
    interp_ct = griddedInterpolant(x, y, 'linear');
    modelos_ct{i} = interp_ct;
end

% Ejemplo de uso de los modelos
J = 0.6; % Valor de ejemplo para el parámetro de avance
angulo_helice = 25; % Ángulo de hélice de ejemplo

% Encuentra la curva correspondiente al ángulo de hélice
indice_curva = find([15, 20, 25, 30, 35, 40, 45] == angulo_helice);
if isempty(indice_curva)
    error('Ángulo de hélice no válido');
end

% Calcula el CP usando el modelo polinómico
p = modelos_cp{indice_curva};
CP = polyval(p, J);

% Calcula el CT usando la interpolación
interp_ct = modelos_ct{indice_curva};
CT = interp_ct(CP);

disp(['Para J = ', num2str(J), ', CP = ', num2str(CP), ' y CT = ', num2str(CT)]);