close all
clear all

% Digitalización de curvas de CT

% Carga de la imagen de la gráfica
img_ct = imread('grafica_ct.jpg');
imshow(img_ct);
hold on;

% Solicitar valores de los ejes
disp('Seleccione el punto de mínimo en los ejes X y Y');
[min_x, min_y] = ginput(1);
min_x_val = input('Ingrese el valor mínimo del eje X: ');
min_y_val = input('Ingrese el valor mínimo del eje Y: ');

disp('Seleccione el punto de máximo en los ejes X y Y');
[max_x, max_y] = ginput(1);
max_x_val = input('Ingrese el valor máximo del eje X: ');
max_y_val = input('Ingrese el valor máximo del eje Y: ');

% Digitalización de las curvas de CT
n_curvas_ct = 7; % Número de curvas de CT
curvas_ct = cell(1, n_curvas_ct);

for i = 1:n_curvas_ct
    disp(['Seleccione los puntos de la curva CT ' num2str(i)]);
    [x, y] = ginput; % Digitalización manual de puntos
    % Escalar los puntos
    x_scaled = min_x_val + (x - min_x) * (max_x_val - min_x_val) / (max_x - min_x);
    y_scaled = min_y_val + (y - min_y) * (max_y_val - min_y_val) / (max_y - min_y);
    curvas_ct{i} = [x_scaled, y_scaled];
    plot(x_scaled, y_scaled, 'o'); % Muestra los puntos seleccionados
end

save('curvas_ct.mat', 'curvas_ct'); % Guardar los datos de las curvas de CT