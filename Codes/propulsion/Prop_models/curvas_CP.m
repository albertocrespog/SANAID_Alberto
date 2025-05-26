close all
clear all

% Digitalización de curvas de CP

% Carga de la imagen de la gráfica
img_cp = imread('grafica_cp.jpg');
imshow(img_cp);
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

% Digitalización de las curvas de CP
n_curvas_cp = 7; % Número de curvas de CP
curvas_cp = cell(1, n_curvas_cp);

for i = 1:n_curvas_cp
    disp(['Seleccione los puntos de la curva CP ' num2str(i)]);
    [x, y] = ginput; % Digitalización manual de puntos
    % Escalar los puntos
    x_scaled = min_x_val + (x - min_x) * (max_x_val - min_x_val) / (max_x - min_x);
    y_scaled = min_y_val + (y - min_y) * (max_y_val - min_y_val) / (max_y - min_y);
    curvas_cp{i} = [x_scaled, y_scaled];
    plot(x_scaled, y_scaled, 'o'); % Muestra los puntos seleccionados
end

save('curvas_cp.mat', 'curvas_cp'); % Guardar los datos de las curvas de CP