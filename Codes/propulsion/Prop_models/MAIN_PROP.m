clear all
close all

% Ejemplo de uso de la función
J = 0.6; % Valor de ejemplo para el parámetro de avance
angulo_helice = 25; % Ángulo de hélice de ejemplo
CP = obtener_cp(J, angulo_helice);
disp(['Para J = ', num2str(J), ', y ángulo de hélice de ', num2str(angulo_helice), ' grados, CP = ', num2str(CP)]);
CT = obtener_ct(J, angulo_helice);
disp(['Para J = ', num2str(J), ', y ángulo de hélice de ', num2str(angulo_helice), ' grados, CT = ', num2str(CT)]);
etha_p = J*CT/CP

% 
% % Cargar los datos de las curvas de CP
% load('curvas_cp.mat', 'curvas_cp');
% load('modelos_cp.mat', 'modelos_cp');
% 
% % Cargar la imagen de la gráfica
% img_cp = imread('grafica_cp.jpg');
% imshow(img_cp);
% hold on;
% title('Curvas CP ajustadas sobre la imagen original');
% xlabel('J');
% ylabel('CP');
% 
% % Solicitar valores de los ejes
% disp('Seleccione el punto de mínimo en los ejes X y Y');
% [min_x, min_y] = ginput(1);
% min_x_val = input('Ingrese el valor mínimo del eje X: ');
% min_y_val = input('Ingrese el valor mínimo del eje Y: ');
% 
% disp('Seleccione el punto de máximo en los ejes X y Y');
% [max_x, max_y] = ginput(1);
% max_x_val = input('Ingrese el valor máximo del eje X: ');
% max_y_val = input('Ingrese el valor máximo del eje Y: ');
% 
% % Interpolar y pintar las curvas ajustadas
% J_values = linspace(0, 2.9, 100); % Rango de J
% 
% for i = 1:length(curvas_cp)
%     % Obtener el modelo polinómico correspondiente
%     p = modelos_cp{i};
% 
%     % Calcular los valores ajustados de CP
%     CP_values = polyval(p, J_values);
% 
%     % Escalar los puntos ajustados de vuelta a las coordenadas de la imagen
%     J_scaled = min_x + (J_values - min_x_val) * (max_x - min_x) / (max_x_val - min_x_val);
%     CP_scaled = min_y + (CP_values - min_y_val) * (max_y - min_y) / (max_y_val - min_y_val);
% 
%     % Pintar las curvas ajustadas
%     plot(J_scaled, CP_scaled, 'LineWidth', 2);
% end
% 
% legend('Curvas ajustadas');
% hold off;