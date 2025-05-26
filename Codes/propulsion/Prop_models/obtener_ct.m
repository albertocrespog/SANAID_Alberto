function CP = obtener_ct(J, angulo_helice)
    % Cargar los modelos polinómicos
    load('modelos_ct.mat', 'modelos_ct');
    
    % Ángulos de hélice disponibles
    angulos_helice = [15, 20, 25, 30, 35, 40, 45];
    
    % Encontrar el índice correspondiente al ángulo de hélice
    indice_curva = find(angulos_helice == angulo_helice);
    if isempty(indice_curva)
        error('Ángulo de hélice no válido. Debe ser uno de los siguientes: [15, 20, 25, 30, 35, 40, 45]');
    end
    
    % Obtener el modelo polinómico correspondiente
    p = modelos_ct{indice_curva};
    
    % Calcular CP usando el modelo polinómico
    CP = polyval(p, J);
end

