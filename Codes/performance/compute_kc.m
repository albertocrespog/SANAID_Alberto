function kc = compute_kc(W_S)
    % Función para modelar k_c en función del wing loading W/S
    
    % Puntos conocidos
    W_S1 = 20;      % Wing loading inicial (psf)
    W_S2 = 100;     % Wing loading final (psf)
    kc1 = 33;       % Valor de k_c en W/S = 20 psf
    kc2 = 28.6;     % Valor de k_c en W/S = 100 psf
    
    % Calcular la pendiente de la recta
    m = (kc2 - kc1) / (W_S2 - W_S1);
    
    % Calcular k_c usando la ecuación de la recta
    kc = m * (W_S - W_S1) + kc1;
    
    % Asegurarse de que k_c no supere el límite de 33 para W/S <= 20
    kc(W_S <= W_S1) = kc1;
    
    % Asegurarse de que k_c se limite a valores inferiores en W/S >= 100
    kc(W_S >= W_S2) = kc2;
end