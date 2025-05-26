function downw = downwash_calc(X_w, X_h, AR_w, lambda_w, Lambda_c4, b_w, h_h)
% Pamadi 195, 3.43
% X_w es la posicion del ala
% X_h es la posicion del estabilizador vertical
% AR_w es el aspect ratio del ala
% lambda_w es el estrechamiento de ala
% Lambda_c4 es la flecha del ala en c/4
% b_w es la envergadura del ala
% h_h es la distancia entre el ala y el estabilizador horizontal, positiva
% cuando el estabilizador se encuentra por encima del ala.
    
    l_h         = (X_h - X_w);
    % Prevents downwash from going to infinite values
    if l_h == 0
        l_h = 0.01;
    end
    Ka          = 1/AR_w - 1/(1+AR_w^1.7);
    Kl          = (10-3*lambda_w)/7;
    Kh          = (1-h_h/b_w)/(2*l_h/b_w)^(1/3);
    downw       = 4.44  * (Ka*Kl*Kh*(cos(Lambda_c4))^1/2)^1.19;
end