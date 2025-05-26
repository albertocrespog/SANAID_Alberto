%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   18 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Proporciona la posicion del punto neutro de la aeronave en unidades de 
% cuerda media del ala.

function [X_na] = pre_aeroCenter(S_ref, S_w, S_h, S_c, X_w, X_h, X_c,...
                                 c_m, CLa_w, CLa_h, CLa_c, eta_h, eta_c,...
                                 lambda_w, AR_w, b_w, h_h, Lambda_c4)
    
    downw   = downwash_calc(X_w, X_h, AR_w, lambda_w, Lambda_c4, b_w, h_h);   
    upw     = upwash_calc(AR_w, X_w, X_c, c_m);
    
    % Aportes de sustentacion de cada una de las superficies    
    CL(1,1) = eta_h*S_h/S_ref*CLa_h*(1-downw);  % Estabilizador horizontal
    CL(1,2) = eta_c*S_c/S_ref*CLa_c*(1+upw);    % Canard
    CL(1,3) = S_w/S_ref*CLa_w;                  % Ala
    X       = [X_c, X_h, X_w];
    X_na = X*CL/sum(CL)/c_m;

end

