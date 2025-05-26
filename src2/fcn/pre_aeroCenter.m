%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   18 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Proporciona la posicion del punto neutro de la aeronave en unidades de 
% cuerda media del ala.

function [X_na] = pre_aeroCenter(model)

%% ADQUISICIÓN DE DATOS 
S_ref       = model.general.Sref;

% Datos del ala
S_w         = model.ala.S;
X_w         = model.ala.Xca;
c_m         = model.ala.MAC;
eta_w       = model.ala.eta;
CLa_w       = model.ala.CLa;
lambda_w    = model.ala.TR;
AR_w        = model.ala.AR;
b_w         = model.ala.b;
Lambda_c4   = model.ala.LAMc4;


switch model.conf
    case 'convencional'
        S_t         = model.horizontal.S;
        if isempty(S_t)
            S_t = 0;
        end
        
        CLa_t       = model.horizontal.CLa;
        if isempty(CLa_t)
            CLa_t = 0;
        end
        
        eta_t       = model.horizontal.eta;
        if isempty(eta_t)
            eta_t = 0;
        end
        
        X_t         = model.horizontal.Xca;
        if isempty(X_t)
            X_t = 0;
        end
        
        if isempty(model.horizontal.Zca) || isempty(model.ala.Zca)
            h_t         = 0;
        else
            h_t = model.horizontal.Zca + model.ala.Zca;
        end
        
        S_c     = 0;
        CLa_c   = 0;
        eta_c   = 0;
        X_c     = 0;
        
    case 'canard'
        S_c         = model.canard.S;
        if isempty(S_c)
            S_c = 0;
        end
        
        CLa_c       = model.canard.CLa;
        if isempty(S_c)
            S_c = 0;
        end
        
        eta_c       = model.canard.eta;
        if isempty(eta_c)
            eta_c = 0;
        end
        
       
        X_c         = model.canard.Xca;
        if isempty(X_c)
            X_c = 0;
        end
        
        S_t     = 0;
        CLa_t   = 0;
        eta_t   = 0;
        X_t     = 0;
        h_t     = 0;
        
    case 'convencional_canard'
        S_t         = model.horizontal.S;
        if isempty(S_t)
            S_t = 0;
        end
        
        CLa_t       = model.horizontal.CLa;
        if isempty(CLa_t)
            CLa_t = 0;
        end
        
        eta_t       = model.horizontal.eta;
        if isempty(eta_t)
            eta_t = 0;
        end
        
        X_t         = model.horizontal.Xca;
        if isempty(X_t)
            X_t = 0;
        end
        
        if isempty(model.horizontal.Zca) || isempty(model.ala.Zca)
            h_t         = 0;
        else
            h_t = model.horizontal.Zca + model.ala.Zca;
        end
        
        S_c         = model.canard.S;
        if isempty(S_c)
            S_c = 0;
        end
        
        CLa_c       = model.canard.CLa;
        if isempty(CLa_c)
            CLa_c = 0;
        end
        
        eta_c       = model.canard.eta;
        if isempty(eta_c)
            eta_c = 0;
        end
        
        X_c         = model.canard.Xca;
        if isempty(X_c)
            X_c = 0;
        end
        
    case 'flyWing'
        S_t     = 0;
        CLa_t   = 0;
        eta_t   = 0;
        X_t     = 0;
        h_t     = 0;
        
        S_c     = 0;
        CLa_c   = 0;
        eta_c   = 0;
        X_c     = 0;
        
end



%% CALCULO DEL CENTRO AERODINAMICO DE LA AERONAVE
    % Obtención del downwash
    downw   = downwash_calc(X_w, X_t, AR_w, lambda_w, Lambda_c4, b_w, h_t);   
    % Obtención del upwash
    upw     = upwash_calc(AR_w, X_w, X_c, c_m);
    
    % Aportes de sustentacion de cada una de las superficies    
    CL(1,1) = eta_t*S_t/S_ref*CLa_t*(1-downw);   % Estabilizador horizontal
    CL(2,1) = eta_c*S_c/S_ref*CLa_c*(1+upw);     % Canard
    CL(3,1) = eta_w*S_w/S_ref*CLa_w;             % Ala
    % Posiciones de las superficies
    X       = [X_t; X_c; X_w];   
    % Centro aerodinámico
    X_na = sum(X.*CL)/sum(CL);

end

