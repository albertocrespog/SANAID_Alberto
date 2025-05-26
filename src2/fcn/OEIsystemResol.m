function [phi, da, dr] = OEIsystemResol(solType, model, derMatrix, rho, beta, V, P)
% Diapositivas de la asignatura Cálculo de Aeronaves. Tema 14_1, diap. 8

% DATOS NECESARIOS PARA LA RESOLUCION DEL SISTEMA

%W           = model.general.mtow*0.98*9.8065;
W           = model.general.mtow*9.8065;
Sref        = model.general.Sref;

b           = model.ala.b;

switch solType
    case 'OEI'
        nV      = length(V);
        nP      = length(P);
        nBeta   = length(beta);
        [n, id] = max([nV,nP]);

    % INICIALIZACION DE LAS MATRICES
        phi = zeros(n,nBeta);
        da  = phi;
        dr  = phi;

        rendProp    = model.propulsion.rendProp;
        Foei        = model.propulsion.F_OEI;
        d           = model.propulsion.Y;

        if id == 1
            par = [V', ones(n,1)*P];
        else
            par = [ones(n,1)*V, P'];
        end
        
    
        for jPar = 1:n
            for jBeta = 1:nBeta
        
                F           = rendProp*par(jPar,2)/par(jPar,1); % rendProp*P/V;
                N           = Foei*F*d;
                q           = 0.5*rho*par(jPar,1)^2; % 0.5*rho*V^2

            % DEFINICION DE LAS MATRICES DEL SISTEMA
                indep       = -derMatrix(:,1)*beta(jBeta) - [0; 0; N/Sref/b/q];
                coefMat     = [W/q/Sref; 0; 0];
                coefMat     = [coefMat, derMatrix(:,2), derMatrix(:,3)];

            % SOLUCION
                sol         = coefMat\indep;
                sol(1)      = asin(sol(1));
        
                phi(jPar,jBeta) = sol(1);
                da(jPar,jBeta)  = sol(2);
                dr(jPar,jBeta)  = sol(3);
        
            end
        end
    case 'Viento'
        n       = length(V);
        nBeta   = length(beta);

    % INICIALIZACION DE LAS MATRICES
        phi = zeros(n,nBeta);
        da  = phi;
        dr  = phi;

        for jV = 1:n
            for jBeta = 1:nBeta
                
                q           = 0.5*rho*V(jV,1)^2; % 0.5*rho*V^2
                
            % DEFINICION DE LAS MATRICES DEL SISTEMA
                indep       = -derMatrix(:,1)*beta(jBeta);
                coefMat     = [W/q/Sref; 0; 0];
                coefMat     = [coefMat, derMatrix(:,2:3)];

            % SOLUCION
                sol         = coefMat\indep;
                sol(1)      = asin(sol(1));
        
                phi(jV,jBeta) = sol(1);
                da(jV,jBeta)  = sol(2);
                dr(jV,jBeta)  = sol(3);
        
            end
        end
        
end

end
