function [fuel_crucero,tiempo_crucero,distancia_crucero,datos] = analisis_crucero_v1(Aero_TH,Aero,Geo_tier,segment,Weight_tier,conv_UNITS)

% function [fuel_crucero,tiempo_crucero,distancia_crucero,datos] = analisis_crucero(propul,aerodinamica,crucero,W_inicial,h_inicial,opcion,i,datos)
g = conv_UNITS.g;

m_TOW = Weight_tier.m_TOW;
W_inicial = Weight_tier.m_TOW*g;
h_inicial = segment{1}.data.h_initial;

% PROPUL
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR
% 4: CONSUMO ESPECIFICO
% 5: DERIVACION(TURBOFANES)
% 6: CAPACIDAD CALORIFICA DEL COMBUSTIBLE
% 7: DIAMETRO DEL MOTOR(TURBOHELICES)
% 8: EFICIENCIA DE LA HELICE (ETA_P)
% 9: AVION CIVIL = 1/MILITAR = 2

% AERODINAMICA
% 1: SUPERFICIE
% 2: CD0
% 3: K1 = K
% 4: K2
% 5: CLMAX LIMPIO

% CRUCERO
% 1: DISTANCIA FINAL
% 2: MACH DE VUELO
% 3: CL DE CRUCERO
% 4: PALANCA DE GASES
% 5: VELOCIDAD INICIAL
% 6: VELOCIDAD FINAL
% 7: COMBUSTIBLE A QUEMAR
% 8: CDO = F(M)
% 9: K1 = F(M)
% 10: K2 = F(M)

% Aero_TH
Geo_tier.S_ref

% AERODINAMICA
% S = aerodinamica(1);
% Cd0 = aerodinamica(2);
% k1 = aerodinamica(3);
% k2 = aerodinamica(4);
% CLmax_limpio = aerodinamica(5);
S = Geo_tier.S_ref;
Cd0 = Aero_TH.CD0;
k2 = -Aero_TH.CD1;
k1 = Aero_TH.CD2;
CLmax_limpio = Aero.C_L_max_w1_CR;

% PROPULSION
% tipo_motor = propul(1);
% num_motores = propul(2);
% T_SL = propul(3);
% c_SL = propul(4);
% civil_mil = propul(5);
% eta_inst = propul(6);
% derivacion = propul(7);
tipo_motor = 2;
num_motores = 1;
T_SL = m_TOW*g;
c_SL = 1;
% civil_mil = propul(5);
eta_inst = 0.7;
derivacion = 2;

a = segment{1}.data.Data_ATM.a;
M = segment{1}.data.V_cr/a;

% CRUCERO
% dist_final = crucero(1);
% Mach = crucero(2);
% C_l = crucero(3);
% delta_T_gases = crucero(4);
% velocidad_ini = crucero(5);
% velocidad_fin = crucero(6);
% fuel_a_quemar = crucero(7);
% Cd0_polar = crucero(8);
% k1_polar = crucero(9);
% k2_polar = crucero(10);
dist_final = segment{1}.data.Range*10;
Mach = M;
% C_l = crucero(3);
% delta_T_gases = crucero(4);
% velocidad_ini = segment{1}.data.V_cr;
% velocidad_fin = segment{1}.data.V_cr;
% fuel_a_quemar = crucero(7);
Cd0_polar = Cd0;
k1_polar = k1;
k2_polar = k2;

% ATMOSFERA ESTANDAR
[Temp_SL,rho_SL,p_SL,a_SL] = atmos_inter_mio(0);
[Temp,rho,p,a] = atmos_inter_mio(h_inicial);

% ------------------------  CONSTANTES --------------------------
g = 9.80665;

W_ini = W_inicial;
h_ini = h_inicial;
x_ini = 0;

iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();   %ESTA LINEA PERMITE USAR CONDICIONES EN LAS FUNCIONES!!!

eta_p = eta_inst;

T_SL = T_SL * num_motores;

%-------------------------------------------------------------------------%
%ECUACIONES DE LA MECANICA DEL VUELO PARA EL CRUCERO

% T(delta_T,V) = D(V,CL) + W/g * dV/dt
% L(V,CL)      = W
% dW/dt        = -C(delta_T,V)*g*T(delta_T,V)

%-------------------------------------------------------------------------%


switch opcion
    
    case 1  %CRUCERO A MACH CONSTANTE Y DISTANCIA DADA
        
        % C_A = 0.5*rho*V^2*S*Cd0;
        % C_B = k2;
        % C_C = k1/(0.5*rho*V^2*S);
        % C_D = sqrt(4*C_A*C_C - C_B^2);
        % C_E = C*g/V;
        % C_F = 2*atan((2*C_C*W_ini-C_B)/C_D)/C_D;
        %
        % % Ecuacion dW/dt = -C*T integrada (BREGUET añadiendole el termino k2)
        %
        % W = @(x) (C_D.*tan(C_D.*(C_F-C_E.*x)./2) + C_B)./(2*C_C);
        
        x_final = dist_final;
        M = Mach;
        V = M * a;
        
        CL = @(W) W./(0.5*rho*V^2*S);
        CD = @(W) Cd0 - k2.*CL(W) + k1.*CL(W).^2;
        D  = @(W) 0.5*rho*V^2*S.*CD(W);
        T  = @(W) D(W);
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        
        switch tipo_motor
            case 1
                delta_T = @(W) T(W)./(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
                correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
                switch derivacion
                    case 1
                        C =  @(W) c_SL.*correccion(W).*(1+1.2*M).*sqrt(Temp/Temp_SL);
                    case 2
                        C =  @(W) c_SL.*correccion(W).*(1+0.33*M).*sqrt(Temp/Temp_SL);
                    case 3
                        C =  @(W) c_SL.*correccion(W).*(1+0.16875*M).*sqrt(Temp/Temp_SL);
                end
            case 2
                delta_T = @(W) T(W)./(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
                correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
                C = @(W) c_SL.*correccion(W).*(1+1.44*M).*sqrt(Temp/Temp_SL).* V./eta_p;
            case 3
                delta_T = @(W) T(W)./(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
                correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
                C = @(W) c_SL.*correccion(W).*V./eta_p;
        end
        
        % Integramos  dW/dt = -C(W)gT(W)
        
        num = 10000;
        X_v = linspace(0,x_final,num);
        
        dwdx = @(x,W) -C(W).*g.*T(W)/V;
        [X,PESO] = ode45(dwdx,X_v,W_ini);
        
        for k = 1:(num-1)
            peso(k) = PESO(k);
            L_D(k) = CL(PESO(k))/CD(PESO(k));
            empuje(k) = T(PESO(k));
            CL2(k) = CL(PESO(k));
            palanca(k) = delta_T(PESO(k));
            tiempo(k) = (X(k+1)-X(k))/V;
        end
        
        fuel_crucero = (peso(1)-peso(end))/g;
        tiempo_crucero = sum(tiempo);
        distancia_crucero = dist_final;
        W = peso(end);
        
        datos(i).crucero.tiempo = cumsum(tiempo);
        datos(i).crucero.fuel   = (PESO(1) - PESO)./g;
        datos(i).crucero.fuel(end) = '';
        datos(i).crucero.distancia = X_v(1:end-1);
        datos(i).crucero.palanca = palanca;
        %datos(i).crucero.empuje = empuje;
        datos(i).crucero.empuje = CL2;
        datos(i).crucero.L_D = L_D;
        datos(i).crucero.peso = peso;
        datos(i).crucero.velocidad = V;%escalar
        
        try
            delta_T = 1;
            W = W_ini;
            
            M = @(V) V./a;
            
            switch tipo_motor
                case 1
                    T = @(V) delta_T.*(T_SL.*(1+0.2.*M(V).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
                case 2
                    T = @(V) delta_T.*(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
                case 3
                    T = @(V) delta_T.*(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
            end
            
            CL = @(V) W./(0.5.*rho.*V^2.*S);
            CD = @(V) Cd0 -k2.*CL(V) + k1.*CL(V).^2;
            D  = @(V) 0.5.*rho.*V.^2.*S.*CD(V);
            
            funcion = @(V) (T(V) - D(V));
            
            syms x positive
            Vmax = double(solve(funcion(x)));
            for kkk = 1:length(Vmax),
                if Vmax(kkk) > 340,
                    Vmax(kkk) = 0;
                end
            end
            Vmax = max(Vmax);
            clear x;
            
            datos(i).crucero.Vmax = Vmax;
        catch
            datos(i).crucero.Vmax = 0;
        end
        %-------------------------------------------------------------------------%
        
    case 2 % CRUCERO CON CL DADO CONSTANTE Y DISTANCIA DADA
        
        x_final = dist_final;
        CL = C_l;
        CD = Cd0 - k2*CL + k1*CL^2;
        
        V = @(W) sqrt(W./(0.5.*rho.*S.*CL));
        dvdw = @(W) 1./sqrt(2*rho.*S.*CL.*W);
        
        M = @(W) V(W)./a;
        D = @(W) 0.5.*rho.*V(W).^2.*S.*CD;
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        C = @(W) c_SL.*(1+1.2.*M(W)).*sqrt(Temp/Temp_SL);
                    case 2
                        C = @(W) c_SL.*(1+0.33.*M(W)).*sqrt(Temp/Temp_SL);
                    case 3
                        C = @(W) c_SL.*(1+0.16875.*M(W)).*sqrt(Temp/Temp_SL);
                end
            case 2
                
                C = @(W) c_SL.*(1+1.44.*M(W)).*sqrt(Temp/Temp_SL).* V(W)./eta_p;
            case 3
                C = @(W) c_SL .*V(W)./eta_p;
        end
        
        dwdx = @(x,W) (-C(W).*g.*D(W))./(V(W).*(1+C(W).*W.*dvdw(W)));
        
        T = @(x,W) D(W) + (W./g).*dvdw(W).*dwdx(x,W).*V(W);
        
        switch tipo_motor
            case 1
                delta_T = @(x,W) T(x,W)./(T_SL.*(1+0.2.*M(W).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(W))).*(rho./rho_SL));
            case 2
                delta_T  = @(x,W) T(x,W)./(T_SL.*eta_p./V(W).*(1+0.2.*M(W).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL.*Temp_SL)));
            case 3
                delta_T = @(x,W) T(x,W)./(T_SL.*eta_p./V(W).*((8.55.*rho./rho_SL - 1)./7.55));
        end
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        correccion = (a_1.*delta_T(0,W_ini).^4+a_2.*delta_T(0,W_ini).^3+a_3.*delta_T(0,W_ini).^2+a_4.*delta_T(0,W_ini)+a_5);
        
        % Simplificamos la compleja operacion de introducir la dependencia de
        % C(delta_T) comprobando que la variacion de la palanca es muy pequeña y
        % por tanto la correccion es bastante constante a lo largo del crucero
        
        V = @(W) sqrt(W./(0.5.*rho.*S.*CL));
        dvdw = @(W) 1./sqrt(2*rho.*S.*CL.*W);
        
        M = @(W) V(W)./a;
        D = @(W) 0.5.*rho.*V(W).^2.*S.*CD;
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        C = @(W) c_SL.*correccion.*(1+1.2.*M(W)).*sqrt(Temp/Temp_SL);
                    case 2
                        C = @(W) c_SL.*correccion.*(1+0.33.*M(W)).*sqrt(Temp/Temp_SL);
                    case 3
                        C = @(W) c_SL.*correccion.*(1+0.16875.*M(W)).*sqrt(Temp/Temp_SL);
                end
            case 2
                
                C = @(W) c_SL.*correccion.*(1+1.44.*M(W)).*sqrt(Temp/Temp_SL).* V(W)./eta_p;
            case 3
                C = @(W) c_SL .*correccion.*V(W)./eta_p;
        end
        
        dwdx = @(x,W) (-C(W).*g.*D(W))./(V(W).*(1+C(W).*W.*dvdw(W)));
        
        T = @(x,W) D(W) + (W./g).*dvdw(W).*dwdx(x,W).*V(W);
        
        switch tipo_motor
            case 1
                delta_T = @(x,W) T(x,W)./(T_SL.*(1+0.2.*M(W).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(W))).*(rho./rho_SL));
            case 2
                delta_T  = @(x,W) T(x,W)./(T_SL.*eta_p./V(W).*(1+0.2.*M(W).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL.*Temp_SL)));
            case 3
                delta_T = @(x,W) T(x,W)./(T_SL.*eta_p./V(W).*((8.55.*rho./rho_SL - 1)./7.55));
        end
        
        n= 10000;
        vector_x = linspace(0,x_final,n);
        
        [X,PESO] = ode45(dwdx,vector_x,W_ini);
        
        for k = 1:n,
            peso(k) = PESO(k);
            L_D(k) = CL/CD;
            empuje(k) = T(X(k),peso(k));
            CL2(k) = CL;
            palanca(k) = delta_T(X(k),peso(k));
            velocidad(k) = V(PESO(k));
            if k < n,
                tiempo(k) = (X(k+1)-X(k))/V(peso(k));
            else
                tiempo(k) = tiempo(k-1);
            end
        end
        
        fuel_crucero = (peso(1)-peso(end))/g;
        tiempo_crucero = sum(tiempo);
        distancia_crucero = dist_final;
        W = peso(end);
        
        datos(i).crucero.tiempo = cumsum(tiempo);
        datos(i).crucero.fuel   = (PESO(1) - PESO)./g;
        datos(i).crucero.distancia = X;
        datos(i).crucero.palanca = palanca;
        %datos(i).crucero.empuje = empuje;
        datos(i).crucero.empuje = CL2;
        datos(i).crucero.L_D = L_D;
        datos(i).crucero.peso = peso;
        datos(i).crucero.velocidad = velocidad;
        
        %-------------------------------------------------------------------------%
        
    case 3   %CRUCERO ACELERADO CON VELOCIDAD INICIAL Y FINAL DADOS Y ALTURA DADA Y PALANCA DE GASES DADA
        
        try
            delta_T = delta_T_gases;
            W = W_ini;
            
            M = @(V) V./a;
            
            a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
            correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
            
            switch tipo_motor
                case 1
                    T = @(V) delta_T.*(T_SL.*(1+0.2.*M(V).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
                case 2
                    T = @(V) delta_T.*(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
                case 3
                    T = @(V) delta_T.*(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
            end
            
            CL = @(V) W./(0.5.*rho.*V^2.*S);
            CD = @(V) Cd0 -k2.*CL(V) + k1.*CL(V).^2;
            D  = @(V) 0.5.*rho.*V.^2.*S.*CD(V);
            
            funcion = @(V) (T(V) - D(V));
            
            syms x positive
            Vmax = double(solve(funcion(x)));
            for kkk = 1:length(Vmax),
                if Vmax(kkk) > 340,
                    Vmax(kkk) = 0;
                end
            end
            Vmax = max(Vmax);
            clear x;
            
            if Vmax < velocidad_fin && Vmax > 0,
                V_fin = Vmax;
            else
                V_fin = velocidad_fin;
            end
            
        catch
            V_fin = velocidad_fin;
            W = W_ini;
        end
        
        delta_T = delta_T_gases;
        V_ini = velocidad_ini;
        %V_fin = velocidad_fin;
        
        M = @(V) V./a;
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        
        switch tipo_motor
            case 1
                T = @(V) delta_T.*(T_SL*(1+0.2.*M(V).^2).^(1.4./0.4)*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
                switch derivacion
                    case 1
                        C = @(V) c_SL.*correccion.*(1+1.2.*M(V)).*sqrt(Temp/Temp_SL);
                    case 2
                        C = @(V) c_SL.*correccion.*(1+0.33.*M(V)).*sqrt(Temp/Temp_SL);
                    case 3
                        C = @(V) c_SL.*correccion.*(1+0.16875.*M(V)).*sqrt(Temp/Temp_SL);
                end
            case 2
                T = @(V) delta_T.*(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
                C = @(V) c_SL.*correccion.*(1+1.44.*M(V)).*sqrt(Temp/Temp_SL).* V./eta_p;
            case 3
                T = @(V) delta_T.*(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
                C = @(V) c_SL .*correccion.*V./eta_p;
        end
        
        CL = @(V,W) W./(0.5.*rho.*V.^2.*S);
        CD = @(V,W) Cd0 -k2.*CL(V,W) + k1.*CL(V,W).^2;
        D  = @(V,W) 0.5.*rho.*V.^2.*S.*CD(V,W);
        
        dwdv = @(V,W) -(W.*C(V).*T(V))./(T(V)-D(V,W));  %Despejando de las ecuaciones del crucero
        
        n= 10000;
        vector_v = linspace(V_ini,V_fin,n);
        
        [VEL,PESO] = ode45(dwdv,vector_v,W_ini);
        
        for k = 1:n,
            peso(k) = PESO(k);
            L_D(k) = CL(VEL(k),PESO(k))/CD(VEL(k),PESO(k));
            CL2(k) = CL(VEL(k),PESO(k));
            empuje(k) = T(VEL(k));
            palanca(k) = delta_T;
            %dt = @(V) PESO(k)./(g.*(T(V)-D(V,PESO(k))));
            %dvdt = @(V) g.*(T(V)-D(V,PESO(k)))./PESO(k);
            %dx = @(V) V./dvdt(V);
            if k < n
                %tiempo(k) = integral(dt,VEL(k),VEL(k+1));
                %distancia(k) = integral(dx,VEL(k),VEL(k+1));
                
                %TARDA 18 SEGUNDOS USANDO INTEGRALES NUMERICAS Y  SEGUNDOS USANDO ESTAS
                %EXPRESIONES Y EL RESULTADO TIENE MENOS DE UN 0.1% DE ERROR
                
                tiempo(k) = PESO(k).*(VEL(k+1)-VEL(k))./(g.*(empuje(k)-D(VEL(k),PESO(k))));
                distancia(k) = ((VEL(k+1)+VEL(k))/2)*tiempo(k);
            elseif k == n
                %tiempo(k) = tiempo(k-1);
                %distancia(k) = distancia(k-1);
                distancia(k) = distancia(k-1);
                tiempo(k) = tiempo(k-1);
            end
        end
        
        fuel_crucero = (peso(1)-peso(end))/g;
        %tiempo_crucero = sum(tiempo);
        tiempo_crucero = sum(tiempo);
        %distancia_crucero = sum(distancia);
        distancia_crucero = sum(distancia);
        W = peso(end);
        
        datos(i).crucero.tiempo = cumsum(tiempo);
        datos(i).crucero.fuel   = (PESO(1) - PESO)./g;
        datos(i).crucero.distancia = cumsum(distancia);
        datos(i).crucero.palanca = palanca;
        %datos(i).crucero.empuje = empuje;
        datos(i).crucero.empuje = CL2;
        datos(i).crucero.L_D = L_D;
        datos(i).crucero.peso = peso;
        datos(i).crucero.velocidad = VEL;
        
        
        %-------------------------------------------------------------------------%
        
        %     case 4 %VELOCIDAD MAXIMA (PALANCA = 1)
        %
        % delta_T = 1;
        % W = W_ini;
        %
        % M = @(V) V./a;
        %
        % switch tipo_motor
        %     case 1
        %         T = @(V) delta_T.*(T_SL.*(1+0.2.*M(V).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
        %     case 2
        %         T = @(V) delta_T.*(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
        %     case 3
        %         T = @(V) delta_T.*(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
        % end
        %
        % CL = @(V) W./(0.5.*rho.*V^2.*S);
        % CD = @(V) Cd0 -k2.*CL(V) + k1.*CL(V).^2;
        % D  = @(V) 0.5.*rho.*V.^2.*S.*CD(V);
        %
        % funcion = @(V) (T(V) - D(V));
        % V_stall = sqrt(2*W/(rho*S*CLmax_limpio));
        %
        % syms x positive
        % Vmax = double(solve(funcion(x)));
        % Vmax = max(Vmax);
        % clear x;
        %
        % datos_crucero = [];
        
        %-------------------------------------------------------------------------%
        
    case 4 % CRUCERO A MACH CONSTANTE CON COEFICIENTES AERODINAMICOS = f(M)
        
        x_final = dist_final;
        M = Mach;
        V = M * a;
        
        % vector_interp = Mach_polar(1):0.001:Mach_polar(end);
        %
        % Cd0_M = interp1(Mach_polar,Cd0_polar,vector_interp,'spline');
        % k2_M = interp1(Mach_polar,k2_polar,vector_interp,'spline');
        % k1_M = interp1(Mach_polar,k1_polar,vector_interp,'spline');
        %
        % Cd0 = Cd0_M(vector_interp == M);
        % k2 = k2_M(vector_interp == M);
        % k1 = k1_M(vector_interp == M);
        
        Cd0 = Cd0_polar;
        k2  = k2_polar;
        k1  = k1_polar;
        
        CL = @(W) W./(0.5*rho*V^2*S);
        CD = @(W) Cd0 - k2.*CL(W) + k1.*CL(W).^2;
        D  = @(W) 0.5*rho*V^2*S.*CD(W);
        T  = @(W) D(W);
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        
        switch tipo_motor
            case 1
                delta_T = @(W) T(W)./(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
                correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
                switch derivacion
                    case 1
                        C =  @(W) c_SL.*correccion(W).*(1+1.2*M).*sqrt(Temp/Temp_SL);
                    case 2
                        C =  @(W) c_SL.*correccion(W).*(1+0.33*M).*sqrt(Temp/Temp_SL);
                    case 3
                        C =  @(W) c_SL.*correccion(W).*(1+0.16875*M).*sqrt(Temp/Temp_SL);
                end
            case 2
                delta_T = @(W) T(W)./(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
                correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
                C = @(W) c_SL.*correccion(W).*(1+1.44*M).*sqrt(Temp/Temp_SL).* V./eta_p;
            case 3
                delta_T = @(W) T(W)./(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
                correccion = @(W)(a_1.*delta_T(W).^4+a_2.*delta_T(W).^3+a_3.*delta_T(W).^2+a_4.*delta_T(W)+a_5);
                C = @(W) c_SL.*correccion(W).*V./eta_p;
        end
        
        % Integramos  dW/dt = -C(W)gT(W)
        
        num = 10000;
        X_v = linspace(0,x_final,num);
        
        dwdx = @(x,W) -C(W).*g.*T(W)/V;
        [X,PESO] = ode45(dwdx,X_v,W_ini);
        
        for k = 1:(num-1)
            peso(k) = PESO(k);
            L_D(k) = CL(PESO(k))/CD(PESO(k));
            empuje(k) = T(PESO(k));
            CL2(k) = CL(PESO(k));
            palanca(k) = delta_T(PESO(k));
            tiempo(k) = (X(k+1)-X(k))/V;
        end
        
        fuel_crucero = (peso(1)-peso(end))/g;
        tiempo_crucero = sum(tiempo);
        distancia_crucero = dist_final;
        W = peso(end);
        
        datos(i).crucero.tiempo = cumsum(tiempo);
        datos(i).crucero.fuel   = (PESO(1) - PESO)./g;
        datos(i).crucero.fuel(end) = '';
        datos(i).crucero.distancia = X_v(1:end-1);
        datos(i).crucero.palanca = palanca;
        %datos(i).crucero.empuje = empuje;
        datos(i).crucero.empuje = CL2;
        datos(i).crucero.L_D = L_D;
        datos(i).crucero.peso = peso;
        datos(i).crucero.velocidad = V;
        
        
        %-------------------------------------------------------------------------%
        
    case 5 %CRUCERO DE MAXIMO ALCANCE DADO WFINAL
        
        W_fin = W_ini - fuel_a_quemar*g;
        
        if W_fin < 0,
            W_fin = 0.01;
        end
        
        % Para hallar el crucero de maximo alcance se va a emplear calculo
        % variacional en la ecuacion integral del alcance, sea esta:
        
        % x = integral entre W0 y W de -V/(C(V)*g*D(V,W)) dW
        
        % De forma que buscamos la ley de velocidades V = V(W) que para cada W
        % proporcione el maximo alcance a la aeronave.
        
        % Derivando el integrando respecto de V e igualando a 0 obtenemos una
        % ecuacion para la obtencion de V(W)
        
        % C*D - V*(dC/dV*D + C*dD/dV) = 0
        
        % Esto implica que segun sea el tipo de consumo (el cual cambia segun el
        % tipo de motor), de una forma u otra sera la expresion V(W) que
        % proporciona el maximo alcance.
        
        % switch tipo_motor
        %     case 1
        %         switch derivacion
        %             case 1
        %                 C = @(V) c_SL.*(1+1.2.*M(V)).*sqrt(Temp/Temp_SL);
        %             case 2
        %                 C = @(V) c_SL.*(1+0.33.*M(V)).*sqrt(Temp/Temp_SL);
        %             case 3
        %                 C = @(V) c_SL.*(1+0.16875.*M(V)).*sqrt(Temp/Temp_SL);
        %         end
        %     case 2
        %         C = @(V) c_SL.*(1+1.44.*M(V)).*sqrt(Temp/Temp_SL).* V./eta_p;
        %     case 3
        %         C = @(V) c_SL .*V./eta_p;
        % end
        %
        % C = A + B*V (JET) ; C = V*(A+B*V) (PROP Y PISTON)
        %
        % D = (0.5 rho S Cd0) V^2 - k2  * W + k1/(0.5 rho S) * (W^2/V^2)
        %   =      A_2        V^2 - B_2 * W +      C_2       * (W^2/V^2)
        %
        
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        A = c_SL * sqrt(Temp/Temp_SL);
                        B = 1.2*A/a;
                    case 2
                        A = c_SL * sqrt(Temp/Temp_SL);
                        B = 0.33*A/a;
                    case 3
                        A = c_SL * sqrt(Temp/Temp_SL);
                        B = 0.16875*A/a;
                end
            case 2
                A = c_SL * sqrt(Temp/Temp_SL)/eta_p;
                B = 1.44*A/a;
            case 3
                A = c_SL/eta_p;
                B = 0;
        end
        
        A_2 = 0.5*rho*S*Cd0;
        B_2 = k2;
        C_2 = k1/(0.5*rho*S);
        
        switch tipo_motor
            case 1
                %         W = @(V) sqrt((2.*A_2.*B.*V.^5 + A.*A_2.*V.^4 + A.*B_2.*V.^2)./(2.*C_2.*B.*V + 3.*A.*C_2));
                %
                %         AA = 2*A_2*B; BB = A*A_2; CC = A*B_2; DD = 2*C_2*B; EE = 3*A*C_2;
                %
                %         dwdv = @(V) (V.*(4.*AA.*DD.*V.^4 + 5.*AA.*EE.*V.^3 + 3.*BB.*DD.*V.^3 + 4.*BB.*EE.*V.^2 + CC.*DD.*V + 2.*CC.*EE))...
                %             ./(2.*(DD.*V+EE).^2.*V.*sqrt((AA.*V.^3+BB.*V.^2+CC)./(DD.*V+EE)));
                
                AA = @(V) 2*B*C_2.*V + 3*A*C_2; BB = @(V) -A*B_2.*V.^2; CC = @(V) -(A*A_2.*V.^4 + 2*A_2*B.*V.^5);
                dAAdv = @(V) 2*B*C_2; dBBdv = @(V) -2*A*B_2.*V; dCCdv = @(V) -(4*A*A_2.*V.^3 + 10*A_2*B.*V.^4);
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
            case 2
                %         W = @(V) sqrt((3.*A_2.*B.*V.^5 + 2.*A.*A_2.*V.^4 - B.*B_2.*V.^3)./(C_2.*B.*V + 2.*A.*C_2));
                %
                %         AA = 3*A_2*B; BB = 2*A*A_2; CC = B*B_2; DD = C_2*B; EE = 2*A*C_2;
                %
                %         dwdv = @(V) (V.^2.*(4.*AA.*DD.*V.^3 + 5.*AA.*EE.*V.^2 + 3.*BB.*DD.*V.^2 + 4.*BB.*EE.*V - 2.*CC.*DD.*V - 3.*CC.*EE))...
                %             ./(2.*(DD.*V+EE).^2.*V.*sqrt(V.*(AA.*V.^2+BB.*V-CC)./(DD.*V+EE)));
                
                AA = @(V) B*C_2.*V + 2*A*C_2; BB = @(V) B*B_2.*V.^3; CC = @(V) -(2*A*A_2.*V.^4 + 3*A_2*B.*V.^5);
                dAAdv = @(V) B*C_2; dBBdv = @(V) 3*B*B_2.*V.^2; dCCdv = @(V) -(8*A*A_2.*V.^3 + 15*A_2*B.*V.^4);
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
            case 3
                W = @(V) sqrt(A_2/C_2).*V.^2;
                
                dwdv = @(V) 2.* sqrt(A_2/C_2).*V;
        end
        
        M = @(V) V./a;
        CL = @(V) W(V)./(0.5.*rho.*V.^2.*S);
        CD = @(V) Cd0 -k2.*CL(V) + k1.*CL(V).^2;
        D = @(V) 0.5.*rho.*V.^2.*S.*CD(V);
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        C = @(V) c_SL.*(1+1.2.*M(V)).*sqrt(Temp/Temp_SL);
                    case 2
                        C = @(V) c_SL.*(1+0.33.*M(V)).*sqrt(Temp/Temp_SL);
                    case 3
                        C = @(V) c_SL.*(1+0.16875.*M(V)).*sqrt(Temp/Temp_SL);
                end
            case 2
                C = @(V) c_SL.*(1+1.44.*M(V)).*sqrt(Temp/Temp_SL).* V./eta_p;
            case 3
                C = @(V) c_SL .*V./eta_p;
        end
        
        dvdt = @(V) -(C(V).*g.*D(V))./(dwdv(V) + C(V).*W(V));
        dt   = @(V) -(dwdv(V) + C(V).*W(V))./(C(V).*g.*D(V));
        dx   = @(V) V.*dt(V);
        
        T    = @(V) D(V) + (W(V)./g).*dvdt(V);
        
        switch tipo_motor
            case 1
                delta_T = @(V) T(V)./(T_SL.*(1+0.2.*M(V).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
            case 2
                delta_T = @(V) T(V)./(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
            case 3
                delta_T = @(V) T(V)./(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
        end
        
        % Con esto ya tenemos todas las funciones dependiendo de V, solo resta
        % hallar las condiciones inicial y final de la velocidad usando la
        % expresion de V=V(W)
        
        syms vel positive
        
        V_ini = double(solve((W(vel)-W_ini),'Real',true));
        
        V_fin = double(solve((W(vel)-W_fin),'Real',true));
        
        distancia_int = @(V) integral(dx,V_ini,V);
        tiempo_int    = @(V) integral(dt,V_ini,V);
        
        n = 1000;
        VEL = linspace(V_ini,V_fin,n);
        
        %------------------------------------------------------------------------%
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        palanca1 = delta_T(VEL(1));
        palanca2 = delta_T(VEL(end));
        palanca_media = (palanca1 + palanca2)/2;
        correccion = (a_1.*palanca_media.^4+a_2.*palanca_media.^3+a_3.*palanca_media.^2+a_4.*palanca_media+a_5);
        %------------------------------------------------------------------------%
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        A = c_SL *correccion * sqrt(Temp/Temp_SL);
                        B = 1.2*A/a;
                    case 2
                        A = c_SL *correccion * sqrt(Temp/Temp_SL);
                        B = 0.33*A/a;
                    case 3
                        A = c_SL *correccion * sqrt(Temp/Temp_SL);
                        B = 0.16875*A/a;
                end
            case 2
                A = c_SL *correccion * sqrt(Temp/Temp_SL)/eta_p;
                B = 1.44*A/a;
            case 3
                A = c_SL*correccion/eta_p;
                B = 0;
        end
        
        A_2 = 0.5*rho*S*Cd0;
        B_2 = k2;
        C_2 = k1/(0.5*rho*S);
        
        switch tipo_motor
            case 1
                %         W = @(V) sqrt((2.*A_2.*B.*V.^5 + A.*A_2.*V.^4 + A.*B_2.*V.^2)./(2.*C_2.*B.*V + 3.*A.*C_2));
                %
                %         AA = 2*A_2*B; BB = A*A_2; CC = A*B_2; DD = 2*C_2*B; EE = 3*A*C_2;
                %
                %         dwdv = @(V) (V.*(4.*AA.*DD.*V.^4 + 5.*AA.*EE.*V.^3 + 3.*BB.*DD.*V.^3 + 4.*BB.*EE.*V.^2 + CC.*DD.*V + 2.*CC.*EE))...
                %             ./(2.*(DD.*V+EE).^2.*V.*sqrt((AA.*V.^3+BB.*V.^2+CC)./(DD.*V+EE)));
                
                AA = @(V) 2*B*C_2.*V + 3*A*C_2; BB = @(V) -A*B_2.*V.^2; CC = @(V) -(A*A_2.*V.^4 + 2*A_2*B.*V.^5);
                dAAdv = @(V) 2*B*C_2; dBBdv = @(V) -2*A*B_2.*V; dCCdv = @(V) -(4*A*A_2.*V.^3 + 10*A_2*B.*V.^4);
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
            case 2
                %         W = @(V) sqrt((3.*A_2.*B.*V.^5 + 2.*A.*A_2.*V.^4 - B.*B_2.*V.^3)./(C_2.*B.*V + 2.*A.*C_2));
                %
                %         AA = 3*A_2*B; BB = 2*A*A_2; CC = B*B_2; DD = C_2*B; EE = 2*A*C_2;
                %
                %         dwdv = @(V) (V.^2.*(4.*AA.*DD.*V.^3 + 5.*AA.*EE.*V.^2 + 3.*BB.*DD.*V.^2 + 4.*BB.*EE.*V - 2.*CC.*DD.*V - 3.*CC.*EE))...
                %             ./(2.*(DD.*V+EE).^2.*V.*sqrt(V.*(AA.*V.^2+BB.*V-CC)./(DD.*V+EE)));
                
                AA = @(V) B*C_2.*V + 2*A*C_2; BB = @(V) B*B_2.*V.^3; CC = @(V) -(2*A*A_2.*V.^4 + 3*A_2*B.*V.^5);
                dAAdv = @(V) B*C_2; dBBdv = @(V) 3*B*B_2.*V.^2; dCCdv = @(V) -(8*A*A_2.*V.^3 + 15*A_2*B.*V.^4);
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
            case 3
                W = @(V) sqrt(A_2/C_2).*V.^2;
                
                dwdv = @(V) 2.* sqrt(A_2/C_2).*V;
        end
        
        M = @(V) V./a;
        CL = @(V) W(V)./(0.5.*rho.*V.^2.*S);
        CD = @(V) Cd0 -k2.*CL(V) + k1.*CL(V).^2;
        D = @(V) 0.5.*rho.*V.^2.*S.*CD(V);
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        C = @(V) c_SL.*correccion .*(1+1.2.*M(V)).*sqrt(Temp/Temp_SL);
                    case 2
                        C = @(V) c_SL.*correccion .*(1+0.33.*M(V)).*sqrt(Temp/Temp_SL);
                    case 3
                        C = @(V) c_SL.*correccion .*(1+0.16875.*M(V)).*sqrt(Temp/Temp_SL);
                end
            case 2
                C = @(V) c_SL.*correccion .*(1+1.44.*M(V)).*sqrt(Temp/Temp_SL).* V./eta_p;
            case 3
                C = @(V) c_SL .*correccion .*V./eta_p;
        end
        
        dvdt = @(V) -(C(V).*g.*D(V))./(dwdv(V) + C(V).*W(V));
        dt   = @(V) -(dwdv(V) + C(V).*W(V))./(C(V).*g.*D(V));
        dx   = @(V) V.*dt(V);
        
        T    = @(V) D(V) + (W(V)./g).*dvdt(V);
        
        switch tipo_motor
            case 1
                delta_T = @(V) T(V)./(T_SL.*(1+0.2.*M(V).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
            case 2
                delta_T = @(V) T(V)./(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
            case 3
                delta_T = @(V) T(V)./(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
        end
        
        % Con esto ya tenemos todas las funciones dependiendo de V, solo resta
        % hallar las condiciones inicial y final de la velocidad usando la
        % expresion de V=V(W)
        
        syms vel positive
        
        V_ini = double(solve((W(vel)-W_ini),'Real',true));
        
        V_fin = double(solve((W(vel)-W_fin),'Real',true));
        
        distancia_int = @(V) integral(dx,V_ini,V);
        tiempo_int    = @(V) integral(dt,V_ini,V);
        
        n = 1000;
        VEL = linspace(V_ini,V_fin,n);
        
        %-------------------------------------------------------------------------%
        
        for k = 1:n,
            peso(k) = W(VEL(k));
            L_D(k) = CL(VEL(k))/CD(VEL(k));
            empuje(k) = T(VEL(k));
            CL2(k) = CL(VEL(k));
            palanca(k) = delta_T(VEL(k));
            tiempo(k) = tiempo_int(VEL(k));
            distancia(k) = distancia_int(VEL(k));
        end
        
        alcance_max = distancia(end);
        tiempo_final= tiempo(end);
        distancia_crucero = distancia(end);
        fuel_crucero = (peso(1)-peso(end))./g;
        tiempo_crucero = tiempo(end);
        
        datos(i).crucero.tiempo = tiempo;
        datos(i).crucero.fuel   = (peso(1) - peso)./g;
        datos(i).crucero.distancia = distancia;
        datos(i).crucero.palanca = palanca;
        %datos(i).crucero.empuje = empuje;
        datos(i).crucero.empuje = CL2;
        datos(i).crucero.L_D = L_D;
        datos(i).crucero.peso = peso;
        datos(i).crucero.velocidad = VEL;
        
        %-------------------------------------------------------------------------%
        
    case 6 %CRUCERO EN CONFIGURACION DE MAX AUTONOMIA DADO WFINAL
        
        W_fin = W_ini - fuel_a_quemar*g;
        
        if W_fin < 0,
            W_fin = 0.01;
        end
        
        % Para hallar el crucero de maxima autonomia se va a emplear calculo
        % variacional en la ecuacion integral de la autonomia, sea esta:
        
        % x = integral entre W0 y W de -1/(C(V)*g*D(V,W)) dW
        
        % De forma que buscamos la ley de velocidades V = V(W) que para cada W
        % proporcione la maxima autonomia a la aeronave.
        
        % Derivando el integrando respecto de V e igualando a 0 obtenemos una
        % ecuacion para la obtencion de V(W)
        
        % (dC/dV*D + C*dD/dV) = 0
        
        % Esto implica que segun sea el tipo de consumo (el cual cambia segun el
        % tipo de motor), de una forma u otra sera la expresion V(W) que
        % proporciona la maxima autonomia.
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        A = c_SL * sqrt(Temp/Temp_SL);
                        B = 1.2*A/a;
                    case 2
                        A = c_SL * sqrt(Temp/Temp_SL);
                        B = 0.33*A/a;
                    case 3
                        A = c_SL * sqrt(Temp/Temp_SL);
                        B = 0.16875*A/a;
                end
            case 2
                A = c_SL * sqrt(Temp/Temp_SL)/eta_p;
                B = 1.44*A/a;
            case 3
                A = c_SL/eta_p;
                B = 0;
        end
        
        A_2 = 0.5*rho*S*Cd0;
        B_2 = k2;
        C_2 = k1/(0.5*rho*S);
        
        switch tipo_motor
            case 1
                %W = @(V) sqrt((3.*A_2.*B.*V.^5 + 2.*A.*A_2.*V.^4 - B.*B_2.*V.^3)./(C_2.*B.*V + 2.*A.*C_2));
                %W = @(V) ( -B.*B_2.*V.^3 + sqrt((B.*B_2).^2.*V.^6 + 4.*(B.*C_2.*V + 2.*A.*C_2).*(2.*A.*A_2.*V.^4 + 3.*A_2.*B.*V.^5)))./(2.*(B.*C_2.*V + 2.*A.*C_2));
                
                AA = @(V) B*C_2.*V + 2*A*C_2; BB = @(V) B*B_2.*V.^3; CC = @(V) -(2*A*A_2.*V.^4 + 3*A_2*B.*V.^5);
                dAAdv = @(V) B*C_2; dBBdv = @(V) 3*B*B_2.*V.^2; dCCdv = @(V) -(8*A*A_2.*V.^3 + 15*A_2*B.*V.^4);
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
                %AA = 3*A_2*B; BB = 2*A*A_2; CC = B*B_2; DD = C_2*B; EE = 2*A*C_2;
                %AA = -B*B_2; BB = B^2*B_2^2 + 12*A_2*B^2*C_2; CC = 32*A*A_2*B*C_2; DD = 16*A^2*A_2*C_2; EE = 2*B*C_2; FF = 4*A*C_2;
                
                %dwdv = @(V) (V.^2.*(4.*AA.*DD.*V.^3 + 5.*AA.*EE.*V.^2 + 3.*BB.*DD.*V.^2 + 4.*BB.*EE.*V - 2.*CC.*DD.*V - 3.*CC.*EE))...
                %./(2.*(DD.*V+EE).^2.*V.*sqrt(V.*(AA.*V.^2+BB.*V-CC)./(DD.*V+EE)));
                %dwdv = @(V) AA.*V.^2.*(2*EE.*V+3*FF)./((EE.*V+FF).^2) - EE.*sqrt(V.^4*(BB.*V.^2+CC.*V+DD))...
                %./((EE.*V+FF).^2) + V.^3.*(6.*BB.*V.^2 + 5.*CC.*V + 4*DD)./(2.*(EE.*V+FF).*sqrt(V.^4.*(BB.*V.^2 + CC.*V + DD)));
                
            case 2
                %W = @(V) sqrt((4.*A_2.*B.*V.^5 + 3.*A.*A_2.*V.^4 - 2.*B.*B_2.*V.^3 - A.*B_2.*V.^2)./(A.*C_2));
                %W2 = @(V) ((V.^2)./(2.*A.*C_2)).*(-(A.*B_2 + 2.*B.*B_2.*V)+sqrt((A.*B_2 + 2.*B.*B_2.*V).^2+4.*A.*C_2.*(3.*A.*A_2 + 4.*A_2.*B.*V)));
                
                AA = @(V) A*C_2; BB = @(V) A*B_2.*V.^2 + 2*B*B_2.*V.^3; CC = @(V) -(3*A*A_2.*V.^4 + 4*A_2*B.*V.^5);
                dAAdv = @(V) 0; dBBdv = @(V) 2*A*B_2.*V + 6*B*B_2.*V.^2; dCCdv = @(V) -(12*A*A_2.*V.^3 + 20*A_2*B.*V.^4);
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
                %AA = 4*A_2*B; BB = 3*A*A_2; CC = 2*B*B_2; DD = B_2*A; EE = A*C_2;
                
                %dwdv = @(V) (V.*(5.*AA(V).*V.^3 + 4.*BB.*V.^2 - 3.*CC.*V - 2.*DD))...
                %./(2.*EE.*sqrt((4.*A_2.*B.*V.^5 + 3.*A.*A_2.*V.^4 - 2.*B.*B_2.*V.^3 - A.*B_2.*V.^2)./(A.*C_2)));
                
            case 3
                %         W = @(V) sqrt((3.*A_2.*V.^4-B_2.*V.^2)./C_2);
                %
                %         dwdv = @(V) (12.*A_2.*V.^3-2.*B_2.*V)./(2.*C_2.*sqrt((3.*A_2.*V.^4-B_2.*V.^2)./C_2));
                
                AA = @(V) C_2; BB = @(V) B_2.*V.^2; CC = @(V) -3*A_2.*V.^4;
                dAAdv = @(V) 0; dBBdv = @(V) 2*B_2.*V; dCCdv = @(V) -12*A_2.*V.^3;
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
        end
        
        M = @(V) V./a;
        CL = @(V) W(V)./(0.5.*rho.*V.^2.*S);
        CD = @(V) Cd0 -k2.*CL(V) + k1.*CL(V).^2;
        D = @(V) 0.5.*rho.*V.^2.*S.*CD(V);
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        C = @(V) c_SL.*(1+1.2.*M(V)).*sqrt(Temp/Temp_SL);
                    case 2
                        C = @(V) c_SL.*(1+0.33.*M(V)).*sqrt(Temp/Temp_SL);
                    case 3
                        C = @(V) c_SL.*(1+0.16875.*M(V)).*sqrt(Temp/Temp_SL);
                end
            case 2
                C = @(V) c_SL.*(1+1.44.*M(V)).*sqrt(Temp/Temp_SL).* V./eta_p;
            case 3
                C = @(V) c_SL .*V./eta_p;
        end
        
        dvdt = @(V) -(C(V).*g.*D(V))./(dwdv(V) + C(V).*W(V));
        dt   = @(V) -(dwdv(V) + C(V).*W(V))./(C(V).*g.*D(V));
        dx   = @(V) V.*dt(V);
        
        T    = @(V) D(V) + (W(V)./g).*dvdt(V);
        
        switch tipo_motor
            case 1
                delta_T = @(V) T(V)./(T_SL*(1+0.2.*M(V).^2).^(1.4./0.4)*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
            case 2
                delta_T = @(V) T(V)./(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
            case 3
                delta_T = @(V) T(V)./(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
        end
        
        % Con esto ya tenemos todas las funciones dependiendo de V, solo resta
        % hallar las condiciones inicial y final de la velocidad usando la
        % expresion de V=V(W)
        
        syms vel positive
        V_ini = double(solve((W(vel)-W_ini),'Real',true));
        V_fin = double(solve((W(vel)-W_fin),'Real',true));
        
        distancia_int = @(V) integral(dx,V_ini,V);
        tiempo_int    = @(V) integral(dt,V_ini,V);
        
        n = 1000;
        VEL = linspace(V_ini,V_fin,n);
        
        %------------------------------------------------------------------------%
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        palanca1 = delta_T(VEL(1));
        palanca2 = delta_T(VEL(end));
        palanca_media = (palanca1 + palanca2)/2;
        correccion = (a_1.*palanca_media.^4+a_2.*palanca_media.^3+a_3.*palanca_media.^2+a_4.*palanca_media+a_5);
        %------------------------------------------------------------------------%
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        A = c_SL *correccion* sqrt(Temp/Temp_SL);
                        B = 1.2*A/a;
                    case 2
                        A = c_SL *correccion* sqrt(Temp/Temp_SL);
                        B = 0.33*A/a;
                    case 3
                        A = c_SL *correccion* sqrt(Temp/Temp_SL);
                        B = 0.16875*A/a;
                end
            case 2
                A = c_SL *correccion* sqrt(Temp/Temp_SL)/eta_p;
                B = 1.44*A/a;
            case 3
                A = c_SL*correccion/eta_p;
                B = 0;
        end
        
        A_2 = 0.5*rho*S*Cd0;
        B_2 = k2;
        C_2 = k1/(0.5*rho*S);
        
        switch tipo_motor
            case 1
                %W = @(V) sqrt((3.*A_2.*B.*V.^5 + 2.*A.*A_2.*V.^4 - B.*B_2.*V.^3)./(C_2.*B.*V + 2.*A.*C_2));
                %W = @(V) ( -B.*B_2.*V.^3 + sqrt((B.*B_2).^2.*V.^6 + 4.*(B.*C_2.*V + 2.*A.*C_2).*(2.*A.*A_2.*V.^4 + 3.*A_2.*B.*V.^5)))./(2.*(B.*C_2.*V + 2.*A.*C_2));
                
                AA = @(V) B*C_2.*V + 2*A*C_2; BB = @(V) B*B_2.*V.^3; CC = @(V) -(2*A*A_2.*V.^4 + 3*A_2*B.*V.^5);
                dAAdv = @(V) B*C_2; dBBdv = @(V) 3*B*B_2.*V.^2; dCCdv = @(V) -(8*A*A_2.*V.^3 + 15*A_2*B.*V.^4);
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
                %AA = 3*A_2*B; BB = 2*A*A_2; CC = B*B_2; DD = C_2*B; EE = 2*A*C_2;
                %AA = -B*B_2; BB = B^2*B_2^2 + 12*A_2*B^2*C_2; CC = 32*A*A_2*B*C_2; DD = 16*A^2*A_2*C_2; EE = 2*B*C_2; FF = 4*A*C_2;
                
                %dwdv = @(V) (V.^2.*(4.*AA.*DD.*V.^3 + 5.*AA.*EE.*V.^2 + 3.*BB.*DD.*V.^2 + 4.*BB.*EE.*V - 2.*CC.*DD.*V - 3.*CC.*EE))...
                %./(2.*(DD.*V+EE).^2.*V.*sqrt(V.*(AA.*V.^2+BB.*V-CC)./(DD.*V+EE)));
                %dwdv = @(V) AA.*V.^2.*(2*EE.*V+3*FF)./((EE.*V+FF).^2) - EE.*sqrt(V.^4*(BB.*V.^2+CC.*V+DD))...
                %./((EE.*V+FF).^2) + V.^3.*(6.*BB.*V.^2 + 5.*CC.*V + 4*DD)./(2.*(EE.*V+FF).*sqrt(V.^4.*(BB.*V.^2 + CC.*V + DD)));
                
            case 2
                %W = @(V) sqrt((4.*A_2.*B.*V.^5 + 3.*A.*A_2.*V.^4 - 2.*B.*B_2.*V.^3 - A.*B_2.*V.^2)./(A.*C_2));
                %W2 = @(V) ((V.^2)./(2.*A.*C_2)).*(-(A.*B_2 + 2.*B.*B_2.*V)+sqrt((A.*B_2 + 2.*B.*B_2.*V).^2+4.*A.*C_2.*(3.*A.*A_2 + 4.*A_2.*B.*V)));
                
                AA = @(V) A*C_2; BB = @(V) A*B_2.*V.^2 + 2*B*B_2.*V.^3; CC = @(V) -(3*A*A_2.*V.^4 + 4*A_2*B.*V.^5);
                dAAdv = @(V) 0; dBBdv = @(V) 2*A*B_2.*V + 6*B*B_2.*V.^2; dCCdv = @(V) -(12*A*A_2.*V.^3 + 20*A_2*B.*V.^4);
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
                %AA = 4*A_2*B; BB = 3*A*A_2; CC = 2*B*B_2; DD = B_2*A; EE = A*C_2;
                
                %dwdv = @(V) (V.*(5.*AA(V).*V.^3 + 4.*BB.*V.^2 - 3.*CC.*V - 2.*DD))...
                %./(2.*EE.*sqrt((4.*A_2.*B.*V.^5 + 3.*A.*A_2.*V.^4 - 2.*B.*B_2.*V.^3 - A.*B_2.*V.^2)./(A.*C_2)));
                
            case 3
                %         W = @(V) sqrt((3.*A_2.*V.^4-B_2.*V.^2)./C_2);
                %
                %         dwdv = @(V) (12.*A_2.*V.^3-2.*B_2.*V)./(2.*C_2.*sqrt((3.*A_2.*V.^4-B_2.*V.^2)./C_2));
                
                AA = @(V) C_2; BB = @(V) B_2.*V.^2; CC = @(V) -3*A_2.*V.^4;
                dAAdv = @(V) 0; dBBdv = @(V) 2*B_2.*V; dCCdv = @(V) -12*A_2.*V.^3;
                
                W = @(V) (-BB(V)+sqrt(BB(V).^2 - 4.*AA(V).*CC(V)))./(2.*AA(V));
                dwdv = @(V) -(dAAdv(V).*W(V).^2 + dBBdv(V).*W(V) + dCCdv(V))./(2.*AA(V).*W(V) + BB(V));
                
        end
        
        M = @(V) V./a;
        CL = @(V) W(V)./(0.5.*rho.*V.^2.*S);
        CD = @(V) Cd0 -k2.*CL(V) + k1.*CL(V).^2;
        D = @(V) 0.5.*rho.*V.^2.*S.*CD(V);
        
        switch tipo_motor
            case 1
                switch derivacion
                    case 1
                        C = @(V) c_SL.*correccion.*(1+1.2.*M(V)).*sqrt(Temp/Temp_SL);
                    case 2
                        C = @(V) c_SL.*correccion.*(1+0.33.*M(V)).*sqrt(Temp/Temp_SL);
                    case 3
                        C = @(V) c_SL.*correccion.*(1+0.16875.*M(V)).*sqrt(Temp/Temp_SL);
                end
            case 2
                C = @(V) c_SL.*correccion.*(1+1.44.*M(V)).*sqrt(Temp/Temp_SL).* V./eta_p;
            case 3
                C = @(V) c_SL .*correccion.*V./eta_p;
        end
        
        dvdt = @(V) -(C(V).*g.*D(V))./(dwdv(V) + C(V).*W(V));
        dt   = @(V) -(dwdv(V) + C(V).*W(V))./(C(V).*g.*D(V));
        dx   = @(V) V.*dt(V);
        
        T    = @(V) D(V) + (W(V)./g).*dvdt(V);
        
        switch tipo_motor
            case 1
                delta_T = @(V) T(V)./(T_SL*(1+0.2.*M(V).^2).^(1.4./0.4)*(1-0.49.*sqrt(M(V))).*(rho./rho_SL));
            case 2
                delta_T = @(V) T(V)./(T_SL.*eta_p./V.*(1+0.2.*M(V).^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
            case 3
                delta_T = @(V) T(V)./(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
        end
        
        % Con esto ya tenemos todas las funciones dependiendo de V, solo resta
        % hallar las condiciones inicial y final de la velocidad usando la
        % expresion de V=V(W)
        
        syms vel positive
        V_ini = double(solve((W(vel)-W_ini),'Real',true));
        V_fin = double(solve((W(vel)-W_fin),'Real',true));
        
        distancia_int = @(V) integral(dx,V_ini,V);
        tiempo_int    = @(V) integral(dt,V_ini,V);
        
        n = 1000;
        VEL = linspace(V_ini,V_fin,n);
        
        %------------------------------------------------------------------------%
        
        for k = 1:n,
            peso(k) = W(VEL(k));
            L_D(k) = CL(VEL(k))/CD(VEL(k));
            empuje(k) = T(VEL(k));
            CL2(k) = CL(VEL(k));
            palanca(k) = delta_T(VEL(k));
            tiempo(k) = tiempo_int(VEL(k));
            distancia(k) = distancia_int(VEL(k));
        end
        
        distancia_fin = distancia(end);
        autonomia_max= tiempo(end);
        distancia_crucero = distancia(end);
        fuel_crucero = (peso(1)-peso(end))./g;
        tiempo_crucero = tiempo(end);
        
        datos(i).crucero.tiempo = tiempo;
        datos(i).crucero.fuel   = (peso(1) - peso)./g;
        datos(i).crucero.distancia = distancia;
        datos(i).crucero.palanca = palanca;
        %datos(i).crucero.empuje = empuje;
        datos(i).crucero.empuje = CL2;
        datos(i).crucero.L_D = L_D;
        datos(i).crucero.peso = peso;
        datos(i).crucero.velocidad = VEL;
        
        
        %-------------------------------------------------------------------------%
        
end

datos(i).nombre = 'Crucero';
datos(i).crucero.altura = h_inicial;

if opcion == 1,
    datos(i).lista_variables = [{''};'Tiempo';'Fuel';'Peso';'Distancia';'Velocidad';'Palanca';'CL';...
        'L/D';'Altura';'Velocidad max'];
else
    datos(i).lista_variables = [{''};'Tiempo';'Fuel';'Peso';'Distancia';'Velocidad';'Palanca';'CL';...
        'L/D';'Altura'];
end

global flag_palanca
if max(datos(i).crucero.palanca) > 1.2,
    flag_palanca = 1;
    return;
else
    flag_palanca = 0;
end


end
