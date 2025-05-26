function [fuel_subida,mbaterias_subida,energia_subida,tiempo_subida,S_subida,datos] = analisis_subida(propul,aerodinamica,subida,W_inicial,h_inicial,opcion,datos,i,Prop_data,data_electric)

%% PROPUL
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR
% 4: CONSUMO ESPECIFICO
% 5: AVION CIVIL =1/MILITAR = 2
% 6: EFICIENCIA DE LA HELICE (ETA_P)
% 7: DERIVACION(TURBOFANES)

% AERODINAMICA
% 1: SUPERFICIE
% 2: CD0
% 3: K1 = K
% 4: K2
% 5: CLMAX LIMPIO

% SUBIDA
% 1: ALTURA FINAL       
% 2: GAMMA DE SUBIDA
% 3: MACH DE VUELO
% 4: VELOCIDAD TAS 
% 5: VELOCIDAD EAS 
% 6: PALANCA DE GASES
% 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)
% 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)

% subida(1) = subida_vec.temp_local;
% subida(2) = subida_vec.h_inicial;
% subida(3) = subida_vec.P_local;
% subida(4) = subida_vec.mu_takeoff;
% subida(5) = subida_vec.h_obstacle;
% subida(6) = subida_vec.gamma_climb;
% subida(7) = subida_vec.delta_T;

% AERODINAMICA
S = aerodinamica(1);
Cd0 = aerodinamica(2);
k1 = aerodinamica(3);
k2 = aerodinamica(4);
CLmax_limpio = aerodinamica(5);

% PROPULSION
tipo_motor = propul(1);
num_motores = propul(2);
T_SL = propul(3);
c_SL = propul(4);
civil_mil = propul(5);
eta_inst = propul(6);
derivacion = propul(7);

% ELECTRIC PROPULSION
e0 = data_electric(1)*3600; %from Wh to J
Dprop = data_electric(2); % in m
SF_prop = data_electric(3);
nps_max0 = data_electric(4);
nps_max = nps_max0*SF_prop;
eta_m = data_electric(7); % Will be replaced by excel
CT_Polyfit = Prop_data.CT_Polyfit;
CP_Polyfit = Prop_data.CP_Polyfit;
etamp_Polyfit = Prop_data.etamp_Polyfit;
tau = data_electric(8)/100;

% SUBIDA 
h_final = subida(1);
angulo = subida(2);
Mach = subida(3);
EAS = subida(4);
TAS = subida(5);
delta_T_gases = subida(6);
velocidad_ini = subida(7);
velocidad_fin = subida(8);

% ATMOSFERA ESTANDAR
[Temp_SL,rho_SL,p_SL,a_SL] = atmos_inter_mio(0);

% ----------  DEFINICION DE LAS FUNCIONES RESPECTO DE H --------------
g = 9.80665;
R = 287.058;                                                               
gamma_atm = 1.4;                                                               
lambda = 6.5*10^-3;                                                                                             

W_ini = W_inicial;
h_ini = h_inicial;

iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();   %ESTA LINEA PERMITE USAR CONDICIONES EN LAS FUNCIONES!!!

Temp = @(h) iif( h<11000, @() Temp_SL - lambda .* h, ...
                 h>=11000,@() 216.65);
rho  = @(h) iif( h<11000, @() rho_SL .*(Temp(h)./Temp_SL).^(g/(R*lambda)-1),...
                 h>=11000,@() 0.3639 .* exp((-g/(R*216.65)).*(h-11000)));
             
a = @(h) sqrt(gamma_atm*R.*Temp(h));

a_prima = @(h) iif (h<11000, @() -sqrt(R*gamma_atm).*(lambda./(2.*sqrt(Temp(h)))),...
                    h>=11000, @() 0);

rho_prima = @(h) iif (h<11000, @() -rho_SL.*(1./Temp_SL).^(g./(R.*lambda)-1).*(g./(R.*lambda)-1).*Temp(h).^(g./(R.*lambda)-2).*lambda,...
                      h>=11000, @() (-g/(R*216.65)).*rho(h));
                
EAS_prima = @(h) EAS.*sqrt(rho_SL).*(-0.5).*(rho(h).^(-1.5)).*rho_prima(h);
                
eta_p = eta_inst;
T_SL = T_SL * num_motores;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

switch opcion
    case 1
        %SUBIDA DADO GAMMA = CONSTANTE, MACH = CONSTANTE, ALTURA FINAL
        M = Mach;
        gamma = angulo;
        V = @(h) M .* a(h);
        
        if tipo_motor==4 %Para el modelo eléctrico no es necesario resolver ecuaciones diferenciales puesto que tendremos
            %un valor de V, CL, CD,y T para cada incremento de altura.
            n = 100;%Numero de intervalos para discretizar la altura
            vector_h = linspace(h_ini,h_final,n);
            Vvect=V(vector_h);
            rhovect=rho(vector_h);
            aprimavect=a_prima(vector_h);
            %Una vez calculada T, se resuelve con un fzero las revoluciones
            %necesarias para dar ese empuje, se calculan palanca, rendimiento y
            %potencia, y se calcula dEdh=P(h)/(V(h)*sin(gamma))
            CL =  W_ini.*cos(gamma)./(0.5.*rhovect.*(Vvect).^2*S);  % L  = W COS(GAMMA)
            CD =  Cd0 - k2.*CL + k1.*(CL.^2);    %POLAR
            D =  0.5.*rhovect.*(Vvect.^2)*S.*CD;     % D
            T =  D + W_ini.*sin(gamma) + (W_ini./g).*Vvect.*sin(gamma).*M.*aprimavect; %T = D + W SIN(GAMMA) + W/G * V *SIN(GAMMA)* DV/DH
        else
            CL = @(h,W) W.*cos(gamma)./(0.5.*rho(h).*(V(h)).^2*S);  % L  = W COS(GAMMA)
            CD = @(h,W) Cd0 - k2.*CL(h,W) + k1.*(CL(h,W).^2);    %POLAR
            D = @(h,W) 0.5.*rho(h).*(V(h).^2)*S.*CD(h,W);     % D
            T = @(h,W) D(h,W) + W.*sin(gamma) + (W./g).*V(h).*sin(gamma).*M.*a_prima(h); %T = D + W SIN(GAMMA) + W/G * V *SIN(GAMMA)* DV/DH
        end
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        
        switch tipo_motor
            case 1
                delta_T = @(h,W) T(h,W)./(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49*sqrt(M)).*(rho(h)./rho_SL));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                switch derivacion
                    case 1
                        C = @(h,W) c_SL.*correccion(h,W).*(1+1.2.*M).*sqrt(Temp(h)./Temp_SL);
                    case 2
                        C = @(h,W) c_SL.*correccion(h,W).*(1+0.33.*M).*sqrt(Temp(h)./Temp_SL);
                    case 3
                        C = @(h,W) c_SL.*correccion(h,W).*(1+0.16875.*M).*sqrt(Temp(h)./Temp_SL);
                end
            case 2
                delta_T  = @(h,W) T(h,W)./(T_SL.*eta_p./V(h).*(1+0.2.*M^2).^(1.4./0.4).*(rho(h).*Temp(h)./(rho_SL*Temp_SL)));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                C = @(h,W) c_SL.*correccion(h,W).*(1+1.44.*M).*sqrt(Temp(h)./Temp_SL).* V(h)./eta_p;
            case 3
                delta_T = @(h,W) T(h,W)./(T_SL.*eta_p./V(h).*((8.55.*rho(h)./rho_SL - 1)./7.55));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                C = @(h,W) c_SL .*correccion(h,W).*V(h)./eta_p;
                
                
            case 4
                
                nps_sol=ones(n,1);
                
                for m=1:n
                    J=@(npsi) Vvect(m)./(npsi*Dprop);
                    CT=@(npsi) 0;
                    
                    for j=1:length(CT_Polyfit)
                        CT_int=@(npsi) CT_Polyfit(j)*J(npsi).^(length(CT_Polyfit)-j);
                        CT=@(npsi) CT_int(npsi)+CT(npsi);%CT_Polyfit(1)+CT_Polyfit(2)*J(nps)+CT_Polyfit(3)*J(nps)^2;
                    end
                    
                    CP=@(npsi) 0;
                    
                    for j=1:length(CP_Polyfit)
                        CP_int=@(npsi) CP_Polyfit(j).*J(npsi).^(length(CP_Polyfit)-j);
                        CP=@(npsi) CP_int(npsi)+CP(npsi);
                    end
                    
                    eta_p=@(npsi) 0;
                    for j=1:length(etamp_Polyfit)
                        eta_p_int=@(npsi) etamp_Polyfit(j).*J(npsi).^(length(etamp_Polyfit)-j);
                        eta_p=@(npsi) eta_p_int(npsi)+eta_p(npsi);
                    end
                    
                    nps_sol(m,1)=  fzero(@(npsi) T(m)-num_motores.*CT(npsi).*rhovect(m).*(npsi^2)*Dprop^4,nps_max);
                    %T=@(h,nps) num_motores*CT(vector_h,nps).*rho(vector_h).*(nps^3).*(Dprop)^5;
                end
                P_sol= num_motores*CP(nps_sol).*rhovect'.*(nps_sol.^3).*(Dprop)^5;
                delta_T=nps_sol/nps_max;
                
        end
        
        % Sabiendo que dh/dt = V*sin(gamma), integrando y con la condicion inicial
        % de t(h = h_ini) = 0, obtenemos esta expresion para el tiempo = tiempo(h)
        
        Constante = 2/(M*sin(gamma)*lambda*sqrt(gamma_atm*R));
        
        tiempo = @(h) iif (h<11000, @() Constante.*(sqrt(Temp_SL - lambda * h_ini)-sqrt(Temp_SL - lambda .* h)),...
            (h-11000).*(h_ini-11000)<0,...
            @() Constante.*(sqrt(Temp_SL - lambda * h_ini)-sqrt(Temp_SL - lambda .* 11000))+(h-11000)./(M.*a(h).*sin(gamma)),...
            (h-11000).*(h_ini-11000)>=0, @() (h-h_ini)./(M.*a(h).*sin(gamma)));
        
        % Usamos dW/dt = dW/dh * dh/dt = dW/dh * V*sin(gamma) = -C(h)T(h,W)
        if tipo_motor==4
            dEdh= P_sol./(Vvect'*sin(gamma).*eta_p(nps_sol)*eta_m); %dEdh=dEdx*dxdh=dEdt*dxdt*1/tan(gamma)=P/(V*tan(gamma))
            ENERGIA=dEdh(:).*(vector_h(:)-h_ini);
            PESO=ENERGIA/e0;
            Newt=PESO*g;
        else
            dwdh = @(h,W) -C(h,W).*g.*T(h,W)./(V(h).*sin(gamma));
            % Condicion inicial
            V_ini = V(h_ini);
            T_ini = T(h_ini,W_ini);
            delta_T_ini = delta_T(h_ini,W_ini);
            n = 500; %Numero de intervalos para discretizar la altura
            
            % Resolvemos la ecuacion diferencial para obtener W(h) y por tanto
            % delta_T(h), aunque para obtener las magnitudes globales solo haria falta
            % integrar una vez con h=h_final
            
            vector_h = linspace(h_ini,h_final,n);
            
            [H,Newt] = ode45(dwdh,vector_h,W_ini);
        end
        k = 1;
        while k <= n
            if tipo_motor==4
                H(k)=vector_h(k);
                palanca(k)=delta_T(k);
                empuje(k)= T(k);
                V_v(k)= Vvect(k)*sin(gamma);
                velocidad(k)=Vvect(k)*cos(gamma);
                CL2(k)= CL(k);
                CD2(k)= CD(k);
                L_D(k)= CL2(k)./CD2(k);
                
            else
                palanca(k) = delta_T(H(k),Newt(k));
                empuje(k) = T(H(k),Newt(k));
                V_v(k) = V(H(k))*sin(gamma);
                velocidad(k) = V(H(k)) * cos(gamma);
                CL2(k) = CL(H(k),Newt(k));
                CD2(k) = CD(H(k),Newt(k));
                L_D(k) = CL2(k)./CD2(k);
                % if palanca(k) > 1,
                % %     techo = H(k);
                % %     msgbox('No es posible alcanzar la altura final definida','Mision','Warn');
                % else
                %     techo = 'No se ha alcanzado el techo';
                % end
            end
            k = k + 1;
        end
        
        tiempo_subida = tiempo(h_final);
        S_subida = (h_final-h_ini)/(tan(gamma));
        W = Newt(end);

        datos(i).segmento.tiempo = (H-H(1))./V_v;
        if tipo_motor==4
            fuel_subida = 0;
            energia_subida = ENERGIA(end);
            mbaterias_subida = (PESO(end));
            datos(i).segmento.fuel = 0; %fuel consumido acumulado
%             datos(i).segmento.mbaterias=(PESO(1) - PESO);
            datos(i).segmento.mbaterias = PESO;
            datos(i).segmento.tiempo = datos(i).segmento.tiempo'; %Para el código de ploteo
%             Newt = W_ini + W*ones(1,length(datos(i).segmento.tiempo)); %Para el código de ploteo
            Newt = W_ini*ones(1,length(datos(i).segmento.tiempo));
        else
            fuel_subida = (Newt(1)-Newt(end))/g;
            energia_subida=0;
            mbaterias_subida=0;
            datos(i).segmento.fuel = (Newt(1) - Newt)./g; %fuel consumido acumulado
            datos(i).segmento.mbaterias=0;
        end
        
        datos(i).segmento.distancia = (H-H(1))./tan(gamma);
        datos(i).segmento.palanca = palanca;
        datos(i).segmento.empuje = empuje;
%         datos(i).segmento.empuje = CL2;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D;
        datos(i).segmento.peso = Newt;
        % datos(i).segmento.techo = techo;
        datos(i).segmento.gamma = gamma; %%% ESCALAR %%%
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.altura = H;
        if tipo_motor==4
            datos(i).segmento.energiav = ENERGIA;
        end
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        
    case 2
        % SUBIDA DADO GAMMA = CONSTANTE, EAS = CONSTANTE, ALTURA FINAL
        
        % TAS = EAS * SQRT(RHO_SL/RHO)
        gamma = angulo;
        V = @(h) EAS .* sqrt(rho_SL./rho(h));
        M = @(h) V(h)./a(h);
        
        CL = @(h,W) W.*cos(gamma)./(0.5.*rho(h).*(V(h)).^2*S);  % L  = W COS(GAMMA)
        CD = @(h,W) Cd0 - k2.*CL(h,W) + k1.*(CL(h,W).^2);    %POLAR
        D = @(h,W) 0.5.*rho(h).*(V(h).^2)*S.*CD(h,W);     % D
        T = @(h,W) D(h,W) + W.*sin(gamma) + (W./g).*V(h).*sin(gamma).*EAS_prima(h); %T = D + W SIN(GAMMA) + W/G * V *SIN(GAMMA)* DV/DH
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        
        switch tipo_motor
            case 1
                delta_T = @(h,W) T(h,W)./(T_SL.*(1+0.2.*M(h).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(h))).*(rho(h)./rho_SL));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                switch derivacion
                    case 1
                        C = @(h,W) c_SL.*correccion(h,W).*(1+1.2.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 2
                        C = @(h,W) c_SL.*correccion(h,W).*(1+0.33.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 3
                        C = @(h,W) c_SL.*correccion(h,W).*(1+0.16875.*M(h)).*sqrt(Temp(h)./Temp_SL);
                end
            case 2
                delta_T  = @(h,W) T(h,W)./(T_SL.*eta_p./V(h).*(1+0.2.*M(h).^2).^(1.4./0.4).*(rho(h).*Temp(h)./(rho_SL*Temp_SL)));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                C = @(h,W) c_SL.*correccion(h,W).*(1+1.44.*M(h)).*sqrt(Temp(h)./Temp_SL).* V(h)./eta_p;
            case 3
                delta_T = @(h,W) T(h,W)./(T_SL.*eta_p./V(h).*((8.55.*rho(h)./rho_SL - 1)./7.55));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                C = @(h,W) c_SL .*correccion(h,W).*V(h)./eta_p;
        end
        
        % Sabiendo que dh/dt = V*sin(gamma), integrando y con la condicion inicial
        % de t(h = h_ini) = 0, obtenemos esta expresion para el tiempo = tiempo(h)
        
        Constante_1 = (1/(EAS*sin(gamma)*0.5*lambda*(g/(R*lambda)+1)))*(1/Temp_SL)^(0.5*(g/(R*lambda)-1));
        Constante_2 = (0.3639*R*216.65/(EAS*g*sin(gamma)*sqrt(rho_SL)));
        
        tiempo = @(h) iif (h<11000, @() Constante_1.*((Temp_SL -lambda.*h_ini).^(0.5*(g/(R*lambda)+1)) - (Temp_SL-lambda.*h).^(0.5*(g/(R*lambda)+1))),...
            (h-11000).*(h_ini-11000)<0,...
            @() Constante_1.*((Temp_SL -lambda.*h_ini).^(0.5*(g/(R*lambda)+1)) - (Temp_SL-lambda.*11000).^(0.5*(g/(R*lambda)+1))) + ...
            Constante_2.*(1-exp(g./(R.*216.65).*(11000-h))),...
            (h-11000).*(h_ini-11000)>=0, @() Constante_2.*(exp(g/(R*216.65)*(11000-h_ini))-exp(g./(R.*216.65).*(11000-h))));
        
        % Usamos dW/dt = dW/dh * dh/dt = dW/dh * V*sin(gamma) = -C(h)T(h,W)
        
        dwdh = @(h,W) -C(h,W).*g.*T(h,W)./(V(h).*sin(gamma));
        
        % Condicion inicial
        V_ini = V(h_ini);
        T_ini = T(h_ini,W_ini);
        delta_T_ini = delta_T(h_ini,W_ini);
        n = 500; %Numero de intervalos para discretizar la altura
        
        % Resolvemos la ecuacion diferencial para obtener W(h) y por tanto
        % delta_T(h), aunque para obtener las magnitudes globales solo haria falta
        % integrar una vez con h=h_final
        
        vector_h = linspace(h_ini,h_final,n);
        
        [H,Newt] = ode45(dwdh,vector_h,W_ini);
        k = 1;
        while k <= n
            palanca(k) = delta_T(H(k),Newt(k));
            empuje(k) = T(H(k),Newt(k));
            V_v(k) = V(H(k))*sin(gamma);
            velocidad(k) = V(H(k)) * cos(gamma);
            L_D(k) = CL(H(k),Newt(k))./CD(H(k),Newt(k));
            CL2(k) = CL(H(k),Newt(k));
            CD2(k) = CD(H(k),Newt(k));
            % if palanca(k) > 1,
            %     techo = H(k);
            %     msgbox('No es posible alcanzar la altura final definida','Mision','Warn');
            % else
            %     techo = 'No se ha alcanzado el techo';
            % end
            k = k + 1;
        end
        
        tiempo_subida = tiempo(h_final);
        S_subida = (h_final-h_ini)/(tan(gamma));
        fuel_subida = (Newt(1)-Newt(end))/g;
        energia_subida=0;
        mbaterias_subida=0;
        W = Newt(end);
        
        datos(i).segmento.tiempo = (H-H(1))./V_v;
        datos(i).segmento.fuel = (Newt(1) - Newt)./g; %fuel consumido acumulado
        datos(i).segmento.distancia = (H-H(1))./tan(gamma);
        datos(i).segmento.palanca = palanca;
        datos(i).segmento.empuje = empuje;
%         datos(i).segmento.empuje = CL2;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D;
        datos(i).segmento.peso = Newt;
        % datos(i).segmento.techo = techo;
        datos(i).segmento.gamma = gamma; %%% ESCALAR %%%
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.altura = H;
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        
    case 3
        % SUBIDA DADO GAMMA = CONSTANTE, V = TAS = CONSTANTE, ALTURA FINAL
        gamma = angulo;
        V = @(h) TAS;
        M = @(h) V(h)./a(h);
        
        if tipo_motor==4 %Para el modelo eléctrico no es necesario resolver ecuaciones diferenciales puesto que tendremos
            %un valor de M, CL, CD,y T para cada incremento de altura.
            n = 100;%Numero de intervalos para discretizar la altura
            vector_h = linspace(h_ini,h_final,n);
            V = TAS;
            Mvect=M(vector_h);
            rhovect=rho(vector_h);
            aprimavect=a_prima(vector_h);
            %Una vez calculada T, se resuelve con un fzero las revoluciones
            %necesarias para dar ese empuje, se calculan palanca, rendimiento y
            %potencia, y se calcula dEdh=P(h)/(V(h)*sin(gamma))
            CL =  W_ini.*cos(gamma)./(0.5*rhovect*(V^2)*S);  % L  = W COS(GAMMA)
            CD =  Cd0 - k2.*CL + k1.*(CL.^2);    %POLAR
            D =  0.5.*rhovect.*(V^2)*S.*CD;     % D
            T =  D + W_ini.*sin(gamma) + (W_ini./g).*V.*sin(gamma).*Mvect.*aprimavect; %T = D + W SIN(GAMMA) + W/G * V *SIN(GAMMA)* DV/DH
        else
            CL = @(h,W) W.*cos(gamma)./(0.5.*rho(h).*(V(h)).^2*S);  % L  = W COS(GAMMA)
            CD = @(h,W) Cd0 - k2.*CL(h,W) + k1.*(CL(h,W).^2);    %POLAR
            D = @(h,W) 0.5.*rho(h).*(V(h).^2)*S.*CD(h,W);     % D
            T = @(h,W) D(h,W) + W.*sin(gamma); %T = D + W SIN(GAMMA) + W/G * V *SIN(GAMMA)* DV/DH
        end
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        
        switch tipo_motor
            case 1
                delta_T = @(h,W) T(h,W)./(T_SL.*(1+0.2.*M(h).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(h))).*(rho(h)./rho_SL));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                switch derivacion
                    case 1
                        C = @(h,W) c_SL.*correccion(h,W).*(1+1.2.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 2
                        C = @(h,W) c_SL.*correccion(h,W).*(1+0.33.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 3
                        C = @(h,W) c_SL.*correccion(h,W).*(1+0.16875.*M(h)).*sqrt(Temp(h)./Temp_SL);
                end
            case 2
                delta_T  = @(h,W) T(h,W)./(T_SL.*eta_p./V(h).*(1+0.2.*M(h).^2).^(1.4./0.4).*(rho(h).*Temp(h)./(rho_SL*Temp_SL)));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                C = @(h,W) c_SL.*correccion(h,W).*(1+1.44.*M(h)).*sqrt(Temp(h)./Temp_SL).* V(h)./eta_p;
            case 3
                delta_T = @(h,W) T(h,W)./(T_SL.*eta_p./V(h).*((8.55.*rho(h)./rho_SL - 1)./7.55));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                C = @(h,W) c_SL .*correccion(h,W).*V(h)./eta_p;
            case 4
                nps_sol=ones(n,1);
                
                
                    J=@(npsi) V./(npsi*Dprop);
                    CT=@(npsi) 0;
                    
                    for j=1:length(CT_Polyfit)
                        CT_int=@(npsi) CT_Polyfit(j)*J(npsi).^(length(CT_Polyfit)-j);
                        CT=@(npsi) CT_int(npsi)+CT(npsi);%CT_Polyfit(1)+CT_Polyfit(2)*J(nps)+CT_Polyfit(3)*J(nps)^2;
                    end
                    
                    CP=@(npsi) 0;
                    
                    for j=1:length(CP_Polyfit)
                        CP_int=@(npsi) CP_Polyfit(j).*J(npsi).^(length(CP_Polyfit)-j);
                        CP=@(npsi) CP_int(npsi)+CP(npsi);
                    end
                    
                    eta_p=@(npsi) 0;
                    for j=1:length(etamp_Polyfit)
                        eta_p_int=@(npsi) etamp_Polyfit(j).*J(npsi).^(length(etamp_Polyfit)-j);
                        eta_p=@(npsi) eta_p_int(npsi)+eta_p(npsi);
                    end
                    
                for m=1:n
                   
                    nps_sol(m,1)=  fzero(@(npsi) T(m)-num_motores.*CT(npsi).*rhovect(m).*(npsi^2)*Dprop^4,nps_max);
                    %T=@(h,nps) num_motores*CT(vector_h,nps).*rho(vector_h).*(nps^3).*(Dprop)^5;
                end
                P_sol= num_motores*CP(nps_sol).*rhovect'.*(nps_sol.^3).*(Dprop)^5;
                delta_T=nps_sol/nps_max;
        end
        
        
        if tipo_motor==4
            dEdh= P_sol./(V'*sin(gamma).*eta_p(nps_sol)*eta_m); %dEdh=dEdx*dxdh=dEdt*dxdt*1/tan(gamma)=P/(V*tan(gamma))
            ENERGIA=dEdh(:).*(vector_h(:)-h_ini);
            PESO=ENERGIA/e0;
            Newt=PESO*g;
            tiempo = (vector_h-h_ini)./(V*sin(gamma));
        else
        % Sabiendo que dh/dt = V*sin(gamma), integrando y con la condicion inicial
        % de t(h = h_ini) = 0, obtenemos esta expresion para el tiempo = tiempo(h)
        
        tiempo = @(h) (h-h_ini)./(V(h)*sin(gamma));
        
        % Usamos dW/dt = dW/dh * dh/dt = dW/dh * V*sin(gamma) = -C(h)T(h,W)
        
        dwdh = @(h,W) -C(h,W).*g.*T(h,W)./(V(h).*sin(gamma));
        
        % Condicion inicial
        V_ini = V(h_ini);
        T_ini = T(h_ini,W_ini);
        delta_T_ini = delta_T(h_ini,W_ini);
        n = 500; %Numero de intervalos para discretizar la altura
        
        % Resolvemos la ecuacion diferencial para obtener W(h) y por tanto
        % delta_T(h), aunque para obtener las magnitudes globales solo haria falta
        % integrar una vez con h=h_final
        
        vector_h = linspace(h_ini,h_final,n);
        
        [H,Newt] = ode45(dwdh,vector_h,W_ini);
        end
        
        k = 1;
        while k <= n
            if tipo_motor==4
                H(k)=vector_h(k);
                palanca(k)=delta_T(k);
                empuje(k)= T(k);
                V_v(k)= V*sin(gamma);
                velocidad(k)=V*cos(gamma);
                CL2(k)= CL(k);
                CD2(k)= CD(k);
                L_D(k)= CL2(k)./CD2(k);
                
            else
            palanca(k) = delta_T(H(k),Newt(k));
            empuje(k) = T(H(k),Newt(k));
            V_v(k) = V(H(k))*sin(gamma);
            velocidad(k) = V(H(k)) * cos(gamma);
            L_D(k) = CL(H(k),Newt(k))./CD(H(k),Newt(k));
            CL2(k) = CL(H(k),Newt(k));
            CD2(k) = CD(H(k),Newt(k));
            % if palanca(k) > 1,
            %     techo = H(k);
            %     msgbox('No es posible alcanzar la altura final definida','Mision','Warn');
            % else
            %     techo = 'No se ha alcanzado el techo';
            % end
            end
            k = k + 1;
        end
        
        tiempo_subida = tiempo(h_final);
        S_subida = (h_final-h_ini)/(tan(gamma));
        datos(i).segmento.tiempo = (H-H(1))./V_v;
        W = Newt(end);
        
        if tipo_motor==4
            fuel_subida = 0;
            energia_subida = ENERGIA(end);
            mbaterias_subida = (PESO(end));
            datos(i).segmento.fuel = 0; %fuel consumido acumulado
%             datos(i).segmento.mbaterias=(PESO(1) - PESO);
            datos(i).segmento.mbaterias = PESO;
            datos(i).segmento.tiempo = datos(i).segmento.tiempo'; %Para el código de ploteo
%             Newt = W_ini + W*ones(1,length(datos(i).segmento.tiempo)); %Para el código de ploteo
            Newt = W_ini*ones(1,length(datos(i).segmento.tiempo));
        else
            datos(i).segmento.fuel = (Newt(1) - Newt)./g; %fuel consumido acumulado
            datos(i).segmento.mbaterias = 0;
            fuel_subida = (Newt(1)-Newt(end))/g;
            energia_subida=0;
            mbaterias_subida=0;
        end
        
        datos(i).segmento.distancia = (H-H(1))./tan(gamma);
        datos(i).segmento.palanca = palanca;
        datos(i).segmento.empuje = empuje;
%         datos(i).segmento.empuje = CL2;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D;
        datos(i).segmento.peso = Newt;
        % datos(i).segmento.techo = techo;
        datos(i).segmento.gamma = gamma; %%% ESCALAR %%%
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.altura = H;
        if tipo_motor==4
            datos(i).segmento.energiav = ENERGIA;
        end
        
        % ------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        
    case 4
        % SUBIDA DADOS PALANCA DE GASES = CONSTANTE, MACH = CONSTANTE, ALTURA FINAL
        M = Mach;
        delta_T = delta_T_gases;
        
        V = @(h) M .* a(h);
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        
        switch tipo_motor
            case 1
                T = @(h) delta_T.*(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49*sqrt(M)).*(rho(h)./rho_SL));
                switch derivacion
                    case 1
                        C = @(h) c_SL.*correccion.*(1+1.2.*M).*sqrt(Temp(h)./Temp_SL);
                    case 2
                        C = @(h) c_SL.*correccion.*(1+0.33.*M).*sqrt(Temp(h)./Temp_SL);
                    case 3
                        C = @(h) c_SL.*correccion.*(1+0.16875.*M).*sqrt(Temp(h)./Temp_SL);
                end
            case 2
                T = @(h) delta_T.*(T_SL.*eta_p./V(h).*(1+0.2.*M^2).^(1.4./0.4).*(rho(h).*Temp(h)./(rho_SL*Temp_SL)));
                C = @(h) c_SL.*correccion.*(1+1.44.*M).*sqrt(Temp(h)./Temp_SL).* V(h)./eta_p;
            case 3
                T = @(h) delta_T.*(T_SL.*eta_p./V(h).*((8.55.*rho(h)./rho_SL - 1)./7.55));
                C = @(h) c_SL .*correccion.*V(h)./eta_p;
        end
        
        % PROBLEMA: No veo forma de implementar las ecuaciones de la mecanica del
        % vuelo directamente: L = W cos(gamma) + W/g V^2 sin(gamma) dgamma/dh,
        % porque despejar gamma implica resolver una ecuacion diferencial y
        % posteriormente otra algebraica de 4 grado, lo cual aunque es posible no
        % veo muy eficiente en principio, y la variacion de gamma es muy muy
        % pequeña de todas formas.
        
        % De modo que se empleará un método iterativo para resolver gamma(h,W)
        
        % Condicion inicial
        V_ini = V(h_ini);
        T_ini = T(h_ini);
        gamma_ini = 0;
        n = 1000; %Numero de intervalos para discretizar la altura
        
        H = linspace(h_ini,h_final,n);
        w(1) = W_ini;
        gamma(1) = gamma_ini;
        
        for k = 1:(n-1)
            contador = 0;
            while contador <= 1000
                CL =  w(k)*cos(gamma(k))/(0.5*rho(H(k))*(V(H(k)))^2*S);  % L = W cos gamma
                CD =  Cd0 - k2*CL + k1*(CL^2);    %POLAR
                acel = 1 + (V(H(k))/g)*M*a_prima(H(k));
                cte_a = 1 + (CL/CD)^2*acel^2;
                cte_b = 2*(T(H(k))/w(k))*acel*(CL/CD)^2;
                cte_c = (CL/CD)^2*(T(H(k))/w(k))^2-1;
                seno_gamma = (cte_b - sqrt(cte_b^2-4*cte_a*cte_c))/(2*cte_a);
                gamma_obt = asin(seno_gamma);
                if abs(gamma(k) - gamma_obt) <= 10^-4,
                    contador = 1001;
                end
                contador = contador + 1;
                gamma(k) = gamma_obt;
            end
            % Aqui ya tendriamos CL,CD y gamma para cada h y w
            L_D(k) = CL/CD;
            palanca(k) = delta_T;
            CL2(k) = CL;
            CD2(k) = CD;
            V_v(k) = V(H(k)) * sin(gamma(k));
            velocidad(k) = V(H(k)) * cos(gamma(k));
            tiempo(k) = (H(k+1)-H(k))/V_v(k);
            distancia(k) = (H(k+1)-H(k))/tan(gamma(k));
            empuje(k) = T(H(k));
            fuel(k) = C(H(k))*T(H(k))*tiempo(k);
            w(k+1) = w(k) - fuel(k) * g;      %Calculamos el nuevo peso que tendra la siguiente h
            gamma(k+1) = gamma(k);
        end
        
        tiempo_subida = sum(tiempo);
        S_subida = sum(distancia);
        fuel_subida = sum(fuel);
        energia_subida=0;
        mbaterias_subida=0;
        W = W_ini - fuel_subida*g;
        
        datos(i).segmento.tiempo = cumsum(tiempo);
        datos(i).segmento.fuel = cumsum(fuel); %fuel consumido acumulado
        datos(i).segmento.distancia = cumsum(distancia);
        datos(i).segmento.palanca = palanca;%%% ESCALAR %%%
        datos(i).segmento.empuje = empuje;
%         datos(i).segmento.empuje = CL2;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D; %%% ESCALAR %%%
        datos(i).segmento.peso = w;
        datos(i).segmento.peso(end) = '';
        datos(i).segmento.gamma = gamma;
        datos(i).segmento.gamma(end) = '';
        %datos(i).segmento.techo = techo;
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.altura = H(1:end-1);
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        
    case 5
        % SUBIDA DADOS PALANCA DE GASES = CONSTANTE, EAS = CONSTANTE, ALTURA FINAL
        delta_T = delta_T_gases;
        
        % TAS = EAS * SQRT(RHO_SL/RHO)
        
        V = @(h) EAS .* sqrt(rho_SL./rho(h));
        M = @(h) V(h)./a(h);
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        
        switch tipo_motor
            case 1
                T = @(h) delta_T.*(T_SL*(1+0.2.*M(h).^2).^(1.4./0.4)*(1-0.49.*sqrt(M(h))).*(rho(h)./rho_SL));
                switch derivacion
                    case 1
                        C = @(h) c_SL.*correccion.*(1+1.2.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 2
                        C = @(h) c_SL.*correccion.*(1+0.33.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 3
                        C = @(h) c_SL.*correccion.*(1+0.16875.*M(h)).*sqrt(Temp(h)./Temp_SL);
                end
            case 2
                T = @(h) delta_T.*(T_SL.*eta_p./V(h).*(1+0.2.*M(h).^2).^(1.4./0.4).*(rho(h).*Temp(h)./(rho_SL*Temp_SL)));
                C = @(h) c_SL.*correccion.*(1+1.44.*M(h)).*sqrt(Temp(h)./Temp_SL).* V(h)./eta_p;
            case 3
                T = @(h) delta_T.*(T_SL.*eta_p./V(h).*((8.55.*rho(h)./rho_SL - 1)./7.55));
                C = @(h) c_SL .*correccion.*V(h)./eta_p;
        end
        
        % PROBLEMA: No veo forma de implementar las ecuaciones de la mecanica del
        % vuelo directamente: L = W cos(gamma) + W/g V^2 sin(gamma) dgamma/dh,
        % porque despejar gamma implica resolver una ecuacion diferencial y
        % posteriormente otra algebraica de 4 grado, lo cual aunque es posible no
        % veo muy eficiente en principio, y la variacion de gamma es muy muy
        % pequeña de todas formas.
        
        % De modo que se empleará un método iterativo para resolver gamma(h,W)
        
        % Condicion inicial
        V_ini = V(h_ini);
        T_ini = T(h_ini);
        gamma_ini = 0;
        n = 1000; %Numero de intervalos para discretizar la altura
        
        H = linspace(h_ini,h_final,n);
        w(1) = W_ini;
        gamma(1) = gamma_ini;
        
        for k = 1:(n-1),
            contador = 0;
            while contador <= 1000,
                CL =  w(k)*cos(gamma(k))/(0.5*rho(H(k))*(V(H(k)))^2*S);  % L = W cos gamma
                CD =  Cd0 - k2*CL + k1*(CL^2);    %POLAR
                acel = 1 + (V(H(k))/g)*EAS_prima(H(k));
                cte_a = 1 + (CL/CD)^2*acel^2;
                cte_b = 2*(T(H(k))/w(k))*acel*(CL/CD)^2;
                cte_c = (CL/CD)^2*(T(H(k))/w(k))^2-1;
                seno_gamma = (cte_b - sqrt(cte_b^2-4*cte_a*cte_c))/(2*cte_a);
                gamma_obt = asin(seno_gamma);
                if abs(gamma(k) - gamma_obt) <= 10^-4,
                    contador = 1001;
                end
                contador = contador + 1;
                gamma(k) = gamma_obt;
            end
            % Aqui ya tendriamos CL,CD y gamma para cada h y w
            L_D(k) = CL/CD;
            palanca(k) = delta_T;
            CL2(k) = CL;
            CD2(k) = CD;
            V_v(k) = V(H(k)) * sin(gamma(k));
            velocidad(k) = V(H(k)) * cos(gamma(k));
            tiempo(k) = (H(k+1)-H(k))/V_v(k);
            distancia(k) = (H(k+1)-H(k))/tan(gamma(k));
            empuje(k) = T(H(k));
            fuel(k) = C(H(k))*T(H(k))*tiempo(k);
            w(k+1) = w(k) - fuel(k) * g;      %Calculamos el nuevo peso que tendra la siguiente h
            gamma(k+1) = gamma(k);
        end
        
        tiempo_subida = sum(tiempo);
        S_subida = sum(distancia);
        fuel_subida = sum(fuel);
        energia_subida=0;
        mbaterias_subida=0;
        W = W_ini - fuel_subida*g;
        
        datos(i).segmento.tiempo = cumsum(tiempo);
        datos(i).segmento.fuel = cumsum(fuel); %fuel consumido acumulado
        datos(i).segmento.distancia = cumsum(distancia);
        datos(i).segmento.palanca = palanca;%%% ESCALAR %%%
        datos(i).segmento.empuje = empuje;
%         datos(i).segmento.empuje = CL2;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D; %%% ESCALAR %%%
        datos(i).segmento.peso = w;
        datos(i).segmento.peso(end) = '';
        datos(i).segmento.gamma = gamma;
        datos(i).segmento.gamma(end) = '';
        %datos(i).segmento.techo = techo;
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.altura = H(1:end-1);
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        
    case 6
        % SUBIDA DADOS PALANCA DE GASES = CONSTANTE, TAS = CONSTANTE, ALTURA FINAL
        delta_T = delta_T_gases;
        
        % TAS = EAS * SQRT(RHO_SL/RHO)
        
        V = @(h) TAS;
        M = @(h) V(h)./a(h);
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        
        switch tipo_motor
            case 1
                T = @(h) delta_T.*(T_SL.*(1+0.2.*M(h).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(h))).*(rho(h)./rho_SL));
                switch derivacion
                    case 1
                        C = @(h) c_SL.*correccion.*(1+1.2.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 2
                        C = @(h) c_SL.*correccion.*(1+0.33.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 3
                        C = @(h) c_SL.*correccion.*(1+0.16875.*M(h)).*sqrt(Temp(h)./Temp_SL);
                end
            case 2
                T = @(h) delta_T.*(T_SL.*eta_p./V(h).*(1+0.2.*M(h).^2).^(1.4./0.4).*(rho(h).*Temp(h)./(rho_SL*Temp_SL)));
                C = @(h) c_SL.*correccion.*(1+1.44.*M(h)).*sqrt(Temp(h)./Temp_SL).* V(h)./eta_p;
            case 3
                T = @(h) delta_T.*(T_SL.*eta_p./V(h).*((8.55.*rho(h)./rho_SL - 1)./7.55));
                C = @(h) c_SL .*correccion.*V(h)./eta_p;
        end
        
        % PROBLEMA: No veo forma de implementar las ecuaciones de la mecanica del
        % vuelo directamente: L = W cos(gamma) + W/g V^2 sin(gamma) dgamma/dh,
        % porque despejar gamma implica resolver una ecuacion diferencial y
        % posteriormente otra algebraica de 4 grado, lo cual aunque es posible no
        % veo muy eficiente en principio, y la variacion de gamma es muy muy
        % pequeña de todas formas.
        
        % De modo que se empleará un método iterativo para resolver gamma(h,W)
        
        % Condicion inicial
        V_ini = V(h_ini);
        T_ini = T(h_ini);
        gamma_ini = 0;
        n = 1000; %Numero de intervalos para discretizar la altura
        
        H = linspace(h_ini,h_final,n);
        w(1) = W_ini;
        gamma(1) = gamma_ini;
        
        for k = 1:(n-1),
            contador = 0;
            while contador <= 1000,
                CL =  w(k)*cos(gamma(k))/(0.5*rho(H(k))*(V(H(k)))^2*S);  % L = W cos gamma
                CD =  Cd0 - k2*CL + k1*(CL^2);    %POLAR
                acel = 1;
                cte_a = 1 + (CL/CD)^2*acel^2;
                cte_b = 2*(T(H(k))/w(k))*acel*(CL/CD)^2;
                cte_c = (CL/CD)^2*(T(H(k))/w(k))^2-1;
                seno_gamma = (cte_b - sqrt(cte_b^2-4*cte_a*cte_c))/(2*cte_a);
                gamma_obt = asin(seno_gamma);
                if abs(gamma(k) - gamma_obt) <= 10^-4,
                    contador = 1001;
                end
                contador = contador + 1;
                gamma(k) = gamma_obt;
            end
            % Aqui ya tendriamos CL,CD y gamma para cada h y w
            L_D(k) = CL/CD;
            palanca(k) = delta_T;
            CL2(k) = CL;
            CD2(k) = CD;
            V_v(k) = V(H(k)) * sin(gamma(k));
            velocidad(k) = V(H(k)) * cos(gamma(k));
            tiempo(k) = (H(k+1)-H(k))/V_v(k);
            distancia(k) = (H(k+1)-H(k))/tan(gamma(k));
            empuje(k) = T(H(k));
            fuel(k) = C(H(k))*T(H(k))*tiempo(k);
            w(k+1) = w(k) - fuel(k) * g;      %Calculamos el nuevo peso que tendra la siguiente h
            gamma(k+1) = gamma(k);
        end
        
        tiempo_subida = sum(tiempo);
        S_subida = sum(distancia);
        fuel_subida = sum(fuel);
        energia_subida=0;
        mbaterias_subida=0;
        W = W_ini - fuel_subida*g;
        
        datos(i).segmento.tiempo = cumsum(tiempo);
        datos(i).segmento.fuel = cumsum(fuel); %fuel consumido acumulado
        datos(i).segmento.distancia = cumsum(distancia);
        datos(i).segmento.palanca = palanca;%%% ESCALAR %%%
        datos(i).segmento.empuje = empuje;
%         datos(i).segmento.empuje = CL2;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D; %%% ESCALAR %%%
        datos(i).segmento.peso = w;
        datos(i).segmento.peso(end) = '';
        datos(i).segmento.gamma = gamma;
        datos(i).segmento.gamma(end) = '';
        %datos(i).segmento.techo = techo;
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.altura = H(1:end-1);
        
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        % SUBIDA ACELERADA LINEALMENTE DADOS GAMMA Y LA VELOCIDAD INICIAL Y FINAL (TAS)
        
    case 7
        Vf = velocidad_fin;
        V0 = velocidad_ini;
        gamma = angulo;
        
        if V0 == 0, V0 = 10^-2; end
        
        V = @(h) ((Vf-V0)/(h_final-h_ini)) .* h + (V0*h_final - Vf*h_ini)/(h_final-h_ini);
        V_prima = ((Vf-V0)/(h_final-h_ini));
        M = @(h) V(h)./a(h);
        
        CL = @(h,W) W.*cos(gamma)./(0.5.*rho(h).*(V(h)).^2*S);  % L  = W COS(GAMMA)
        CD = @(h,W) Cd0 - k2.*CL(h,W) + k1.*(CL(h,W).^2);    %POLAR
        D = @(h,W) 0.5.*rho(h).*(V(h).^2)*S.*CD(h,W);     % D
        T = @(h,W) D(h,W) + W.*sin(gamma) + (W./g).*V(h).*sin(gamma).*V_prima; %T = D + W SIN(GAMMA) + W/G * V *SIN(GAMMA)* DV/DH
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        
        switch tipo_motor
            case 1
                delta_T = @(h,W) T(h,W)./(T_SL.*(1+0.2.*M(h).^2).^(1.4./0.4).*(1-0.49.*sqrt(M(h))).*(rho(h)./rho_SL));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                switch derivacion
                    case 1
                        C = @(h,W) c_SL.*correccion(h,W).*(1+1.2.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 2
                        C = @(h,W) c_SL.*correccion(h,W).*(1+0.33.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 3
                        C = @(h,W) c_SL.*correccion(h,W).*(1+0.16875.*M(h)).*sqrt(Temp(h)./Temp_SL);
                end
            case 2
                delta_T  = @(h,W) T(h,W)./(T_SL.*eta_p./V(h).*(1+0.2.*M(h).^2).^(1.4./0.4).*(rho(h).*Temp(h)./(rho_SL*Temp_SL)));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                C = @(h,W) c_SL.*correccion(h,W).*(1+1.44.*M(h)).*sqrt(Temp(h)./Temp_SL).* V(h)./eta_p;
            case 3
                delta_T = @(h,W) T(h,W)./(T_SL.*eta_p./V(h).*((8.55.*rho(h)./rho_SL - 1)./7.55));
                correccion = @(h,W) (a_1.*delta_T(h,W).^4+a_2.*delta_T(h,W).^3+a_3.*delta_T(h,W).^2+a_4.*delta_T(h,W)+a_5);
                C = @(h,W) c_SL .*correccion(h,W).*V(h)./eta_p;
        end
        
        % Sabiendo que dh/dt = V*sin(gamma), integrando y con la condicion inicial
        % de t(h = h_ini) = 0, obtenemos esta expresion para el tiempo = tiempo(h)
        
        C_A = ((Vf-V0)/(h_final-h_ini));
        C_B = (V0*h_final - Vf*h_ini)/(h_final-h_ini);
        
        tiempo = @(h) log((C_A.*h + C_B)./(C_A*h_ini + C_B))./(C_A*sin(gamma));
        
        % Usamos dW/dt = dW/dh * dh/dt = dW/dh * V*sin(gamma) = -C(h)T(h,W)
        
        dwdh = @(h,W) -C(h,W).*g.*T(h,W)./(V(h).*sin(gamma));
        
        % Condicion inicial
        V_ini = V0;
        T_ini = T(h_ini,W_ini);
        delta_T_ini = delta_T(h_ini,W_ini);
        n = 1000; %Numero de intervalos para discretizar la altura
        
        % Resolvemos la ecuacion diferencial para obtener W(h) y por tanto
        % delta_T(h), aunque para obtener las magnitudes globales solo haria falta
        % integrar una vez con h=h_final
        
        vector_h = linspace(h_ini,h_final,n);
        
        [H,Newt] = ode45(dwdh,vector_h,W_ini);
        k = 1;
        while k <= n,
            palanca(k) = delta_T(H(k),Newt(k));
            empuje(k) = T(H(k),Newt(k));
            velocidad(k) = V(H(k))* cos(gamma);
            V_v(k) = V(H(k))*sin(gamma);
            L_D(k) = CL(H(k),Newt(k))./CD(H(k),Newt(k));
            CL2(k) = CL(H(k),Newt(k));
            CD2(k) = CD(H(k),Newt(k));
            % if palanca(k) > 1,
            %     techo = H(k);
            %     msgbox('No es posible alcanzar la altura final definida','Mision','Warn');
            % else
            %     techo = 'No se ha alcanzado el techo';
            % end
            k = k + 1;
        end
        
        tiempo_subida = tiempo(h_final);
        S_subida = (h_final-h_ini)/(tan(gamma));
        fuel_subida = (Newt(1)-Newt(end))/g;
        energia_subida=0;
        mbaterias_subida=0;
        W = Newt(end);
        
        datos(i).segmento.tiempo = (H-H(1))./V_v;
        datos(i).segmento.fuel = (Newt(1) - Newt)./g; %fuel consumido acumulado
        datos(i).segmento.distancia = (H-H(1))./tan(gamma);
        datos(i).segmento.palanca = palanca;
        datos(i).segmento.empuje = empuje;
%         datos(i).segmento.empuje = CL2;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D;
        datos(i).segmento.peso = Newt;
        % datos(i).segmento.techo = techo;
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.gamma = gamma; %%% ESCALAR %%%
        datos(i).segmento.altura = H;
        
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        % SUBIDA STEPPEST CLIMB DADO PALANCA DE GASES Y ALTURA FINAL
        
    case 8
        delta_T = delta_T_gases;
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        
        % Condicion inicial
        V_Stall = sqrt((2*W_ini)/(rho(h_ini)*S*CLmax_limpio));
        V_ini = 1.2*V_Stall;
        V_final = 3 * V_ini;
        V_media = (V_ini + V_final)/2;
        gamma_ini = 0;
        
        n = 500; %Numero de intervalos para discretizar la altura
        num = 1000; %Numero de posibles velocidades para cada altura para hallar gamma max
        H = linspace(h_ini,h_final,n);
        %V = linspace(V_ini,V_final,num);
        
        V(1) = V_ini;
        epsilon = 10^-4;
        V_aux(1) = V(1) + epsilon;
        w(1) = W_ini;
        gamma(1) = gamma_ini;
        gamma_aux(1) = gamma_ini + epsilon;
        indice = 1;
        
        for k = 1:(n-1),    %PARA CADA ALTURA DESDE H_INICIAL A H_FINAL
            memo_1 = V(indice);
            memo_2 = V_aux(indice);         %Con esto inicializamos las variables que estamos optimizando en cada bucle
            memo_3 = gamma(indice);
            memo_4 = gamma_aux(indice);
            clear V;
            clear V_aux;
            clear gamma;
            clear gamma_aux;
            ii = 1;
            giros = 0;
            V(1) = memo_1;
            V_aux(1) = memo_2;
            gamma(1) = memo_3;
            gamma_aux(1) = memo_4;
            
            if k >= 2, V_media = V(1)/10; end;
            %V_media define el ancho de los saltos que pega el algoritmo optimizador
            %Para cada altura, la siguiente velocidad optima se sabe que se encuentra
            %muy cerca del optimo anterior.
            
            % El metodo emplea el siguiente algoritmo: se definen 2 velocidades, V y
            % V_aux, la cual se encuentra a una distancia epsilon de V.
            
            % Se evaluan ambos puntos y se obtiene la tendencia que siguen, de forma
            % que se decide si buscar el maximo aumentando o disminuyendo V.
            
            % Una vez encontrado el primer optimo, el valor inicial para la busqueda
            % del siguiente optimo es el valor hallado previamente, de forma que
            % agiliza enormemente la busqueda del optimo
            
            % Evaluando las 1000 posibles velocidades para 500 alturas, tardaría 9
            % minutos y 10 segundos. Empleando este algoritmo que he inventado, 2.4 segundos.
            
            while ii < num,  %PARA CADA VELOCIDAD
                if ii >=2, gamma(ii) = gamma(ii-1); gamma_aux(ii) = gamma_aux(ii-1); end;
                M(k,ii) = V(ii)/a(H(k));
                M_aux(k,ii) = V_aux(ii)/a(H(k));
                
                switch tipo_motor
                    case 1
                        T(k,ii) = delta_T*(T_SL*(1+0.2*M(k,ii)^2)^(1.4./0.4)*(1-0.49*sqrt(M(k,ii)))*(rho(H(k))/rho_SL));
                        T_aux(k,ii) = delta_T*(T_SL*(1+0.2*M_aux(k,ii)^2)^(1.4./0.4)*(1-0.49*sqrt(M_aux(k,ii)))*(rho(H(k))/rho_SL));
                        switch derivacion
                            case 1
                                C(k,ii) = c_SL*correccion*(1+1.2*M(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                                C_aux(k,ii) = c_SL*correccion*(1+1.2*M_aux(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                            case 2
                                C(k,ii) = c_SL*correccion*(1+0.33*M(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                                C_aux(k,ii) = c_SL*correccion*(1+0.33*M_aux(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                            case 3
                                C(k,ii) = c_SL*correccion*(1+0.16875*M(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                                C_aux(k,ii) = c_SL*correccion*(1+0.16875*M_aux(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                        end
                    case 2
                        T(k,ii) = delta_T*(T_SL.*eta_p/V(ii)*(1+0.2.*M(k,ii)^2).^(1.4./0.4)*(rho(H(k))*Temp(H(k))/(rho_SL*Temp_SL)));
                        C(k,ii) = c_SL*correccion*(1+1.44*M(k,ii))*sqrt(Temp(H(k))/Temp_SL)* V(ii)/eta_p;
                        
                        T_aux(k,ii) = delta_T*(T_SL.*eta_p/V_aux(ii)*(1+0.2.*M_aux(k,ii)^2).^(1.4./0.4)*(rho(H(k))*Temp(H(k))/(rho_SL*Temp_SL)));
                        C_aux(k,ii) = c_SL*correccion*(1+1.44*M_aux(k,ii))*sqrt(Temp(H(k))/Temp_SL)* V_aux(ii)/eta_p;
                    case 3
                        T(k,ii) = delta_T*(T_SL*eta_p/V(ii)*((8.55*rho(H(k))/rho_SL - 1)/7.55));
                        C(k,ii) = c_SL*correccion*V(ii)/eta_p;
                        
                        T_aux(k,ii) = delta_T*(T_SL*eta_p/V_aux(ii)*((8.55*rho(H(k))/rho_SL - 1)/7.55));
                        C_aux(k,ii) = c_SL*correccion*V_aux(ii)/eta_p;
                end
                
                contador = 0;
                while contador <= 1000,
                    CL =  w(k)*cos(gamma(ii))/(0.5*rho(H(k))*(V(ii))^2*S);  % L = W cos gamma
                    CD =  Cd0 - k2*CL + k1*(CL^2);    %POLAR
                    D  =  0.5*rho(H(k))*V(ii)^2*S*CD;
                    seno_gamma = (T(k,ii)-D)/(w(k));
                    gamma_obt = asin(seno_gamma);
                    
                    CL_aux =  w(k)*cos(gamma_aux(ii))/(0.5*rho(H(k))*(V_aux(ii))^2*S);  % L = W cos gamma
                    CD_aux =  Cd0 - k2*CL_aux + k1*(CL_aux^2);    %POLAR
                    D_aux  =  0.5*rho(H(k))*V_aux(ii)^2*S*CD_aux;
                    seno_gamma_aux = (T_aux(k,ii)-D_aux)/(w(k));
                    gamma_obt_aux = asin(seno_gamma_aux);
                    
                    if abs(gamma(ii) - gamma_obt) <= 10^-4,
                        contador = 1001;
                    end
                    contador = contador + 1;
                    gamma(ii) = gamma_obt;
                    gamma_aux(ii) = gamma_obt_aux;
                end
                
                if gamma_aux(ii) > gamma(ii),
                    if ii >=2,
                        if gamma_aux(ii-1) < gamma(ii-1)
                            giros = giros + 1;
                        end
                    end
                    V(ii+1) = V(ii) + V_media/(giros+1);
                    V_aux(ii+1) = V(ii+1) + epsilon;
                else
                    if ii >=2,
                        if gamma_aux(ii-1) > gamma(ii-1)
                            giros = giros + 1;
                        end
                    end
                    V(ii+1) = V(ii) - V_media/(giros+1);
                    V_aux(ii+1) = V(ii+1) + epsilon;
                end
                ii = ii + 1;
                if ii > 2,
                    if abs(gamma(ii-1)-gamma(ii-2)) <= 10^-4, break; end;
                end
                
            end
            
            [gam(k),indice] = max(gamma);   %Nos quedamos para cada altura con el gamma maximo y su velocidad horizontal correspondiente
            velocidad(k) = V(indice)*cos(gam(k));
            empuje(k) = T(k,indice);
            consumo(k) = C(k,indice);
            
            % Aqui ya tendriamos CL,CD, la palanca que hace el gamma max para el h y
            % w dado
            L_D(k) = CL/CD;
            CL2(k) = CL;
            CD2(k) = CD;
            V_v(k) = velocidad(k) * sin(gam(k));
            tiempo(k) = (H(k+1)-H(k))/V_v(k);
            distancia(k) = (H(k+1)-H(k))/tan(gam(k));
            fuel(k) = consumo(k)*empuje(k)*tiempo(k);
            w(k+1) = w(k) - fuel(k) * g;      %Calculamos el nuevo peso que tendra la siguiente h
        end
        
        tiempo_subida = sum(tiempo);
        S_subida = sum(distancia);
        fuel_subida = sum(fuel);
        W = W_ini - fuel_subida*g;
        
        datos(i).segmento.tiempo = cumsum(tiempo);
        datos(i).segmento.fuel = cumsum(fuel);
        datos(i).segmento.distancia = cumsum(distancia);
        datos(i).segmento.palanca = delta_T;
        datos(i).segmento.empuje = empuje;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D;
        datos(i).segmento.peso = w(1:end-1);
        %datos(i).segmento.techo = techo;
        datos(i).segmento.gamma = gam;
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.altura = H(1:end-1);
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        
        
        % SUBIDA FASTEST CLIMB DADO PALANCA DE GASES Y ALTURA FINAL
        
    case 9
        delta_T = delta_T_gases;
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        
        % Condicion inicial
        V_Stall = sqrt((2*W_ini)/(rho(h_ini)*S*CLmax_limpio));
        V_ini = 0.5*V_Stall;
        V_final = 4 * V_ini;
        V_media = (V_ini + V_final)/2;
        gamma_ini = 0;
        
        n = 500; %Numero de intervalos para discretizar la altura
        num = 500; %Numero de posibles velocidades para cada altura para hallar gamma max
        H = linspace(h_ini,h_final,n);
        %V = linspace(V_ini,V_final,num);
        
        V(1) = V_ini;
        epsilon = 10^-4;
        V_aux(1) = V(1) + epsilon;
        w(1) = W_ini;
        gamma(1) = gamma_ini;
        gamma_aux(1) = gamma_ini + epsilon;
        indice = 1;
        
        for k = 1:(n-1),    %PARA CADA ALTURA DESDE H_INICIAL A H_FINAL
            memo_1 = V(indice);
            memo_2 = V_aux(indice);         %Con esto inicializamos las variables que estamos optimizando en cada bucle
            memo_3 = gamma(indice);
            memo_4 = gamma_aux(indice);
            clear V;
            clear V_aux;
            clear gamma;
            clear gamma_aux;
            clear veloc_vert;
            clear veloc_vert_aux;
            ii = 1;
            giros = 1;
            V(1) = memo_1;
            V_aux(1) = memo_2;
            gamma(1) = memo_3;
            gamma_aux(1) = memo_4;
            
            if k >= 2, V_media = V(1)/10; end;
            %V_media define el ancho de los saltos que pega el algoritmo optimizador
            %Para cada altura, la siguiente velocidad optima se sabe que se encuentra
            %muy cerca del optimo anterior.
            
            % El metodo emplea el siguiente algoritmo: se definen 2 velocidades, V y
            % V_aux, la cual se encuentra a una distancia epsilon de V.
            
            % Se evaluan ambos puntos y se obtiene la tendencia que siguen, de forma
            % que se decide si buscar el maximo aumentando o disminuyendo V.
            
            % Una vez encontrado el primer optimo, el valor inicial para la busqueda
            % del siguiente optimo es el valor hallado previamente, de forma que
            % agiliza enormemente la busqueda del optimo
            
            % Evaluando las 1000 posibles velocidades para 500 alturas, tardaría 9
            % minutos y 10 segundos. Empleando este algoritmo que he inventado, 2.4 segundos.
            
            while ii < num,  %PARA CADA VELOCIDAD
                if ii >=2, gamma(ii) = gamma(ii-1); gamma_aux(ii) = gamma_aux(ii-1); end;
                M(k,ii) = V(ii)/a(H(k));
                M_aux(k,ii) = V_aux(ii)/a(H(k));
                
                switch tipo_motor
                    case 1
                        T(k,ii) = delta_T*(T_SL*(1+0.2*M(k,ii)^2)^(1.4./0.4)*(1-0.49*sqrt(M(k,ii)))*(rho(H(k))/rho_SL));
                        T_aux(k,ii) = delta_T*(T_SL*(1+0.2*M_aux(k,ii)^2)^(1.4./0.4)*(1-0.49*sqrt(M_aux(k,ii)))*(rho(H(k))/rho_SL));
                        switch derivacion
                            case 1
                                C(k,ii) = c_SL*(1+1.2*M(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                                C_aux(k,ii) = c_SL*correccion*(1+1.2*M_aux(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                            case 2
                                C(k,ii) = c_SL*(1+0.33*M(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                                C_aux(k,ii) = c_SL*correccion*(1+0.33*M_aux(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                            case 3
                                C(k,ii) = c_SL*(1+0.16875*M(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                                C_aux(k,ii) = c_SL*correccion*(1+0.16875*M_aux(k,ii))*sqrt(Temp(H(k))/Temp_SL);
                        end
                    case 2
                        T(k,ii) = delta_T*(T_SL.*eta_p/V(ii)*(1+0.2.*M(k,ii)^2).^(1.4./0.4)*(rho(H(k))*Temp(H(k))/(rho_SL*Temp_SL)));
                        C(k,ii) = c_SL*correccion*(1+1.44*M(k,ii))*sqrt(Temp(H(k))/Temp_SL)* V(ii)/eta_p;
                        
                        T_aux(k,ii) = delta_T*(T_SL.*eta_p/V_aux(ii)*(1+0.2.*M_aux(k,ii)^2).^(1.4./0.4)*(rho(H(k))*Temp(H(k))/(rho_SL*Temp_SL)));
                        C_aux(k,ii) = c_SL*correccion*(1+1.44*M_aux(k,ii))*sqrt(Temp(H(k))/Temp_SL)* V_aux(ii)/eta_p;
                    case 3
                        T(k,ii) = delta_T*(T_SL*eta_p/V(ii)*((8.55*rho(H(k))/rho_SL - 1)/7.55));
                        C(k,ii) = c_SL *correccion*V(ii)/eta_p;
                        
                        T_aux(k,ii) = delta_T*(T_SL*eta_p/V_aux(ii)*((8.55*rho(H(k))/rho_SL - 1)/7.55));
                        C_aux(k,ii) = c_SL *correccion*V_aux(ii)/eta_p;
                end
                
                contador = 0;
                while contador <= 1000,
                    CL =  w(k)*cos(gamma(ii))/(0.5*rho(H(k))*(V(ii))^2*S);  % L = W cos gamma
                    CD =  Cd0 - k2*CL + k1*(CL^2);    %POLAR
                    D  =  0.5*rho(H(k))*V(ii)^2*S*CD;
                    seno_gamma = (T(k,ii)-D)/(w(k));
                    gamma_obt = asin(seno_gamma);
                    
                    CL_aux =  w(k)*cos(gamma_aux(ii))/(0.5*rho(H(k))*(V_aux(ii))^2*S);  % L = W cos gamma
                    CD_aux =  Cd0 - k2*CL_aux + k1*(CL_aux^2);    %POLAR
                    D_aux  =  0.5*rho(H(k))*V_aux(ii)^2*S*CD_aux;
                    seno_gamma_aux = (T_aux(k,ii)-D_aux)/(w(k));
                    gamma_obt_aux = asin(seno_gamma_aux);
                    
                    if abs(gamma(ii) - gamma_obt) <= 10^-4,
                        contador = 1001;
                    end
                    contador = contador + 1;
                    gamma(ii) = gamma_obt;
                    gamma_aux(ii) = gamma_obt_aux;
                end
                
                veloc_vert(ii) = V(ii)*sin(gamma(ii));
                veloc_vert_aux(ii) = V_aux(ii)*sin(gamma_aux(ii));
                
                if veloc_vert_aux(ii) > veloc_vert(ii),
                    if ii >=2,
                        if veloc_vert_aux(ii-1) < veloc_vert(ii-1)
                            giros = giros + 1;
                        end
                    end
                    V(ii+1) = V(ii) + V_media/(giros+1);
                    V_aux(ii+1) = V(ii+1) + epsilon;
                else
                    if ii >=2,
                        if veloc_vert_aux(ii-1) > veloc_vert(ii-1)
                            giros = giros + 1;
                        end
                    end
                    V(ii+1) = V(ii) - V_media/(giros+1);
                    V_aux(ii+1) = V(ii+1) + epsilon;
                end
                ii = ii + 1;
                if ii > 2,
                    if abs(veloc_vert(ii-1)-veloc_vert(ii-2)) <= 10^-2, break; end;
                end
                
            end
            
            [max_Vv(k),indice] = max(veloc_vert);   %Nos quedamos para cada altura con la Velocidad vertical maxima y su velocidad correspondiente
            gam(k) = gamma(indice);
            velocidad(k) = V(indice)*cos(gam(k));
            empuje(k) = T(k,indice);
            consumo(k) = C(k,indice);
            L_D(k) = CL/CD;
            CL2(k) = CL;
            CD2(k) = CD;
            % Aqui ya tendriamos CL,CD, la palanca que hace el V_vertical max para el h y
            % w dado
            V_v(k) = max_Vv(k);
            tiempo(k) = (H(k+1)-H(k))/V_v(k);
            distancia(k) = (H(k+1)-H(k))/tan(gam(k));
            fuel(k) = consumo(k)*empuje(k)*tiempo(k);
            w(k+1) = w(k) - fuel(k) * g;      %Calculamos el nuevo peso que tendra la siguiente h
        end
        
        tiempo_subida = sum(tiempo);
        S_subida = sum(distancia);
        fuel_subida = sum(fuel);
        W = W_ini - fuel_subida*g;
        
        datos(i).segmento.tiempo = cumsum(tiempo);
        datos(i).segmento.fuel = cumsum(fuel);
        datos(i).segmento.distancia = cumsum(distancia);
        datos(i).segmento.palanca = delta_T;
        datos(i).segmento.empuje = empuje;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D;
        datos(i).segmento.peso = w(1:end-1);
        %datos(i).segmento.techo = techo;
        datos(i).segmento.gamma = gam;
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.altura = H(1:end-1);
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        
    case 10
        % SUBIDA DADOS PALANCA DE GASES = CONSTANTE, VELOCIDAD INICIAL Y FINAL
        Vf = velocidad_fin;
        V0 = velocidad_ini;
        delta_T = delta_T_gases;
        
        a_1=3.559957437510763; a_2=-10.739698199171459; a_3= 11.989635150373475; a_4=-5.869876557884609; a_5=2.059994459180667;
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        
        V = @(h) ((Vf-V0)/(h_final-h_ini)) .* h + (V0*h_final - Vf*h_ini)/(h_final-h_ini);
        V_prima = ((Vf-V0)/(h_final-h_ini));
        M = @(h) V(h)./a(h);
        
        switch tipo_motor
            case 1
                T = @(h) delta_T.*(T_SL*(1+0.2*M(h)^2)^(1.4./0.4)*(1-0.49*sqrt(M(h))).*(rho(h)./rho_SL));
                switch derivacion
                    case 1
                        C = @(h) c_SL.*correccion.*(1+1.2.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 2
                        C = @(h) c_SL.*correccion.*(1+0.33.*M(h)).*sqrt(Temp(h)./Temp_SL);
                    case 3
                        C = @(h) c_SL.*correccion.*(1+0.16875.*M(h)).*sqrt(Temp(h)./Temp_SL);
                end
            case 2
                T = @(h) delta_T.*(T_SL.*eta_p./V(h).*(1+0.2.*M(h)^2).^(1.4./0.4).*(rho(h).*Temp(h)./(rho_SL*Temp_SL)));
                C = @(h) c_SL.*correccion.*(1+1.44.*M(h)).*sqrt(Temp(h)./Temp_SL).* V(h)./eta_p;
            case 3
                T = @(h) delta_T.*(T_SL.*eta_p./V(h).*((8.55.*rho(h)./rho_SL - 1)./7.55));
                C = @(h) c_SL .*correccion.*V(h)./eta_p;
        end
        
        % PROBLEMA: No veo forma de implementar las ecuaciones de la mecanica del
        % vuelo directamente: L = W cos(gamma) + W/g V^2 sin(gamma) dgamma/dh,
        % porque despejar gamma implica resolver una ecuacion diferencial y
        % posteriormente otra algebraica de 4 grado, lo cual aunque es posible no
        % veo muy eficiente en principio, y la variacion de gamma es muy muy
        % pequeña de todas formas.
        
        % De modo que se empleará un método iterativo para resolver gamma(h,W)
        
        % Condicion inicial
        V_ini = V(h_ini);
        T_ini = T(h_ini);
        gamma_ini = 0;
        n = 1000; %Numero de intervalos para discretizar la altura
        
        H = linspace(h_ini,h_final,n);
        w(1) = W_ini;
        gamma(1) = gamma_ini;
        
        for k = 1:(n-1),
            contador = 0;
            while contador <= 1000,
                CL =  w(k)*cos(gamma(k))/(0.5*rho(H(k))*(V(H(k)))^2*S);  % L = W cos gamma
                CD =  Cd0 - k2*CL + k1*(CL^2);    %POLAR
                acel = 1 + (V(H(k))/g)*V_prima;
                cte_a = 1 + (CL/CD)^2*acel^2;
                cte_b = 2*(T(H(k))/w(k))*acel*(CL/CD)^2;
                cte_c = (CL/CD)^2*(T(H(k))/w(k))^2-1;
                seno_gamma = real((cte_b - sqrt(cte_b^2-4*cte_a*cte_c))/(2*cte_a));
                gamma_obt = asin(seno_gamma);
                if abs(gamma(k) - gamma_obt) <= 10^-4,
                    contador = 1001;
                end
                contador = contador + 1;
                gamma(k) = gamma_obt;
            end
            % Aqui ya tendriamos CL,CD y gamma para cada h y w
            L_D(k) = CL/CD;
            palanca(k) = delta_T;
            CL2(k) = CL;
            CD2(k) = CD;
            V_v(k) = V(H(k)) * sin(gamma(k));
            velocidad(k) = V(H(k)) * cos(gamma(k));
            tiempo(k) = (H(k+1)-H(k))/V_v(k);
            distancia(k) = (H(k+1)-H(k))/tan(gamma(k));
            empuje(k) = T(H(k));
            fuel(k) = C(H(k))*T(H(k))*tiempo(k);
            w(k+1) = w(k) - fuel(k) * g;      %Calculamos el nuevo peso que tendra la siguiente h
            gamma(k+1) = gamma(k);
        end
        
        tiempo_subida = sum(tiempo);
        distancia_subida = sum(distancia);
        S_subida = distancia_subida;
        fuel_subida = sum(fuel);
        W = W_ini - fuel_subida*g;
        
        datos(i).segmento.tiempo = cumsum(tiempo);
        datos(i).segmento.fuel = cumsum(fuel);
        datos(i).segmento.distancia = cumsum(distancia);
        datos(i).segmento.palanca = palanca;
        datos(i).segmento.empuje = empuje;
%         datos(i).segmento.empuje = CL2;
        datos(i).segmento.veloc_vertical = V_v;
        datos(i).segmento.CL = CL2;
        datos(i).segmento.CD = CD2;
        datos(i).segmento.L_D = L_D;
        datos(i).segmento.peso = w;
        datos(i).segmento.peso(end) = '';
        %datos(i).segmento.techo = techo;
        datos(i).segmento.gamma = gamma;
        datos(i).segmento.gamma(end) = '';
        datos(i).segmento.velocidad = velocidad;
        datos(i).segmento.altura = H(1:end-1);
        
        
        
end
% datos(i).nombre = 'Subida';
% datos(i).lista_variables = [{''};'Tiempo';'Fuel';'Peso';'Distancia';'Angulo de subida';'Velocidad';'Velocidad vertical';'Palanca';'CL';...
%     'L/D';'Altura'];

global flag_palanca
if max(datos(i).segmento.palanca) > 1.2
    flag_palanca = 1;
else
    flag_palanca = 0;
end

end
