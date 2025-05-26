function [fuel_vtolclimb,mbaterias_vtolclimb,energia_vtolclimb,tiempo_vtolclimb,S_climb,datos] =  analisis_vtolclimb(propul,aerodinamica,vtol_climb,W_inicial,h_inicial,opcion,datos,i,handles,Prop_data,data_electric)

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


% VTOL Climb

% 1: Altura inicial
% 2: Altura final
% 3: T hovering
% 4: Palanca
% 5: Masa de baterías bruta

%       vtolclimb(1) Altura inicial
%       vtolclimb(2) %Altura final
%       vtolclimb(3) %Tiempo de espera en hovering
%       vtolclimb(4) %Palanca de gases
%       vtolclimb(5) %Masa de baterías bruta

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

% PESOS

empty_weight = handles.pesos(1);
payload =   handles.pesos(2);

%Datos VTOL climb

% h_inicial = vtol_climb(1);
h_final = vtol_climb(2);
t_espera = vtol_climb(3);
delta_T = vtol_climb(4);
vclimb = vtol_climb(5);
mbat_nominal = vtol_climb(6); %Masa EN BRUTO de baterías


% ATMOSFERA ESTANDAR
[Temp_SL,rho_SL,p_SL,a_SL] = atmos_inter_mio(0);

% ----------  DEFINICION DE LAS FUNCIONES RESPECTO DE H --------------
g = 9.80665;
R = 287.058;
gamma_atm = 1.4;
lambda = 6.5*10^-3;

W_ini = W_inicial;
h_ini = h_inicial;
E_ini=0;
Vv_ini=0;

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
        % Subida entre dos alturas a palanca constante
        
        n=100; %numero de intervalos
        vector_h=linspace(h_ini,h_final,n);
        rhovect=rho(vector_h);
        
        if tipo_motor==4 %Para el modelo eléctrico no es necesario resolver ecuaciones diferenciales puesto que tendremos
            %un valor de V, CL, CD,y T para cada incremento de altura.
            
            %Una vez calculada T, se resuelve con un fzero las revoluciones
            %necesarias para dar ese empuje, se calculan palanca, rendimiento y
            %potencia, y se calcula dEdh=P(h)/(V(h)*sin(gamma))
            
            CD =  0.5;    %POLAR
            D = @(Vv,h) 0.5.*rho(h).*(Vv^2)*S.*CD;     % D
            nps= delta_T*nps_max;
            
            J=@(Vv) Vv/(Dprop*nps);
            CT=@(Vv) 0;
            
            for j=1:length(CT_Polyfit)
                CT_int=@(Vv) CT_Polyfit(j)*J(Vv).^(length(CT_Polyfit)-j);
                CT=@(Vv) CT_int(Vv)+CT(Vv);%CT_Polyfit(1)+CT_Polyfit(2)*J(nps)+CT_Polyfit(3)*J(nps)^2;
            end
            
            T = @(Vv,h) num_motores*CT(Vv).*rho(h).*nps^2*Dprop^4;
        end
        
        %%% INTEGRAMOS
        if tipo_motor==4
            
            dVdh= @(Vv,h) g*(T(Vv,h)-D(Vv,h)-W_ini)/(W_ini*Vv);
            [H,VVERTICAL]= ode45(dVdh,vector_h,Vv_ini);
            
        end
        
        J_sol=VVERTICAL./(nps*Dprop);
        
        CP= ones(length(n));
        
        for j=1:length(CP_Polyfit)
            CP_int= CP_Polyfit(j).*J_sol.^(length(CP_Polyfit)-j);
            CP= CP_int+CP;
        end
        
        eta_mp_sol= 0;
        for j=1:length(etamp_Polyfit)
            eta_mp_int = etamp_Polyfit(j).*J_sol.^(length(etamp_Polyfit)-j);
            eta_mp_sol= eta_mp_int+eta_mp_sol;
        end
        
        P_sol=(num_motores.*CP'.*rho(H)'.*(nps^3).*(Dprop)^5)'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tiempo = ones(length(n));
        ENERGIA = ones(length(n));
        tiempo(1)=0;
        ENERGIA(1)=E_ini;
        for k=1:n-1
            tiempo(k+1) = (H(k+1)-H(k))/((VVERTICAL(k+1)+VVERTICAL(k))/2)+tiempo(k);
            ENERGIA(k+1) = ENERGIA(k)+ (P_sol(k+1)+P_sol(k))*(tiempo(k+1)-tiempo(k))/((eta_mp_sol(k+1)+eta_mp_sol(k))*eta_m) ;
        end
        
        palanca=ones(length(n))*delta_T;
        T=num_motores*CT(VVERTICAL).*rho(H).*nps^2*Dprop^4;
%         ENERGIAtot=sum(ENERGIA);
        ENERGIAtot=ENERGIA(end);
        mbaterias=ENERGIAtot/e0;
        S_climb=0; %NO RECORRE DISTANCIA HORIZONTAL!!
        
        
        %k=1;
        % while k <= n,
        %
        %     if tipo_motor==4
        %         palanca(k)=delta_T*ones();
        %         empuje(k)= T(k);
        %         V_v(k)= Vvect(k)*sin(gamma);
        %         velocidad(k)=Vvect(k)*cos(gamma);
        %         L_D(k)= CL(k)./CD(k);
        %         CL2(k)= CL(k);
        %
        %     else
        % palanca(k) = delta_T(H(k),Newt(k));
        % empuje(k) = T(H(k),Newt(k));
        % V_v(k) = V(H(k))*sin(gamma);
        % velocidad(k) = V(H(k)) * cos(gamma);
        % L_D(k) = CL(H(k),Newt(k))./CD(H(k),Newt(k));
        % CL2(k) = CL(H(k),Newt(k));
        % if palanca(k) > 1,
        % %     techo = H(k);
        % %     msgbox('No es posible alcanzar la altura final definida','Mision','Warn');
        % else
        %     techo = 'No se ha alcanzado el techo';
        % end
        
        
        
        
        if tipo_motor==4
            fuel_vtolclimb = 0;
            energia_vtolclimb=ENERGIA(end);
            mbaterias_vtolclimb = mbaterias;
%             tiempo_vtolclimb=sum(tiempo);
            tiempo_vtolclimb = tiempo(end);
            
        else
            % fuel_subida = (Newt(1)-Newt(end))/g;
            %  energia_subida=0;
            % mbaterias_subida=0;
            % datos(i).subida.fuel = (Newt(1) - Newt)./g; %fuel consumido acumulado
            % datos(i).subida.mbaterias=0;
            
        end
        datos(i).segmento.fuel = fuel_vtolclimb; %fuel consumido acumulado
        datos(i).segmento.mbaterias = ENERGIA/e0;
        datos(i).segmento.tiempo = tiempo;
        datos(i).segmento.palanca = palanca;
        datos(i).segmento.empuje = T;
        datos(i).segmento.peso = W_ini;
        datos(i).segmento.veloc_vertical = VVERTICAL;
        datos(i).segmento.altura = H;
        datos(i).segmento.distancia=0; %No recorre distancia horizontal!!
        datos(i).segmento.CD = CD;
        datos(i).segmento.energiav = ENERGIA;
        
        
        
        
    case 2
        % Subida entre dos alturas a velocidad constante
        
        n=100; %numero de intervalos
        vector_h=linspace(h_ini,h_final,n);
        rhovect=rho(vector_h);
        
        if tipo_motor==4 %Para el modelo eléctrico no es necesario resolver ecuaciones diferenciales puesto que tendremos
            %un valor de V, CL, CD,y T para cada incremento de altura.
            
            
            
            CD = 2.97;    %POLAR
            D = 0.5.*rhovect.*(vclimb.^2)*S.*CD;     % D
            T = W_ini + D;
            for m=1:n
                J=@(npsi) vclimb./(Dprop*npsi);
                
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
                
                eta_mp=@(npsi) 0;
                for j=1:length(etamp_Polyfit)
                    eta_mp_int=@(npsi) etamp_Polyfit(j).*J(npsi).^(length(etamp_Polyfit)-j);
                    eta_mp=@(npsi) eta_mp_int(npsi)+eta_mp(npsi);
                end
                
                
                nps_sol(m,1)=  fzero(@(npsi) T(m)-num_motores.*CT(npsi).*rhovect(m).*(npsi^2)*Dprop^4,nps_max);
                
            end
        end
        
        
        delta_T = nps_sol/nps_max;
        P_sol= num_motores*CP(nps_sol).*rhovect'.*(nps_sol.^3).*(Dprop)^5;
        eta_mp_sol= eta_mp(nps_sol);
        
        %%% INTEGRAMOS
        VVERTICAL = ones(n,1)*vclimb;
        tiempo=ones(n,1);
        ENERGIA=ones(n,1);
        tiempo(1)=0;
        ENERGIA(1)=E_ini;
        for k=1:n-1
            tiempo(k+1) = (vector_h(k+1)-vector_h(k))/((VVERTICAL(k+1)+VVERTICAL(k))/2)+tiempo(k);
            ENERGIA(k+1) = ENERGIA(k)+ (P_sol(k+1)+P_sol(k))*(tiempo(k+1)-tiempo(k))/((eta_mp_sol(k+1)+eta_mp_sol(k))*eta_m) ;
        end
        tiempo_vtolclimb=tiempo(end);
        fuel_vtolclimb = 0;
        
        fprintf('La energía gastada es ' )
        energia_vtolclimb=ENERGIA(end);
        mbaterias_vtolclimb=energia_vtolclimb/e0;
        S_climb=0; %NO RECORRE DISTANCIA HORIZONTAL!!
        
        
        %k=1;
        % while k <= n,
        %
        %     if tipo_motor==4
        %         palanca(k)=delta_T*ones();
        %         empuje(k)= T(k);
        %         V_v(k)= Vvect(k)*sin(gamma);
        %         velocidad(k)=Vvect(k)*cos(gamma);
        %         L_D(k)= CL(k)./CD(k);
        %         CL2(k)= CL(k);
        %
        %     else
        % palanca(k) = delta_T(H(k),Newt(k));
        % empuje(k) = T(H(k),Newt(k));
        % V_v(k) = V(H(k))*sin(gamma);
        % velocidad(k) = V(H(k)) * cos(gamma);
        % L_D(k) = CL(H(k),Newt(k))./CD(H(k),Newt(k));
        % CL2(k) = CL(H(k),Newt(k));
        % if palanca(k) > 1,
        % %     techo = H(k);
        % %     msgbox('No es posible alcanzar la altura final definida','Mision','Warn');
        % else
        %     techo = 'No se ha alcanzado el techo';
        % end
        
        
        
        
        
        datos(i).segmento.fuel = fuel_vtolclimb; %fuel consumido acumulado
        datos(i).segmento.mbaterias = ENERGIA/e0;
        datos(i).segmento.energia = energia_vtolclimb;
        datos(i).segmento.energiav = ENERGIA;
        datos(i).segmento.tiempo = tiempo;
        datos(i).segmento.palanca = delta_T;
        datos(i).segmento.empuje = T;
        datos(i).segmento.peso = W_ini;
        datos(i).segmento.veloc_vertical = VVERTICAL;
        datos(i).segmento.altura = vector_h;
        datos(i).segmento.distancia=0; %No recorre distancia horizontal!!
        
    case 3
        
        % Ascenso entre dos alturas a velocidad constante
        
        n=100; %numero de intervalos
        vector_h=linspace(h_ini,h_final,n);
        rhovect=rho(vector_h);
       
        if tipo_motor==4 %Para el modelo eléctrico no es necesario resolver ecuaciones diferenciales puesto que tendremos 
    %un valor de V, CL, CD,y T para cada incremento de altura. 
    

    
    CD =  0.5;    %POLAR 
    D = 0.5.*rhovect.*(vclimb.^2)*S.*CD;     % D 
    T = W_ini +D;
    for m=1:n
        J=@(npsi) vclimb./(Dprop*npsi);

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
        
        eta_mp=@(npsi) 0;
        for j=1:length(etamp_Polyfit)
            eta_mp_int=@(npsi) etamp_Polyfit(j).*J(npsi).^(length(etamp_Polyfit)-j);
            eta_mp=@(npsi) eta_mp_int(npsi)+eta_mp(npsi);
        end
        
       
        nps_sol(m,1)=  fzero(@(npsi) T(m)-num_motores.*CT(npsi).*rhovect(m).*(npsi^2)*Dprop^4,nps_max);
        
    end
        end
        
        
        delta_T=nps_sol/nps_max;
        P_sol= num_motores*CP(nps_sol).*rhovect'.*(nps_sol.^3).*(Dprop)^5;
        eta_mp_sol= eta_mp(nps_sol);
        
%%% INTEGRAMOS 
VVERTICAL = ones(n,1)*vclimb;
tiempo=ones(n,1);
tiempo(1)=0;
ENERGIA(1)=E_ini;
for k=1:n-1
   tiempo(k+1) = (vector_h(k+1)-vector_h(k))/((VVERTICAL(k+1)+VVERTICAL(k))/2)+tiempo(k);
   ENERGIA(k+1) = ENERGIA(k)+ (P_sol(k+1)+P_sol(k))*(tiempo(k+1)-tiempo(k))/((eta_mp_sol(k+1)+eta_mp_sol(k))*eta_m);  
end
    

tiempo_vtolclimb=tiempo(end);
fuel_vtolclimb = 0;
energia_vtolclimb=ENERGIA(end);
mbaterias_vtolclimb=energia_vtolclimb/e0;
S_climb=0; %NO RECORRE DISTANCIA HORIZONTAL!!


%k=1;
% while k <= n,
%     
%     if tipo_motor==4
%         palanca(k)=delta_T*ones();
%         empuje(k)= T(k);
%         V_v(k)= Vvect(k)*sin(gamma);
%         velocidad(k)=Vvect(k)*cos(gamma);
%         L_D(k)= CL(k)./CD(k);
%         CL2(k)= CL(k);
%         
%     else
% palanca(k) = delta_T(H(k),Newt(k));
% empuje(k) = T(H(k),Newt(k));
% V_v(k) = V(H(k))*sin(gamma);
% velocidad(k) = V(H(k)) * cos(gamma);
% L_D(k) = CL(H(k),Newt(k))./CD(H(k),Newt(k));
% CL2(k) = CL(H(k),Newt(k));
% if palanca(k) > 1,
% %     techo = H(k);
% %     msgbox('No es posible alcanzar la altura final definida','Mision','Warn');
% else
%     techo = 'No se ha alcanzado el techo';
% end





datos(i).segmento.fuel = fuel_vtolclimb; %fuel consumido acumulado
datos(i).segmento.mbaterias = ENERGIA/e0;
datos(i).segmento.energia = energia_vtolclimb;
datos(i).segmento.tiempo = tiempo;
datos(i).segmento.palanca = delta_T;
datos(i).segmento.empuje = T;
datos(i).segmento.peso = W_ini;
datos(i).segmento.veloc_vertical = VVERTICAL;
datos(i).segmento.altura = vector_h;
datos(i).segmento.distancia=0; %No recorre distancia horizontal!!
datos(i).segmento.CD = CD;
datos(i).segmento.energiav = ENERGIA;


        
    case 4 % Hovering dada masa de baterías
        
        n=100; %numero de intervalos
        rho_ini=rho(h_ini);
        
        %Carga de pago
        %Masa de baterías
        
        %Masa de baterías actualizada = %Masa de bat estandar + carga de pago estandar - carga de pago usada
        %LOS VALORES ESTANDAR SON LOS COMPLICADOS DE SACAR, NECESITAN DEL
        %ESTUDIO DE SENSIBILIDAD
        
        
        
        if tipo_motor==4 %Para el modelo eléctrico no es necesario resolver ecuaciones diferenciales
            nps_max = 150000/(Dprop/0.0254)/60;
            %Una vez calculada T, se resuelve con un fzero las revoluciones
            %necesarias para dar ese empuje, se calculan palanca, rendimiento y
            %potencia, y se calcula dEdh=P(h)/(V(h)*sin(gamma))
            
            
            mbat_util = mbat_nominal*(1-tau);
            energia_util = mbat_util*e0;
            
            W_ini = (mbat_nominal + payload + empty_weight)*g;
            
            %%%%%%%%%%%%% J=0; %%%%%%%%%%
            
            CT = CT_Polyfit(end);
            CP = CP_Polyfit(end);
            eta_mp = etamp_Polyfit(end);
            
            T = @(nps)  num_motores*CT.*rho_ini.*nps^2*Dprop^4;
            nps_sol=  fzero(@(nps) T(nps)-W_ini,nps_max);
            
            FM =0.3369; %FIGURA DE MÉRITO, CALCULADA CON CT Y CP
            P_sol = (num_motores.*CP.*rho_ini'.*(nps_sol^3).*(Dprop)^5)/(FM); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        fprintf ('Los valores que quieres son: ')
        tiempo_vtolclimb = energia_util/(P_sol/(eta_mp*eta_m));
        delta_T = nps_sol/nps_max*ones(n,1);
        fprintf ('Delta vale:')
        delta_T(1)
        fprintf ('La potencia necesaria es: ')
        format long
        P_sol
        format short
        fuel_vtolclimb = 0;
        mbaterias_vtolclimb = mbat_util;
        energia_vtolclimb = energia_util;
        S_climb= 0;
        
        fprintf ('Aquí acaban los valores. ')
        
        datos(i).segmento.fuel = fuel_vtolclimb; %fuel consumido acumulado
        datos(i).segmento.mbaterias = mbaterias_vtolclimb;
        datos(i).segmento.tiempo = tiempo_vtolclimb;
        datos(i).segmento.palanca = delta_T;
        datos(i).segmento.empuje = T(nps_sol);
        datos(i).segmento.veloc_vertical = 0;
        datos(i).segmento.altura = h_ini;
        datos(i).segmento.distancia=0; %No recorre distancia horizontal!!
        
        
end




end



