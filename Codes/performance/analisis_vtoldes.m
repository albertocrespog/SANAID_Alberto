 function [fuel_vtoldes,mbaterias_vtoldes,energia_vtoldes,tiempo_vtoldes,S_des,datos] =  analisis_vtoldes(propul,aerodinamica,vtol_des,W_inicial,h_inicial,opcion,datos,i,Prop_data,data_electric)

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

%       vtoldes(1) Altura inicial
%       vtoldes(2) %Altura final
%       vtoldes(3) %Tiempo de espera en hovering
%       vtoldes(4) %Palanca de gases

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

%Datos VTOL climb

h_inicial = vtol_des(1);
h_final = vtol_des(2);
vdes = vtol_des(3);

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
        % Descenso entre dos alturas a velocidad constante
        
        n=100; %numero de intervalos
        vector_h=linspace(h_ini,h_final,n);
        rhovect=rho(vector_h);
       
        if tipo_motor==4 %Para el modelo eléctrico no es necesario resolver ecuaciones diferenciales puesto que tendremos 
    %un valor de V, CL, CD,y T para cada incremento de altura. 
    

    
    CD =  0.5;    %POLAR 
    D = 0.5.*rhovect.*(vdes.^2)*S.*CD;     % D 
    T = W_ini -D;
    for m=1:n
        J=@(npsi) vdes./(Dprop*npsi);

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
VVERTICAL = ones(n,1)*vdes;
tiempo=ones(n,1);
tiempo(1)=0;
ENERGIA(1)=E_ini;
for k=1:n-1
   tiempo(k+1) = (-vector_h(k+1)+vector_h(k))/((VVERTICAL(k+1)+VVERTICAL(k))/2)+tiempo(k);
   ENERGIA(k+1) = ENERGIA(k)+ (P_sol(k+1)+P_sol(k))*(tiempo(k+1)-tiempo(k))/((eta_mp_sol(k+1)+eta_mp_sol(k))*eta_m) ;  
end
    
end
tiempo_vtoldes=tiempo(end);
fuel_vtoldes = 0;
energia_vtoldes=ENERGIA(end);
mbaterias_vtoldes=energia_vtoldes/e0;
S_des=0; %NO RECORRE DISTANCIA HORIZONTAL!!


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





datos(i).segmento.fuel = fuel_vtoldes; %fuel consumido acumulado
datos(i).segmento.mbaterias = ENERGIA/e0;
datos(i).segmento.energia = energia_vtoldes;
datos(i).segmento.tiempo = tiempo;
datos(i).segmento.palanca = delta_T;
datos(i).segmento.empuje = T;
datos(i).segmento.peso = W_ini;
datos(i).segmento.veloc_vertical = -VVERTICAL;
datos(i).segmento.altura = vector_h;
datos(i).segmento.distancia=0; %No recorre distancia horizontal!!
datos(i).segmento.CD = CD;
datos(i).segmento.energiav = ENERGIA;


end


