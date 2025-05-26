function actuaciones_diagrama_carga_de_pago_autonomia
%Funci�n que representa el diagrama de carga de pago-autonom�a

%Valores operativos
h_CR      = 3000;   %[m] Altitud de crucero
h_f       = 100;    %[m] Altitud a la que se acaba la misi�n y se abre el paraca�das

%Datos de entrada del UAV:
S          = 2.2;    %[m^2] Superficie alar
PL_max     = 12;     %[kg] Carga de pago m�xima
PL_min     = 4.075+0.4;   %[kg] Carga de pago m�nima (imprescindible)
mF_trapped = 0.69; %[kg] Combustible inutlizable
ZFM        = 57.831 + mF_trapped;  %[kg] Zero fuel mass
MTOM       = 75;%[kg] Masa m�xima al despegue
m_TO       = [ceil((ZFM+PL_min)/5)*5:5:MTOM];
mF_max     = 17.95*.75/1.05 - mF_trapped; %[kg] Masa m�xima de combustible utilizable = capacidad del tanque expresada en kg - combustible inutilizable

%Coef aerodin�micos despegando desde catapulta y sin misiles
C_L_max_CR = 1.637378*2.05/2.2;  %[-] Coef de sustentaci�n m�ximo en configuraci�n de crucero
C_D0_CR    = 0.031429185864671; %[-] Coef de resistencia en configuraci�n de crucero
C_D1_CR    = -0.047465554112273; %[-] Coef de resistencia en configuraci�n de crucero
C_D2_CR    = 0.069695334554813;  %[-] Coef de resistencia en configuraci�n de crucero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagrama carga de pago - autonom�a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Datos de misi�n
fprintf('******************************************************************\n')
disp('Diagrama carga de pago - autonom�a, sin tren y con 0 misiles')
diagrama_carga_de_pago(PL_max,PL_min,ZFM,m_TO,mF_max,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,h_CR,h_f);
fprintf('******************************************************************\n')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function diagrama_carga_de_pago(PL_max,PL_min,ZFM,m_TO,mF_max,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,h_CR,h_f)
% 
% figure, grid on, hold on
% xlabel('Autonom�a [h]')
% ylabel('Carga de pago [kg]')
% axis([0 13 PL_min PL_max+2])
% h = text(10,10,['Masa al despegue ' num2str(m_TO(end)) ' kg']);
% set(h,'rotation',-40);
% 
% misops_fzero  = optimset('TolX',1e-7); %Opciones para fzero
% misops_fsolve = optimset('Display','off','TolX',1e-7,'TolFun',1e-7,'MaxIter',10,'DiffMinChange',1e-8,'FinDiffType','central'); %Opciones para fsolve
% 
% %El punto A se corresponde al caso de no fuel
% PL_A  = PL_max;
% End_A = 0;
% 
% %El punto B se corresponde al caso de alcanzar la masa m�xima al despegue
% %(debe coincidir con las 8 horas de vuelo)
% PL_B  = PL_max;
% End_B = fzero(@(End) ZFM + PL_B - combustible_analisis(m_TO(end),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[1 15],misops_fzero);
% 
% plot([End_A End_B],[PL_A PL_B],'linewidth',1)
% 
% %El punto C se corresponde con el corte entre la curva de masa m�xima al
% %despegue y m�xima capacidad del tanque (el combustible consumido debe
% %coincidir con la capacidad m�xima del tanque)
% flag = -1;
% fval = 1;
% x_0 = [PL_min 0.9*mF_max];
% while flag <= 0 || norm(fval) > 1e-3 %Si la bandera es mnor que 0, entonces es que fsolve ha tenido problemas. Reanudamos la resuloci�n partiendo del �ltimo punto propuesto
%     %Tambi�n continuamos repitiendo el proceso si la norma de la
%     %funci�n objetivo es superior a 1e-3 (aproximadamente 10 gramos)
%     [x_C,fval,flag] = fsolve(@(x) [(ZFM + x(1) - combustible_analisis(m_TO(end),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,x(2)*3600,h_CR,h_f))/ZFM
%         (ZFM + x(1) - combustible_analisis(ZFM + x(1) + mF_max,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,x(2)*3600,h_CR,h_f))/ZFM],...
%         x_0,misops_fsolve); %La primera estimaci�n est� basada en que se consume aproximadamente 0.9 kg por hora
%     x_0 = x_C;
% end
% PL_C  = x_C(1);
% End_C = x_C(2);
% 
% if PL_C > PL_min
%     PL = linspace(PL_B,PL_C,5);
%     End(1) = End_B; End(5) = End_C;
%     for j = 2:length(PL)-1
%         End(j) = fzero(@(End) ZFM + PL(j) - combustible_analisis(m_TO(end),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[End_B End_C],misops_fzero);
%     end
% else
%     PL = linspace(PL_B,PL_min,5);
%     End(1) = End_B;
%     for j = 2:length(PL)
%         End(j) = fzero(@(End) ZFM + PL(j) - combustible_analisis(m_TO(end),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[End_B End_C],misops_fzero);
%     end
% end
% plot(End,PL,'linewidth',1)
% 
% for i = (length(m_TO)-1):-1:1
%     m_TO(i)
%     %Primero determinamos cu�l ser�a la carga de pago que se corresponder�a
%     %con autonom�a 0
%     PL_End0 = combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,0,h_CR,h_f) - ZFM;
%     
%     %Si la carga de pago que se corresponde con autonom�a 0 es mayor que
%     %PL_max, entonces el punto B tiene el mismo significado que antes
%     if PL_End0 >= PL_max
%         PL_B  = PL_max;
%         End_B = fzero(@(End) ZFM + PL_B - combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[0 8],misops_fzero);
%     else
%         %Si no, ahora el punto B es el punto en el que se tiene autonom�a 0
%         PL_B = PL_End0;
%         End_B = 0;        
%     end
%     
%     %El punto C se corresponde con el corte entre la curva de masa m�xima al
%     %despegue y m�xima capacidad del tanque
%     flag = -1;
%     fval = 1;
%     x_0 = x_C; %Usamos como punto de arranque el punto C obtenido en la masa anterior
%     while flag <= 0 || norm(fval) > 1e-3 %Si la bandera es menor que 0, entonces es que fsolve ha tenido problemas. Reanudamos la resuloci�n partiendo del �ltimo punto propuesto
%         %Tambi�n continuamos repitiendo el proceso si la norma de la
%         %funci�n objetivo es superior a 1e-3 (aproximadamente 10 gramos)
%         [x_C,fval,flag] = fsolve(@(x) [(ZFM + x(1) - combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,x(2)*3600,h_CR,h_f))/ZFM
%             (ZFM + x(1) - combustible_analisis(ZFM + x(1) + mF_max,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,x(2)*3600,h_CR,h_f))/ZFM],...
%             x_0,misops_fsolve); 
%         x_0 = x_C;
%     end
%     PL_C  = x_C(1);
%     End_C = x_C(2);
%     
%     %Si PL_C es negativo, entonces hay que cortar en el punto en el que PL = PL_min
%     if PL_C < 0
%         PL_C = PL_min;
%         End_C = fzero(@(End) ZFM + PL_C - combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[End_B End_C],misops_fzero);
%     end
%     
%     %Representamos enre PL_B y PL_C
%     PL = linspace(PL_B,PL_C,5);
%     End(1) = End_B; End(5) = End_C;
%     for j = 2:length(PL)-1
%         End(j) = fzero(@(End) ZFM + PL(j) - combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[End_B End_C],misops_fzero);
%     end
%     plot(End,PL,'linewidth',1)
% end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function m_end = combustible_analisis(m0,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,tf,h_CR,h_f)
% rho_SL = densityISA(0);   %[kg/m^3] Densidad a nivel del mar
% g      = 9.80665;         %[m/s^2] Gravedad
% misops = odeset('RelTol',1e-9,'AbsTol',1e-9); %Tolerancias para el integrador
% misops_opt = optimset('Display','off','TolX',1e-9,'TolFun',1e-18); %Tolerancias para el optimizador. Ponemos una 
% %TolFun exageradamente peque�a para que el criterio de parada sea siempre TolX
% 
% %Comparamos las velocidades m�nima y de m�xima autonom�a a h_CR
% rho = densityISA(h_CR);
% V_stall = sqrt(m0*g./(0.5*rho*C_L_max_CR*S));
% V_min = 1.3*V_stall;
% %Utilizamos como primera estimaci�n a la velocidad de m�xima autonom�a la 
% %que se obtiene suponiendo que el consumo espec�fico y el rendimiento 
% %propulsivo no dependen de la velocidad de vuelo:
% V_E0 = sqrt(m0*g/.5/rho/S)*sqrt(-(C_D1_CR/C_D0_CR/6) + sqrt((C_D1_CR/C_D0_CR/6)^2+C_D2_CR/C_D0_CR/3)); 
% V_E = fmincon(@(V) 1e4*consumo_crucero(V,h_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,m0,g),V_E0,[],[],[],[],V_stall,2*V_stall,[],misops_opt);
% V_EAS = max(V_min,V_E)*sqrt(rho/rho_SL);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Subida
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Ascendemos hasta la altitud de crucero a velocidad EAS contante
% h_A = 0;
% h_B = h_CR;
% m = m0;
% 
% [h,y] = ode45(@(h,y) subida(h,y,V_EAS,S,C_D0_CR,C_D1_CR,C_D2_CR,g,rho_SL),[h_A h_B],m,misops);
% 
% mF_subida = y(1) - y(end);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Crucero
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h   = h_CR;
% rho = densityISA(h);
% m   = m - mF_subida;
% 
% V = V_EAS *sqrt(rho_SL/rho);
% if tf == 0;
%     t = 0;
%     y = m;
% else
%     [t,y] = ode45(@(t,y) crucero(t,y,V,S,C_D0_CR,C_D1_CR,C_D2_CR,g,rho,h),[0 tf],m,misops);
% end
% mF_crucero = y(1,1) - y(end,1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Descenso
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Descendemos hasta 100 m
% h_A = h_CR;
% h_B = h_f;
% m = m - mF_crucero;
% 
% [h,y] = ode45(@(h,y) descenso(h,y,V_EAS,S,C_D0_CR,C_D1_CR,C_D2_CR,g,rho_SL),[h_A h_B],m,misops);
% 
% mF_descenso = y(1) - y(end);
% 
% m_end = m - mF_descenso

  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = subida(h,y,V_EAS,S,C_D0,C_D1,C_D2,g,rho_SL)
% m = y(1);
% 
% [rho,drho_dh] = densityISA(h); 
% 
% V     = V_EAS*sqrt(rho_SL/rho);
% [T,c] = propulsive_model(V,h,1);
% 
% gamma = (T-.5*rho*V^2*S*(C_D0+C_D1*(m*g/(.5*rho*S*V^2))+C_D2*(m*g/(.5*rho*S*V^2))^2))/(m*g-.5*m*V^2/rho*drho_dh);
% 
% f = -c/(V*gamma);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = crucero(t,y,V,S,C_D0,C_D1,C_D2,g,rho,h)
% m = y(1);
% 
% persistent delta_p
% if isempty(delta_p)
%     delta_p = 0.5;
% end
% delta_p = fzero(@(delta_p) .5*rho*V^2*S*(C_D0+C_D1*(m*g/(.5*rho*S*V^2))+C_D2*(m*g/(.5*rho*S*V^2))^2) - propulsive_model(V,h,delta_p),delta_p);
% 
% [T,c] = propulsive_model(V,h,delta_p);
% 
% f = -c;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = descenso(h,y,V_EAS,S,C_D0,C_D1,C_D2,g,rho_SL)
% m = y(1);
% 
% [rho,drho_dh] = densityISA(h); 
% 
% V     = V_EAS*sqrt(rho_SL/rho);
% T = 0;
% c  = .01/60; %[kg/s] %Consumo con el motor en ralent�
% 
% gamma = (T-.5*rho*V^2*S*(C_D0+C_D1*(m*g/(.5*rho*S*V^2))+C_D2*(m*g/(.5*rho*S*V^2))^2))/(m*g-.5*m*V^2/rho*drho_dh);
% 
% f = -c/(V*gamma);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = consumo_crucero(V,h,S,C_D0,C_D1,C_D2,m,g)
% %Funci�n que devuelve el consumo en crucero (es necesario determinar la 
% %posici�n de palanca)
% rho = densityISA(h);
% 
% persistent delta_p
% if isempty(delta_p)
%     delta_p = 0.5;
% end
% 
% delta_p = fzero(@(delta_p) .5*rho*V^2*S*(C_D0+C_D1*(m*g/(.5*rho*S*V^2))+C_D2*(m*g/(.5*rho*S*V^2))^2) - propulsive_model(V,h,delta_p),delta_p);
% [T,c] = propulsive_model(V,h,delta_p);
% f = c;