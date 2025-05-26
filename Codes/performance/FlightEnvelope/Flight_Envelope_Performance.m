function Flight_Envelope_Performance

%% Rate of Climb in Service Ceiling - Excess of vertical velocity
ROC_SC    = 100*0.3048/60; %[m/s] 

%% aricraft Data:
S        = 2.2;  %[m^2] Superficie alar
MTOM     = 75;   %[kg] Masa máxima al despegue

%% Aerodynamic model 
C_L_max_CR = 1.637378*2.05/2.2;  %[-] Coef de sustentación máximo en configuración de crucero
C_D0_CR    = 0.031429185864671; %[-] Coef de resistencia en configuración de crucero
C_D1_CR    = -0.047465554112273; %[-] Coef de resistencia en configuración de crucero
C_D2_CR    = 0.069695334554813;  %[-] Coef de resistencia en configuración de crucero

SF = 1.3; % Safety Factor
h_1 = 3000; % altitude to check flight envelope (m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Envolvente de vuelo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('******************************************************************\n')
disp('Flight Envelope')
disp('Cruise Configuration')
flight_envelope(MTOM,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,ROC_SC,SF,h_1)
fprintf('******************************************************************\n')
fprintf('\n')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function flight_envelope(m,C_L_max,S,C_D0,C_D1,C_D2,ROC_SC,SF,h_1)
% 
% rho_SL     = densityISA(0);   %[kg/m^3] Densidad a nivel del mar
% g          = 9.80665;         %[m/s^2] Gravedad
% misops     = optimset('Display','off','TolX',1e-9,'TolFun',1e-12); %Opciones para fsolve y fmincon
% warning('off','all'); %Para deshabilitir el warning que devuelve fmincon
% 
% figure, grid on, hold on
% xlabel('V [m/s]')
% ylabel('h [m]')
% 
% axis([15 55 0 9000])
% 
% %% Define corners of the Fliught Envelope Diagram
% 
% % Stall speed @ sea level
% V_stall_SL = sqrt(m*g./(0.5*rho_SL*C_L_max*S));
% V_A = V_stall_SL;
% h_A = 0;
% x_A = [V_A rho_SL];
% 
% % Cross between stall speed and absolute ceiling 
% x_B = fsolve(@(x) [x(1) - sqrt(m*g/(0.5*x(2)*C_L_max*S))
%      propulsive_model(x(1),densityISA2ALT(x(2)),1) - .5*x(2)*x(1)^2*S*(C_D0+C_D1*(m*g/(0.5*x(2)*x(1)^2*S))+C_D2*(m*g/(0.5*x(2)*x(1)^2*S))^2)],....
%     [V_stall_SL,rho_SL],misops);
% V_B = x_B(1);
% h_B = densityISA2ALT(x_B(2));
% 
% % Max speed @ sea level without climb capabilities. 
% % 2 possibilities need to be tested: if max speed has been achieved or not
% V_max_SL = fzero(@(V) propulsive_model(V,0,1) - .5*rho_SL*V^2*S*(C_D0+C_D1*(m*g/(0.5*rho_SL*V^2*S))+C_D2*(m*g/(0.5*rho_SL*V^2*S))^2),...
%     V_stall_SL*2);
% V_C = V_max_SL;
% h_C = 0;
% x_C = [V_C rho_SL];
% 
% % Cross between stall speed and service ceiling (ROC)
% x_D = fsolve(@(x) [x(1) - sqrt(m*g/(0.5*x(2)*C_L_max*S))
%     propulsive_model(x(1),densityISA2ALT(x(2)),1) - (0.5*x(2)*x(1)^2*S*(C_D0+C_D1*(m*g/(0.5*x(2)*x(1)^2*S))+C_D2*(m*g/(0.5*x(2)*x(1)^2*S))^2)+m*g*ROC_SC/x(1))],...
%     [V_stall_SL,rho_SL],misops);
% V_D = x_D(1);
% h_D = densityISA2ALT(x_D(2));
% 
% % Max speed @ sea level with climb capabilities
% V_max_SL_SC = fzero(@(V) propulsive_model(V,0,1) - (0.5*rho_SL*V^2*S*(C_D0+C_D1*(m*g/(0.5*rho_SL*V^2*S))+C_D2*(m*g/(0.5*rho_SL*V^2*S))^2)+m*g*ROC_SC/V),...
%     V_stall_SL*2);
% V_E = V_max_SL_SC;
% h_E = 0;
% 
% % Min speed @ sea level 
% V_min_SL = 1.3*V_stall_SL;
% V_F = V_min_SL;
% h_F = 0;
% 
% % Cross between minimum speed and absolute ceiling 
% x_G = fsolve(@(x) [x(1) - SF*sqrt(m*g/(0.5*x(2)*C_L_max*S))
%     propulsive_model(x(1),densityISA2ALT(x(2)),1) - (.5*x(2)*x(1)^2*S*(C_D0+C_D1*(m*g/(0.5*x(2)*x(1)^2*S))+C_D2*(m*g/(0.5*x(2)*x(1)^2*S))^2))],....
%     [V_stall_SL,rho_SL],misops);
% V_G = x_G(1);
% h_G = densityISA2ALT(x_G(2));
% 
% %% Plotting %%
% 
% %% Stall Speed 
% h = linspace(h_A,h_B,50);
% rho = densityISA(h);
% V_stall = sqrt(m*g./(0.5*rho*C_L_max*S));
% plot(V_stall,h,'linewidth',1)
% 
% 
% %% Minimum Velocit associated to Safety factor
% h = linspace(h_F,h_G,25);
% rho = densityISA(h);
% V_min = SF*sqrt(m*g./(0.5*rho*C_L_max*S));
% plot(V_min,h,'--','linewidth',1)
% 
% %% Service Ceiling - Associated to ROC
% clear rho
% V = linspace(V_D,V_E,100);
% for j = 1:length(V)
%     if j == 1,
%         x0 = x_D(2);
%     elseif j == length(V)
%         x0 = rho_SL;
%     else
%         x0 = rho(j-1);
%     end
%     rho(j) = fzero(@(rho)  propulsive_model(V(j),densityISA2ALT(rho),1) -...
%        (.5*rho*V(j)^2*S*(C_D0+C_D1*(m*g/(0.5*rho*V(j)^2*S))+C_D2*(m*g/(0.5*rho*V(j)^2*S))^2)+m*g*ROC_SC/V(j)),x0);
% end
% h_SC = densityISA2ALT(min(rho,rho_SL));
% plot(V,h_SC,':','linewidth',1)
% 
% %Techo absoluto
% V = linspace(V_B,V_C,100);
% clear rho
% for j = 1:length(V)
%     if j == 1,
%         x0 = x_B(2);
%     elseif j == length(V)
%         x0 = rho_SL;
%     else
%         x0 = rho(j-1);
%     end
%     rho(j) = fzero(@(rho) propulsive_model(V(j),densityISA2ALT(rho),1) -...
%         (.5*rho*V(j)^2*S*(C_D0+C_D1*(m*g/(0.5*rho*V(j)^2*S))+C_D2*(m*g/(0.5*rho*V(j)^2*S))^2)),x0);
% end
% h_AC = densityISA2ALT(min(rho,rho_SL));
% plot(V,h_AC,'linewidth',1)
% 
% %Buscamos los valores de los techos
% %Techo absoluto
% x_J = fmincon(@(x) x(2),x_B,[],[],[],[],[x_A(1) 0],[x_C(1) rho_SL],@(x) nonlcon(x,S,C_D0,C_D1,C_D2,m,g,0),misops);
% V_J = x_J(1);
% h_J = densityISA2ALT(x_J(2));
% %Techo de servicio
% x_K = fmincon(@(x) x(2),x_B,[],[],[],[],[x_A(1) 0],[x_C(1) rho_SL],@(x) nonlcon(x,S,C_D0,C_D1,C_D2,m,g,ROC_SC),misops);
% V_K = x_K(1);
% h_K = densityISA2ALT(x_K(2));
% 
% 
% %Velocidad de máxima autonomía
% %Representamos la velocidad de máxima autonomía hasta que crucemos con el
% %techo absoluto
% h = linspace(0,h_J,25);
% for i = 1:length(h)
%     if i == 1
%         x0 = 1.3*V_stall_SL;
%     else
%         x0 = V_max_endurance(i-1);
%     end
%     try
%         V_max_endurance(i) = fmincon(@(V) 1e4*consumo_crucero(V,h(i),S,C_D0,C_D1,C_D2,m,g),x0,[],[],[],[],x_A(1),x_C(1),[],misops);
%     catch
%         h(i:end)=[];
%         break
%     end
% end
% plot(V_max_endurance,h,':','linewidth',1)
% legend('V_{Stall}','V_{min}=1.3V_{Stall}','Service Ceiling','Absolute Ceiling','V_{Max-Endurance}')
% 
% %Representamos la temperatura a la cual no es seguro operar la electrónica
% Temp_min = -20; %[ºC] Temperatura mínima marcada por la electrónica del motor
% %Según el modelo ISA, se corresponde con una altitud de
% h_Temp_min = (15 - Temp_min)/(6.5e-3);
% rho = densityISA(h_Temp_min);
% clear V
% V(1) = sqrt(m*g./(0.5*rho*C_L_max*S));
% V(2) = fzero(@(V) propulsive_model(V,h_Temp_min,1) -...
%         (.5*rho*V^2*S*(C_D0+C_D1*(m*g/(0.5*rho*V^2*S))+C_D2*(m*g/(0.5*rho*V^2*S))^2)),V_C);
% % plot(V,h_Temp_min*ones(1,2),'--','linewidth',1)
% 
% %% Prints on Screen
% disp(['Stall Speed @ 0 m: ' num2str(V_A) ' m/s'])
% disp(['Minimum Speed @ 0 m: ' num2str(V_F) ' m/s'])
% disp(['Maximum Endurance Speed @ 0 m: ' num2str(V_max_endurance(1)) ' m/s'])
% disp(['Maximum Speed @ 0 m (with climb capabilities): ' num2str(V_E) ' m/s'])
% disp(['Maximum Speed @ 0 m (without climb capabilities): ' num2str(V_C) ' m/s'])
% disp(['Service Ceiling: ' num2str(h_K) ' m (' num2str(V_K) ' m/s)'])
% disp(['Absolute Ceiling: ' num2str(h_J) ' m (' num2str(V_J) ' m/s)'])
% fprintf('------------------------------------------------------------------\n')
% 
% %% Study for a given altitude
% rho = densityISA(h_1);
% V_stall = sqrt(m*g./(0.5*rho*C_L_max*S));
% V_min = 1.3*V_stall;
% V_max_endurance = fmincon(@(V) 1e4*consumo_crucero(V,h_1,S,C_D0,C_D1,C_D2,m,g),V_max_endurance(end),[],[],[],[],x_A(1),x_C(1),[],misops);
% V_max_SC = fzero(@(V) propulsive_model(V,h_1,1) -...
%         (.5*rho*V^2*S*(C_D0+C_D1*(m*g/(0.5*rho*V^2*S))+C_D2*(m*g/(0.5*rho*V^2*S))^2)+m*g*ROC_SC/V),1.5*V_min);
% V_max = fzero(@(V) propulsive_model(V,h_1,1) -...
%         (.5*rho*V^2*S*(C_D0+C_D1*(m*g/(0.5*rho*V^2*S))+C_D2*(m*g/(0.5*rho*V^2*S))^2)),1.5*V_min);
% fprintf('** Flight Envelope at a given altitude**\n')
% disp(['Stall Speed @ ' num2str(h_1) ' m: ' num2str(V_stall) ' m/s'])
% disp(['Minimum Speed @ ' num2str(h_1) ' m: ' num2str(V_min) ' m/s'])
% disp(['Maximum Endurance Speed @ ' num2str(h_1) ' m: ' num2str(V_max_endurance) ' m/s'])
% disp(['Maximum Speed @ ' num2str(h_1) ' m: (with climb capabilities): ' num2str(V_max_SC) ' m/s'])
% disp(['Maximum Speed @ ' num2str(h_1) ' m: (without climb capabilities): ' num2str(V_max) ' m/s'])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [c,ceq] = nonlcon(x,S,C_D0,C_D1,C_D2,m,g,ROC_SC)
% c = [];
% ceq = propulsive_model(x(1),densityISA2ALT(x(2)),1) - (0.5*x(2)*x(1)^2*S*(C_D0+C_D1*(m*g/(0.5*x(2)*x(1)^2*S))+C_D2*(m*g/(0.5*x(2)*x(1)^2*S))^2)+m*g*ROC_SC/x(1));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = consumo_crucero(V,h,S,C_D0,C_D1,C_D2,m,g)
% %Función que devuelve el consumo en crucero (es necesario determinar la 
% %posición de palanca)
% rho = densityISA(h);
% 
% [T,c,Vnmax,delta_p_max] = propulsive_model(V,h,1);
% 
% delta_p = fzero(@(delta_p) .5*rho*V^2*S*(C_D0+C_D1*(m*g/(.5*rho*S*V^2))+C_D2*(m*g/(.5*rho*S*V^2))^2) - propulsive_model(V,h,delta_p),[0 delta_p_max]);
% [T,c] = propulsive_model(V,h,delta_p);
% f = c;
