function [seg, handles] = Generation_Mission_Segments(conv_UNITS,num_missions,type_missions,cruise_mode,climb_mode,turn_mode,Performance,case_AC)

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
ft2m = conv_UNITS.ft2m;
m2ft = conv_UNITS.m2ft;

% Asks user for the input data        
% Enter type of mission segments betwee brackets being
% - 1 Taxy
% - 2 TakeOff
% - 3 = Climb
% - 4 = VTOL Climb
% - 5 = Cruise
% - 6 = Load Deployment
% - 7 = Turn
% - 8 = Descent
% - 9 = Descent (VTOL)
% - 10 = Alternative Airport climb to 3000ft
% - 11 = Turn loitter 45 min
% - 12 = Landing
% - 13 = Dummy to complete the 3 segment requirement for AP

for i=1:num_missions
    type_mission = type_missions(i);
    switch type_mission
        case 1 % case 1 Taxy
            mission_tex = 'Taxy';
            temp_local_taxy = 15; % (Celsius)
            h_inicial_taxy = 0; % (m)
            P_local_taxy = 1.013268093075000e+05; % (Pa)
            delta_T_taxy = 0.05; % 4: PALANCA DE RALENTI EN TAXI = 0.05
            V_taxy = 10; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            t_tazy = 5.60; % 6: TIEMPO DE ESPERA EN TAXI (s)
            
            seg(i).nombre = mission_tex;
            seg(i).datos.mision = type_mission; 
            seg(i).datos.temp_local = temp_local_taxy;% (Celsius)
            seg(i).datos.h_inicial = h_inicial_taxy; % (m)
            seg(i).datos.P_local = P_local_taxy;% (Pa)
            seg(i).datos.delta_T = delta_T_taxy;% (%)
            seg(i).datos.V_taxy = V_taxy;% (%)
            seg(i).datos.t_tazy = t_tazy;% (%)
            
            taxi(1) = temp_local_taxy + 273.15; % 1: TEMPERATURA LOCAL (K)
            taxi(2) = h_inicial_taxy; % 2: ALTURA LOCAL (m)
            taxi(3) = P_local_taxy; % 3: PRESION LOCAL (Pa)
            taxi(4) = delta_T_taxy; % 4: PALANCA DE RALENTI EN TAXI = 0.05
            taxi(5) = V_taxy; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            taxi(6) = t_tazy; % 6: TIEMPO DE ESPERA EN TAXI (s)
            
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
            seg(i).datos.Data_ATM = Data_ATM;
            seg(i).datos.Performance = Performance;

        case 2 % case 2 TakeOff
            mission_tex = 'TakeOff';
            temp_local_TO = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_TO = 0; % 2: ALTURA LOCAL (m)
            P_local_TO = 1.013268093075000e+05; % 3: PRESION LOCAL (Pa)
            mu_TO = 0.02; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)  
            h_obstacle_TO = 10; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            gamma_climb_TO = 2*D2R; % 6: GAMMA DE SUBIDA MINIMO
            delta_T_TO = 1; % 7: PALANCA DE GASES PARA DESPEGUE
            
            seg(i).nombre = mission_tex;
            seg(i).datos.mision = type_mission; 
            seg(i).datos.temp_local = temp_local_TO;% (Celsius)
            seg(i).datos.h_inicial = h_inicial_TO; % (m)
            seg(i).datos.P_local = P_local_TO;% (Pa)
            seg(i).datos.mu_takeoff = mu_TO;% (-)
            seg(i).datos.h_obstacle = h_obstacle_TO;% (m)
            seg(i).datos.gamma_climb = gamma_climb_TO;% (rads)
            seg(i).datos.delta_T = delta_T_TO;% (%)
            
            despegue(1) = temp_local_TO + 273.15; % 1: TEMPERATURA LOCAL (K)
            despegue(2) = h_inicial_TO; % 2: ALTURA LOCAL (m)
            despegue(3) = P_local_TO; % 3: PRESION LOCAL (Pa)
            despegue(4) = mu_TO; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            despegue(5) = h_obstacle_TO; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            despegue(6) = gamma_climb_TO; % 6: GAMMA DE SUBIDA MINIMO
            despegue(7) = delta_T_TO  ; % 7: PALANCA DE GASES PARA DESPEGUE
        
            handles.despegue(1) = despegue(1); % 1: TEMPERATURA LOCAL (K)
            handles.despegue(2) = despegue(2); % 2: ALTURA LOCAL (m)
            handles.despegue(3) = despegue(3); % 3: PRESION LOCAL (Pa)
            handles.despegue(4) = despegue(4); % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            handles.despegue(5) = despegue(5); % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            handles.despegue(6) = despegue(6); % 6: GAMMA DE SUBIDA MINIMO
            handles.despegue(7) = despegue(7); % 7: PALANCA DE GASES PARA DESPEGUE
            
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;
            
        case 3 % case 3 Climb
            mission_tex = 'Climb';
            seg(i).datos.mision = type_mission; 
            seg(i).nombre = mission_tex;
            h_inicial_cl = 0; % 1: ALTURA INICIAL - [m]
            h_final_cl = 3000; % 1: ALTURA FINAL - [m]
            
            subida.h_inicial = h_inicial_cl; % 1: ALTURA INICIAL - [m]
            subida.h_final = h_final_cl; % 1: ALTURA FINAL - [m]
            handles.subida.h_inicial = subida.h_inicial; % 1: ALTURA INICIAL - [m]
            handles.subida.h_final = subida.h_final; % 1: ALTURA FINAL - [m]
            
            seg(i).datos.h_inicial = h_inicial_cl;% (%)
            seg(i).datos.h_final = h_final_cl;% (%)
            
            switch climb_mode
                case 2 % 'Subida dados M y gamma';
                    Mach_cl = 0.14;                   
                    gamma_cl = 3*D2R;                   
                    subida.Mach = Mach_cl; % 3: MACH DE VUELO - [-]
                    subida.gamma = gamma_cl; % 2: GAMMA DE SUBIDA
                    handles.subida.Mach = subida.Mach; % 3: MACH DE VUELO - [-]
                    handles.subida.gamma = subida.gamma; % 2: GAMMA DE SUBIDA
                case 3 % 'Subida dados EAS y gamma';
%                     EAS_cl =
%                     gamma_cl = 3*D2R;                   
%                     subida.EAS = ; % 5: VELOCIDAD EAS  - [m/s]
%                     subida.gamma = ; % 2: GAMMA DE SUBIDA  - [-]
%                     handles.subida.EAS = str2double(get(handles.edit12,'String')); % 5: VELOCIDAD EAS  - [m/s]
%                     handles.subida.gamma = str2double(get(handles.edit13,'String')); % 2: GAMMA DE SUBIDA  - [-]
                case 4 % 'Subida dados TAS y gamma';
%                     subida.TAS = ; % 4: VELOCIDAD TAS  - [m/s]
%                     subida.gamma = ; % 2: GAMMA DE SUBIDA  - [-]
%                     handles.subida.TAS = str2double(get(handles.edit12,'String')); % 4: VELOCIDAD TAS  - [m/s]
%                     handles.subida.gamma = str2double(get(handles.edit13,'String')); % 2: GAMMA DE SUBIDA  - [-]
                case 5 % 'Subida dados M y palanca';
%                     subida.Mach = ; % 3: MACH DE VUELO  - [-]
%                     subida.palanca = ; % 6: PALANCA DE GASES  - [-]
%                     handles.subida.Mach = str2double(get(handles.edit12,'String')); % 3: MACH DE VUELO  - [-]
%                     handles.subida.palanca = str2double(get(handles.edit13,'String')); % 6: PALANCA DE GASES  - [-]
                case 6 % 'Subida dados EAS y palanca';
%                     subida.EAS = ; % 5: VELOCIDAD EAS  - [m/s]
%                     subida.palanca = ; % 6: PALANCA DE GASES  - [-]
%                     handles.subida.EAS = str2double(get(handles.edit12,'String')); % 5: VELOCIDAD EAS  - [m/s]
%                     handles.subida.palanca = str2double(get(handles.edit13,'String')); % 6: PALANCA DE GASES  - [-]
                case 7 % 'Subida dados TAS y palanca';
%                     subida.TAS = ; % 4: VELOCIDAD TAS  - [m/s]
%                     subida.palanca = ; % 6: PALANCA DE GASES  - [-]
%                     handles.subida.TAS = str2double(get(handles.edit12,'String')); % 4: VELOCIDAD TAS  - [m/s]
%                     handles.subida.palanca = str2double(get(handles.edit13,'String')); % 6: PALANCA DE GASES  - [-]
                case 8 % 'Subida dados V inicial,final y gamma';
%                     subida.V_ini = ; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
%                     subida.V_fin = ; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
%                     subida.gamma = ; % 2: GAMMA DE SUBIDA  - [-]
%                     handles.subida.V_ini = str2double(get(handles.edit12,'String')); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
%                     handles.subida.V_fin = str2double(get(handles.edit13,'String')); % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
%                     handles.subida.gamma = str2double(get(handles.edit14,'String')); % 2: GAMMA DE SUBIDA  - [-]
                case 9 % 'Subida steppest climb';
%                     subida.palanca = ; % 6: PALANCA DE GASES  - [-]
%                     handles.subida.palanca = str2double(get(handles.edit12,'String')); % 6: PALANCA DE GASES  - [-]
                case 10 % 'Subida fastest climb';
%                     subida.palanca = ; % 6: PALANCA DE GASES  - [-]
%                     handles.subida.palanca = str2double(get(handles.edit12,'String')); % 6: PALANCA DE GASES  - [-]
                case 11 % 'Subida dados V inicial,final y palanca'
%                     subida.V_ini = ; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA) - [m/s]
%                     subida.V_fin = ;% 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA) - [m/s]
%                     subida.palanca = ; % 6: PALANCA DE GASES - [-]
%                     handles.subida.V_ini = str2double(get(handles.edit12,'String')); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA) - [m/s]
%                     handles.subida.V_fin = str2double(get(handles.edit13,'String'));% 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA) - [m/s]
%                     handles.subida.palanca = str2double(get(handles.edit14,'String')); % 6: PALANCA DE GASES - [-]
            end
                               
            seg(i).datos.h_inicial = h_inicial_cl;% (%)
            seg(i).datos.Mach = Mach_cl;
            seg(i).opcion = climb_mode - 1;
            
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;

        case 4 % case 4 VTOL Climb
            mission_tex = 'VTOL_Climb';
            seg(i).datos.mision = type_mission; 
        case 5 % case 5 Cruise
            mission_tex = 'Cruise';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 
            switch cruise_mode
                case 2 % 'Crucero dado M y distancia'
                    h_inicial_cr = 3000;% - [m]
                    dist_final_cr = 40*1000;% - [m]
                    V_cr = 30;
                    [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_cr);
                    seg(i).datos.Data_ATM = Data_ATM;
                    seg(i).datos.Performance = Performance;
                    Mach_cr = V_cr/Data_ATM.a;% - [-]
                    crucero.h_inicial = h_inicial_cr;% - [m]
                    crucero.dist_final = dist_final_cr;% - [m]
                    crucero.Mach = Mach_cr;% - [-]
                    handles.crucero.h_inicial = crucero.h_inicial;% - [m]
                    handles.crucero.dist_final = dist_final_cr;% - [m]
                    handles.crucero.Mach = crucero.Mach;% - [-]
                    
                    seg(i).datos.h_inicial = h_inicial_cr;% (%)
                    seg(i).datos.dist_final = dist_final_cr;% (%)
                    seg(i).datos.Mach = Mach_cr;
                    seg(i).opcion = cruise_mode - 1;

                case 3 % ;'Crucero dado CL y distancia'
%                     crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     crucero.dist_final = str2double(get(handles.edit11,'String'));% - [m]
%                     crucero.CL = str2double(get(handles.edit12,'String'));% - [-]
%                     handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     handles.crucero.dist_final = str2double(get(handles.edit11,'String'));% - [m]
%                     handles.crucero.CL = str2double(get(handles.edit12,'String'));% - [-]
                case 4 % ;'Crucero dados V inicial,final y palanca'
%                     crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     crucero.dist_final = str2double(get(handles.edit11,'String'));% - [m]
%                     crucero.V_ini = str2double(get(handles.edit12,'String'));% - [m/s]
%                     crucero.V_fin = str2double(get(handles.edit13,'String'));% - [m/s]
%                     crucero.palanca = str2double(get(handles.edit14,'String'));% - []
%                     handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     handles.crucero.dist_final = str2double(get(handles.edit11,'String'));% - [m]
%                     handles.crucero.V_ini = str2double(get(handles.edit12,'String'));% - [m/s]
%                     handles.crucero.V_fin = str2double(get(handles.edit13,'String'));% - [m/s]
%                     handles.crucero.palanca = str2double(get(handles.edit14,'String'));% - []
                case 5 % 'Crucero con polar en funcion de M';
%                     crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     crucero.dist_final = str2double(get(handles.edit11,'String'));% - [m]
%                     crucero.Mach = str2double(get(handles.edit12,'String'));% - []
%                     crucero.Cd0 = str2double(get(handles.edit13,'String'));% - []
%                     crucero.k1 = str2double(get(handles.edit14,'String'));% - []
%                     crucero.k2 = str2double(get(handles.edit15,'String'));% - []
%                     handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     handles.crucero.dist_final = str2double(get(handles.edit11,'String'));% - [m]
%                     handles.crucero.Mach = str2double(get(handles.edit12,'String'));% - []
%                     handles.crucero.Cd0 = str2double(get(handles.edit13,'String'));% - []
%                     handles.crucero.k1 = str2double(get(handles.edit14,'String'));% - []
%                     handles.crucero.k2 = str2double(get(handles.edit15,'String'));% - []
                case 6 % 'Crucero de max alcance dado peso final'
%                     crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     crucero.fuel = str2double(get(handles.edit11,'String'));% - [kg]
%                     handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     handles.crucero.fuel = str2double(get(handles.edit11,'String'));% - [kg]
                case 7 % 'Crucero de max autonomia dado peso final'
%                     crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     crucero.fuel = str2double(get(handles.edit11,'String'));% - [kg]
%                     handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));% - [m]
%                     handles.crucero.fuel = str2double(get(handles.edit11,'String'));% - [kg]
            end
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_cr);
            seg(i).datos.Data_ATM = Data_ATM;
            seg(i).datos.Performance = Performance;
 
        case 6 % case 6 Load Deployment
            mission_tex = 'Load_Deployment';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 
        case 7 % case 7 - Turn
            mission_tex = 'Turn';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 
            
            h_inicial_turn = 3000;% - [m]
            t_final_turn = 45*60;% - [seg]
            V_turn = 30;
            Mach_turn = V_turn/340;% - [-]
            
            
            seg(i).datos.h_inicial = h_inicial_turn;% (%)
            seg(i).datos.t_final = t_final_turn;% (%)
            seg(i).datos.Mach = Mach_cr;
            seg(i).opcion = turn_mode - 1;

            viraje.h_inicial = h_inicial_turn;
            viraje.tiempo_final = t_final_turn;
            handles.viraje.h_inicial = viraje.h_inicial;
            handles.viraje.tiempo_final = viraje.tiempo_final;
            
            switch turn_mode
                case 2 % 'Viraje horizontal dado V y palanca'
%                     viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     viraje.palanca = str2double(get(handles.edit13,'String'));% - []
%                     handles.viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     handles.viraje.palanca = str2double(get(handles.edit13,'String'));% - []
                case 3 % 'Viraje horizontal dado V y CL'
%                     viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     viraje.CL = str2double(get(handles.edit13,'String'));% - []
%                     handles.viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     handles.viraje.CL = str2double(get(handles.edit13,'String'));% - []
                case 4 % 'Viraje horizontal dado V y balance'
%                     viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     viraje.balance = str2double(get(handles.edit13,'String'));% - []
%                     handles.viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     handles.viraje.balance = str2double(get(handles.edit13,'String'));% - []
                case 5 % 'Viraje horizontal dado V y n'
                    V_turn = V_turn;
                    n_turn = 1.25;
                    viraje.velocidad = V_turn;% - [m/s]
                    viraje.n = n_turn;% - []
                    handles.viraje.velocidad = viraje.velocidad;% - [m/s]
                    handles.viraje.n = viraje.n;% - []
                    seg(i).datos.n = n_turn;
                case 6 % 'V.H dado V y radio de giro '
%                     viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     viraje.radio = str2double(get(handles.edit13,'String'));% - [m]
%                     handles.viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     handles.viraje.radio = str2double(get(handles.edit13,'String'));% - [m]
                case 7 % 'V.H dado V y velocidad de guiñada';...
%                     viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     viraje.vel_guiniada = str2double(get(handles.edit13,'String'));% - [rad/s]
%                     handles.viraje.velocidad = str2double(get(handles.edit12,'String'));% - [m/s]
%                     handles.viraje.vel_guiniada = str2double(get(handles.edit13,'String'));% - [rad/s]
                case 8 % 'V.H dado palanca y a factor de carga max'
%                     viraje.palanca = str2double(get(handles.edit12,'String'));% - []
%                     handles.viraje.palanca = str2double(get(handles.edit12,'String'));% - []
                case 9 % 'V.H dado palanca y a v de guiñada max'
%                     viraje.palanca = str2double(get(handles.edit12,'String'));% - []
%                     handles.viraje.palanca = str2double(get(handles.edit12,'String'));% - []
                case 10 % 'V.H dado palanca y a radio min'
%                     viraje.palanca = str2double(get(handles.edit12,'String'));% - []
%                     handles.viraje.palanca = str2double(get(handles.edit12,'String'));% - []
            end
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_turn,V_turn);
            seg(i).datos.Data_ATM = Data_ATM;
            seg(i).datos.Performance = Performance;
                    
        case 8 % case 8 Descent
            mission_tex = 'Descent';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 

        case 9 % case 9 Descent (VTOL)
            mission_tex = 'Descent_VTOL';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission;
            
        case 10 % case 10 Climb Waiting Area to 3000ft
            mission_tex = 'Climb';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission;

            h_inicial_cl = 0; % 1: ALTURA INICIAL - [m]
            h_final_cl = 3000*ft2m; % 1: ALTURA FINAL - [m]
            
            subida.h_inicial = h_inicial_cl; % 1: ALTURA INICIAL - [m]
            subida.h_final = h_final_cl; % 1: ALTURA FINAL - [m]
            handles.subida.h_inicial = subida.h_inicial; % 1: ALTURA INICIAL - [m]
            handles.subida.h_final = subida.h_final; % 1: ALTURA FINAL - [m]
            
            seg(i).datos.h_inicial = h_inicial_cl;% (%)
            seg(i).datos.h_final = h_final_cl;% (%)

            switch climb_mode
                case 2 % 'Subida dados M y gamma';
                    Mach_cl = 0.14;                   
                    gamma_cl = 3*D2R;                   
                    subida.Mach = Mach_cl; % 3: MACH DE VUELO - [-]
                    subida.gamma = gamma_cl; % 2: GAMMA DE SUBIDA
                    handles.subida.Mach = subida.Mach; % 3: MACH DE VUELO - [-]
                    handles.subida.gamma = subida.gamma; % 2: GAMMA DE SUBIDA
            end

        case 11 % case 11 Turn loitter 45 min
            mission_tex = 'Turn';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 
            
            h_inicial_turn = 3000*ft2m;% - [m]
            t_final_turn = 45*60;% - [seg]
            V_turn = 30;
            Mach_turn = V_turn/340;% - [-]

            seg(i).datos.h_inicial = h_inicial_turn;% (%)
            seg(i).datos.t_final = t_final_turn;% (%)
            seg(i).datos.Mach = Mach_cr;
            seg(i).opcion = turn_mode - 1;

            viraje.h_inicial = h_inicial_turn;
            viraje.tiempo_final = t_final_turn;
            handles.viraje.h_inicial = viraje.h_inicial;
            handles.viraje.tiempo_final = viraje.tiempo_final;

         case 12 % case 10 Landing
            mission_tex = 'Landing';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 
        case 13 % case 10 Dummy to account for the 3 segments
            mission_tex = 'Cruise';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 

            h_inicial_cr = 3000;% - [m]
            dist_final_cr = 1;% - [m]
            V_cr = 30;
            Mach_cr = V_cr/340;% - [-]
            crucero.h_inicial = h_inicial_cr;% - [m]
            crucero.dist_final = dist_final_cr;% - [m]
            crucero.Mach = Mach_cr;% - [-]
            handles.crucero.h_inicial = crucero.h_inicial;% - [m]
            handles.crucero.dist_final = dist_final_cr;% - [m]
            handles.crucero.Mach = crucero.Mach;% - [-]
            seg(i).datos.h_inicial = h_inicial_cr;% (%)
            seg(i).datos.dist_final = dist_final_cr;% (%)
            seg(i).datos.Mach = Mach_cr;
            seg(i).opcion = cruise_mode - 1;
    end
    
%     switch type_mission
%         case 1 % case 1 TakeOff
%             % Enter the initial altitude (m):');
%             h_inicial = 0;
%             seg(i).data.h_initial = h_inicial;
%             % Enter the final altitude (m):');
%             h_final = 0;
%             seg(i).data.h_final = h_final;
%             % Enter Take Off Speed (m/s) (It is Just an Estimae):');
%             V_TO = 15;
%             seg(i).data.mision = type_mission; 
%             seg(i).data.V_TO = V_TO;
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(seg(i).data.h_initial,seg(i).data.V_TO);
%             seg(i).data.Data_ATM = Data_ATM;
%             seg(i).data.Performance = Performance;
%         case 2 % case 2 Climb
%             % Enter the initial altitude (m)
%             h_inicial = 0;
%             seg(i).data.h_initial = h_inicial;
%             % Enter the final altitude (m)
%             h_final = 0;
%             seg(i).data.h_final = h_final;
%             % Enter the climb speed (m/s)
%             V_cl = 22;
%             seg(i).data.V_cl = h_final;
%             % Enter the flight ascent path angle (deg):
%             gamma_cl = 3*D2R;
%             seg(i).data.gamma_cl = gamma_cl;
%             
%             seg(i).data.mision = type_mission; 
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(seg(i).data.h_initial,seg(i).data.V_cl);
%             seg(i).data.V_V = seg(i).data.V_cl*sin(seg(i).data.gamma_cl);
%             seg(i).data.Data_ATM = Data_ATM;
%             seg(i).data.Performance = Performance;
%         case 3 % case 3 VTOL Climb
%             % Enter the initial altitude (m) - VTOL Flight
%             h_inicial = 0;
%             seg(i).data.h_initial = h_inicial;
%             % Enter the final altitude (m) - VTOL Flight
%             h_final = 200;
%             seg(i).data.h_final = h_final;
%             % Enter the climb speed (m/s) - VTOL Flight
%             V_VTOL = 5;
%             seg(i).data.V_VTOL = V_VTOL;
%             
%             seg(i).data.mision = type_mission; 
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(seg(i).data.h_initial,seg(i).data.V_VTOL);
%             seg(i).data.Data_ATM = Data_ATM;
%             seg(i).data.Performance = Performance;
%         case 4 % case 4 Cruise(Range)
%             % Enter the initial altitude (m)
%             h_inicial = 200;
%             seg(i).data.h_initial = h_inicial;
%             % Enter the final altitude (m)
%             h_final = 200;
%             seg(i).data.h_final = h_final;
%             % Enter the cruise speed (m/s)
%             V_cr = 26;
%             seg(i).data.V_cr = V_cr;
%             % Enter the Range (Km)
%             Range = 20;
%             seg(i).data.Range = Range*1000;
%             
%             seg(i).data.mision = type_mission;
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(seg(i).data.h_initial,seg(i).data.V_cr);
%             seg(i).data.Data_ATM = Data_ATM;
%             seg(i).data.Performance = Performance;
%         case 5 % case 5 Cruise (Endurance)
%             % Enter the initial altitude (m)
%             h_inicial = 200;
%             seg(i).data.h_initial = h_inicial;
%             % Enter the final altitude (m)
%             h_final = 200;
%             seg(i).data.h_final = h_final;
%             % Enter the cruise speed (m/s)
%             V_cr = 25;
%             seg(i).data.V_cr = V_cr;
%             % Enter the Endurance Time (min)
%             Endurance = 10;
%             seg(i).data.Endurance = Endurance*60;
%             
%             seg(i).data.mision = type_mission; 
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(seg(i).data.h_initial,seg(i).data.V_cr);
%             seg(i).data.Data_ATM = Data_ATM;
%             seg(i).data.Performance = Performance;
%         case 6 % case 6 Descent
%             % Enter the initial altitude (m)
%             h_inicial = 200;
%             seg(i).data.h_initial = h_inicial;
%             % Enter the final altitude (m)
%             h_final = 0;
%             seg(i).data.h_final = h_final;
%             % Enter the descent speed (m/s)
%             V_dsc = 5;
%             seg(i).data.V_dsc = V_dsc;
%             % Enter the flight path angle (min)
%             gamma_dsc = 5*D2R;
%             seg(i).data.gamma_dsc = gamma_dsc;
%             
%             seg(i).data.mision = type_mission; 
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(seg(i).data.h_initial,seg(i).data.V_dsc);
%             seg(i).data.V_V = seg(i).data.V_dsc*sin(seg(i).data.gamma_dsc);
%             seg(i).data.Data_ATM = Data_ATM;
%             seg(i).data.Performance = Performance;
%         case 7 % case 7 Descent (VTOL)
%             % Enter the initial altitude (m)
%             h_inicial = 200;
%             seg(i).data.h_initial = h_inicial;
%             % Enter the final altitude (m)
%             h_final = 0;
%             seg(i).data.h_final = h_final;
%             % Enter the descent speed (m/s) - VTOL
%             V_VTOL = -2;
%             seg(i).data.V_VTOL = V_VTOL;
% 
%             seg(i).data.mision = type_mission; 
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(seg(i).data.h_initial,seg(i).data.V_VTOL);
%             seg(i).data.Data_ATM = Data_ATM;
%             seg(i).data.Performance = Performance;
%     end
end

