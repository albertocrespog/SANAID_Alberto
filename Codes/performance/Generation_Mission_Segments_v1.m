function [seg] = Generation_Mission_Segments_v1(conv_UNITS,num_missions,type_missions,Performance,case_AC,Segments)

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

cruise_mode = Segments{1}.MODES.cruise_mode;
climb_mode = Segments{1}.MODES.climb_mode;
turn_mode = Segments{1}.MODES.turn_mode;
descent_mode = Segments{1}.MODES.descent_mode;

taxy = Segments{1}.taxy;
despegue = Segments{1}.despegue;
subida = Segments{1}.subida;
crucero = Segments{1}.crucero;
viraje = Segments{1}.viraje;
viraje_wt = Segments{1}.viraje_wt;
descenso = Segments{1}.descenso;
aterrizaje = Segments{1}.aterrizaje;
dummy = Segments{1}.dummy;

for i=1:num_missions
    type_mission = type_missions(i);
    switch type_mission
        case 1 % case 1 Taxy
            mission_tex = 'Taxy';            
            seg(i).nombre = mission_tex;
            seg(i).datos.mision = type_mission; 
            seg(i).datos.temp_local = taxy(1);%  % 1: TEMPERATURA LOCAL (K)
            seg(i).datos.h_inicial = taxy(2); % 2: ALTURA LOCAL (m)
            seg(i).datos.P_local = taxy(3); % 3: PRESION LOCAL (Pa)
            seg(i).datos.delta_T = taxy(4); % 4: PALANCA DE RALENTI EN TAXI = 0.05
            seg(i).datos.V_taxy = taxy(5); % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            seg(i).datos.t_tazy = taxy(6); % 6: TIEMPO DE ESPERA EN TAXI (s)
            
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;

        case 2 % case 2 TakeOff
            mission_tex = 'TakeOff';
            seg(i).nombre = mission_tex;
            seg(i).datos.mision = type_mission; 
            seg(i).datos.temp_local = despegue(1);% 1: TEMPERATURA LOCAL (K)
            seg(i).datos.h_inicial = despegue(2); % 2: ALTURA LOCAL (m)
            seg(i).datos.P_local = despegue(3); % 3: PRESION LOCAL (Pa)
            seg(i).datos.mu_takeoff = despegue(4); % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            seg(i).datos.h_obstacle = despegue(5); % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            seg(i).datos.gamma_climb = despegue(6); % 6: GAMMA DE SUBIDA MINIMO
            seg(i).datos.delta_T = despegue(7); % 7: PALANCA DE GASES PARA DESPEGUE
            
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;
            
        case 3 % case 3 Climb
            mission_tex = 'Climb';
            seg(i).datos.mision = type_mission;
            seg(i).nombre = mission_tex;
            seg(i).datos.h_inicial = subida(1);% (%)
            seg(i).datos.h_final = subida(2);% (%)
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;
                        
            switch climb_mode
                case 2 % 'Subida dados M y gamma';
                    seg(i).datos.Mach = subida(4);% (%)
                    seg(i).datos.gamma = subida(3);% (%)
                case 3 % 'Subida dados EAS y gamma';
                    seg(i).datos.EAS = subida(6); % 5: VELOCIDAD EAS  - [m/s]
                    seg(i).datos.gamma = subida(3); % 2: GAMMA DE SUBIDA  - [-]
                case 4 % 'Subida dados TAS y gamma';
                    seg(i).datos.TAS = subida(5); % 4: VELOCIDAD TAS  - [m/s]
                    seg(i).datos.gamma = subida(3); % 2: GAMMA DE SUBIDA  - [-]
                case 5 % 'Subida dados M y palanca';
                    seg(i).datos.Mach = subda(4); % 3: MACH DE VUELO  - [-]
                    handles.subida.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 6 % 'Subida dados EAS y palanca';
                    seg(i).datos.EAS = subida(6); % 5: VELOCIDAD EAS  - [m/s]
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 7 % 'Subida dados TAS y palanca';
                    seg(i).datos.TAS = subida(5); % 4: VELOCIDAD TAS  - [m/s]
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 8 % 'Subida dados V inicial,final y gamma';
                    seg(i).datos.V_ini = subida(8); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
                    seg(i).datos.V_fin = subida(9); % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
                    seg(i).datos.gamma = subida(3); % 2: GAMMA DE SUBIDA  - [-]
                case 9 % 'Subida steppest climb';
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 10 % 'Subida fastest climb';
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 11 % 'Subida dados V inicial,final y palanca'
                    seg(i).datos.V_ini = subida(8); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA) - [m/s]
                    seg(i).datos.V_fin = subida(9);% 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA) - [m/s]
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES - [-]
            end
            seg(i).opcion = climb_mode - 1;
        case 4 % case 4 VTOL Climb
            mission_tex = 'VTOL_Climb';
            seg(i).datos.mision = type_mission; 
        case 5 % case 5 Cruise
            mission_tex = 'Cruise';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 
%             crucero(1) = h_inicial_cr;% - [m] % Altura inicial
%             crucero(2) = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
%             crucero(3) = V_cr; % Velocidad de crucero
%             crucero(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
%             crucero(5) = CL_cr;% - [-]% 3: CL DE CRUCERO
%             crucero(6) = delta_T_cr;% - []% 4: PALANCA DE GASES
%             crucero(7) = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
%             crucero(8) = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
%             crucero(9) = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
%             crucero(10) = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
%             crucero(11) = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
%             crucero(12) = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL

            switch cruise_mode
                case 2 % 'Crucero dado M y distancia'
                    seg(i).datos.h_inicial = crucero(1); % - [m] % Altura inicial
                    seg(i).datos.dist_final = crucero(2); % - [m] % 1: DISTANCIA FINAL
                    seg(i).datos.V = crucero(3); % Velocidad de crucero
                    seg(i).datos.Mach = crucero(4); % 2: MACH DE VUELO
                case 3 % ;'Crucero dado CL y distancia'
                    seg(i).datos.h_inicial = crucero(1); % - [m] % Altura inicial
                    seg(i).datos.dist_final = crucero(2); % - [m] % 1: DISTANCIA FINAL
                    seg(i).datos.CL = crucero(5); % 3: CL DE CRUCERO
                case 4 % ;'Crucero dados V inicial,final y palanca'
                    seg(i).datos.h_inicial = crucero(1); % - [m] % Altura inicial
                    seg(i).datos.dist_final = crucero(2); % - [m] % 1: DISTANCIA FINAL
                    seg(i).datos.palanca = crucero(6); % 4: PALANCA DE GASES
                    seg(i).datos.V_ini = crucero(7); % 5: VELOCIDAD INICIAL
                    seg(i).datos.V_fin = crucero(8); % 6: VELOCIDAD FINAL
                case 5 % 'Crucero con polar en funcion de M';
                    seg(i).datos.h_inicial = crucero(1); % - [m] % Altura inicial
                    seg(i).datos.dist_final = crucero(2); % - [m] % 1: DISTANCIA FINAL
                    seg(i).datos.V = crucero(3); % Velocidad de crucero
                    seg(i).datos.Mach = crucero(4); % 2: MACH DE VUELO
                    seg(i).datos.Cd0 = crucero(10); % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
                    seg(i).datos.k1 = crucero(11); % 9: K1 = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
                    seg(i).datos.k2 = crucero(12); % 10: K2 = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
                case 6 % 'Crucero de max alcance dado peso final'
                    seg(i).datos.h_inicial = crucero(1); % - [m] % Altura inicial
                    seg(i).datos.fuel = crucero(9); % 7: COMBUSTIBLE A QUEMAR
                case 7 % 'Crucero de max autonomia dado peso final'
                    seg(i).datos.h_inicial = crucero(1); % - [m] % Altura inicial
                    seg(i).datos.fuel = crucero(9); % 7: COMBUSTIBLE A QUEMAR
            end
            seg(i).opcion = cruise_mode - 1;
                
            [Data_ATM Performance] = Flight_Conditions_2020_v1(crucero(1),crucero(3));
            seg(i).datos.Data_ATM = Data_ATM;
            seg(i).datos.Performance = Performance;
 
        case 6 % case 6 Load Deployment
            mission_tex = 'Load_Deployment';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 
        case 7 % case 7 - Turn
             
%             viraje(1) = h_inicial_tr;% - [m]
%             viraje(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
%             viraje(3) = V_turn; % turn velocity
%             viraje(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
%             viraje(5) = delta_T_tr; % 3: PALANCA DE GASES
%             viraje(6) = CL_tr;% 4: CL DE VIRAJE
%             viraje(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
%             viraje(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
%             viraje(9) = n_tr;% 7: FACTOR DE CARGA
%             viraje(10) = R_tr;% 8: RADIO DE GIRO (m)

            mission_tex = 'Turn';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 
            seg(i).datos.h_inicial = viraje(1);% (%)
            seg(i).datos.t_final = viraje(2);% (%)
            switch turn_mode
                case 2 % 'Viraje horizontal dado V y palanca'
                    seg(i).datos.velocidad = viraje(3);% - [m/s]
                    seg(i).datos.palanca = viraje(5);% - []
                case 3 % 'Viraje horizontal dado V y CL'
                    seg(i).datos.velocidad = viraje(3);% - [m/s]
                    seg(i).datos.CL = viraje(6);% - []
                case 4 % 'Viraje horizontal dado V y balance'
                    seg(i).datos.velocidad = viraje(3);
                    seg(i).datos.balance = viraje(7);
                case 5 % 'Viraje horizontal dado V y n'
                    seg(i).datos.velocidad = viraje(3);
                    seg(i).datos.n = viraje(9);
                case 6 % 'V.H dado V y radio de giro '
                    seg(i).datos.velocidad = viraje(3);
                    seg(i).datos.radio = viraje(10);
                case 7 % 'V.H dado V y velocidad de guiñada';...
                    seg(i).datos.velocidad = viraje(3);
                    seg(i).datos.vel_guiniada = viraje(8);
                case 8 % 'V.H dado palanca y a factor de carga max'
                    seg(i).datos.palanca = viraje(5);
                case 9 % 'V.H dado palanca y a v de guiñada max'
                    seg(i).datos.palanca = viraje(5);
                case 10 % 'V.H dado palanca y a radio min'
                    seg(i).datos.palanca = viraje(5);
            end
            seg(i).opcion = turn_mode - 1;

            [Data_ATM Performance] = Flight_Conditions_2020_v1(viraje(1),viraje(3));
            seg(i).datos.Data_ATM = Data_ATM;
            seg(i).datos.Performance = Performance;
                    
        case 8 % case 8 Descent
            mission_tex = 'Descent';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission;
            %             descenso(1) = h_inicial_d;
            %             descenso(2) = h_final_d;
            %             descenso(3) = gamma_d;
            %             descenso(4) = V_d;
            %             descenso(5) = Mach_d;
            %             descenso(6) = EAS_d;
            %             descenso(7) = TAS_d;
            %             descenso(8) = delta_T_d;
            %             descenso(9) = V_ini_d;
            %             descenso(10) = V_fin_d;
            seg(i).datos.h_inicial = descenso(1);
            seg(i).datos.h_final = descenso(2);
            switch descent_mode
                case 2 % 'Descenso dados M y gamma'
                    seg(i).datos.Mach = descenso(5);
                    seg(i).datos.gamma = descenso(3);
                case 3 % 'Descenso dados EAS y gamma';
                    seg(i).datos.EAS = descenso(6);
                    seg(i).datos.gamma = descenso(3);
                case 4 % 'Descenso dados TAS y gamma'
                    seg(i).datos.TAS = descenso(7);
                    seg(i).datos.gamma = descenso(3);
                case 5 % 'Descenso dados M y palanca'
                    seg(i).datos.Mach = descenso(5);
                    seg(i).datos.palanca = descenso(8);
                case 6 % 'Descenso dados EAS y palanca'
                    seg(i).datos.EAS = descenso(6);
                    seg(i).datos.palanca = descenso(8);
                case 7 % 'Descenso dados TAS y palanca'
                    seg(i).datos.TAS = descenso(7);
                    seg(i).datos.palanca = descenso(8);
                case 8 % 'Descenso dados V inicial,final y gamma'
                    seg(i).datos.V_ini = descenso(9);
                    seg(i).datos.V_fin = descenso(10);
                    seg(i).datos.gamma = descenso(3);
                case 9 % 'Descenso a minimo gamma'
                    seg(i).datos.palanca = descenso(8);
                case 10 % 'Descenso "slowest sink"'
                    seg(i).datos.palanca = descenso(8);
                case 11 % 'Descenso dados V inicial,final y palanca'
                    seg(i).datos.V_ini = descenso(9);
                    seg(i).datos.V_fin = descenso(10);
                    seg(i).datos.palanca = descenso(8);
            end
            seg(i).opcion = descent_mode - 1;
        case 9 % case 9 Descent (VTOL)
            mission_tex = 'Descent_VTOL';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission;
            
        case 10 % case 10 Climb Waiting Area to 3000ft
            mission_tex = 'Climb_waiting';
            seg(i).datos.mision = type_mission;
            seg(i).nombre = mission_tex;
            seg(i).datos.h_inicial = subida(1);% (%)
            seg(i).datos.h_final = subida(2);% (%)
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;
            
            switch climb_mode
                case 2 % 'Subida dados M y gamma';
                    seg(i).datos.Mach = subida(4);% (%)
                    seg(i).datos.gamma = subida(3);% (%)
                case 3 % 'Subida dados EAS y gamma';
                    seg(i).datos.EAS = subida(6); % 5: VELOCIDAD EAS  - [m/s]
                    seg(i).datos.gamma = subida(3); % 2: GAMMA DE SUBIDA  - [-]
                case 4 % 'Subida dados TAS y gamma';
                    seg(i).datos.TAS = subida(5); % 4: VELOCIDAD TAS  - [m/s]
                    seg(i).datos.gamma = subida(3); % 2: GAMMA DE SUBIDA  - [-]
                case 5 % 'Subida dados M y palanca';
                    seg(i).datos.Mach = subda(4); % 3: MACH DE VUELO  - [-]
                    handles.subida.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 6 % 'Subida dados EAS y palanca';
                    seg(i).datos.EAS = subida(6); % 5: VELOCIDAD EAS  - [m/s]
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 7 % 'Subida dados TAS y palanca';
                    seg(i).datos.TAS = subida(5); % 4: VELOCIDAD TAS  - [m/s]
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 8 % 'Subida dados V inicial,final y gamma';
                    seg(i).datos.V_ini = subida(8); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
                    seg(i).datos.V_fin = subida(9); % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
                    seg(i).datos.gamma = subida(3); % 2: GAMMA DE SUBIDA  - [-]
                case 9 % 'Subida steppest climb';
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 10 % 'Subida fastest climb';
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES  - [-]
                case 11 % 'Subida dados V inicial,final y palanca'
                    seg(i).datos.V_ini = subida(8); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA) - [m/s]
                    seg(i).datos.V_fin = subida(9);% 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA) - [m/s]
                    seg(i).datos.palanca = subida(7); % 6: PALANCA DE GASES - [-]
            end
            seg(i).opcion = climb_mode - 1;

        case 11 % case 11 Turn loitter 45 min
            mission_tex = 'Turn_waiting';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 
            %             viraje(1) = h_inicial_tr;% - [m]
%             viraje(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
%             viraje(3) = V_turn; % turn velocity
%             viraje(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
%             viraje(5) = delta_T_tr; % 3: PALANCA DE GASES
%             viraje(6) = CL_tr;% 4: CL DE VIRAJE
%             viraje(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
%             viraje(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
%             viraje(9) = n_tr;% 7: FACTOR DE CARGA
%             viraje(10) = R_tr;% 8: RADIO DE GIRO (m)

            seg(i).datos.h_inicial = viraje_wt(1);% (%)
            seg(i).datos.t_final = viraje_wt(2);% (%)
            switch turn_mode
                case 2 % 'Viraje horizontal dado V y palanca'
                    seg(i).datos.velocidad = viraje_wt(3);% - [m/s]
                    seg(i).datos.palanca = viraje_wt(5);% - []
                case 3 % 'Viraje horizontal dado V y CL'
                    seg(i).datos.velocidad = viraje_wt(3);% - [m/s]
                    seg(i).datos.CL = viraje_wt(6);% - []
                case 4 % 'Viraje horizontal dado V y balance'
                    seg(i).datos.velocidad = viraje_wt(3);
                    seg(i).datos.balance = viraje_wt(7);
                case 5 % 'Viraje horizontal dado V y n'
                    seg(i).datos.velocidad = viraje_wt(3);
                    seg(i).datos.n = viraje_wt(9);
                case 6 % 'V.H dado V y radio de giro '
                    seg(i).datos.velocidad = viraje_wt(3);
                    seg(i).datos.radio = viraje_wt(10);
                case 7 % 'V.H dado V y velocidad de guiñada';...
                    seg(i).datos.velocidad = viraje_wt(3);
                    seg(i).datos.vel_guiniada = viraje_wt(8);
                case 8 % 'V.H dado palanca y a factor de carga max'
                    seg(i).datos.palanca = viraje_wt(5);
                case 9 % 'V.H dado palanca y a v de guiñada max'
                    seg(i).datos.palanca = viraje_wt(5);
                case 10 % 'V.H dado palanca y a radio min'
                    seg(i).datos.palanca = viraje_wt(5);
            end
            seg(i).opcion = turn_mode - 1;

            [Data_ATM Performance] = Flight_Conditions_2020_v1(viraje_wt(1),viraje_wt(3));
            seg(i).datos.Data_ATM = Data_ATM;
            seg(i).datos.Performance = Performance;
%             
%             
%             
%             h_inicial_turn = 3000*ft2m;% - [m]
%             t_final_turn = 45*60;% - [seg]
%             V_turn = 30;
%             Mach_turn = V_turn/340;% - [-]
%             
%             seg(i).datos.h_inicial = h_inicial_turn;% (%)
%             seg(i).datos.t_final = t_final_turn;% (%)
%             seg(i).datos.Mach = Mach_turn;
%             seg(i).opcion = turn_mode - 1;
            
        case 12 % case 12 Landing
            mission_tex = 'Landing';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission;
                        
            mission_tex = 'Landing';
            seg(i).nombre = mission_tex;
            seg(i).datos.mision = type_mission;
            
            seg(i).datos.temp_local = aterrizaje(1);% 1: TEMPERATURA LOCAL (K)
            seg(i).datos.h_inicial = aterrizaje(2); % 2: ALTURA LOCAL (m)
            seg(i).datos.P_local = aterrizaje(3); % 3: PRESION LOCAL (Pa)
            seg(i).datos.mu_landing = aterrizaje(4); % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            seg(i).datos.palanca = aterrizaje(5); % 5: PALANCA DE GASES PARA DESPEGUE
            seg(i).datos.t_brake = aterrizaje(6); % 'Tiempo en activar frenos' - [s]
            
            %             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
            %             seg(i).datos.Data_ATM = Data_ATM;
            %             seg(i).datos.Performance = Performance;
            
        case 13 % case 10 Dummy to account for the 3 segments
            mission_tex = 'Cruise';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 

            seg(i).datos.h_inicial = dummy(1);% (%)
            seg(i).datos.dist_final = dummy(2);% (%)
            seg(i).datos.V = dummy(3);
            seg(i).datos.Mach = dummy(4);
            seg(i).opcion = cruise_mode - 1;
    end
end