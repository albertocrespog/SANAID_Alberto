function [seg] = Generation_Mission_Segments_v3(conv_UNITS,num_missions,type_missions,case_AC,Segments,FlightMODE)

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
ft2m = conv_UNITS.ft2m;
m2ft = conv_UNITS.m2ft;

% Asks user for the input data        
% Enter type of mission segments between brackets being
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
            % Assigns data for type of mission
            taxy = Segments{1}.taxy;
            mission_tex = 'Taxy';            
            seg(i).nombre = mission_tex;
            seg(i).datos.mision = type_mission; 
            seg(i).datos.temp_local = taxy.T(i);%  % 1: TEMPERATURA LOCAL (K)
            seg(i).datos.h_inicial = taxy.h(i); % 2: ALTURA LOCAL (m)
            seg(i).datos.P_local = taxy.P(i); % 3: PRESION LOCAL (Pa)
            seg(i).datos.delta_T = taxy.dT(i); % 4: PALANCA DE RALENTI EN TAXI = 0.05
            seg(i).datos.V_taxy = taxy.V(i); % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            seg(i).datos.t_taxy = taxy.t(i); % 6: TIEMPO DE ESPERA EN TAXI (s)
            
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;

        case 2 % case 2 TakeOff
            % Assigns data for type of mission
            despegue = Segments{1}.despegue;
            mission_tex = 'TakeOff';
            seg(i).nombre = mission_tex;
            seg(i).datos.mision = type_mission; 
            seg(i).datos.temp_local = despegue.T(i);% 1: TEMPERATURA LOCAL (K)
            seg(i).datos.h_inicial = despegue.h(i); % 2: ALTURA LOCAL (m)
            seg(i).datos.P_local = despegue.P(i); % 3: PRESION LOCAL (Pa)
            seg(i).datos.mu_takeoff = despegue.mu(i); % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            seg(i).datos.h_obstacle = despegue.ho(i); % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            seg(i).datos.gamma_climb = despegue.g(i); % 6: GAMMA DE SUBIDA MINIMO
            seg(i).datos.delta_T = despegue.dT(i); % 7: PALANCA DE GASES PARA DESPEGUE
            
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;
            
        case 3 % case 3 Climb
            % Assigns data for type of mission
            subida = Segments{1}.subida;
            mission_tex = 'Climb';
            seg(i).datos.mision = type_mission;
            seg(i).nombre = mission_tex;
            seg(i).datos.h_inicial = subida.hi(i);% (%)
            seg(i).datos.h_final = subida.hf(i);% (%)
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;
            climb_mode = FlightMODE(i)+1;
            switch climb_mode
                case 2 % 'Subida dados M y gamma';
                    seg(i).datos.Mach = subida.M(i);% (%)
                    seg(i).datos.gamma = subida.g(i);% (%)
                case 3 % 'Subida dados EAS y gamma';
                    seg(i).datos.EAS = subida.E(i); % 5: VELOCIDAD EAS  - [m/s]
                    seg(i).datos.gamma = subida.g(i); % 2: GAMMA DE SUBIDA  - [-]
                case 4 % 'Subida dados TAS y gamma';
                    seg(i).datos.TAS = subida.T(i); % 4: VELOCIDAD TAS  - [m/s]
                    seg(i).datos.gamma = subida.g(i); % 2: GAMMA DE SUBIDA  - [-]
                case 5 % 'Subida dados M y palanca';
                    seg(i).datos.Mach = subida.M(i); % 3: MACH DE VUELO  - [-]
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES  - [-]
                case 6 % 'Subida dados EAS y palanca';
                    seg(i).datos.EAS = subida.E(i); % 5: VELOCIDAD EAS  - [m/s]
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES  - [-]
                case 7 % 'Subida dados TAS y palanca';
                    seg(i).datos.TAS = subida.T(i); % 4: VELOCIDAD TAS  - [m/s]
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES  - [-]
                case 8 % 'Subida dados V inicial,final y gamma';
                    seg(i).datos.V_ini = subida.vi(i); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
                    seg(i).datos.V_fin = subida.vf(i); % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
                    seg(i).datos.gamma = subida.g(i); % 2: GAMMA DE SUBIDA  - [-]
                case 9 % 'Subida steppest climb';
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES  - [-]
                case 10 % 'Subida fastest climb';
                    seg(i).datos.palanca = subidadT(i); % 6: PALANCA DE GASES  - [-]
                case 11 % 'Subida dados V inicial,final y palanca'
                    seg(i).datos.V_ini = subida.vi(i); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA) - [m/s]
                    seg(i).datos.V_fin = subida.vf(i);% 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA) - [m/s]
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES - [-]
            end
            seg(i).opcion = climb_mode - 1;
        case 4 % case 4 VTOL Climb
            mission_tex = 'VTOL_Climb';
            seg(i).datos.mision = type_mission;% Assigns data for type of mission
            subidavt = Segments{1}.subidavt;
            seg(i).nombre = mission_tex;
            seg(i).datos.h_inicial = subidavt.hi(i);% 
            climbvt_mode = FlightMODE(i)+1;
            switch climbvt_mode
                case 2
                    seg(i).datos.h_final = subidavt.hf(i);
                    seg(i).datos.palanca = subidavt.dT(i);
                case 3
                    seg(i).datos.thover = subidavt.thover(i);
                case 4
                    seg(i).datos.h_final = subidavt.hf(i);
                    seg(i).datos.vclimb = subidavt.vc(i);
                case 5
                    seg(i).datos.mbathover = subidavt.mbathover(i);
            end
            seg(i).opcion = climbvt_mode - 1;
        case 5 % case 5 Cruise
            
            % Assigns data for type of mission
            crucero = Segments{1}.crucero;
            mission_tex = 'Cruise';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission;
            %             crucero.hi(i) = h_inicial_cr;% - [m] % Altura inicial
            %             crucero.d(i) = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
            %             crucero.vcr(i) = V_cr; % Velocidad de crucero
            %             crucero.M(i) = Mach_cr;% - [-] % 2: MACH DE VUELO
            %             crucero.CL(i) = CL_cr;% - [-]% 3: CL DE CRUCERO
            %             crucero.dT(i) = delta_T_cr;% - []% 4: PALANCA DE GASES
            %             crucero.vi(i) = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
            %             crucero.vf(i) = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
            %             crucero.f(i) = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            %             crucero.cd0(i) = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            %             crucero.k1(i) = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            %             crucero.k2(i) = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            
            cruise_mode = FlightMODE(i)+1;
            
            switch cruise_mode
                case 2 % 'Crucero dado M y distancia'
                    seg(i).datos.h_inicial = crucero.hi(i); % - [m] % Altura inicial
                    seg(i).datos.dist_final = crucero.d(i); % - [m] % 1: DISTANCIA FINAL
                    seg(i).datos.V = crucero.vcr(i); % Velocidad de crucero
                    seg(i).datos.Mach = crucero.M(i); % 2: MACH DE VUELO
                case 3 % ;'Crucero dado CL y distancia'
                    seg(i).datos.h_inicial = crucero.hi(i); % - [m] % Altura inicial
                    seg(i).datos.dist_final = crucero.d(i); % - [m] % 1: DISTANCIA FINAL
                    seg(i).datos.CL = crucero.CL;
                    % 3: CL DE CRUCERO
                case 4 % ;'Crucero dados V inicial,final y palanca'
                    seg(i).datos.h_inicial = crucero.hi(i); % - [m] % Altura inicial
                    seg(i).datos.dist_final = crucero.d(i); % - [m] % 1: DISTANCIA FINAL
                    seg(i).datos.palanca = crucero.dT(i); % 4: PALANCA DE GASES
                    seg(i).datos.V_ini = crucero.vi(i); % 5: VELOCIDAD INICIAL
                    seg(i).datos.V_fin = crucero.vf(i); % 6: VELOCIDAD FINAL
                case 5 % 'Crucero con polar en funcion de M';
                    seg(i).datos.h_inicial = crucero.hi(i); % - [m] % Altura inicial
                    seg(i).datos.dist_final = crucero.d(i); % - [m] % 1: DISTANCIA FINAL
                    seg(i).datos.V = crucero.vcr(i); % Velocidad de crucero
                    seg(i).datos.Mach = crucero.M(i); % 2: MACH DE VUELO
                    seg(i).datos.Cd0 = crucero.cd0(i); % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
                    seg(i).datos.k1 = crucero.k1(i); % 9: K1 = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
                    seg(i).datos.k2 = crucero.k2(i); % 10: K2 = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
                case 6 % 'Crucero de max alcance dado peso final'
                    seg(i).datos.h_inicial = crucero.hi(i); % - [m] % Altura inicial
                    seg(i).datos.fuel = crucero.f(i); % 7: COMBUSTIBLE A QUEMAR
                case 7 % 'Crucero de max autonomia dado peso final'
                    seg(i).datos.h_inicial = crucero.hi(i); % - [m] % Altura inicial
                    seg(i).datos.fuel = crucero.f(i); % 7: COMBUSTIBLE A QUEMAR
                case 8 % 'Crucero dado Mach y masa en bruto de baterías disponibles'
                    seg(i).datos.h_inicial = crucero.hi(i); % - [m] % Altura inicial
                    seg(i).datos.Mach = crucero.M(i); % 2: MACH DE VUELO
                    seg(i).datos.m_bat = crucero.mbat(i); % 11: MASA DE BAT
            end
            seg(i).opcion = cruise_mode - 1;
                
            [Data_ATM Performance] = Flight_Conditions_2020_v1(crucero.hi(i),crucero.vcr(i));
            seg(i).datos.Data_ATM = Data_ATM;
            seg(i).datos.Performance = Performance;
 
        case 6 % case 6 Load Deployment
            % Assigns data for type of mission
%             crucero = Segments{1}.crucero;
            loadep = Segments{1}.loadep;
            mission_tex = 'Load_Deployment';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission;
            seg(i).datos.carga = loadep.carga(i);
        case 7 % case 7 - Turn
            % Assigns data for type of mission
            viraje = Segments{1}.viraje;
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
            seg(i).datos.h_inicial = viraje.hi(i);% (%)
            seg(i).datos.t_final = viraje.tf(i);% (%)
            
            turn_mode = FlightMODE(i) + 1;
            switch turn_mode
                case 2 % 'Viraje horizontal dado V y palanca'
                    seg(i).datos.velocidad = viraje.vtr(i);% - [m/s]
                    seg(i).datos.palanca = viraje.dT(i);% - []
                case 3 % 'Viraje horizontal dado V y CL'
                    seg(i).datos.velocidad = viraje.vtr(i);% - [m/s]
                    seg(i).datos.CL = viraje.CL(i);% - []
                case 4 % 'Viraje horizontal dado V y balance'
                    seg(i).datos.velocidad = viraje.vtr(i);
                    seg(i).datos.balance = viraje.phi(i);
                case 5 % 'Viraje horizontal dado V y n'
                    seg(i).datos.velocidad = viraje.vtr(i);
                    seg(i).datos.n = viraje.n(i);
                case 6 % 'V.H dado V y radio de giro '
                    seg(i).datos.velocidad = viraje.vtr(i);
                    seg(i).datos.radio = viraje.R(i);
                case 7 % 'V.H dado V y velocidad de guiñada';...
                    seg(i).datos.velocidad = viraje.vtr(i);
                    seg(i).datos.vel_guiniada = viraje.vpsi(i);
                case 8 % 'V.H dado palanca y a factor de carga max'
                    seg(i).datos.palanca = viraje.dT(i);
                case 9 % 'V.H dado palanca y a v de guiñada max'
                    seg(i).datos.palanca = viraje.dT(i);
                case 10 % 'V.H dado palanca y a radio min'
                    seg(i).datos.palanca = viraje.dT(i);
            end
            seg(i).opcion = turn_mode - 1;

            [Data_ATM Performance] = Flight_Conditions_2020_v1(viraje.hi(i),viraje.vtr(i));
            seg(i).datos.Data_ATM = Data_ATM;
            seg(i).datos.Performance = Performance;
                    
        case 8 % case 8 Descent
            % Assigns data for type of mission
            descenso = Segments{1}.descenso;
            mission_tex = 'Descent';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission;
            %             descenso.hi(i) = h_inicial_d;
            %             descenso.hf(i) = h_final_d;
            %             descenso.g(i) = gamma_d;
            %             descenso(4) = V_d;
            %             descenso.M(i) = Mach_d;
            %             descenso.E(i) = EAS_d;
            %             descenso.T(i) = TAS_d;
            %             descenso.dT(i) = delta_T_d;
            %             descenso.vi(i) = V_ini_d;
            %             descenso.vf(i) = V_fin_d;
            seg(i).datos.h_inicial = descenso.hi(i);
            seg(i).datos.h_final = descenso.hf(i);
    
            descent_mode = FlightMODE(i)+1;
            switch descent_mode
                case 2 % 'Descenso dados M y gamma'
                    seg(i).datos.Mach = descenso.M(i);
                    seg(i).datos.gamma = descenso.g(i);
                case 3 % 'Descenso dados EAS y gamma';
                    seg(i).datos.EAS = descenso.E(i);
                    seg(i).datos.gamma = descenso.g(i);
                case 4 % 'Descenso dados TAS y gamma'
                    seg(i).datos.TAS = descenso.T(i);
                    seg(i).datos.gamma = descenso.g(i);
                case 5 % 'Descenso dados M y palanca'
                    seg(i).datos.Mach = descenso.M(i);
                    seg(i).datos.palanca = descenso.dT(i);
                case 6 % 'Descenso dados EAS y palanca'
                    seg(i).datos.EAS = descenso.E(i);
                    seg(i).datos.palanca = descenso.dT(i);
                case 7 % 'Descenso dados TAS y palanca'
                    seg(i).datos.TAS = descenso.T(i);
                    seg(i).datos.palanca = descenso.dT(i);
                case 8 % 'Descenso dados V inicial,final y gamma'
                    seg(i).datos.V_ini = descenso.vi(i);
                    seg(i).datos.V_fin = descenso.vf(i);
                    seg(i).datos.gamma = descenso.g(i);
                case 9 % 'Descenso a minimo gamma'
                    seg(i).datos.palanca = descenso.dT(i);
                case 10 % 'Descenso "slowest sink"'
                    seg(i).datos.palanca = descenso.dT(i);
                case 11 % 'Descenso dados V inicial,final y palanca'
                    seg(i).datos.V_ini = descenso.vi(i);
                    seg(i).datos.V_fin = descenso.vf(i);
                    seg(i).datos.palanca = descenso.dT(i);
            end
            seg(i).opcion = descent_mode - 1;
        case 9 % case 9 Descent (VTOL)
            mission_tex = 'Descent_VTOL';
            seg(i).datos.mision = type_mission;% Assigns data for type of mission
            descensovt = Segments{1}.descensovt;
            seg(i).nombre = mission_tex;
            seg(i).datos.h_inicial = descensovt.hi(i);% 
            descentvt_mode = FlightMODE(i)+1;
            switch descentvt_mode
                case 2
                    seg(i).datos.h_final = descensovt.hf(i);
                    seg(i).datos.vdes = descensovt.vd(i);
            end
            seg(i).opcion = descentvt_mode - 1;
        case 10 % case 10 Climb Waiting Area to 3000ft
            % Assigns data for type of mission
            subida = Segments{1}.subida;
            mission_tex = 'Climb_waiting';
            seg(i).datos.mision = type_mission;
            seg(i).nombre = mission_tex;
            seg(i).datos.h_inicial = subida.hi(i);% (%)
            seg(i).datos.h_final = subida.hf(i);% (%)
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
%             seg(i).datos.Data_ATM = Data_ATM;
%             seg(i).datos.Performance = Performance;
            climb_mode = FlightMODE(i)+1;
            switch climb_mode
                case 2 % 'Subida dados M y gamma';
                    seg(i).datos.Mach = subida.M(i);% (%)
                    seg(i).datos.gamma = subida.g(i);% (%)
                case 3 % 'Subida dados EAS y gamma';
                    seg(i).datos.EAS = subida.E(i); % 5: VELOCIDAD EAS  - [m/s]
                    seg(i).datos.gamma = subida.g(i); % 2: GAMMA DE SUBIDA  - [-]
                case 4 % 'Subida dados TAS y gamma';
                    seg(i).datos.TAS = subida.T(i); % 4: VELOCIDAD TAS  - [m/s]
                    seg(i).datos.gamma = subida.g(i); % 2: GAMMA DE SUBIDA  - [-]
                case 5 % 'Subida dados M y palanca';
                    seg(i).datos.Mach = subida.M(i); % 3: MACH DE VUELO  - [-]
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES  - [-]
                case 6 % 'Subida dados EAS y palanca';
                    seg(i).datos.EAS = subida.E(i); % 5: VELOCIDAD EAS  - [m/s]
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES  - [-]
                case 7 % 'Subida dados TAS y palanca';
                    seg(i).datos.TAS = subida.T(i); % 4: VELOCIDAD TAS  - [m/s]
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES  - [-]
                case 8 % 'Subida dados V inicial,final y gamma';
                    seg(i).datos.V_ini = subida.vi(i); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
                    seg(i).datos.V_fin = subida.vf(i); % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
                    seg(i).datos.gamma = subida.g(i); % 2: GAMMA DE SUBIDA  - [-]
                case 9 % 'Subida steppest climb';
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES  - [-]
                case 10 % 'Subida fastest climb';
                    seg(i).datos.palanca = subidadT(i); % 6: PALANCA DE GASES  - [-]
                case 11 % 'Subida dados V inicial,final y palanca'
                    seg(i).datos.V_ini = subida.vi(i); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA) - [m/s]
                    seg(i).datos.V_fin = subida.vf(i);% 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA) - [m/s]
                    seg(i).datos.palanca = subida.dT(i); % 6: PALANCA DE GASES - [-]
            end
            seg(i).opcion = climb_mode - 1;
        case 11 % case 11 Turn loitter 45 min
            % Assigns data for type of mission
            viraje_wt = Segments{1}.viraje_wt;

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
            
            % Assigns data for type of mission
            aterrizaje = Segments{1}.aterrizaje;
            
            mission_tex = 'Landing';
            seg(i).nombre = mission_tex;
            seg(i).datos.mision = type_mission;
            
            seg(i).datos.temp_local = aterrizaje.T(i);% 1: TEMPERATURA LOCAL (K)
            seg(i).datos.h_inicial = aterrizaje.h(i); % 2: ALTURA LOCAL (m)
            seg(i).datos.P_local = aterrizaje.P(i); % 3: PRESION LOCAL (Pa)
            seg(i).datos.mu_landing = aterrizaje.mu(i); % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            seg(i).datos.palanca = aterrizaje.dT(i); % 5: PALANCA DE GASES PARA DESPEGUE
            seg(i).datos.t_brake = aterrizaje.tb(i); % 'Tiempo en activar frenos' - [s]
            
            %             [Data_ATM Performance] = Flight_Conditions_2020_v1(h_initial,V_taxy);
            %             seg(i).datos.Data_ATM = Data_ATM;
            %             seg(i).datos.Performance = Performance;
            
        case 13 % case 10 Dummy to account for the 3 segments
            
            % Assigns data for type of mission
            dummy = Segments{1}.dummy;
            
            mission_tex = 'Dummy';
            seg(i).nombre = mission_tex;%
            seg(i).datos.mision = type_mission; 

            seg(i).datos.h_inicial = dummy(1);% (%)
            seg(i).datos.dist_final = dummy(2);% (%)
            seg(i).datos.V = dummy(3);
            seg(i).datos.Mach = dummy(4);
            seg(i).opcion = 1;
    end
end