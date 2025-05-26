function Segments = Define_Segments(Aero,Aero_TH,case_AC,conv_UNITS,FlightMODE,num_missions)

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
ft2m = conv_UNITS.ft2m;
m2ft = conv_UNITS.m2ft;

climb_mode = FlightMODE.climb_mode;
turn_mode = FlightMODE.turn_mode;
descent_mode = FlightMODE.descent_mode;
cruise_mode = FlightMODE.cruise_mode;

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

%% Climb options
% case 2 % 'Subida dados M y gamma';
% case 3 % 'Subida dados EAS y gamma';
% case 4 % 'Subida dados TAS y gamma';
% case 5 % 'Subida dados M y palanca';
% case 6 % 'Subida dados EAS y palanca';
% case 7 % 'Subida dados TAS y palanca';
% case 8 % 'Subida dados V inicial,final y gamma';
% case 9 % 'Subida steppest climb';
% case 10 % 'Subida fastest climb';
% case 11 % 'Subida dados V inicial,final y palanca'
%% Cruise options
% case 2 % 'Crucero dado M y distancia'
% case 3 % ;'Crucero dado CL y distancia'
% case 4 % ;'Crucero dados V inicial,final y palanca'
% case 5 % 'Crucero con polar en funcion de M';
% case 6 % 'Crucero de max alcance dado peso final'
% case 7 % 'Crucero de max autonomia dado peso final'
%% Turn options
% case 2 % 'Viraje horizontal dado V y palanca'
% case 3 % 'Viraje horizontal dado V y CL'
% case 4 % 'Viraje horizontal dado V y balance'
% case 5 % 'Viraje horizontal dado V y n'
% case 6 % 'V.H dado V y radio de giro '
% case 7 % 'V.H dado V y velocidad de guiñada';...
% case 8 % 'V.H dado palanca y a factor de carga max'
% case 9 % 'V.H dado palanca y a v de guiñada max'
% case 10 % 'V.H dado palanca y a radio min'
%% Descent options
% case 2 % 'Descenso dados M y gamma'
% case 3 % 'Descenso dados EAS y gamma';
% case 4 % 'Descenso dados TAS y gamma'
% case 5 % 'Descenso dados M y palanca'
% case 6 % 'Descenso dados EAS y palanca'
% case 7 % 'Descenso dados TAS y palanca'
% case 8 % 'Descenso dados V inicial,final y gamma'
% case 9 % 'Descenso a minimo gamma'
% case 10 % 'Descenso "slowest sink"'
% case 11 % 'Descenso dados V inicial,final y palanca'
                                           

%% Define segment conditions for each AC
% introduce the value for each input condition, select the climb, cruise
% and turn options and if a segment is repeated, introduce as a vector with
% as many elemts as segments.
% For Those inputs that are not 
% for i=1:num_missions
    switch case_AC
        case 1 % case_AC = 1 - EMERGENTIA 1:1
            
            % Climb
            h_inicial_cl = 0; % 1: ALTURA INICIAL - [m]
            h_final_cl = 5; % 1: ALTURA FINAL - [m]
            gamma_cl = 3*D2R; % 2: GAMMA DE SUBIDA  - [-]
            Mach_cl = 0.14; % 3: MACH DE VUELO  - [-]
            TAS_cl = 20; % 4: VELOCIDAD TAS  - [m/s]
            EAS_cl = 20; % 5: VELOCIDAD EAS  - [m/s]
            delta_T_cl = 0.95; % 6: PALANCA DE GASES  - [-]
            V_ini_cl = 20; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            V_fin_cl = 20; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            subida(1) = h_inicial_cl; % 1: ALTURA INICIAL - [m]
            subida(2) = h_final_cl; % 1: ALTURA FINAL - [m]
            subida(3) = gamma_cl; % 2: GAMMA DE SUBIDA  - [-]
            subida(4) = Mach_cl; % 3: MACH DE VUELO  - [-]
            subida(5) = TAS_cl; % 4: VELOCIDAD TAS  - [m/s]
            subida(6) = EAS_cl; % 5: VELOCIDAD EAS  - [m/s]
            subida(7) = delta_T_cl; % 6: PALANCA DE GASES  - [-]
            subida(8) = V_ini_cl; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            subida(9) = V_fin_cl; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            % Cruise
            h_inicial_cr = 5;% - [m] % Altura inicial
            dist_final_cr = 20*1000;% - [m] % 1: DISTANCIA FINAL
            V_cr = 25;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_cr);
            Mach_cr = V_cr/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            CL_cr = sqrt(Aero_TH.CD0/Aero_TH.CD2);% - [-]% 3: CL DE CRUCERO
            delta_T_cr = 0.9;% - []% 4: PALANCA DE GASES
            V_ini_cr = 25;% - [m/s] % 5: VELOCIDAD INICIAL
            V_fin_cr = 30;% - [m/s] % 6: VELOCIDAD FINAL
            fuel_cr = 1;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            Cd0_cr = -1;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            k1_cr = -1;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            k2_cr = -1;% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
            
            crucero(1) = h_inicial_cr;% - [m] % Altura inicial
            crucero(2) = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
            crucero(3) = V_cr; % Velocidad de crucero
            crucero(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            crucero(5) = CL_cr;% - [-]% 3: CL DE CRUCERO
            crucero(6) = delta_T_cr;% - []% 4: PALANCA DE GASES
            crucero(7) = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
            crucero(8) = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
            crucero(9) = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            crucero(10) = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            crucero(11) = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            crucero(12) = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            
            % Turn
            h_inicial_tr = 10*ft2m;% - [m]
            t_final_tr = 45*60;% 1: TIEMPO FINAL (seg)
            V_turn = 30; % turn velocity
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
            phi_tr = 25*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.05;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje(1) = h_inicial_tr;% - [m]
            viraje(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje(3) = V_turn; % turn velocity
            viraje(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje(6) = CL_tr;% 4: CL DE VIRAJE
            viraje(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje(9) = n_tr;% 7: FACTOR DE CARGA
            viraje(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Descent
            h_inicial_d = 10;
            h_final_d = 0;
            gamma_d = 3*D2R;
            V_d = 120;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_d,V_d);
            Mach_d = V_d/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            EAS_d = V_d;
            TAS_d = V_d;
            delta_T_d = 0.10;
            V_ini_d = V_d;
            V_fin_d = V_d*1.05;
            
            descenso(1) = h_inicial_d;
            descenso(2) = h_final_d;
            descenso(3) = gamma_d;
            descenso(4) = V_d;
            descenso(5) = Mach_d;
            descenso(6) = EAS_d;
            descenso(7) = TAS_d;
            descenso(8) = delta_T_d;
            descenso(9) = V_ini_d;
            descenso(10) = V_fin_d;
            
            % Turn waiting area
            h_inicial_tr = 10*ft2m;% - [m]
            t_final_tr = 45*60;% 1: TIEMPO FINAL (seg)
            V_turn = 90; % turn velocity
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
            phi_tr = 10*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.05;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje_wt(1) = h_inicial_tr;% - [m]
            viraje_wt(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje_wt(3) = V_turn; % turn velocity
            viraje_wt(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje_wt(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje_wt(6) = CL_tr;% 4: CL DE VIRAJE
            viraje_wt(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje_wt(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje_wt(9) = n_tr;% 7: FACTOR DE CARGA
            viraje_wt(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Dummy segment to account for 3 segments
            dummy(1) = h_inicial_cr;% - [m] % Altura inicial
            dummy(2) = 1;% - [m] % 1: DISTANCIA FINAL
            dummy(3) = V_cr; % Velocidad de crucero
            dummy(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            
            % Flight Modes
            MODES.cruise_mode = cruise_mode;
            MODES.climb_mode = climb_mode;
            MODES.turn_mode = turn_mode;
            MODES.descent_mode = descent_mode;
            
            Segments{1}.subida = subida;
            Segments{1}.crucero = crucero;
            Segments{1}.viraje = viraje;
            Segments{1}.viraje_wt = viraje_wt;
            Segments{1}.descenso = descenso;
            Segments{1}.dummy = dummy;
            Segments{1}.MODES = MODES;
                       
            Segments{1}.MODES = MODES;
            
        case 2 % case_AC = 2 - EMERGENTIA 1:2
            MODES.cruise_mode = cruise_mode;
            MODES.climb_mode = climb_mode;
            MODES.turn_mode = turn_mode;
            MODES.descent_mode = descent_mode;

            %             Segments{1}.taxy = taxy;
%             Segments{1}.despegue = despegue;
%             Segments{1}.subida = subida;
%             Segments{1}.crucero = crucero;
%             Segments{1}.viraje = viraje;
%             Segments{1}.viraje_wt = viraje_wt;
%             Segments{1}.descenso = descenso;
%             Segments{1}.aterrizaje = aterrizaje;
%             Segments{1}.dummy = dummy;
            Segments{1}.MODES = MODES;

        case 3 % case_AC = 3 - PEPIÑO XXL
            % Taxy
            temp_local_taxy = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_taxy = 0; % (m)
            P_local_taxy = 1.013268093075000e+05; % (Pa)
            delta_T_taxy = 0.05; % 4: PALANCA DE RALENTI EN TAXI = 0.05
            V_taxy = 10; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            t_taxy = 5*60; % 6: TIEMPO DE ESPERA EN TAXI (s)
            
            taxy(1) = temp_local_taxy;  % 1: TEMPERATURA LOCAL (Celsius)
            taxy(2) = h_inicial_taxy; % (m)
            taxy(3) = P_local_taxy; % (Pa)
            taxy(4) = delta_T_taxy; % 4: PALANCA DE RALENTI EN TAXI = 0.05
            taxy(5) = V_taxy; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            taxy(6) = t_taxy; % 6: TIEMPO DE ESPERA EN TAXI (s)
            
            % Takeoff
            temp_local_TO = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_TO = 0; % 2: ALTURA LOCAL (m)
            P_local_TO = 1.013268093075000e+05; % 3: PRESION LOCAL (Pa)
            mu_TO = 0.02; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            h_obstacle_TO = 10; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            gamma_climb_TO = 2*D2R; % 6: GAMMA DE SUBIDA MINIMO
            delta_T_TO = 1; % 7: PALANCA DE GASES PARA DESPEGUE
            
            despegue(1) = temp_local_TO; % 1: TEMPERATURA LOCAL (Celsius)
            despegue(2) = h_inicial_TO; % 2: ALTURA LOCAL (m)
            despegue(3) = P_local_TO; % 3: PRESION LOCAL (Pa)
            despegue(4) = mu_TO; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            despegue(5) = h_obstacle_TO; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            despegue(6) = gamma_climb_TO; % 6: GAMMA DE SUBIDA MINIMO
            despegue(7) = delta_T_TO; % 7: PALANCA DE GASES PARA DESPEGUE
            
            % Climb
            h_inicial_cl = 0; % 1: ALTURA INICIAL - [m]
            h_final_cl = 3000; % 1: ALTURA FINAL - [m]
            gamma_cl = 3*D2R; % 2: GAMMA DE SUBIDA  - [-]
            Mach_cl = 0.17; % 3: MACH DE VUELO  - [-]
            TAS_cl = 70; % 4: VELOCIDAD TAS  - [m/s]
            EAS_cl = 70; % 5: VELOCIDAD EAS  - [m/s]
            delta_T_cl = 0.95; % 6: PALANCA DE GASES  - [-]
            V_ini_cl = 70; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            V_fin_cl = 70; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            subida(1) = h_inicial_cl; % 1: ALTURA INICIAL - [m]
            subida(2) = h_final_cl; % 1: ALTURA FINAL - [m]
            subida(3) = gamma_cl; % 2: GAMMA DE SUBIDA  - [-]
            subida(4) = Mach_cl; % 3: MACH DE VUELO  - [-]
            subida(5) = TAS_cl; % 4: VELOCIDAD TAS  - [m/s]
            subida(6) = EAS_cl; % 5: VELOCIDAD EAS  - [m/s]
            subida(7) = delta_T_cl; % 6: PALANCA DE GASES  - [-]
            subida(8) = V_ini_cl; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            subida(9) = V_fin_cl; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            % Cruise
            h_inicial_cr = 3000;% - [m] % Altura inicial
            dist_final_cr = 200*1000;% - [m] % 1: DISTANCIA FINAL
            V_cr = 80;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_cr);
            Mach_cr = V_cr/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            CL_cr = sqrt(Aero_TH.CD0/Aero_TH.CD2);% - [-]% 3: CL DE CRUCERO
            delta_T_cr = 0.9;% - []% 4: PALANCA DE GASES
            V_ini_cr = 80;% - [m/s] % 5: VELOCIDAD INICIAL
            V_fin_cr = 80;% - [m/s] % 6: VELOCIDAD FINAL
            fuel_cr = 25;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            Cd0_cr = -1;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            k1_cr = -1;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            k2_cr = -1;% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
            
            crucero(1) = h_inicial_cr;% - [m] % Altura inicial
            crucero(2) = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
            crucero(3) = V_cr; % Velocidad de crucero
            crucero(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            crucero(5) = CL_cr;% - [-]% 3: CL DE CRUCERO
            crucero(6) = delta_T_cr;% - []% 4: PALANCA DE GASES
            crucero(7) = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
            crucero(8) = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
            crucero(9) = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            crucero(10) = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            crucero(11) = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            crucero(12) = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            
            % Turn
            h_inicial_tr = 3000;% - [m]
            t_final_tr = 45*60;% 1: TIEMPO FINAL (seg)
            V_turn = 70; % turn velocity
            [Data_ATM ,Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
            phi_tr = 30*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.05;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje(1) = h_inicial_tr;% - [m]
            viraje(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje(3) = V_turn; % turn velocity
            viraje(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje(6) = CL_tr;% 4: CL DE VIRAJE
            viraje(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje(9) = n_tr;% 7: FACTOR DE CARGA
            viraje(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Descent
            h_inicial_d = 3000;
            h_final_d = 0;
            gamma_d = 3*D2R;
            V_d = 60;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_d,V_d);
            Mach_d = V_d/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            EAS_d = V_d;
            TAS_d = V_d;
            delta_T_d = 0.10;
            V_ini_d = V_d;
            V_fin_d = V_d*1.05;
            
            descenso(1) = h_inicial_d;
            descenso(2) = h_final_d;
            descenso(3) = gamma_d;
            descenso(4) = V_d;
            descenso(5) = Mach_d;
            descenso(6) = EAS_d;
            descenso(7) = TAS_d;
            descenso(8) = delta_T_d;
            descenso(9) = V_ini_d;
            descenso(10) = V_fin_d;
            
            % Turn waiting area
            h_inicial_tr = 1500;% - [m]
            t_final_tr = 20*60;% 1: TIEMPO FINAL (seg)
            V_turn = 70; % turn velocity
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
            phi_tr = 10*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.05;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje_wt(1) = h_inicial_tr;% - [m]
            viraje_wt(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje_wt(3) = V_turn; % turn velocity
            viraje_wt(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje_wt(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje_wt(6) = CL_tr;% 4: CL DE VIRAJE
            viraje_wt(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje_wt(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje_wt(9) = n_tr;% 7: FACTOR DE CARGA
            viraje_wt(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Landing
            temp_local_LND = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_LND = 0; % 2: ALTURA LOCAL (m)
            P_local_LND = 1.013268093075000e+05; % 3: PRESION LOCAL (Pa)
            mu_LND = 0.2; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            delta_T_LND = 0.5; % PALANCA DE REVERSA
            t_brake = 4; % 'Tiempo en activar frenos' - [s]
            
            aterrizaje(1) = temp_local_LND;  % 1: TEMPERATURA LOCAL (Celsius)
            aterrizaje(2) = h_inicial_LND; % 2: ALTURA LOCAL (m)
            aterrizaje(3) = P_local_LND; % 3: PRESION LOCAL (Pa)
            aterrizaje(4) = mu_LND; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            aterrizaje(5) = delta_T_LND; % PALANCA DE REVERSA
            aterrizaje(6) = t_brake; % 'Tiempo en activar frenos' - [s]
            
            % Dummy segment to account for 3 segments
            dummy(1) = h_inicial_cr;% - [m] % Altura inicial
            dummy(2) = 1;% - [m] % 1: DISTANCIA FINAL
            dummy(3) = V_cr; % Velocidad de crucero
            dummy(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            
            % Flight Modes
            MODES.cruise_mode = cruise_mode;
            MODES.climb_mode = climb_mode;
            MODES.turn_mode = turn_mode;
            MODES.descent_mode = descent_mode;
            
            Segments{1}.taxy = taxy;
            Segments{1}.despegue = despegue;
            Segments{1}.subida = subida;
            Segments{1}.crucero = crucero;
            Segments{1}.viraje = viraje;
            Segments{1}.viraje_wt = viraje_wt;
            Segments{1}.descenso = descenso;
            Segments{1}.aterrizaje = aterrizaje;
            Segments{1}.dummy = dummy;
            Segments{1}.MODES = MODES;

        case 4 % Commercial example
            % Taxy
            temp_local_taxy = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_taxy = 0; % (m)
            P_local_taxy = 1.013268093075000e+05; % (Pa)
            delta_T_taxy = 0.05; % 4: PALANCA DE RALENTI EN TAXI = 0.05
            V_taxy = 10; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            t_taxy = 5*60; % 6: TIEMPO DE ESPERA EN TAXI (s)
            
            taxy(1) = temp_local_taxy;  % 1: TEMPERATURA LOCAL (Celsius)
            taxy(2) = h_inicial_taxy; % (m)
            taxy(3) = P_local_taxy; % (Pa)
            taxy(4) = delta_T_taxy; % 4: PALANCA DE RALENTI EN TAXI = 0.05
            taxy(5) = V_taxy; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            taxy(6) = t_taxy; % 6: TIEMPO DE ESPERA EN TAXI (s)
            
            % Takeoff
            temp_local_TO = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_TO = 0; % 2: ALTURA LOCAL (m)
            P_local_TO = 1.013268093075000e+05; % 3: PRESION LOCAL (Pa)
            mu_TO = 0.02; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            h_obstacle_TO = 10; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            gamma_climb_TO = 2*D2R; % 6: GAMMA DE SUBIDA MINIMO
            delta_T_TO = 1; % 7: PALANCA DE GASES PARA DESPEGUE
            
            despegue(1) = temp_local_TO; % 1: TEMPERATURA LOCAL (Celsius)
            despegue(2) = h_inicial_TO; % 2: ALTURA LOCAL (m)
            despegue(3) = P_local_TO; % 3: PRESION LOCAL (Pa)
            despegue(4) = mu_TO; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            despegue(5) = h_obstacle_TO; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            despegue(6) = gamma_climb_TO; % 6: GAMMA DE SUBIDA MINIMO
            despegue(7) = delta_T_TO; % 7: PALANCA DE GASES PARA DESPEGUE
            
            % Climb
            h_inicial_cl = 0; % 1: ALTURA INICIAL - [m]
            h_final_cl = 3000; % 1: ALTURA FINAL - [m]
            gamma_cl = 3*D2R; % 2: GAMMA DE SUBIDA  - [-]
            Mach_cl = 0.14; % 3: MACH DE VUELO  - [-]
            TAS_cl = 60; % 4: VELOCIDAD TAS  - [m/s]
            EAS_cl = 60; % 5: VELOCIDAD EAS  - [m/s]
            delta_T_cl = 0.95; % 6: PALANCA DE GASES  - [-]
            V_ini_cl = 60; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            V_fin_cl = 70; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            subida(1) = h_inicial_cl; % 1: ALTURA INICIAL - [m]
            subida(2) = h_final_cl; % 1: ALTURA FINAL - [m]
            subida(3) = gamma_cl; % 2: GAMMA DE SUBIDA  - [-]
            subida(4) = Mach_cl; % 3: MACH DE VUELO  - [-]
            subida(5) = TAS_cl; % 4: VELOCIDAD TAS  - [m/s]
            subida(6) = EAS_cl; % 5: VELOCIDAD EAS  - [m/s]
            subida(7) = delta_T_cl; % 6: PALANCA DE GASES  - [-]
            subida(8) = V_ini_cl; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            subida(9) = V_fin_cl; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            % Cruise
            h_inicial_cr = 3000;% - [m] % Altura inicial
            dist_final_cr = 250*1000;% - [m] % 1: DISTANCIA FINAL
            V_cr = 92.6;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_cr);
            Mach_cr = V_cr/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            CL_cr = 0.7;% - [-]% 3: CL DE CRUCERO
            delta_T_cr = 0.8;% - []% 4: PALANCA DE GASES
            V_ini_cr = 80;% - [m/s] % 5: VELOCIDAD INICIAL
            V_fin_cr = 90;% - [m/s] % 6: VELOCIDAD FINAL
            fuel_cr = 300;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            Cd0_cr = -1;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            k1_cr = -1;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            k2_cr = -1;% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
            
            crucero(1) = h_inicial_cr;% - [m] % Altura inicial
            crucero(2) = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
            crucero(3) = V_cr; % Velocidad de crucero
            crucero(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            crucero(5) = CL_cr;% - [-]% 3: CL DE CRUCERO
            crucero(6) = delta_T_cr;% - []% 4: PALANCA DE GASES
            crucero(7) = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
            crucero(8) = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
            crucero(9) = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            crucero(10) = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            crucero(11) = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            crucero(12) = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            
            % Turn
            h_inicial_tr = 3000*ft2m;% - [m]
            t_final_tr = 45*60;% 1: TIEMPO FINAL (seg)
            V_turn = 50; % turn velocity
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = 0.7;% 4: CL DE VIRAJE
            phi_tr = 30*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.20;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje(1) = h_inicial_tr;% - [m]
            viraje(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje(3) = V_turn; % turn velocity
            viraje(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje(6) = CL_tr;% 4: CL DE VIRAJE
            viraje(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje(9) = n_tr;% 7: FACTOR DE CARGA
            viraje(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Descent
            h_inicial_d = 3000;
            h_final_d = 10;
            gamma_d = 3*D2R;
            V_d = 40;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_d,V_d);
            Mach_d = V_d/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            EAS_d = V_d;
            TAS_d = V_d;
            delta_T_d = 0.10;
            V_ini_d = V_d;
            V_fin_d = V_d*1.05;
            
            descenso(1) = h_inicial_d;
            descenso(2) = h_final_d;
            descenso(3) = gamma_d;
            descenso(4) = V_d;
            descenso(5) = Mach_d;
            descenso(6) = EAS_d;
            descenso(7) = TAS_d;
            descenso(8) = delta_T_d;
            descenso(9) = V_ini_d;
            descenso(10) = V_fin_d;
            
            % Turn waiting area
            h_inicial_tr = 10*ft2m;% - [m]
            t_final_tr = 45*60;% 1: TIEMPO FINAL (seg)
            V_turn = 50; % turn velocity
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = 0.7;% 4: CL DE VIRAJE
            phi_tr = 30*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.20;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje_wt(1) = h_inicial_tr;% - [m]
            viraje_wt(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje_wt(3) = V_turn; % turn velocity
            viraje_wt(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje_wt(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje_wt(6) = CL_tr;% 4: CL DE VIRAJE
            viraje_wt(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje_wt(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje_wt(9) = n_tr;% 7: FACTOR DE CARGA
            viraje_wt(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Landing
            temp_local_LND = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_LND = 0; % 2: ALTURA LOCAL (m)
            P_local_LND = 1.013268093075000e+05; % 3: PRESION LOCAL (Pa)
            mu_LND = 0.2; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            delta_T_LND = 0.5; % PALANCA DE REVERSA
            t_brake = 4; % 'Tiempo en activar frenos' - [s]
            
            aterrizaje(1) = temp_local_LND;  % 1: TEMPERATURA LOCAL (Celsius)
            aterrizaje(2) = h_inicial_LND; % 2: ALTURA LOCAL (m)
            aterrizaje(3) = P_local_LND; % 3: PRESION LOCAL (Pa)
            aterrizaje(4) = mu_LND; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            aterrizaje(5) = delta_T_LND; % PALANCA DE REVERSA
            aterrizaje(6) = t_brake; % 'Tiempo en activar frenos' - [s]
            
            % Dummy segment to account for 3 segments
            dummy(1) = h_inicial_cr;% - [m] % Altura inicial
            dummy(2) = 1;% - [m] % 1: DISTANCIA FINAL
            dummy(3) = V_cr; % Velocidad de crucero
            dummy(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            
            % Flight Modes
            MODES.cruise_mode = cruise_mode;
            MODES.climb_mode = climb_mode;
            MODES.turn_mode = turn_mode;
            MODES.descent_mode = descent_mode;
            
            Segments{1}.taxy = taxy;
            Segments{1}.despegue = despegue;
            Segments{1}.subida = subida;
            Segments{1}.crucero = crucero;
            Segments{1}.viraje = viraje;
            Segments{1}.viraje_wt = viraje_wt;
            Segments{1}.descenso = descenso;
            Segments{1}.aterrizaje = aterrizaje;
            Segments{1}.dummy = dummy;
            Segments{1}.MODES = MODES;
            
        case 5 % WIG
            % Taxy
            temp_local_taxy = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_taxy = 0; % (m)
            P_local_taxy = 1.013268093075000e+05; % (Pa)
            delta_T_taxy = 0.05; % 4: PALANCA DE RALENTI EN TAXI = 0.05
            V_taxy = 10; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            t_taxy = 5*60; % 6: TIEMPO DE ESPERA EN TAXI (s)
            
            taxy(1) = temp_local_taxy;  % 1: TEMPERATURA LOCAL (Celsius)
            taxy(2) = h_inicial_taxy; % (m)
            taxy(3) = P_local_taxy; % (Pa)
            taxy(4) = delta_T_taxy; % 4: PALANCA DE RALENTI EN TAXI = 0.05
            taxy(5) = V_taxy; % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
            taxy(6) = t_taxy; % 6: TIEMPO DE ESPERA EN TAXI (s)
            
            % Takeoff
            temp_local_TO = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_TO = 0; % 2: ALTURA LOCAL (m)
            P_local_TO = 1.013268093075000e+05; % 3: PRESION LOCAL (Pa)
            mu_TO = 0.02; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            h_obstacle_TO = 10; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            gamma_climb_TO = 2*D2R; % 6: GAMMA DE SUBIDA MINIMO
            delta_T_TO = 1; % 7: PALANCA DE GASES PARA DESPEGUE
            
            despegue(1) = temp_local_TO; % 1: TEMPERATURA LOCAL (Celsius)
            despegue(2) = h_inicial_TO; % 2: ALTURA LOCAL (m)
            despegue(3) = P_local_TO; % 3: PRESION LOCAL (Pa)
            despegue(4) = mu_TO; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            despegue(5) = h_obstacle_TO; % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
            despegue(6) = gamma_climb_TO; % 6: GAMMA DE SUBIDA MINIMO
            despegue(7) = delta_T_TO; % 7: PALANCA DE GASES PARA DESPEGUE
            
            % Climb
            h_inicial_cl = 0; % 1: ALTURA INICIAL - [m]
            h_final_cl = 5; % 1: ALTURA FINAL - [m]
            gamma_cl = 3*D2R; % 2: GAMMA DE SUBIDA  - [-]
            Mach_cl = 0.14; % 3: MACH DE VUELO  - [-]
            TAS_cl = 60; % 4: VELOCIDAD TAS  - [m/s]
            EAS_cl = 80; % 5: VELOCIDAD EAS  - [m/s]
            delta_T_cl = 0.95; % 6: PALANCA DE GASES  - [-]
            V_ini_cl = 125; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            V_fin_cl = 125; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            subida(1) = h_inicial_cl; % 1: ALTURA INICIAL - [m]
            subida(2) = h_final_cl; % 1: ALTURA FINAL - [m]
            subida(3) = gamma_cl; % 2: GAMMA DE SUBIDA  - [-]
            subida(4) = Mach_cl; % 3: MACH DE VUELO  - [-]
            subida(5) = TAS_cl; % 4: VELOCIDAD TAS  - [m/s]
            subida(6) = EAS_cl; % 5: VELOCIDAD EAS  - [m/s]
            subida(7) = delta_T_cl; % 6: PALANCA DE GASES  - [-]
            subida(8) = V_ini_cl; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            subida(9) = V_fin_cl; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            % Cruise
            h_inicial_cr = 5;% - [m] % Altura inicial
            dist_final_cr = 2000*1000;% - [m] % 1: DISTANCIA FINAL
            V_cr = 90;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_cr);
            Mach_cr = V_cr/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            CL_cr = sqrt(Aero_TH.CD0/Aero_TH.CD2);% - [-]% 3: CL DE CRUCERO
            delta_T_cr = 0.9;% - []% 4: PALANCA DE GASES
            V_ini_cr = 90;% - [m/s] % 5: VELOCIDAD INICIAL
            V_fin_cr = 125;% - [m/s] % 6: VELOCIDAD FINAL
            fuel_cr = 13500;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            Cd0_cr = -1;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            k1_cr = -1;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            k2_cr = -1;% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
            
            crucero(1) = h_inicial_cr;% - [m] % Altura inicial
            crucero(2) = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
            crucero(3) = V_cr; % Velocidad de crucero
            crucero(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            crucero(5) = CL_cr;% - [-]% 3: CL DE CRUCERO
            crucero(6) = delta_T_cr;% - []% 4: PALANCA DE GASES
            crucero(7) = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
            crucero(8) = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
            crucero(9) = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            crucero(10) = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            crucero(11) = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            crucero(12) = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            
            % Turn
            h_inicial_tr = 10*ft2m;% - [m]
            t_final_tr = 45*60;% 1: TIEMPO FINAL (seg)
            V_turn = 90; % turn velocity
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
            phi_tr = 30*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.05;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje(1) = h_inicial_tr;% - [m]
            viraje(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje(3) = V_turn; % turn velocity
            viraje(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje(6) = CL_tr;% 4: CL DE VIRAJE
            viraje(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje(9) = n_tr;% 7: FACTOR DE CARGA
            viraje(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Descent
            h_inicial_d = 10;
            h_final_d = 0;
            gamma_d = 3*D2R;
            V_d = 120;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_d,V_d);
            Mach_d = V_d/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            EAS_d = V_d;
            TAS_d = V_d;
            delta_T_d = 0.10;
            V_ini_d = V_d;
            V_fin_d = V_d*1.05;
            
            descenso(1) = h_inicial_d;
            descenso(2) = h_final_d;
            descenso(3) = gamma_d;
            descenso(4) = V_d;
            descenso(5) = Mach_d;
            descenso(6) = EAS_d;
            descenso(7) = TAS_d;
            descenso(8) = delta_T_d;
            descenso(9) = V_ini_d;
            descenso(10) = V_fin_d;
            
            % Turn waiting area
            h_inicial_tr = 10*ft2m;% - [m]
            t_final_tr = 45*60;% 1: TIEMPO FINAL (seg)
            V_turn = 90; % turn velocity
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
            phi_tr = 10*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.05;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje_wt(1) = h_inicial_tr;% - [m]
            viraje_wt(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje_wt(3) = V_turn; % turn velocity
            viraje_wt(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje_wt(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje_wt(6) = CL_tr;% 4: CL DE VIRAJE
            viraje_wt(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje_wt(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje_wt(9) = n_tr;% 7: FACTOR DE CARGA
            viraje_wt(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Landing
            temp_local_LND = 15;  % 1: TEMPERATURA LOCAL (Celsius)
            h_inicial_LND = 0; % 2: ALTURA LOCAL (m)
            P_local_LND = 1.013268093075000e+05; % 3: PRESION LOCAL (Pa)
            mu_LND = 0.2; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            delta_T_LND = 0.5; % PALANCA DE REVERSA
            t_brake = 4; % 'Tiempo en activar frenos' - [s]
            
            aterrizaje(1) = temp_local_LND;  % 1: TEMPERATURA LOCAL (Celsius)
            aterrizaje(2) = h_inicial_LND; % 2: ALTURA LOCAL (m)
            aterrizaje(3) = P_local_LND; % 3: PRESION LOCAL (Pa)
            aterrizaje(4) = mu_LND; % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
            aterrizaje(5) = delta_T_LND; % PALANCA DE REVERSA
            aterrizaje(6) = t_brake; % 'Tiempo en activar frenos' - [s]
            
            % Dummy segment to account for 3 segments
            dummy(1) = h_inicial_cr;% - [m] % Altura inicial
            dummy(2) = 1;% - [m] % 1: DISTANCIA FINAL
            dummy(3) = V_cr; % Velocidad de crucero
            dummy(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            
            % Flight Modes
            MODES.cruise_mode = cruise_mode;
            MODES.climb_mode = climb_mode;
            MODES.turn_mode = turn_mode;
            MODES.descent_mode = descent_mode;
            
            Segments{1}.taxy = taxy;
            Segments{1}.despegue = despegue;
            Segments{1}.subida = subida;
            Segments{1}.crucero = crucero;
            Segments{1}.viraje = viraje;
            Segments{1}.viraje_wt = viraje_wt;
            Segments{1}.descenso = descenso;
            Segments{1}.aterrizaje = aterrizaje;
            Segments{1}.dummy = dummy;
            Segments{1}.MODES = MODES;
            
        case 6 % case_AC = 6 - CERVERA
            
            % Climb
            h_inicial_cl = 0; % 1: ALTURA INICIAL - [m]
            h_final_cl = 5; % 1: ALTURA FINAL - [m]
            gamma_cl = 3*D2R; % 2: GAMMA DE SUBIDA  - [-]
            Mach_cl = 0.14; % 3: MACH DE VUELO  - [-]
            TAS_cl = 20; % 4: VELOCIDAD TAS  - [m/s]
            EAS_cl = 20; % 5: VELOCIDAD EAS  - [m/s]
            delta_T_cl = 0.95; % 6: PALANCA DE GASES  - [-]
            V_ini_cl = 20; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            V_fin_cl = 20; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            subida(1) = h_inicial_cl; % 1: ALTURA INICIAL - [m]
            subida(2) = h_final_cl; % 1: ALTURA FINAL - [m]
            subida(3) = gamma_cl; % 2: GAMMA DE SUBIDA  - [-]
            subida(4) = Mach_cl; % 3: MACH DE VUELO  - [-]
            subida(5) = TAS_cl; % 4: VELOCIDAD TAS  - [m/s]
            subida(6) = EAS_cl; % 5: VELOCIDAD EAS  - [m/s]
            subida(7) = delta_T_cl; % 6: PALANCA DE GASES  - [-]
            subida(8) = V_ini_cl; % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)  - [m/s]
            subida(9) = V_fin_cl; % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)  - [m/s]
            
            % Cruise
            h_inicial_cr = 5;% - [m] % Altura inicial
            dist_final_cr = 20*1000;% - [m] % 1: DISTANCIA FINAL
            V_cr = 25;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_cr);
            Mach_cr = V_cr/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            CL_cr = sqrt(Aero_TH.CD0/Aero_TH.CD2);% - [-]% 3: CL DE CRUCERO
            delta_T_cr = 0.9;% - []% 4: PALANCA DE GASES
            V_ini_cr = 25;% - [m/s] % 5: VELOCIDAD INICIAL
            V_fin_cr = 30;% - [m/s] % 6: VELOCIDAD FINAL
            fuel_cr = 1;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            Cd0_cr = -1;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            k1_cr = -1;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            k2_cr = -1;% - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
            
            crucero(1) = h_inicial_cr;% - [m] % Altura inicial
            crucero(2) = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
            crucero(3) = V_cr; % Velocidad de crucero
            crucero(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            crucero(5) = CL_cr;% - [-]% 3: CL DE CRUCERO
            crucero(6) = delta_T_cr;% - []% 4: PALANCA DE GASES
            crucero(7) = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
            crucero(8) = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
            crucero(9) = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            crucero(10) = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            crucero(11) = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            crucero(12) = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            
            % Turn
            h_inicial_tr = 10*ft2m;% - [m]
            t_final_tr = 45*60;% 1: TIEMPO FINAL (seg)
            V_turn = 30; % turn velocity
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
            phi_tr = 25*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.05;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje(1) = h_inicial_tr;% - [m]
            viraje(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje(3) = V_turn; % turn velocity
            viraje(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje(6) = CL_tr;% 4: CL DE VIRAJE
            viraje(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje(9) = n_tr;% 7: FACTOR DE CARGA
            viraje(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Descent
            h_inicial_d = 10;
            h_final_d = 0;
            gamma_d = 3*D2R;
            V_d = 120;
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_d,V_d);
            Mach_d = V_d/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            EAS_d = V_d;
            TAS_d = V_d;
            delta_T_d = 0.10;
            V_ini_d = V_d;
            V_fin_d = V_d*1.05;
            
            descenso(1) = h_inicial_d;
            descenso(2) = h_final_d;
            descenso(3) = gamma_d;
            descenso(4) = V_d;
            descenso(5) = Mach_d;
            descenso(6) = EAS_d;
            descenso(7) = TAS_d;
            descenso(8) = delta_T_d;
            descenso(9) = V_ini_d;
            descenso(10) = V_fin_d;
            
            % Turn waiting area
            h_inicial_tr = 10*ft2m;% - [m]
            t_final_tr = 45*60;% 1: TIEMPO FINAL (seg)
            V_turn = 90; % turn velocity
            [Data_ATM Performance] = Flight_Conditions_2020_v1(h_inicial_cr,V_turn);
            Mach_tr = V_turn/Data_ATM.a;% - [-] % 2: MACH DE VUELO
            delta_T_tr = 1; % 3: PALANCA DE GASES
            CL_tr = sqrt(3*Aero_TH.CD0/Aero_TH.CD2);% 4: CL DE VIRAJE
            phi_tr = 10*D2R;% 5: ANGULO DE ALABEO (rads)
            V_psi = 10*D2R;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            n_tr = 1.05;% 7: FACTOR DE CARGA
            R_tr = 300;% 8: RADIO DE GIRO (m)
            
            viraje_wt(1) = h_inicial_tr;% - [m]
            viraje_wt(2) = t_final_tr;% 1: TIEMPO FINAL (seg)
            viraje_wt(3) = V_turn; % turn velocity
            viraje_wt(4) = Mach_tr;% - [-] % 2: MACH DE VUELO
            viraje_wt(5) = delta_T_tr; % 3: PALANCA DE GASES
            viraje_wt(6) = CL_tr;% 4: CL DE VIRAJE
            viraje_wt(7) = phi_tr;% 5: ANGULO DE ALABEO (rads)
            viraje_wt(8) = V_psi;% 6: VELOCIDAD DE GUIÑADA (rads/seg)
            viraje_wt(9) = n_tr;% 7: FACTOR DE CARGA
            viraje_wt(10) = R_tr;% 8: RADIO DE GIRO (m)
            
            % Dummy segment to account for 3 segments
            dummy(1) = h_inicial_cr;% - [m] % Altura inicial
            dummy(2) = 1;% - [m] % 1: DISTANCIA FINAL
            dummy(3) = V_cr; % Velocidad de crucero
            dummy(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            
            % Flight Modes
            MODES.cruise_mode = cruise_mode;
            MODES.climb_mode = climb_mode;
            MODES.turn_mode = turn_mode;
            MODES.descent_mode = descent_mode;
            
            Segments{1}.subida = subida;
            Segments{1}.crucero = crucero;
            Segments{1}.viraje = viraje;
            Segments{1}.viraje_wt = viraje_wt;
            Segments{1}.descenso = descenso;
            Segments{1}.dummy = dummy;
            Segments{1}.MODES = MODES;
                       
            Segments{1}.MODES = MODES;
  
    end
