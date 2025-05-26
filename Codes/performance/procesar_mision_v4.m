function [Weights_AP,Total_Datos,datos] = procesar_mision_v4(seg,tramos,handles,FlightMODE,Geo_tier,Prop_data,AC_CONFIGURATION)

propul = handles.propul;
data_electric = handles.electric_propul;
aerodinamica = handles.aerodinamica;
aero_despegue = handles.aero_despegue;
aero_aterrizaje = handles.aero_aterrizaje;
pesos = handles.pesos;

% cruise_mode = Segments{1}.MODES.cruise_mode;
% climb_mode = Segments{1}.MODES.climb_mode;
% turn_mode = Segments{1}.MODES.turn_mode;
% descent_mode = Segments{1}.MODES.descent_mode;

g = 9.80665;
fuel_inicial = 0;
mbaterias_inicial = 0;
contador = 0;
saltar = 0;
datos.inicializar = 0;
W(1) = 0;

while contador < 50
    clear W
    
    W(1) = (pesos(1) + pesos(2) + pesos(3) + fuel_inicial + mbaterias_inicial)*g;
    for i=1:tramos
%         climb_mode = FlightMODE.climb_mode(i);
%         turn_mode = FlightMODE.turn_mode(i);
%         descent_mode = FlightMODE.descent_mode(i);
%         cruise_mode = FlightMODE.cruise_mode(i);
        
        %% Case 1: TAXY
        if strcmp(seg(i).nombre,'Taxy') == 1
            [fuel(i),tiempo(i),distancia(i)] = analisis_taxi(seg(i).datos,propul,i,Prop_data,data_electric); %kg
            datos(i).segmento.fuel = fuel(i);
            datos(i).segmento.tiempo = tiempo(i);
            datos(i).segmento.distancia = distancia(i);
            datos(i).nombre = 'Taxy';
%             datos(i).lista_variables = [{''};{'Fuel'}];
        end
        
        %% Case 2: TAKE OFF
        if strcmp(seg(i).nombre,'TakeOff') == 1
            [fuel(i),tiempo(i),distancia(i),datos] = analisis_despegue(seg(i).datos,propul,aerodinamica,aero_despegue,W(i),datos,i,Prop_data,data_electric);
            datos(i).segmento.fuel = fuel(i);
            datos(i).nombre = 'TakeOff';
%             datos(i).lista_variables = [{''};{'Fuel'}];
        end
        
        %% Case 3: CLIMB
        if strcmp(seg(i).nombre,'Climb') == 1
            h_inicial = seg(i).datos.h_inicial;
            climb_mode    = seg(i).opcion;
            subida(1) = seg(i).datos.h_final;
            %             subida(2) = seg(i).datos.gamma;
            %             subida(3) = seg(i).datos.Mach;
            %             subida(4) = seg(i).datos.EAS;
            %             subida(5) = seg(i).datos.TAS;
            %             subida(6) = seg(i).datos.palanca;
            %             subida(7) = seg(i).datos.V_ini;
            %             subida(8) = seg(i).datos.V_fin;
            
            switch climb_mode
                case 1
                    subida(2) = seg(i).datos.gamma;
                    subida(3) = seg(i).datos.Mach;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = -1;
                    subida(7) = -1;
                    subida(8) = -1;
                case 2
                    subida(2) = seg(i).datos.gamma;
                    subida(3) = -1;
                    subida(4) = seg(i).datos.EAS;
                    subida(5) = -1;
                    subida(6) = -1;
                    subida(7) = -1;
                    subida(8) = -1;
                case 3
                    subida(2) = seg(i).datos.gamma;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = seg(i).datos.TAS;
                    subida(6) = -1;
                    subida(7) = -1;
                    subida(8) = -1;
                case 4
                    subida(2) = -1;
                    subida(3) = seg(i).datos.Mach;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 5
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = seg(i).datos.EAS;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 6
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = seg(i).datos.TAS;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 7
                    subida(2) = seg(i).datos.gamma;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = -1;
                    subida(7) = seg(i).datos.V_ini;
                    subida(8) = seg(i).datos.V_fin;
                case 8
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 9
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 10
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = seg(i).datos.V_ini;
                    subida(8) = seg(i).datos.V_fin;
            end
            [fuel(i),mbaterias(i),energia(i),tiempo(i),distancia(i),datos] = analisis_subida(propul,aerodinamica,subida,W(i),h_inicial,climb_mode,datos,i,Prop_data,data_electric);
            datos(i).nombre = 'Climb';
        end
        
        %% Case 4: VTOL CLIMB
        if strcmp(seg(i).nombre,'VTOL_Climb') == 1
            %             carga_soltada = seg(i).datos.carga;
            %             fuel(i) = 0;
            %             datos(i).soltar_carga.carga_soltada = carga_soltada;
            %             datos(i).soltar_carga.fuel = 0;
            %             datos(i).nombre = 'Soltar carga';
            %             datos(i).lista_variables = [{''};{'Carga soltada'}];
            %             W(i+1) = W(i) - carga_soltada * g;
            %             saltar = 1;
            vtolclimb_mode = seg(i).opcion;
            h_inicial = seg(i).datos.h_inicial;
            vtolclimb(1) = h_inicial;
            switch vtolclimb_mode
                
                case 1
                    vtolclimb(2) = seg(i).datos.h_final;
                    vtolclimb(3) = -1;
                    vtolclimb(4) = seg(i).datos.palanca;
                    vtolclimb(5) = -1;
                    vtolclimb(6) = -1;
                case 2
                    vtolclimb(2) = -1;
                    vtolclimb(3) = seg(i).datos.thover;
                    vtolclimb(4) = -1;
                    vtolclimb(5) = -1;
                    vtolclimb(6) = -1;

                case 3 
                    vtolclimb(2) = seg(i).datos.h_final;
                    vtolclimb(3) = -1;
                    vtolclimb(4) = -1;
                    vtolclimb(5) = seg(i).datos.vclimb;
                    vtolclimb(6) = -1;

                case 4
                    vtolclimb(2) = -1;
                    vtolclimb(3) = -1;
                    vtolclimb(4) = -1;
                    vtolclimb(5) = -1;
                    vtolclimb(6) = seg(i).datos.mbathover;
            end
            
           [fuel(i),mbaterias(i),energia(i),tiempo(i),distancia(i),datos] = analisis_vtolclimb(propul,aerodinamica,vtolclimb,W(i),h_inicial,vtolclimb_mode,datos,i,handles,Prop_data,data_electric);
           datos(i).nombre = 'VTOL Climb';
        end
        
        %% Case 5: CRUISE
        if strcmp(seg(i).nombre,'Cruise') == 1
            h_inicial = seg(i).datos.h_inicial;
            cruise_mode    = seg(i).opcion;
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
            
            % crucero(1) = h_inicial_cr;% - [m] % Altura inicial
            % crucero(2) = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
            % crucero(3) = V_cr; % Velocidad de crucero
            % crucero(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            % crucero(5) = CL_cr;% - [-]% 3: CL DE CRUCERO
            % crucero(6) = delta_T_cr;% - []% 4: PALANCA DE GASES
            % crucero(7) = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
            % crucero(8) = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
            % crucero(9) = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            % crucero(10) = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            % crucero(11) = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            % crucero(12) = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            
            % variables 11 y 12 obsoletas. Se incorpora la masa de baterías
            % disponible para el case 7
            
            % crucero(11) = mbat_disp % - [kg] % 11: MASA DE BATERÍAS DISPONIBLE
            
            switch cruise_mode
                case 1
                    crucero(1) = seg(i).datos.dist_final;
                    crucero(2) = seg(i).datos.Mach;
                    crucero(3) = -1;
                    crucero(4) = -1;
                    crucero(5) = -1;
                    crucero(6) = -1;
                    crucero(7) = -1;
                    crucero(8) = -1;
                    crucero(9) = -1;
                    crucero(10) = -1;
                    crucero(11) = -1;
                case 2
                    crucero(1) = seg(i).datos.dist_final;
                    crucero(2) = -1;
                    crucero(3) = seg(i).datos.CL;
                    crucero(4) = -1;
                    crucero(5) = -1;
                    crucero(6) = -1;
                    crucero(7) = -1;
                    crucero(8) = -1;
                    crucero(9) = -1;
                    crucero(10) = -1;
                    crucero(11) = -1;
                case 3
                    crucero(1) = seg(i).datos.dist_final;
                    crucero(2) = -1;
                    crucero(3) = -1;
                    crucero(4) = seg(i).datos.palanca;
                    crucero(5) = seg(i).datos.V_ini;
                    crucero(6) = seg(i).datos.V_fin;
                    crucero(7) = -1;
                    crucero(8) = -1;
                    crucero(9) = -1;
                    crucero(10) = -1;
                    crucero(11) = -1;
                case 4
                    crucero(1) = seg(i).datos.dist_final;
                    crucero(2) = seg(i).datos.Mach;
                    crucero(3) = -1;
                    crucero(4) = -1;
                    crucero(5) = -1;
                    crucero(6) = -1;
                    crucero(7) = -1;
                    crucero(8) = seg(i).datos.Cd0;
                    crucero(9) = seg(i).datos.k1;
                    crucero(10) = seg(i).datos.k2;
                    crucero(11) = -1;
                case 5
                    crucero(1) = -1;
                    crucero(2) = -1;
                    crucero(3) = -1;
                    crucero(4) = -1;
                    crucero(5) = -1;
                    crucero(6) = -1;
                    crucero(7) = seg(i).datos.fuel;
                    crucero(8) = -1;
                    crucero(9) = -1;
                    crucero(10) = -1;
                    crucero(11) = -1;
                case 6
                    crucero(1) = -1;
                    crucero(2) = -1;
                    crucero(3) = -1;
                    crucero(4) = -1;
                    crucero(5) = -1;
                    crucero(6) = -1;
                    crucero(7) = seg(i).datos.fuel;
                    crucero(8) = -1;
                    crucero(9) = -1;
                    crucero(10) = -1;
                    crucero(11) = -1;
                case 7
                    crucero(1) = -1;
                    crucero(2) = seg(i).datos.Mach;
                    crucero(3) = -1;
                    crucero(4) = -1;
                    crucero(5) = -1;
                    crucero(6) = -1;
                    crucero(7) = -1;
                    crucero(8) = -1;
                    crucero(9) = -1;
                    crucero(10) = -1;
                    crucero(11) = seg(i).datos.m_bat;
            end
%             [fuel(i),tiempo(i),distancia(i),datos] = analisis_crucero_v2(propul,aerodinamica,crucero,W(i),h_inicial,cruise_mode,i,datos,Geo_tier,Prop_data,AC_CONFIGURATION);
            [fuel(i),mbaterias(i),energia(i),tiempo(i),distancia(i),datos] = analisis_crucero_v3(propul,aerodinamica,crucero,W(i),h_inicial,cruise_mode,datos,i,Geo_tier,AC_CONFIGURATION,Prop_data,data_electric,pesos);
            datos(i).nombre = 'Cruise';
        end
        
        %% Case 6: LOAD DEPLOYMENT
        if strcmp(seg(i).nombre,'Load_Deployment') == 1
            carga_soltada = seg(i).datos.carga;
            fuel(i) = 0;
            tiempo(i) = 0;
            distancia(i) = 0;
            mbaterias(i) = 0;
            datos(i).segmento.tiempo = 0;
            datos(i).segmento.distancia = 0;
            datos(i).segmento.carga_soltada = carga_soltada;
            datos(i).segmento.fuel = 0;
            datos(i).segmento.mbaterias = 0;
            datos(i).segmento.energiav = 0;
            datos(i).nombre = 'Load Deployment';
%             datos(i).lista_variables = [{''};{'Load Deployment'}];
            W(i+1) = W(i) - carga_soltada * g;
            saltar = 1;
        end
        
        %% Case 7: TURN
        if strcmp(seg(i).nombre,'Turn') == 1
            % VIRAJE
            % 1: TIEMPO FINAL
            % 2: MACH DE VUELO
            % 3: PALANCA DE GASES
            % 4: CL DE VIRAJE
            % 5: ANGULO DE ALABEO
            % 6: VELOCIDAD DE GUIÑADA
            % 7: FACTOR DE CARGA
            % 8: RADIO DE GIRO
            h_inicial = seg(i).datos.h_inicial;
            viraje(1) = seg(i).datos.t_final;
            turn_mode    = seg(i).opcion;
            switch turn_mode
                case 1
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = seg(i).datos.palanca;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 2
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = seg(i).datos.CL;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 3
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = -1;
                    viraje(5) = seg(i).datos.balance;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 4
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = seg(i).datos.n;
                    viraje(8) = -1;
                case 5
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = seg(i).datos.radio;
                case 6
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = seg(i).datos.V_psi;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 7
                    viraje(2) = -1;
                    viraje(3) = seg(i).datos.palanca;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 8
                    viraje(2) = -1;
                    viraje(3) = seg(i).datos.palanca;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 9
                    viraje(2) = -1;
                    viraje(3) = seg(i).datos.palanca;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
            end
            [fuel(i),tiempo(i),distancia(i),datos] = analisis_viraje(propul,aerodinamica,viraje,W(i),h_inicial,turn_mode,datos,i,Prop_data,data_electric);
            tiempo(i) = viraje(1);
            datos(i).nombre = 'Turn';
        end
        
        %% Case 8: DESCENT
        if strcmp(seg(i).nombre,'Descent') == 1
            h_inicial = seg(i).datos.h_inicial;
            descent_mode    = seg(i).opcion;
            descenso(1) = seg(i).datos.h_final;
            
            switch descent_mode
                case 1
                    descenso(2) = seg(i).datos.gamma;
                    descenso(3) = seg(i).datos.Mach;
                    descenso(4) = -1;
                    descenso(5) = -1;
                    descenso(6) = -1;
                    descenso(7) = -1;
                    descenso(8) = -1;
                case 2
                    descenso(2) = seg(i).datos.gamma;
                    descenso(3) = -1;
                    descenso(4) = seg(i).datos.EAS;
                    descenso(5) = -1;
                    descenso(6) = -1;
                    descenso(7) = -1;
                    descenso(8) = -1;
                case 3
                    descenso(2) = seg(i).datos.gamma;
                    descenso(3) = -1;
                    descenso(4) = -1;
                    descenso(5) = seg(i).datos.TAS;
                    descenso(6) = -1;
                    descenso(7) = -1;
                    descenso(8) = -1;
                case 4
                    descenso(2) = -1;
                    descenso(3) = seg(i).datos.Mach;
                    descenso(4) = -1;
                    descenso(5) = -1;
                    descenso(6) = seg(i).datos.palanca;
                    descenso(7) = -1;
                    descenso(8) = -1;
                case 5
                    descenso(2) = -1;
                    descenso(3) = -1;
                    descenso(4) = seg(i).datos.EAS;
                    descenso(5) = -1;
                    descenso(6) = seg(i).datos.palanca;
                    descenso(7) = -1;
                    descenso(8) = -1;
                case 6
                    descenso(2) = -1;
                    descenso(3) = -1;
                    descenso(4) = -1;
                    descenso(5) = seg(i).datos.TAS;
                    descenso(6) = seg(i).datos.palanca;
                    descenso(7) = -1;
                    descenso(8) = -1;
                case 7
                    descenso(2) = seg(i).datos.gamma;
                    descenso(3) = -1;
                    descenso(4) = -1;
                    descenso(5) = -1;
                    descenso(6) = -1;
                    descenso(7) = seg(i).datos.V_ini;
                    descenso(8) = seg(i).datos.V_fin;
                case 8
                    descenso(2) = -1;
                    descenso(3) = -1;
                    descenso(4) = -1;
                    descenso(5) = -1;
                    descenso(6) = seg(i).datos.palanca;
                    descenso(7) = -1;
                    descenso(8) = -1;
                case 9
                    descenso(2) = -1;
                    descenso(3) = -1;
                    descenso(4) = -1;
                    descenso(5) = -1;
                    descenso(6) = seg(i).datos.palanca;
                    descenso(7) = -1;
                    descenso(8) = -1;
                case 10
                    descenso(2) = -1;
                    descenso(3) = -1;
                    descenso(4) = -1;
                    descenso(5) = -1;
                    descenso(6) = seg(i).datos.palanca;
                    descenso(7) = seg(i).datos.V_ini;
                    descenso(8) = seg(i).datos.V_fin;
            end
            [fuel(i),mbaterias(i),energia(i),tiempo(i),distancia(i),datos] = analisis_descenso(propul,aerodinamica,descenso,W(i),h_inicial,descent_mode,datos,i,Prop_data,data_electric);
            datos(i).nombre = 'Descent';
        end
        
        %% Case 9: VTOL DESCENT
        if strcmp(seg(i).nombre,'Descent_VTOL') == 1
            vtoldes_mode = seg(i).opcion;
            h_inicial    = seg(i).datos.h_inicial;
            vtoldes(1)   = seg(i).datos.h_inicial;

            switch vtoldes_mode
                case 1
                    vtoldes(2) = seg(i).datos.h_final;
                    vtoldes(3) = seg(i).datos.vdes;
            end
           [fuel(i),mbaterias(i),energia(i),tiempo(i),distancia(i),datos] = analisis_vtoldes(propul,aerodinamica,vtoldes,W(i),h_inicial,vtoldes_mode,datos,i,Prop_data,data_electric);
           datos(i).nombre = 'VTOL Descent';
        end
        
        %% Case 10: Climb Waiting Area to 3000ft
        if strcmp(seg(i).nombre,'Climb_waiting') == 1
            h_inicial = seg(i).datos.h_inicial;
            climb_mode    = seg(i).opcion;
            subida(1) = seg(i).datos.h_final;
            %             subida(2) = seg(i).datos.gamma;
            %             subida(3) = seg(i).datos.Mach;
            %             subida(4) = seg(i).datos.EAS;
            %             subida(5) = seg(i).datos.TAS;
            %             subida(6) = seg(i).datos.palanca;
            %             subida(7) = seg(i).datos.V_ini;
            %             subida(8) = seg(i).datos.V_fin;
            
            switch climb_mode
                case 1
                    subida(2) = seg(i).datos.gamma;
                    subida(3) = seg(i).datos.Mach;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = -1;
                    subida(7) = -1;
                    subida(8) = -1;
                case 2
                    subida(2) = seg(i).datos.gamma;
                    subida(3) = -1;
                    subida(4) = seg(i).datos.EAS;
                    subida(5) = -1;
                    subida(6) = -1;
                    subida(7) = -1;
                    subida(8) = -1;
                case 3
                    subida(2) = seg(i).datos.gamma;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = seg(i).datos.TAS;
                    subida(6) = -1;
                    subida(7) = -1;
                    subida(8) = -1;
                case 4
                    subida(2) = -1;
                    subida(3) = seg(i).datos.Mach;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 5
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = seg(i).datos.EAS;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 6
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = seg(i).datos.TAS;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 7
                    subida(2) = seg(i).datos.gamma;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = -1;
                    subida(7) = seg(i).datos.V_ini;
                    subida(8) = seg(i).datos.V_fin;
                case 8
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 9
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = -1;
                    subida(8) = -1;
                case 10
                    subida(2) = -1;
                    subida(3) = -1;
                    subida(4) = -1;
                    subida(5) = -1;
                    subida(6) = seg(i).datos.palanca;
                    subida(7) = seg(i).datos.V_ini;
                    subida(8) = seg(i).datos.V_fin;
            end
            [fuel(i),tiempo(i),distancia(i),datos] = analisis_subida(propul,aerodinamica,subida,W(i),h_inicial,climb_mode,datos,i,Prop_data,data_electric);
            datos(i).nombre = 'Climb Waiting';
        end
        
        %% Case 11: Turn Waiting Area to 3000ft
        if strcmp(seg(i).nombre,'Turn_waiting') == 1
            % VIRAJE
            % 1: TIEMPO FINAL
            % 2: MACH DE VUELO
            % 3: PALANCA DE GASES
            % 4: CL DE VIRAJE
            % 5: ANGULO DE ALABEO
            % 6: VELOCIDAD DE GUIÑADA
            % 7: FACTOR DE CARGA
            % 8: RADIO DE GIRO
            h_inicial = seg(i).datos.h_inicial;
            viraje(1) = seg(i).datos.t_final;
            turn_mode    = seg(i).opcion;
            switch turn_mode
                case 1
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = seg(i).datos.palanca;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 2
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = seg(i).datos.CL;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 3
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = -1;
                    viraje(5) = seg(i).datos.balance;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 4
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = seg(i).datos.n;
                    viraje(8) = -1;
                case 5
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = seg(i).datos.radio;
                case 6
                    viraje(2) = seg(i).datos.velocidad; %Mach;
                    viraje(3) = -1;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = seg(i).datos.V_psi;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 7
                    viraje(2) = -1;
                    viraje(3) = seg(i).datos.palanca;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 8
                    viraje(2) = -1;
                    viraje(3) = seg(i).datos.palanca;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
                case 9
                    viraje(2) = -1;
                    viraje(3) = seg(i).datos.palanca;
                    viraje(4) = -1;
                    viraje(5) = -1;
                    viraje(6) = -1;
                    viraje(7) = -1;
                    viraje(8) = -1;
            end
            [fuel(i),tiempo(i),distancia(i),datos] = analisis_viraje(propul,aerodinamica,viraje,W(i),h_inicial,turn_mode,datos,i,Prop_data,data_electric);
            tiempo(i) = viraje(1);
            datos(i).nombre = 'Turn Waiting';
        end
        
        %% Case 12: LANDING
        if strcmp(seg(i).nombre,'Landing') == 1
            [fuel(i),tiempo(i),distancia(i),datos] = analisis_aterrizaje(seg(i).datos,propul,aerodinamica,aero_aterrizaje,W(i),datos,i,Prop_data,data_electric);
            datos(i).nombre = 'Landing';
        end
        
%         progressbar([],i/tramos);
        
%         fh=findall(0,'type','figure');
%         if length(fh) == 1,
%             fuel_total = -777; tiempo_total = -777; distancia_total = -777; W = -777; datos=-777;
%             return; end;

        %% Case 13: DUMMY
        if strcmp(seg(i).nombre,'Dummy') == 1
            
            h_inicial = seg(i).datos.h_inicial;
            cruise_mode    = seg(i).opcion;
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
            
            % crucero(1) = h_inicial_cr;% - [m] % Altura inicial
            % crucero(2) = dist_final_cr;% - [m] % 1: DISTANCIA FINAL
            % crucero(3) = V_cr; % Velocidad de crucero
            % crucero(4) = Mach_cr;% - [-] % 2: MACH DE VUELO
            % crucero(5) = CL_cr;% - [-]% 3: CL DE CRUCERO
            % crucero(6) = delta_T_cr;% - []% 4: PALANCA DE GASES
            % crucero(7) = V_ini_cr;% - [m/s] % 5: VELOCIDAD INICIAL
            % crucero(8) = V_fin_cr;% - [m/s] % 6: VELOCIDAD FINAL
            % crucero(9) = fuel_cr;% - [kg] % 7: COMBUSTIBLE A QUEMAR
            % crucero(10) = Cd0_cr;% - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
            % crucero(11) = k1_cr;% - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
            % crucero(12) = k2_cr;% - []% 10: K2 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL

                    crucero(1) = seg(i).datos.dist_final;
                    crucero(2) = seg(i).datos.Mach;
                    crucero(3) = -1;
                    crucero(4) = -1;
                    crucero(5) = -1;
                    crucero(6) = -1;
                    crucero(7) = -1;
                    crucero(8) = -1;
                    crucero(9) = -1;
                    crucero(10) = -1;

            [fuel(i),tiempo(i),distancia(i),datos] = analisis_crucero_v3(propul,aerodinamica,crucero,W(i),h_inicial,cruise_mode,datos,i,Geo_tier,AC_CONFIGURATION,Prop_data,data_electric);
            datos(i).nombre = 'Dummy';
        end
        
        
        %% Fuel total
        if saltar == 0
            W(i+1) = W(i)-fuel(i)*g;
        else
            saltar = 0;
        end
        
        global flag_palanca
        
        if flag_palanca == 1
            fuel_total = -777; tiempo_total = -777; distancia_total = -777; W = -777; datos=-777;
            return;
        end
    end
    
    fuel_total = sum(fuel);
    tiempo_total = sum(tiempo);
    distancia_total = sum(distancia);
    mbaterias_total = sum(mbaterias);

    if propul(1) == 4
%         mbaterias_inicial = mbaterias_total + (1-data_electric(8))*mbaterias_inicial/100;
        mbaterias_inicial = mbaterias_total;
        contador = contador + 1;
    else
        if fuel_inicial - fuel_total - (pesos(4)/100)*fuel_inicial < -10
            %prog = fuel_inicial - fuel_total - (pesos(4)/100)*fuel_inicial
            fuel_inicial = fuel_total + pesos(4)*fuel_inicial/100;
            contador = contador + 1;
        else
            contador = 51;
        end
    end
    
%     progressbar(contador/10)
end


%% PESOS
% 1: PESO EN VACIO
% 2: CARGA DE PAGO INICIAL
% 3: PESO TRIPULACION
% 4: FUEL RESTANTE AL ACABAR

m_f = fuel_total;
m_f_res = fuel_inicial;
m_bat_res = mbaterias_inicial;
m_TOW = pesos(1) + pesos(2) + pesos(3) + m_f_res + m_bat_res;
m_e = pesos(1);
m_p = pesos(2) + pesos(3);

m_f_W0 = m_f_res/m_TOW;
m_e_W0 = m_e/m_TOW;

Weights_AP.m_bat = m_bat_res;
Weights_AP.m_TOW = m_TOW;
Weights_AP.m_e = m_e;
Weights_AP.m_f = m_f;
Weights_AP.m_F = m_TOW - m_f;
for i=1:tramos
    if strcmp(seg(i).nombre,'Load_Deployment') == 1
        Weights_AP.m_F = Weights_AP.m_F - seg(i).datos.carga;
    end
end
Weights_AP.m_p = m_p;
Weights_AP.m_f_W0 = m_f_W0;
Weights_AP.m_e_W0 = m_e_W0;

Total_Datos.fuel_total = fuel_total;
Total_Datos.tiempo_total = tiempo_total;
Total_Datos.distancia_total = distancia_total;
Total_Datos.W = W;

%CALCULO DEL CASM
% Cost Index
CI = 0.3; % Kg/seg
cost_fuel_USD = 2.948; % $/gall
cost_fuel_cent = cost_fuel_USD*100; % cent/gall
Gallon2liters = 3.785412;
density_fuel = 6.70; %lb/gallon
lb2kg = 0.4535924; % libras/kg
nm2m = 1852; % nautical mile/meter
C_fuel = cost_fuel_cent/density_fuel/lb2kg; % coste combustible centavos/kg
DOC = (tiempo_total*CI + fuel_total)*C_fuel;
ASM = pesos(2)/100 * distancia_total * (1/1852);
CASM = DOC/ASM;

datos(end+1).TOTAL.CASM = CASM;
datos(end).TOTAL.fuel = fuel;
datos(end).TOTAL.distancia = distancia;
datos(end).TOTAL.tiempo = tiempo;
% datos(end).lista_variables = [{''};{'Fuel total'};{'Distancia total'};{'Tiempo total'};{'CASM'}];
%-------------
% progressbar(1);
end
