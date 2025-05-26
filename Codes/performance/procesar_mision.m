function [fuel_total,tiempo_total,distancia_total,W,datos] = procesar_mision(seg,tramos,propul,aerodinamica,pesos,aero_despegue,aero_aterrizaje)

%% PROPUL
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR
% 4: CONSUMO ESPECIFICO
% 5: AVION CIVIL =1/MILITAR = 2
% 6: EFICIENCIA DE LA HELICE (ETA_P)
% 7: DERIVACION(TURBOFANES)

%% PESOS
% 1: PESO EN VACIO
% 2: CARGA DE PAGO INICIAL
% 3: PESO TRIPULACION
% 4: %FUEL RESTANTE AL ACABAR

%% AERODINAMICA
% 1: SUPERFICIE
% 2: CD0
% 3: K1 = K
% 4: K2
% 5: CLMAX LIMPIO

%% AERODINAMICA_DESPEGUE
%       % 1: CD0 EN DESPEGUE (CONSIDERANDO TREN + FLAPS)
%       % 2: CL EN DESPEGUE
%       % 3: K CON EFECTO SUELO (O ALTURA DEL ALA SOBRE EL SUELO Y
%       ENVERGADURA DEL AVION)
%       % 4: CLMAX SUCIO DESPEGUE (FLAPS = 20º)

%% AERODINAMICA_ATERRIZAJE
%       % 1: CD0 EN ATERRIZAJE (CONSIDERANDO TREN + FLAPS + SPOILERS)
%       % 2: CL EN ATERRIZAJE
%       % 3: K CON EFECTO SUELO (O ALTURA DEL ALA SOBRE EL SUELO Y
%       ENVERGADURA DEL AVION)
%       % 4: CLMAX SUCIO ATERRIZAJE (FLAPS = 30º

%% MISION
% --> %% ARRANQUE DE LOS MOTORES + WARM-UP = TAXI
%       % 1: TEMPERATURA LOCAL
%       % 2: ALTURA LOCAL
%       % 3: PRESION LOCAL
%               % 4: PALANCA DE RALENTI EN TAXI = 0.05
%       % 5: VELOCIDAD A LA QUE HACE EL TAXI
%       % 6: TIEMPO DE ESPERA EN TAXI
%
% --> %% DESPEGUE
%       % 1: TEMPERATURA LOCAL
%       % 2: ALTURA LOCAL       
%       % 3: PRESION LOCAL
%       % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
%               % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
%               % 6: GAMMA DE SUBIDA (0.024,0.027,0.03 2,3,4 MOTORES)
%       % 7: PALANCA DE GASES PARA DESPEGUE
%
% --> %% SUBIDA
%               % 1: ALTURA FINAL       
%               % 2: GAMMA DE SUBIDA
%               % 3: MACH DE VUELO
%               % 4: VELOCIDAD TAS = 
%               % 5: VELOCIDAD EAS = 
%               % 6: PALANCA DE GASES = 
%               % 7: VELOCIDAD INICIAL = 
%               % 8: VELOCIDAD FINAL = 
% --> %% CRUCERO
%
% 1: DISTANCIA FINAL       
% 2: MACH DE VUELO
% 3: CL DE CRUCERO
% 4: PALANCA DE GASES 
% 5: VELOCIDAD INICIAL 
% 6: VELOCIDAD FINAL
% 7: COMBUSTIBLE A QUEMAR
% --> % DATOS POLAR
% 1: NUMEROS DE MACH PARA LOS QUE SE VA A DAR LA POLAR = f(M)
% 2: Cd0 = Cd0(M)
% 3: k1 = k1(M)
% 4: k2 = k2(M)

% DESCENSO
% 1: ALTURA FINAL       
% 2: GAMMA DE DESCENSO
% 3: MACH DE VUELO
% 4: VELOCIDAD TAS 
% 5: VELOCIDAD EAS 
% 6: PALANCA DE GASES
% 7: VELOCIDAD INICIAL (PARA DESCENSO ACELERADO)
% 8: VELOCIDAD FINAL (PARA DESCENSO ACELERADO)

% --> %% ATERRIZAJE
%       % 1: TEMPERATURA LOCAL
%       % 2: ALTURA LOCAL       
%       % 3: PRESION LOCAL
%       % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
%       % 5: PALANCA DE GAS PARA REVERSA
%       % 6: TIEMPO QUE SE TARDA EN ACTIVAR LOS FRENOS AL ATERRIZAR

%Peso inicial = 29545.6 kg
%Peso vacio = 16167 kg
%Pago = 7104 kg
g = 9.80665;
fuel_inicial = 0;
contador = 0;
saltar = 0;
datos.inicializar = 0;
W(1) = 0;

progressbar('Iteracion','Segmento');

while contador < 10,
    
clear W
W(1) = (pesos(1) + pesos(2) + pesos(3) + fuel_inicial)*g;

for i=1:tramos,
    progressbar([],0)

    if strcmp(seg(i).nombre,'Taxi') == 1,
        [fuel(i)] = analisis_taxi(seg(i).datos,propul,i); %kg
        datos(i).taxi.fuel = fuel(i);
        datos(i).nombre = 'Taxi';
        datos(i).lista_variables = [{''};{'Fuel'}];
    end
    
    if strcmp(seg(i).nombre,'Despegue') == 1,
        [fuel(i),tiempo(i),distancia(i),datos] = analisis_despegue(seg(i).datos,propul,aerodinamica,aero_despegue,W(i),i,datos);    
    end
    
    if strcmp(seg(i).nombre,'Subida') == 1,
    h_inicial = seg(i).datos.h_inicial;
    opcion    = seg(i).opcion;
    subida(1) = seg(i).datos.h_final;
%             subida(2) = seg(i).datos.gamma;
%             subida(3) = seg(i).datos.Mach;
%             subida(4) = seg(i).datos.EAS;
%             subida(5) = seg(i).datos.TAS;
%             subida(6) = seg(i).datos.palanca;
%             subida(7) = seg(i).datos.V_ini;
%             subida(8) = seg(i).datos.V_fin;
    switch opcion
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
    [fuel(i),tiempo(i),distancia(i),datos] = analisis_subida(propul,aerodinamica,subida,W(i),h_inicial,opcion,i,datos);
    end
    
    if strcmp(seg(i).nombre,'Crucero') == 1,
    h_inicial = seg(i).datos.h_inicial;
    opcion    = seg(i).opcion;
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
    switch opcion
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
    end
    [fuel(i),tiempo(i),distancia(i),datos] = analisis_crucero(propul,aerodinamica,crucero,W(i),h_inicial,opcion,i,datos);
    end
    
    if strcmp(seg(i).nombre,'Soltar carga') == 1,
        carga_soltada = seg(i).datos.carga;
        fuel(i) = 0;
        datos(i).soltar_carga.carga_soltada = carga_soltada;
        datos(i).soltar_carga.fuel = 0;
        datos(i).nombre = 'Soltar carga';
        datos(i).lista_variables = [{''};{'Carga soltada'}];
        W(i+1) = W(i) - carga_soltada * g;
        saltar = 1;
    end
    
    if strcmp(seg(i).nombre,'Viraje') == 1,
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
    viraje(1) = seg(i).datos.tiempo_final;
    opcion    = seg(i).opcion;
    
    switch opcion
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
        viraje(6) = seg(i).datos.vel_guiniada;
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
    [fuel(i),tiempo(i),distancia(i),datos] = analisis_viraje(propul,aerodinamica,viraje,W(i),h_inicial,opcion,i,datos);
    tiempo(i) = viraje(1);    
    end
    
    if strcmp(seg(i).nombre,'Descenso') == 1,
    h_inicial = seg(i).datos.h_inicial;
    opcion    = seg(i).opcion;
    descenso(1) = seg(i).datos.h_final;
    
    switch opcion
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
    [fuel(i),tiempo(i),distancia(i),datos] = analisis_descenso(propul,aerodinamica,descenso,W(i),h_inicial,opcion,i,datos);
    end
    
    if strcmp(seg(i).nombre,'Aterrizaje') == 1,
        [fuel(i),tiempo(i),distancia(i),datos] = analisis_aterrizaje(seg(i).datos,propul,aerodinamica,aero_aterrizaje,W(i),i,datos);    
    end

    if saltar == 0,
        W(i+1) = W(i)-fuel(i)*g;
    else
        saltar = 0;
    end
    
    global flag_palanca
    if flag_palanca == 1,
        fuel_total = -777; tiempo_total = -777; distancia_total = -777; W = -777; datos=-777;
        return;
    end
    
    progressbar([],i/tramos);
    
    fh=findall(0,'type','figure');
    if length(fh) == 1,
    fuel_total = -777; tiempo_total = -777; distancia_total = -777; W = -777; datos=-777;
    return; end; 
end

fuel_total = sum(fuel);
tiempo_total = sum(tiempo);
distancia_total = sum(distancia);


if fuel_inicial - fuel_total - (pesos(4)/100)*fuel_inicial < -10,
    %prog = fuel_inicial - fuel_total - (pesos(4)/100)*fuel_inicial
    fuel_inicial = fuel_total + pesos(4)*fuel_inicial/100;
    contador = contador + 1;   
else
    contador = 11;
end

progressbar(contador/10)
end

%CALCULO DEL CASM
DOC = (tiempo_total + fuel_total)*97.003;
ASM = pesos(2)/100 * distancia_total * (1/1852);
CASM = DOC/ASM;
datos(end+1).TOTAL.CASM = CASM;
datos(end).TOTAL.fuel = fuel;
datos(end).TOTAL.distancia = distancia;
datos(end).TOTAL.tiempo = tiempo;
datos(end).lista_variables = [{''};{'Fuel total'};{'Distancia total'};{'Tiempo total'};{'CASM'}];
%-------------
progressbar(1);
end