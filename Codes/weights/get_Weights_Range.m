function Weight_Range = get_Weights_Range(case_AC,conditions,OUTPUT_read_XLSX)

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
    case 2 % case_AC = 2 - EMERGENTIA 1:2
    case 3 % case_AC = 3 - PEPIÑO XXL
    case 4 % case_AC = 4 - COMERCIAL
    case 5 % case_AC = 5 - WIGL
    case 6 % case_AC = 6 - CERVERA
    case 7 % Existing Aircraft = 7
    case 8 % case_AC = 8 - TAMIZ
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
    case 11 % case_AC = 11 - EMERGENTIA Manufacturing
    case 12 % case_AC = 12 - ALO
    case 13 % case_AC = 13 - ALO Fuel Cell
    case 14 % case_AC = 14 - A320
    case 15 % case_AC = 15 - EMERGENTIA Tandem Wing
    case 16 % case_AC = 16 - Tarsis T75
    case 17 % case_AC = 17 - Tarsis T120
        % Weight Configuration from AERTEC
        MOTOR = conditions.MOTOR; % 1 para SP210, 2 para DA215
        RACK = conditions.RACK; % 1 == rack, 0 == sin rack
        n_MSL = conditions.n_MSL; %0,1,2,3,4
        DEPOSITO = conditions.DEPOSITO; % 1 == deposito original, 2 == depósito húmedo
        W_PL = conditions.W_PL;

        % full fuel weight
        CATIAT120 = OUTPUT_read_XLSX.Weights_flags.CATIAT120;
        FUEL = 1;
        % valor del fuel disponible dependiendo de la máxima cantidad de fuel que se puede almacenar
        porcentaje_actual = 1;
        
        if conditions.T120 == 1
            OUTPUT_f = CG_TARSIS120(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            OUTPUT_f = CG_TARSIS120_KSA(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            CG_DEP = OUTPUT_read_XLSX.Weights_flags.CG_DEP;
            OUTPUT_f = CG_TARSIS120_KSA_V5(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120,CG_DEP);
        elseif conditions.T120 == 4
           % Lectura datos T120
           Tabla_cg_inercias_T120 = OUTPUT_read_XLSX.Weights_flags.Tabla_cg_inercias_T120;
           % Lectura datos CG depósito
           % CATIAT120 = OUTPUT_read_XLSX.Weights_flags.CATIAT120;
           [percentage_FUEL OUTPUT_f] = CG_TARSIS120_KSA_V6B(W_PL, Tabla_cg_inercias_T120, CATIAT120, n_MSL,porcentaje_actual);
           conditions.FUEL = percentage_FUEL; % overwrites predefine value
        end

        %         OUTPUT_f = CG_TARSIS120(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        
        % half fuel weight
        FUEL = 0.5;
        % valor del fuel disponible dependiendo de la máxima cantidad de fuel que se puede almacenar
        porcentaje_actual = 0.5;

        if conditions.T120 == 1
            OUTPUT_2f = CG_TARSIS120(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            OUTPUT_2f = CG_TARSIS120_KSA(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            CG_DEP = OUTPUT_read_XLSX.Weights_flags.CG_DEP;
            OUTPUT_2f = CG_TARSIS120_KSA_V5(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120,CG_DEP);
        elseif conditions.T120 == 4
           % Lectura datos T120
           Tabla_cg_inercias_T120 = OUTPUT_read_XLSX.Weights_flags.Tabla_cg_inercias_T120;
           % Lectura datos CG depósito
           %FUEL_XLSX = OUTPUT_read_XLSX.Weights_flags.FUEL;
           [percentage_FUEL OUTPUT_2f] = CG_TARSIS120_KSA_V6B(W_PL, Tabla_cg_inercias_T120, CATIAT120, n_MSL,porcentaje_actual);
           conditions.FUEL = percentage_FUEL; % overwrites predefine value
        end

        %         OUTPUT_2f = CG_TARSIS120(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        
        % Empty fuel weight
        FUEL = 0;
        % valor del fuel disponible dependiendo de la máxima cantidad de fuel que se puede almacenar
        porcentaje_actual = 0.0;

        if conditions.T120 == 1
            OUTPUT_0f = CG_TARSIS120(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            OUTPUT_0f = CG_TARSIS120_KSA(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            CG_DEP = OUTPUT_read_XLSX.Weights_flags.CG_DEP;
            OUTPUT_0f = CG_TARSIS120_KSA_V5(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120,CG_DEP);
        elseif conditions.T120 == 4
           % Lectura datos T120
           Tabla_cg_inercias_T120 = OUTPUT_read_XLSX.Weights_flags.Tabla_cg_inercias_T120;
           % Lectura datos CG depósito
           % FUEL_XLSX = OUTPUT_read_XLSX.Weights_flags.FUEL;
           [percentage_FUEL OUTPUT_0f] = CG_TARSIS120_KSA_V6B(W_PL, Tabla_cg_inercias_T120, CATIAT120, n_MSL,porcentaje_actual);
           conditions.FUEL = percentage_FUEL; % overwrites predefine value
        end

        %         OUTPUT_0f = CG_TARSIS120(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        
        Weight_Range = [OUTPUT_f(1);OUTPUT_2f(1);OUTPUT_0f(1)];
        
%      Weight_Range = Weight_Range_Estimation_TARSIS(conf_missiles);
%      x_XCG = get_Variable_x_XCG(m_TOW,case_AC,conf_missiles) + Dxcg;% - 0.078381885757103;
%  
%      m_TOW = max(Weight_Range);
%      x_XCG = get_Variable_x_XCG(m_TOW,case_AC,conf_missiles) + Dxcg;% - 0.078381885757103;
%      conditions.m_TOW = m_TOW;
%      conditions.x_XCG = x_XCG;  
    case 18 % case_AC = 18 - BAT
    case 19 % case_AC = 19 - FALCON200
    case 20 % case_AC = 20 - SOLARTII
    case 21 % case_AC = 21 - EVOL3
        
end

