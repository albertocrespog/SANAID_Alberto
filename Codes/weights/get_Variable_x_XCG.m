function [conditions,OUTPUT_read_XLSX] = get_Variable_x_XCG(m_TOW,case_AC,conditions,OUTPUT_read_XLSX,Modifications)

Dxcg = Modifications.Dxcg;

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
        %% Rutina para determina el procentaje de combustible disponible para poder determinar el peso de entrada y calcular el XCG y las inercias
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
        porcentaje_actual_f = 1;

        % Estimación peso máximo con máxima disponibilidad de fuel hasta los 120kg
        FUEL = 1;
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
           %FUEL_XLSX = OUTPUT_read_XLSX.Weights_flags.FUEL;
           [percentage_FUEL OUTPUT_f] = CG_TARSIS120_KSA_V6B(W_PL, Tabla_cg_inercias_T120, CATIAT120, n_MSL,porcentaje_actual_f);
        end

        % almacena pero máximo posible
        m_TOW_f = OUTPUT_f(1);

        % Estimación peso máximo con sin combustible
        FUEL = 0;
        % valor del fuel disponible dependiendo de la máxima cantidad de fuel que se puede almacenar
        porcentaje_actual_0f = 0;
        
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
           %FUEL_XLSX = OUTPUT_read_XLSX.Weights_flags.FUEL;
           [percentage_FUEL OUTPUT_0f] = CG_TARSIS120_KSA_V6B(W_PL, Tabla_cg_inercias_T120, CATIAT120, n_MSL,porcentaje_actual_0f);
        end
        
        % almacena pero posible sin fuel
        m_TOW_0f = OUTPUT_0f(1);
        % valor del fuel disponible dependiendo de la máxima cantidad de fuel que se puede almacenar
        porcentaje_actual = 1;
        
        % calcula porcentaje de combustible disponible cantidad de combustible disponible
        m_fuel_remaining = (m_TOW - m_TOW_0f)/(m_TOW_f - m_TOW_0f);
        
        %% Calcula ya las condicviones exactas para el peso asociado
        if conditions.T120 == 1
            OUTPUT = CG_TARSIS120(m_fuel_remaining,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        elseif conditions.T120 == 2
            OUTPUT = CG_TARSIS120_KSA(m_fuel_remaining,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120);
        elseif conditions.T120 == 3
            CG_DEP = OUTPUT_read_XLSX.Weights_flags.CG_DEP;
            OUTPUT = CG_TARSIS120_KSA_V5(m_fuel_remaining,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120,CG_DEP);
        elseif conditions.T120 == 4
           % Lectura datos T120
           Tabla_cg_inercias_T120 = OUTPUT_read_XLSX.Weights_flags.Tabla_cg_inercias_T120;
           % Lectura datos CG depósito
           %FUEL_XLSX = OUTPUT_read_XLSX.Weights_flags.FUEL;
           [percentage_FUEL OUTPUT] = CG_TARSIS120_KSA_V6B(W_PL, Tabla_cg_inercias_T120, CATIAT120, n_MSL,m_fuel_remaining);
        end


        %% Almacena los valores de masa, posición CDG x,y,z e inercias
        % Conversion factor depending on the units of the CG_TARSIS120 (1 for m 1000 for mm)
        % conv_factor = 1 old version
        % conv_factor = 1000 new version
        conv_factor = OUTPUT_read_XLSX.Weights_flags.conv_factor;
                
        conditions.m_TOW = OUTPUT(1);
        conditions.x_XCG = OUTPUT(2)/conv_factor  + Dxcg;
        
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = OUTPUT(2)/conv_factor + Dxcg;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = OUTPUT(3)/conv_factor;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = OUTPUT(4)/conv_factor;

        OUTPUT_read_XLSX.Weights_flags.I_xx = OUTPUT(5); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = OUTPUT(6); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = OUTPUT(7); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = OUTPUT(8); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = OUTPUT(9); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = OUTPUT(10);
end

