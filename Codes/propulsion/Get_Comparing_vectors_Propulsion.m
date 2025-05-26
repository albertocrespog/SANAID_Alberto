function VECTOR_Prop = Get_Comparing_vectors_Propulsion(case_AC,OUTPUT_read_XLSX)

compare_prop = OUTPUT_read_XLSX.Propulsive_flags.compare_prop;

%% Data for analysis of Prop DATA
% Determines the plots that want to show
% Compares Propeller properties for different sources
% The VECTOR in order to analyze the results of Props consist of 2 grous of
% vector each prop model
% compare_prop = 1 compares APC Models alone
% compare_prop = 2 compares APC with wind tunnel models
% VECTOR_Prop.compare_props = compare_props;

%% Propulsive Model Wind Tunnel Model 1
% Number of Prop used
% 1 - APC 20x8
% 2 - APC 22x10
% 3 - APC 22x12
% 4 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
% 5 - APC 21x13W
% 6 - APC 21x14

%% Propulsive Model Wind Tunnel Model 2 (Rai models with alpha)
% Number of Prop used
% 1 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/

%% Generates the file for Aerodynamic Data
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 3 % case_AC = 3 - PEPIÑO XXL
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 4 % case_AC = 4 COMERCIAL
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 5 % case_AC = 5 WIG
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 6 % case_AC = 1 - CERVERA
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 7 % case_AC = 7 - Existing Aircraft
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 8 % case_AC = 8 - TAMIZ
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 11 % case_AC = 11 - EMERGENTIA Manufactured
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 12 % case_AC = 12 - ALO 
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 13 % case_AC = 13 - ALO Fuel Cell
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 14 % case_AC = 11 - EMERGENTIA Manufactured
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 15 % case_AC = 11 - EMERGENTIA Manufactured
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 16 % case_AC = 17 - TARSIS 75
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 17 % case_AC = 17 - TARSIS 120
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 18 % case_AC = 18 - BAT
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 19 % case_AC = 19 - FALCON2000L
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 20 % case_AC = 20 - SOLARTII
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 21 % case_AC = 21 - VANTUS
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 22 % case_AC = 22 - Cessna 208
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 23 % case_AC = 23 - King Air 305
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 24 % case_AC = 24 - Future 1
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 25 % case_AC = 25 - Future 2
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 26 % case_AC = 26 - Future 3
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 27 % case_AC = 27 - Future 4
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 28 % case_AC = 28 - Future 5
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 29 % case_AC =  29 - Future 6
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
    case 30 % case_AC = 30 - Future 7
        switch compare_prop
            case 1 % compare_prop = 1 compares APC Models alone
                VECTOR_Prop.v1 = [1];
            case 2 % compare_prop = 2 compares APC with wind tunnel models
                VECTOR_Prop.v1 = [1];
                VECTOR_Prop.v2 = [1,2];
        end
end
