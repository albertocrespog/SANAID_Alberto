function  Aero_approx = get_approximate_polar_values(case_AC,Aero,Aero_TH,Geo_tier)

    %% Selects the file to red aero files
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
    case 2 % case_AC = 2 - EMERGENTIA 1:2
    case 3 % case_AC = 3 - PEPIÑO XXL
    case 4 % case_AC = 4 - COMERCIAL
    case 5 % case_AC = 5 - WIGL
    case 6 % case_AC = 6 - MILVUS
        % Aero_approx = get_Aero_MILVUS(Aero,Aero_TH,Geo_tier);
        Aero_approx = get_Aero_MILVUS_POLAR(Aero,Aero_TH,Geo_tier);
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
        Aero_approx = get_Aero_T75(Aero,Aero_TH,Geo_tier);
    case 17 % case_AC = 17 - Tarsis T120
        Aero_approx = get_Aero_T75(Aero,Aero_TH,Geo_tier);
    case 18 % case_AC = 18 - BAT
    case 19 % case_AC = 19 - FALCON2000
    case 20 % case_AC = 20 - SOLARTII

    case 21 % case_AC = 22
    case 22 % case_AC = 22
    case 23 % case_AC = 22
    case 24 % case_AC = 22
    case 25 % case_AC = 22
    case 26 % case_AC = 22
    case 27 % case_AC = 22
    case 28 % case_AC = 22
    case 29 % case_AC = 22
        Aero_approx = get_Aero_MILVUS_POLAR(Aero,Aero_TH,Geo_tier);
    case 30 % case_AC = 22
   
end

