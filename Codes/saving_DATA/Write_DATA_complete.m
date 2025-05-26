%% Function that defines where to store the datra in teh excel depending onm the aircraft selected
function Write_DATA_complete(Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_AERO_DATA_1,Storing_PERFORMANCE_DATA_1,Storing_PERFORMANCE_DATA_2,Storing_STABILITY_DATA_1,...
        Storing_STABILITY_DATA_2,Storing_STABILITY_DATA_2B,Storing_STABILITY_DATA_3,Storing_STABILITY_DATA_4A,Storing_STABILITY_DATA_4B,...
        Storing_STABILITY_DATA_4C,Storing_STABILITY_DATA_4D,Storing_STABILITY_DATA_5,...
        Storing_PROPULSION_DATA,M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX)
    
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        filename = '../Results/EMERGENTIA/Results_AC1.xlsx';
        Sheet_AC = 'AC1';
        Sheet_AC2 = 'Stab_der_Parts';
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        filename = '../Results/EMERGENTIA/Results_AC2.xlsx';
        Sheet_AC = 'AC2';
    case 3 % case_AC = 3 - PEPIÑO XXL
        filename = '../Results/PEPINO/Results_AC3.xlsx';
        Sheet_AC = 'AC3';
    case 4 % case_AC = 4 - A400
        filename = '../Results/A400/Results_AC4.xlsx';
        Sheet_AC = 'AC4';
    case 5 % case_AC = 5 - MISC
        filename = '../Results/MISC/Results_AC5.xlsx';
        Sheet_AC = 'AC5';        
    case 6 % case_AC = 6 - MILVUS
        filename = '../Results/MILVUS/Results_AC6.xlsx';
        Sheet_AC = 'AC6';        
    case 7 % Existing Aircraft = 7
        filename = '../Results/MISC/Results_AC7.xlsx';
        Sheet_AC = 'AC7';        
    case 8 % case_AC = 8 - TAMIZ
        filename = '../Results/E26TAMIZ/Results_AC8.xlsx';
        Sheet_AC = 'AC8';        
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
        filename = '../Results/EMERGENTIA/Results_AC9.xlsx';
        Sheet_AC = 'AC9';
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
        filename = '../Results/EMERGENTIA/Results_AC10.xlsx';
        Sheet_AC = 'AC10';
    case 11 % case_AC = 11 - EMERGENTIA Manufacturing
        filename = '../Results/EMERGENTIA_Manufactured/Results_AC11.xlsx';
        Sheet_AC = 'AC11';
    case 12 % case_AC = 12 - ALO
        filename = '../Results/ALO/Results_AC12.xlsx';
        Sheet_AC = 'AC12';Write_DATA_complete_AC_v2
    case 13 % case_AC = 13 - ALO Fuel Cell
        filename = '../Results/ALO_Fuel_Cell/Results_AC13.xlsx';
        Sheet_AC = 'AC13';
    case 14 % case_AC = 14 - A320
        filename = '../Results/ONEiRE/Results_AC14.xlsx';
        Sheet_AC = 'AC14';
    case 15 % case_AC = 15 - EMERGENTIA Tandem Wing
        filename = '../Results/EMERGENTIA_v2/Results_AC15.xlsx';
        Sheet_AC = 'AC15';
    case 16 % case_AC = 16 - Tarsis T75
        filename = '../Results/TARSIS_75/Results_AC16.xlsx';
        Sheet_AC = 'AC16';
    case 17 % case_AC = 17 - Tarsis T120
        filename = '../Results/TARSIS_120/Results_AC17.xlsx';
        Sheet_AC = 'AC17';
    case 18 % case_AC = 18 - BAT
        filename = '../Results/BAT/Results_AC18.xlsx';
        Sheet_AC = 'AC18';
    case 19 % case_AC = 19 - FALCON2000
        filename = '../Results/FALCON2000/Results_AC19.xlsx';
        Sheet_AC = 'AC19';
    case 20 % case_AC = 20 - SOLAR TII
        filename = '../Results/SOLARTII/Results_AC20.xlsx';
        Sheet_AC = 'AC20';
    case 21 % case_AC = 21 - VANTUS 24
        filename = '../Results/VANTUS/Results_AC21.xlsx';
        Sheet_AC = 'AC21';
    case 22 % case_AC = 22 - Cessna 208
        filename = '../Results/Cessna208/Results_AC22.xlsx';
        Sheet_AC = 'AC22';
    case 23 % case_AC = 23 - King Air 350i
        filename = '../Results/KingAir350/Results_AC23.xlsx';
        Sheet_AC = 'AC23';
    case 24 % case_AC = 24 - Future 1
        filename = '../Results/MK84/Results_AC24.xlsx';
        Sheet_AC = 'AC24';
    case 25 % case_AC = 25 - Future 2
        filename = '../Results/MK84B/Results_AC25.xlsx';
        Sheet_AC = 'AC25';
    case 26 % case_AC = 26 - Future 3
        filename = '../Results/Future3/Results_AC26.xlsx';
        Sheet_AC = 'AC26';
    case 27 % case_AC = 27 - Future 4
        filename = '../Results/Future4/Results_AC27.xlsx';
        Sheet_AC = 'AC27';
    case 28 % case_AC = 28 - Future 5
        filename = '../Results/Future5/Results_AC28.xlsx';
        Sheet_AC = 'AC28';
    case 29 % case_AC = 29 - Future 6
        filename = '../Results/Future6/Results_AC29.xlsx';
        Sheet_AC = 'AC29';
    case 30 % case_AC = 30 - Future 7
        filename = '../Results/Future7/Results_AC30.xlsx';
        Sheet_AC = 'AC30';
end


Sheet_AC = 'Results';
Sheet_AC2 = 'Stab_der_Parts';


Write_DATA_complete_AC_v2(Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_AERO_DATA_1,Storing_PERFORMANCE_DATA_1,Storing_PERFORMANCE_DATA_2,Storing_STABILITY_DATA_1,...
    Storing_STABILITY_DATA_2,Storing_STABILITY_DATA_2B,Storing_STABILITY_DATA_3,Storing_STABILITY_DATA_4A,Storing_STABILITY_DATA_4B,...
    Storing_STABILITY_DATA_4C,Storing_STABILITY_DATA_4D,Storing_STABILITY_DATA_5,...
    Storing_PROPULSION_DATA,M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filename,Sheet_AC,Sheet_AC2);

