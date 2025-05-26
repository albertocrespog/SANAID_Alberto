function Saving_DATA_complete_v2(Plot_Options,Storing_DATA,M_alpha_cero,V_alpha_cero,case_AC)


Storing_TRIM_VAR_DATA.M_alpha_cero = M_alpha_cero;
Storing_TRIM_VAR_DATA.V_alpha_cero = V_alpha_cero;

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        filename = '../Results/EMERGENTIA';
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        filename = '../Results/EMERGENTIA';
    case 3 % case_AC = 3 - PEPIÑO XXL
        filename = '../Results/PEPINO';
    case 4 % case_AC = 4 - A400
        filename = '../Results/A400';
    case 5 % case_AC = 5 - MISC
        filename = '../Results/MISC';
    case 6 % case_AC = 6 - MILVUS
        filename = '../Results/MILVUS';
    case 7 % Existing Aircraft = 7
        filename = '../Results/MISC';
    case 8 % case_AC = 8 - TAMIZ
        filename = '../Results/E26TAMIZ';
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
        filename = '../Results/EMERGENTIA';
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
        filename = '../Results/EMERGENTIA';
    case 11 % case_AC = 11 - EMERGENTIA Manufacturing
        filename = '../Results/EMERGENTIA_Manufactured';
    case 12 % case_AC = 12 - ALO
        filename = '../Results/ALO';
    case 13 % case_AC = 13 - ALO Fuel Cell
        filename = '../Results/ALO_Fuel_Cell';
    case 14 % case_AC = 14 - A320
        filename = '../Results/ONEiRE';
    case 15 % case_AC = 15 - EMERGENTIA Tandem Wing
        filename = '../Results/EMERGENTIA_v2';
    case 16 % case_AC = 16 - Tarsis T75
        filename = '../Results/TARSIS_75';
    case 17 % case_AC = 17 - Tarsis T120
        filename = '../Results/TARSIS_120';
    case 18 % case_AC = 18 - BAT
        filename = '../Results/BAT';
    case 19 % case_AC = 19 - FALCON2000
        filename = '../Results/FALCON2000';
    case 20 % case_AC = 20 - SOLAR TII
        filename = '../Results/SOLARTII';
    case 21 % case_AC = 21 - VANTUS 24
        filename = '../Results/VANTUS24';
    case 22 % case_AC = 22 - Cessna 208
        filename = '../Results/Cessna208';
    case 23 % case_AC = 23 - King Air 350i
        filename = '../Results/KingAir350';
    case 24 % case_AC = 24 - Future1
        filename = '../Results/MK84';
    case 25 % case_AC = 25 - Future2
        filename = '../Results/MK84B';
    case 26 % case_AC = 26 - Future3
        filename = '../Results/Future3';
    case 27 % case_AC = 27 - Future4
        filename = '../Results/Future4';
    case 28 % case_AC = 28 - Future5
        filename = '../Results/Future5';
    case 29 % case_AC = 29 - Future6
        filename = '../Results/Future6';
    case 30 % case_AC = 30 - Future7
        filename = '../Results/Future7';
end

% Writes the data in the selected floder according to aircraft selected
save1 = '/Plot_Options.mat';
name1A   = strcat(filename,save1);
save2 = '/Storing_GEO_DATA_1.mat';
name2A   = strcat(filename,save2);
save3 = '/Storing_WEIGHT_DATA_1.mat';
name3A   = strcat(filename,save3);
save4 = '/Storing_AERO_DATA_1.mat';
name4A   = strcat(filename,save4);
save5 = '/Storing_PERFORMANCE_DATA_1.mat';
name5A   = strcat(filename,save5);
save5b = '/Storing_PERFORMANCE_DATA_2.mat';
name5B   = strcat(filename,save5b);
save5c = '/Storing_PERFORMANCE_DATA_3.mat';
name5C   = strcat(filename,save5c);
save5d = '/Storing_PROPULSION_DATA_1.mat';
name5D   = strcat(filename,save5d);
save6 = '/Storing_STABILITY_DATA_1.mat';
name6A   = strcat(filename,save6);
save7 = '/Storing_STABILITY_DATA_2.mat';
name7A   = strcat(filename,save7);
save8 = '/Storing_STABILITY_DATA_2B.mat';
name8A   = strcat(filename,save8);
save9 = '/Storing_STABILITY_DATA_3.mat';
name9A   = strcat(filename,save9);
save10 = '/Storing_STABILITY_DATA_4A.mat';
name10A   = strcat(filename,save10);
save11 = '/Storing_STABILITY_DATA_4B.mat';
name11A   = strcat(filename,save11);
save12 = '/Storing_STABILITY_DATA_4C.mat';
name12A   = strcat(filename,save12);
save12b = '/Storing_STABILITY_DATA_4D.mat';
name12B   = strcat(filename,save12b);
save13 = '/Storing_STABILITY_DATA_5.mat';
name13A   = strcat(filename,save13);
save14 = '/Storing_PROPULSION_DATA.mat';
name14A   = strcat(filename,save14);
save15 = '/Storing_TRIM_VAR_DATA.mat';
name15A = strcat(filename,save15);


Storing_GEO_DATA_1 = Storing_DATA.Storing_GEO_DATA_1;
Storing_WEIGHT_DATA_1 = Storing_DATA.Storing_WEIGHT_DATA_1;
Storing_AERO_DATA_1 = Storing_DATA.Storing_AERO_DATA_1;
Storing_PROPULSION_DATA_1 = Storing_DATA.Storing_PROPULSION_DATA_1;
Storing_PERFORMANCE_DATA_1 = Storing_DATA.Storing_PERFORMANCE_DATA_1;
Storing_PERFORMANCE_DATA_2 = Storing_DATA.Storing_PERFORMANCE_DATA_2;
Storing_PERFORMANCE_DATA_3 = Storing_DATA.Storing_PERFORMANCE_DATA_3;
Storing_STABILITY_DATA_1 = Storing_DATA.Storing_STABILITY_DATA_1;
Storing_STABILITY_DATA_2 = Storing_DATA.Storing_STABILITY_DATA_2;
Storing_STABILITY_DATA_2B = Storing_DATA.Storing_STABILITY_DATA_2B;
Storing_STABILITY_DATA_3 = Storing_DATA.Storing_STABILITY_DATA_3;
Storing_STABILITY_DATA_4A = Storing_DATA.Storing_STABILITY_DATA_4A;
Storing_STABILITY_DATA_4B = Storing_DATA.Storing_STABILITY_DATA_4B;
Storing_STABILITY_DATA_4C = Storing_DATA.Storing_STABILITY_DATA_4C;
Storing_STABILITY_DATA_4D = Storing_DATA.Storing_STABILITY_DATA_4D;
Storing_STABILITY_DATA_5 = Storing_DATA.Storing_STABILITY_DATA_5;

save(name1A, 'Plot_Options')
save(name2A, 'Storing_GEO_DATA_1')
save(name3A, 'Storing_WEIGHT_DATA_1')
save(name4A, 'Storing_AERO_DATA_1')
save(name5A, 'Storing_PERFORMANCE_DATA_1')
save(name5B, 'Storing_PERFORMANCE_DATA_2')
save(name5C, 'Storing_PERFORMANCE_DATA_3')
save(name5D, 'Storing_PROPULSION_DATA_1')
save(name6A, 'Storing_STABILITY_DATA_1')
save(name7A, 'Storing_STABILITY_DATA_2')
save(name8A, 'Storing_STABILITY_DATA_2B')
save(name9A, 'Storing_STABILITY_DATA_3')
save(name10A, 'Storing_STABILITY_DATA_4A')
save(name11A, 'Storing_STABILITY_DATA_4B')
save(name12A, 'Storing_STABILITY_DATA_4C')
save(name12B, 'Storing_STABILITY_DATA_4D')
save(name13A, 'Storing_STABILITY_DATA_5')
save(name14A, 'Storing_PROPULSION_DATA')
save(name15A, 'Storing_TRIM_VAR_DATA')

% save('data/GLOBAL_DATA/Plot_Options.mat', 'Plot_Options')
% save('data/GLOBAL_DATA/Storing_GEO_DATA_1.mat', 'Storing_GEO_DATA_1')
% save('data/GLOBAL_DATA/Storing_WEIGHT_DATA_1.mat', 'Storing_WEIGHT_DATA_1')
% save('data/GLOBAL_DATA/Storing_AERO_DATA_1.mat', 'Storing_AERO_DATA_1')
% save('data/GLOBAL_DATA/Storing_PERFORMANCE_DATA_1.mat', 'Storing_PERFORMANCE_DATA_1')
% save('data/GLOBAL_DATA/Storing_STABILITY_DATA_1.mat', 'Storing_STABILITY_DATA_1')
% save('data/GLOBAL_DATA/Storing_STABILITY_DATA_2.mat', 'Storing_STABILITY_DATA_2')
% save('data/GLOBAL_DATA/Storing_STABILITY_DATA_2B.mat', 'Storing_STABILITY_DATA_2B')
% save('data/GLOBAL_DATA/Storing_STABILITY_DATA_3.mat', 'Storing_STABILITY_DATA_3')
% save('data/GLOBAL_DATA/Storing_STABILITY_DATA_4.mat', 'Storing_STABILITY_DATA_4')
% save('data/GLOBAL_DATA/Storing_STABILITY_DATA_5.mat', 'Storing_STABILITY_DATA_5')
% save('data/GLOBAL_DATA/Storing_PROPULSION_DATA.mat', 'Storing_PROPULSION_DATA')
% save('data/GLOBAL_DATA/Storing_TRIM_VAR_DATA.mat', 'Storing_TRIM_VAR_DATA')
