function Saving_DATA_complete(Plot_Options,Storing_DATA,M_alpha_cero,V_alpha_cero,case_AC,OUTPUT_read_XLSX,filenameS,Sheet_AC_var)
                             
% Stores the names of folders zand address
filename = filenameS.filenameB;

Sheet_AC1 = Sheet_AC_var.Sheet_AC1;
Sheet_AC2 = Sheet_AC_var.Sheet_AC2;
Sheet_AC3 = Sheet_AC_var.Sheet_AC3;
Sheet_AC4 = Sheet_AC_var.Sheet_AC4;

Storing_TRIM_VAR_DATA.M_alpha_cer = M_alpha_cero;
Storing_TRIM_VAR_DATA.V_alpha_cero = V_alpha_cero;
Storing_GEO_DATA_1 = Storing_DATA.Storing_GEO_DATA_1;
Storing_WEIGHT_DATA_1 = Storing_DATA.Storing_WEIGHT_DATA_1;
Storing_AERO_DATA_1 = Storing_DATA.Storing_AERO_DATA_1;
% Performance
Storing_PROPULSION_DATA_1 = Storing_DATA.Storing_PROPULSION_DATA_1;
Storing_PERFORMANCE_DATA_1 = Storing_DATA.Storing_PERFORMANCE_DATA_1;
Storing_PERFORMANCE_DATA_21 = Storing_DATA.Storing_PERFORMANCE_DATA_21;
Storing_PERFORMANCE_DATA_22 = Storing_DATA.Storing_PERFORMANCE_DATA_22;
Storing_PERFORMANCE_DATA_23 = Storing_DATA.Storing_PERFORMANCE_DATA_23;
Storing_PERFORMANCE_DATA_24 = Storing_DATA.Storing_PERFORMANCE_DATA_24;
Storing_PERFORMANCE_DATA_25 = Storing_DATA.Storing_PERFORMANCE_DATA_25;
Storing_PERFORMANCE_DATA_26 = Storing_DATA.Storing_PERFORMANCE_DATA_26;
Storing_PERFORMANCE_DATA_27 = Storing_DATA.Storing_PERFORMANCE_DATA_27;
% Stability
Storing_STABILITY_DATA_1 = Storing_DATA.Storing_STABILITY_DATA_1;
Storing_STABILITY_DATA_2 = Storing_DATA.Storing_STABILITY_DATA_2;
Storing_STABILITY_DATA_2B = Storing_DATA.Storing_STABILITY_DATA_2B;
Storing_STABILITY_DATA_3 = Storing_DATA.Storing_STABILITY_DATA_3;
Storing_STABILITY_DATA_4A = Storing_DATA.Storing_STABILITY_DATA_4A;
Storing_STABILITY_DATA_4B = Storing_DATA.Storing_STABILITY_DATA_4B;
Storing_STABILITY_DATA_4C = Storing_DATA.Storing_STABILITY_DATA_4C;
Storing_STABILITY_DATA_4D = Storing_DATA.Storing_STABILITY_DATA_4D;
Storing_STABILITY_DATA_5 = Storing_DATA.Storing_STABILITY_DATA_5;

% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filename_DATA = strcat(filename,filenamed);
filename = filenameS.filename_DATA;

%         filenameb = 'Results_AC2.xlsx';
%         filenamec = 'Figs';
%         filenameB = strcat(filename,filenameb);
%         filename_Plots = strcat(filename,filenamec);


% Writes the data in the selected floder according to aircraft selected
save1 = '/Plot_Options.mat';
name1A   = strcat(filename,save1);
save2 = '/Storing_GEO_DATA_1.mat';
name2A   = strcat(filename,save2);
save3 = '/Storing_WEIGHT_DATA_1.mat';
name3A   = strcat(filename,save3);
save4 = '/Storing_AERO_DATA_1.mat';
name4A   = strcat(filename,save4);
% Performance
save5 = '/Storing_PERFORMANCE_DATA_1.mat';
name5A   = strcat(filename,save5);
save5b = '/Storing_PERFORMANCE_DATA_21.mat';
name5B   = strcat(filename,save5b);
save5c = '/Storing_PERFORMANCE_DATA_22.mat';
name5C   = strcat(filename,save5c);
save5d = '/Storing_PERFORMANCE_DATA_23.mat';
name5D   = strcat(filename,save5d);
save5e = '/Storing_PERFORMANCE_DATA_24.mat';
name5E   = strcat(filename,save5e);
save5f = '/Storing_PERFORMANCE_DATA_25.mat';
name5F   = strcat(filename,save5f);
save5g = '/Storing_PERFORMANCE_DATA_26.mat';
name5G   = strcat(filename,save5g);
save5h = '/Storing_PERFORMANCE_DATA_27.mat';
name5H   = strcat(filename,save5h);

% Stability
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

save12b = '/Storing_STABILITY_DATA_4C.mat';
name12B   = strcat(filename,save12b);

save13 = '/Storing_STABILITY_DATA_5.mat';
name13A   = strcat(filename,save13);
save14 = '/Storing_PROPULSION_DATA.mat';
name14A   = strcat(filename,save14);
save15 = '/Storing_TRIM_VAR_DATA.mat';
name15A = strcat(filename,save15);

save(name1A, 'Plot_Options')
save(name2A, 'Storing_GEO_DATA_1')
save(name3A, 'Storing_WEIGHT_DATA_1')
save(name4A, 'Storing_AERO_DATA_1')
% Performance
save(name5A, 'Storing_PERFORMANCE_DATA_1')
save(name5B, 'Storing_PERFORMANCE_DATA_21')
save(name5C, 'Storing_PERFORMANCE_DATA_22')
save(name5D, 'Storing_PERFORMANCE_DATA_23')
save(name5E, 'Storing_PERFORMANCE_DATA_24')
save(name5F, 'Storing_PERFORMANCE_DATA_25')
save(name5G, 'Storing_PERFORMANCE_DATA_26')
save(name5H, 'Storing_PERFORMANCE_DATA_27')
%Stability
save(name6A, 'Storing_STABILITY_DATA_1')
save(name7A, 'Storing_STABILITY_DATA_2')
save(name8A, 'Storing_STABILITY_DATA_2B')
save(name9A, 'Storing_STABILITY_DATA_3')
save(name10A, 'Storing_STABILITY_DATA_4A')
save(name11A, 'Storing_STABILITY_DATA_4B')
save(name12A, 'Storing_STABILITY_DATA_4C')
save(name12B, 'Storing_STABILITY_DATA_4D')
save(name13A, 'Storing_STABILITY_DATA_5')
save(name14A, 'Storing_PROPULSION_DATA_1')
save(name15A, 'Storing_TRIM_VAR_DATA')
