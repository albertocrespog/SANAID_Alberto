%% Function that writes into the excel the selected structure of data
function Write_DATA_complete_AC_v4(Storing_DATA,M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filenameS,Sheet_AC_var)
                                  
% Changes to the different columns where to store the results
switch mission_actual
    case 0
        prefix = 'D';
    case 1
        prefix = 'E';
    case 2
        prefix = 'F';
    case 3
        prefix = 'G';
    case 4
        prefix = 'H';
    case 5
        prefix = 'I';
    case 6
        prefix = 'J';
    case 7
        prefix = 'K';
    case 8
        prefix = 'L';
    case 9
        prefix = 'M';
    case 10
        prefix = 'N';
    case 11
        prefix = 'O';
    case 12
        prefix = 'P';
    case 13
        prefix = 'Q';
    case 14
        prefix = 'R';
    case 15
        prefix = 'S';
    case 16
        prefix = 'T';
    case 17
        prefix = 'U';
    case 18
        prefix = 'V';
    case 19
        prefix = 'W';
    case 20
        prefix = 'X';
    case 21
        prefix = 'Y';
    case 22
        prefix = 'z';
end

% Stores the names of folders zand address
filename = filenameS.filenameB;

Sheet_AC1 = Sheet_AC_var.Sheet_AC1;
Sheet_AC2 = Sheet_AC_var.Sheet_AC2;
Sheet_AC3 = Sheet_AC_var.Sheet_AC3;
Sheet_AC4 = Sheet_AC_var.Sheet_AC4;
Sheet_AC5 = Sheet_AC_var.Sheet_AC5;

%% First Sheet of the Excel DATA
[Geo,Aero,Stab,Perfo,Misc,name_geo,name_aero,name_stab,name_perfo,name_misc] = Write_DATA_complete_AC_sheet1(Storing_DATA,...
    M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filename,Sheet_AC1,prefix);

writematrix(Geo,filename,'Sheet',Sheet_AC1,'Range',name_geo);
writematrix(Aero,filename,'Sheet',Sheet_AC1,'Range',name_aero);
if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
    writematrix(Stab,filename,'Sheet',Sheet_AC1,'Range',name_stab)

    %% Second Sheet of the Excel DATA
    [Aero_comparativa,name_aero_comprativa] = Write_DATA_complete_AC_sheet2(Storing_DATA,...
        M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filename,Sheet_AC2,prefix);
    writematrix(Aero_comparativa,filename,'Sheet',Sheet_AC2,'Range',name_aero_comprativa);
    
    %% Third Sheet of the Excel DATA
    [Stab_der_parts_longitudinal name_stab_der_parts_longitudinal] = Write_DATA_complete_AC_sheet3(Storing_DATA,...
        M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filename,Sheet_AC3,prefix);
    writematrix(Stab_der_parts_longitudinal,filename,'Sheet',Sheet_AC3,'Range',name_stab_der_parts_longitudinal);

    %% Fourth Sheet of the Excel DATA
    [Stab_der_parts_lateraldirectional name_stab_der_parts_lateraldirectional] = Write_DATA_complete_AC_sheet4(Storing_DATA,...
        M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filename,Sheet_AC4,prefix);
    writematrix(Stab_der_parts_lateraldirectional,filename,'Sheet',Sheet_AC4,'Range',name_stab_der_parts_lateraldirectional);

    %% Fifth Sheet of the Excel DATA
    [data_STUDY_dyn_DIMENSION name_stab_der_parts_DIMENSION] = Write_DATA_complete_AC_sheet5(Storing_DATA,...
        M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filename,Sheet_AC4,prefix);
    writematrix(data_STUDY_dyn_DIMENSION,filename,'Sheet',Sheet_AC5,'Range',name_stab_der_parts_DIMENSION);

end

if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
    writematrix(Perfo,filename,'Sheet',Sheet_AC1,'Range',name_perfo)
end

if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat_Trim_TAB == 1
    writematrix(Misc,filename,'Sheet',Sheet_AC1,'Range',name_misc)
end


