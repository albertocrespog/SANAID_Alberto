function [case_AC OUTPUT_read_XLSX] = Initial_prompt

%% Selects Automatization process
prompt1 = sprintf('Would you like to Conduct:\n Automatic Study According to Excel Selection EMERGENTIA 100%% (Press Yes) \n or \n Select Aircraft form list (Press No)?');
answer_excel =questdlg (prompt1,'Automatic Selection','Yes','No','No');
%Handle response
switch answer_excel
    case 'Yes'
        automatic_excel=1;
        load OUTPUT_read_XLSX_EMERGENTIA_nom.mat
        OUTPUT_read_XLSX = OUTPUT_read_XLSX_EMERGENTIA;
        case_AC = 1;
    case 'No'
        automatic_excel=0;
        % Prompt aircraft Selection
        % Defines de aircraft to be analyzed
        % case_AC = 1 - EMERGENTIA 1:1 case_AC = 2 - EMERGENTIA 1:2 case_AC = 3 - PEPIÑOXXL case_AC = 4 - Comercial AC case_AC = 5 - WIG case_AC = 6 - CERVERA case_AC = 7 - Existing AC - B747       % case_AC = 8 - Tamiz
        prompt = sprintf('ENTER AICRAFT SELECTION:\n1 - EMERGENTIA 100%% SCALE\n2 - EMERGENTIA MANUFACTURED%% SCALE\n3 - PEPIÑOXXL\n4 - Comertial AC\n5 - WIG\n6 - CERVERA\n7 - Existing AC,\n8 - Tamiz\n9 - EMERGENTIA Wind Tunnel (wing no sweep)\n10 - EMERGENTIA Wind Tunnel (wing with sweep)\n11 - EMERGENTIA Manufactured\n12 - ALO INTA\n13 - ALO INTA with Fuel Cell');
        dlgtitle = 'Aircraft Design Selection';
        dims = [1 35];
        definput = {'1'};
        opts.Interpreter = 'tex';
        answer_CASE_AC = inputdlg(prompt,dlgtitle,dims,definput);
        case_AC = str2num(answer_CASE_AC{1});
        
        %% Pre-load variables stored in Excel
        % Function that reads all the data for the different aircraft
        prompt2 = sprintf('Would you like to use:\n Stored Aircraft Data from a previous read of "SANAID_AIRCRAFT_DATA.xlsx" (Press Yes)\n or \n Read Excel file "SANAID_AIRCRAFT_DATA.xlsx" again (Press No)?');
        answer_stored = questdlg(prompt2,'Stored Aircraft Data','Yes','No','Yes');
        % Handle response
        switch answer_stored
            case 'Yes'
                switch case_AC
                    case 1 % case_AC = 1 - EMERGENTIA 1:1
                        load OUTPUT_read_XLSX_EMERGENTIA_100.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_EMERGENTIA;
                    case 2 % case_AC = 2 - EMERGENTIA 1:2
                        load OUTPUT_read_XLSX_EMERGENTIA_Mft.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_EMERGENTIA;
                    case 3 % case_AC = 3 - PEPIÑO XXL
                        load OUTPUT_read_XLSX_PEPINO.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_PEPINO;
                    case 4 % Commercial example
                        load OUTPUT_read_XLSX_Comertial.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_Comertial;
                    case 5 % WIG Aircraft
                        load OUTPUT_read_XLSX_WIG.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_WIG;
                    case 6 % case_AC = 6 CERVERA
                        load OUTPUT_read_XLSX_CERVERA.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_CERVERA;
                    case 7 % case_AC = 7 - Existing AC - B747
                        load OUTPUT_read_XLSX_EAC.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_EAC;
                    case 8 % case_AC = 8 - E26 Tamiz
                        load OUTPUT_read_XLSX_Tamiz.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_Tamiz;
                    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model w2 (no sweep)
                        load OUTPUT_read_XLSX_EMERGENTIA_WT2.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_EMERGENTIA_WT2;
                    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model w1 (sweep)
                        load OUTPUT_read_XLSX_EMERGENTIA_WT1.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_EMERGENTIA_WT1;
                    case 11 % case_AC = 11 - EMERGENTIA Manufactured
                        load OUTPUT_read_XLSX_EMERGENTIA_MFT.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_EMERGENTIA_MFT;
                    case 12 % case_AC = 12 - ALO
                        load OUTPUT_read_XLSX_EMERGENTIA_ALO.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_EMERGENTIA_ALO;
                    case 13 % case_AC = 13 - ALO Fuel Cell
                        load OUTPUT_read_XLSX_EMERGENTIA_ALO_FC.mat
                        OUTPUT_read_XLSX = OUTPUT_read_XLSX_EMERGENTIA_ALO_FC;
                end
            case 'No'
                switch case_AC
                    case 1 % case_AC = 1 - EMERGENTIA 1:1
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_EMERGENTIA = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_EMERGENTIA_100.mat
                    case 2 % case_AC = 2 - EMERGENTIA Manufactured
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_EMERGENTIA = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_EMERGENTIA_Mft.mat
                    case 3 % case_AC = 3 - PEPIÑO XXL
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_PEPINO = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_PEPINO.mat
                    case 4 % Commercial example
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_Comertial = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_Comertial.mat
                    case 5 % WIG Aircraft
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_WIG = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_WIG.mat
                    case 6 % case_AC = 6 CERVERA
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_CERVERA = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_CERVERA.mat
                    case 7 % case_AC = 7 - Existing AC - B747
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_EAC = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_EAC.mat
                    case 8 % case_AC = 8 - E26 Tamiz
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_Tamiz = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_Tamiz.mat
                    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - (no sweep)
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_EMERGENTIA_WT2 = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_EMERGENTIA_WT2.mat
                    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - (sweep)
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_EMERGENTIA_WT1 = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_EMERGENTIA_WT1.mat
                    case 11 % case_AC = 11 - EMERGENTIA Manufactured
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_EMERGENTIA_MFT = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_EMERGENTIA_MFT.mat
                    case 12 % case_AC = 12 - ALO
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_EMERGENTIA_ALO = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_EMERGENTIA_ALO.mat
                    case 13 % case_AC = 13 - ALO Fuel Cell
                        OUTPUT_read_XLSX = read_DATA_XLSX(case_AC);
                        OUTPUT_read_XLSX_EMERGENTIA_ALO_FC = OUTPUT_read_XLSX;
                        save OUTPUT_read_XLSX_EMERGENTIA_ALO_FC.mat
                end
        end
end

