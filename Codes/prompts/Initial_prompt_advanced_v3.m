function [case_AC OUTPUT_read_XLSX] = Initial_prompt_advanced_v3

%% Selects Automatization process

% Prompt aircraft Selection
% Defines de aircraft to be analyzed
% case_AC = 1 - EMERGENTIA 1:1 case_AC = 2 - EMERGENTIA 1:2 case_AC = 3 - PEPIÑOXXL case_AC = 4 - A400m 
% case_AC = 5 - WIG case_AC = 6 - MILVUS case_AC = 7 - Existing AC - B747  % case_AC = 8 - Tamiz
% case_AC = 9 - EMERGENTIA WINND TUNNEL case_AC = 10 - EMERGENTIA Wind Tunnel (wing with sweep) case_AC = 11 - EMERGENTIA Manufactured
% case_AC = 12 - ALO INTA case_AC = 13 - ALO INTA Fuel Cells  case_AC = 14 - A320-200  case_AC = 15 - EMERGENTIA tandem-wing
%  case_AC = 16 - TARSIS 75  case_AC = 17 - TARSIS 120 case_AC = 18 - BAT case_AC = 19 - EVOL1 case_AC = 20 - EVOL2 case_AC = 21 - EVOL3);

prompt = sprintf(['CIERRA EL EXCEL!!!\n \n ENTER AICRAFT SELECTION:\n1 - EMERGENTIA 100%% SCALE\n2 - EMERGENTIA MANUFACTURED%% SCALE...' ...
    '\n3 - PEPIÑOXXL\n4 - A400M \n5 - WIG\n6 - MILVUS\n7 - Existing AC,\n8 - Tamiz\n9 - EMERGENTIA Wind Tunnel (wing no sweep)...' ...
    '\n10 - EMERGENTIA Wind Tunnel (wing with sweep)\n11 - EMERGENTIA Manufactured\n12 - ALO INTA\n13 - ALO INTA with Fuel Cell...' ...
    '\n14 - A320-200\n15 - EMERGENTIA tandem-wing\n16 - TARSIS 75\n17 - TARSIS 120\n18 - BAT\n19 - FALCON2000\n20 - SOLARTII\n21 - VANTUS-24'...
    '\n22 - Cessna 208\n23 - Kink Air 350i \n24 - MK84\n25 - MK84B\n26 - HA\n27 - Fast Response\n28 - H2\n29 - MILVUS2\n30 - Future 7'...
    ]);
dlgtitle = 'Aircraft Design Selection';
dims = [1 35];
definput = {'1'};
opts.Interpreter = 'tex';
answer_CASE_AC = inputdlg(prompt,dlgtitle,dims,definput);
case_AC = str2num(answer_CASE_AC{1});

% Defines path where the MAT files are stored
pathname = fileparts ('Codes/ac_models');

OUTPUT_read_XLSX = read_DATA_XLSX_advanced_v3(case_AC);

       
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        % OUTPUT_read_XLSX_EMERGENTIA = OUTPUT_read_XLSX;
        % save('Codes/ac_models/01_OUTPUT_read_XLSX_EMERGENTIA_100.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/01_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 2 % case_AC = 2 - EMERGENTIA Manufactured
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        % OUTPUT_read_XLSX_EMERGENTIA = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_EMERGENTIA_Mft.mat
        % save('Codes/Codes/ac_models/02_OUTPUT_read_XLSX_EMERGENTIA_Mft.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/02_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 3 % case_AC = 3 - PEPIÑO XXL
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_PEPINO = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_PEPINO.mat
        % save('Codes/ac_models/03_OUTPUT_read_XLSX_PEPINO.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/03_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 4 % Commercial example
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_Comertial = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_A400M.mat
        % save('Codes/ac_models/04_OUTPUT_read_XLSX_A400M.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/04_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 5 % WIG Aircraft
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_WIG = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_WIG.mat
        % save('Codes/ac_models/05_OUTPUT_read_XLSX_WIG.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/05_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 6 % case_AC = 6 MILVUS
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_MILVUS = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_MILVUS.mat
        % save('Codes/ac_models/06_OUTPUT_read_XLSX_MILVUS.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/06_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 7 % case_AC = 7 - Existing AC - B747
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EAC = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_EAC.mat
        % save('Codes/ac_models/07_OUTPUT_read_XLSX_EAC.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/07_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 8 % case_AC = 8 - E26 Tamiz
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_Tamiz = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_Tamiz.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Tamiz.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/08_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - (no sweep)
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_WT2 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_EMERGENTIA_WT2.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_EMERGENTIA_WT2.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/09_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - (sweep)
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        % OUTPUT_read_XLSX_EMERGENTIA_WT1 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_EMERGENTIA_WT1.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_EMERGENTIA_WT1.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/10_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 11 % case_AC = 11 - EMERGENTIA Manufactured
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        % OUTPUT_read_XLSX_EMERGENTIA_MFT = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_EMERGENTIA_MFT.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_EMERGENTIA_MFT.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/11_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 12 % case_AC = 12 - ALO
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        % OUTPUT_read_XLSX_EMERGENTIA_ALO = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_ALO.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_ALO.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/12_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 13 % case_AC = 13 - ALO Fuel Cell
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        % OUTPUT_read_XLSX_EMERGENTIA_ALO_FC = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_ALO_FC.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_ALO_FC.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/13_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 14 % case_AC = 14 - A320-200
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        % OUTPUT_read_XLSX_EMERGENTIA_A320_FC = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_A320_200.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_A320_200.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/14_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 15 % case_AC = 15 - EMERGENTIA tandem-wing
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        % OUTPUT_read_XLSX_EMERGENTIA_TW_FC = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_EMERGENTIA_Tandem_wing.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_EMERGENTIA_Tandem_wing.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/15_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 16 % case_AC = 16 - Tarsis 75
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        % OUTPUT_read_XLSX_EMERGENTIA_T95_FC = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_Tarsis75.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Tarsis95.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/16_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 17 % case_AC = 17 - Tarsis 120
%         OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        % OUTPUT_read_XLSX_EMERGENTIA_T120 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_Tarsis120.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Tarsis120.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/17_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 18 % case_AC = 18 - BAT
        % OUTPUT_read_XLSX_BAT = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_BAT.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_BAT.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/18_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 19 % case_AC = 19 - FALCON2000
        % OUTPUT_read_XLSX_FALCON2000 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_FALCON2000.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_FALCON2000.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/19_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 20 % case_AC = 20 - SOLARTII
        % OUTPUT_read_XLSX_SOLARTII = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_EVOL2.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_SOLARTII.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/20_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 21 % case_AC = 21 - VANTUS24
        % OUTPUT_read_XLSX_VANTUS24 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_EVOL3.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_VANTUS24.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/21_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 22 % case_AC = 22
        % OUTPUT_read_XLSX_Cessna208 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_EVOL3.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Cessna208.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/22_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 23 % case_AC = 23
        % OUTPUT_read_XLSX_KingAir350 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_KingAir350.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_KingAir350.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/23_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 24 % case_AC = 24
        % OUTPUT_read_XLSX_Future1 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_Future1.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_MK84.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/24_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 25 % case_AC = 25
        % OUTPUT_read_XLSX_Future1 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_Future2.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Future2.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/25_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 26 % case_AC = 26
        % OUTPUT_read_XLSX_Future1 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_Future3.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Future3.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/26_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 27 % case_AC = 27
        % OUTPUT_read_XLSX_Future1 = OUTPUT_read_XLSX;
        % %         save OUTPUT_read_XLSX_Future4.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Future4.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/27_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 28 % case_AC = 28
        % OUTPUT_read_XLSX_Future1 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_Future5.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Future5.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/28_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 29 % case_AC = 29
        % OUTPUT_read_XLSX_Future1 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_Future6.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Future6.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/29_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
    case 30 % case_AC = 30
        % OUTPUT_read_XLSX_Future1 = OUTPUT_read_XLSX;
        %         save OUTPUT_read_XLSX_Future7.mat
        % save('Codes/ac_models/OUTPUT_read_XLSX_Future8.mat', 'OUTPUT_read_XLSX')
        save('Codes/ac_models/30_OUTPUT_read_XLSX.mat', 'OUTPUT_read_XLSX')
end

% Prompt for performance studies and getting inputs from users
[OUTPUT_read_XLSX] = performance_prompt_advanced_v2(OUTPUT_read_XLSX);
