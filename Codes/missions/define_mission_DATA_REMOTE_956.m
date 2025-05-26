% Defines the geometru changes associated to mission and geometry changes
function missions = define_mission_DATA(case_AC,OUTPUT_read_XLSX)

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        missions = [1];
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        missions = [1];
    case 3 % case_AC = 3 - PEPIÑO XXL
        missions = [1];
    case 4 % case_AC = 4 - A400
        missions = [1];
    case 5 % case_AC = 5 - WIGL
        missions = [1];
    case 6 % case_AC = 6 - MILVUS
        missions = [1];
    case 7 % Existing Aircraft = 7
        missions = [1];
    case 8 % case_AC = 8 - TAMIZ
        missions = [1];
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
        missions = [1];
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
        missions = [1];
    case 11 % case_AC = 11 - EMERGENTIA Manufacturing
        missions = [1];
    case 12 % case_AC = 12 - ALO
        missions = [1];
    case 13 % case_AC = 13 - ALO Fuel Cell
        missions = [1];
    case 14 % case_AC = 14 - A320
        missions = [0];
%         missions = [0 1 2 3 4 5 6 7 8 9];
%         missions = [0];
    case 15 % case_AC = 15 - EMERGENTIA Tandem Wing
        missions = [1];
    case 16 % case_AC = 16 - Tarsis T75
        missions = [1];
    case 17 % case_AC = 17 - Tarsis T120
%         missions = [0 1 2 3 4 5 6 7 8 9];
        missions = [5 6 7];
end