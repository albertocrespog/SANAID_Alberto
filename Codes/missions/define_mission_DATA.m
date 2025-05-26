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
        % Stability
        % missions = [13 14];
        missions = [1];
        % Performence
        % missions = [1 2 3 4];
    case 7 % Existing Aircraft = 7
        missions = [1];
    case 8 % case_AC = 8 - TAMIZ
        missions = [1];
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
        missions = [1];
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
        missions = [1];
    case 11 % case_AC = 11 - EMERGENTIA Manufacturing
        %% Modificación temporal
        %         var_i_w1_emerg = [0:3];
        %         var_i_vee_emerg = [0:3];
        %         num_casos = length(var_i_w1_emerg)*length(var_i_vee_emerg);
        %         missions = [1:num_casos];
%         missions = [3];
         %missions = [1 2 3 4 5 ];
         missions = [5];
%         missions = [4 7 8 9];
%         missions = [1 4 9];
%         missions = [4];
% missions = [1 2 3 4 5 6 7 8 9];
        %         missions = [0];
    case 12 % case_AC = 12 - ALO
        missions = [1];
    case 13 % case_AC = 13 - ALO Fuel Cell
        missions = [1];
    case 14 % case_AC = 14 - A320
%         missions = [0];
        % missions = [0 1 2 3 4 5 6 7 8 9];
        missions = [1];
    case 15 % case_AC = 15 - EMERGENTIA Tandem Wing
        missions = [1];
    case 16 % case_AC = 16 - Tarsis T75
        missions = [1];
    case 17 % case_AC = 17 - Tarsis T120
%         missions = [0 1 2 3 4 5 6 7 8 9];
%          missions = [1];
%         missions = [5 6 7 8 9];
        % Var Missiles
%         missions = [1 2 3 4 5];
        % Case 5 Final Studies
%         missions = [1 2 3 4 5 6 7 8 9 10 11 12];
        % Case 5 Final Studies KSA
         % missions = [3];
        missions = [1 2 3 4 5 6 7 8 9];
    case 18 % case_AC = 18 - BAT
        missions = [1];
    case 19 % case_AC = 19 - FALCON2000
        missions = [1];
    case 20 % case_AC = 20 - SOLARTII
        % missions = [0];
        % missions = [3];
        missions = [2];
        % missions = [1 2 3 4 5 6 7 8 9 10 11];
        % missions = [1 2 3 4 5 6 7 8 9];
    case 21 % case_AC = 21 - VANTUS
        missions = [1];
    case 22 % case_AC = 22 - Cessna 208
        % missions = [1 2 3 ];
        % missions = [1 2 3 4 5 6];
        % missions = [1 2 3 4 5 6 7 8 9 10 11];
        missions = [1];
    case 23 % case_AC = 23 - King Air 350
        missions = [1];
    case 24 % case_AC = 24 - MK84
        missions = [1];
    case 25 % case_AC = 25 - MK84B
        missions = [3 1 2];
    case 26 % case_AC = 26 - HA - Hunter Aero
        missions = [2];
    case 27 % case_AC = 27 - FR - Fast Response
        missions = [1];
    case 28 % case_AC = 28 - Future 5
        missions = [1];
    case 29 % case_AC = 29 - MILVUS 2
        % missions = [1]
        missions = [1];
    case 30 % case_AC = 30 - Future 7
        missions = [1];
end