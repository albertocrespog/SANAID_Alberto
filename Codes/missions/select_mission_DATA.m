% Defines the geometru changes associated to mission and geometry changes
function [conditions OUTPUT_read_XLSX] = select_mission_DATA(case_AC,OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        % Selects mission according to Stability of Performance variation
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
        [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
        end
        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

    case 2 % case_AC = 2 - EMERGENTIA 1:2
        [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
    case 3 % case_AC = 3 - PEPIÑO XXL
        % Selects mission according to Stability of Performance variation
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_PEPINO_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end
        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_PEPINO_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end
    case 4 % case_AC = 4 - A400
        [conditions OUTPUT_read_XLSX] = mission_changes_Pepino(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);

    case 5 % case_AC = 5 - WIGL

    case 6 % case_AC = 6 - MILVUS
        % Selects mission according to Stability of Performance variation
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_MILVUS_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end
        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_MILVUS_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

    case 7 % Existing Aircraft = 7

    case 8 % case_AC = 8 - TAMIZ

    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
        [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
        [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
    case 11 % case_AC = 11 - EMERGENTIA Manufacturing
        %         [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
        % [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_vee_dihedral_var(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
        % [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_vee_inc_var(OUTPUT_read_XLSX,mission_actual,conv_UNITS);

        %Vtail geometry variation cases
        [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_vee_tails_case1(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
        %               [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_vee_tails_case2a(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
        %               [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_vee_tails_case2b(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
    case 12 % case_AC = 12 - ALO
        [conditions OUTPUT_read_XLSX] = mission_changes_ALO_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);

    case 13 % case_AC = 13 - ALO Fuel Cell
        [conditions OUTPUT_read_XLSX] = mission_changes_ALO_Fuel_Cell_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);

    case 14 % case_AC = 14 - A320
        % [conditions OUTPUT_read_XLSX] = mission_changes_A320_tails_v0(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
        [conditions OUTPUT_read_XLSX] = mission_changes_A320_Performance_study_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
    case 15 % case_AC = 15 - EMERGENTIA Tandem Wing

    case 16 % case_AC = 16 - Tarsis T75
        [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS(OUTPUT_read_XLSX,mission_actual,conv_UNITS);
    case 17 % case_AC = 17 - Tarsis T120
        %          [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        %          [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS_CASE5(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        %          [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS_CASE5_varmissiles(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        % Final results T120
        %           [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS_CASE5_Final_Study_Dic2022(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        % Final results T120 - Fase 2
        %            [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS_CASE5_Final_Study_July2023(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        % Studies TARSIS KSA September 2023
        %            [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS_KSA_June2023(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
        % Studies TARSIS KSA November  2023
        %           [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS_KSA_Nov2023(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
        % Studies TARSIS KSA December 2023
        [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS_KSA_Dic2023(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

    case 18 % case_AC = 18 - BAT
        [conditions OUTPUT_read_XLSX] = mission_changes_BAT_July2023(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
    case 19 % case_AC = 19 - FALCON2000
        [conditions OUTPUT_read_XLSX] = mission_changes_FALCON2000(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
    case 20 % case_AC = 20 - SOLARTII
        % [conditions OUTPUT_read_XLSX] = mission_changes_SOLARTII(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        % [conditions OUTPUT_read_XLSX] = mission_changes_SOLARTII_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        %[conditions OUTPUT_read_XLSX] = mission_changes_SOLARTII_v2(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        [conditions OUTPUT_read_XLSX] = mission_changes_SOLARTII_v3(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
    case 21 % case_AC = 21 - VANTUS
        % [conditions OUTPUT_read_XLSX] = mission_changes_VANTUS(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        [conditions OUTPUT_read_XLSX] = mission_changes_VANTUS_24_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);

    case 22 % case_AC = 22 - Cessna 208
        % [conditions OUTPUT_read_XLSX] = mission_changes_Cessna208_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        % Second set of mission checks with the max allowable aileron
        % stick and rudder pedal forces as inputs to solve the 5 equations
        % Selects mission according to Stability of Performance variation
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_CESSNA208_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);

            mission_changes_CESSNA208_STABILITY_v1
        end
        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_CESSNA208_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end


    case 23 % case_AC = 23 - King Air 350
        % Selects mission according to Stability of Performance variation
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_KingAir_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end
        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_KingAir_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

    case 24 % case_AC = 24 - MK84
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_MK84_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end
        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_MK84_PERFORMANCE_v2(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

    case 25 % case_AC = 25 - Future 1
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_MK84B_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end
        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_MK84B_PERFORMANCE_v2(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

    case 26 % case_AC = 26 - Hunter Aero
        % Selects mission according to Stability of Performance variation
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_HA_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            % [conditions OUTPUT_read_XLSX] = mission_changes_HA_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
            [conditions OUTPUT_read_XLSX] = mission_changes_HA_PERFORMANCE_v2(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

    case 27 % case_AC = Fast Response - Future 1
        % Selects mission according to Stability of Performance variation
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_FR_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_FR_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

    case 28 % case_AC = 24 - Future 1
        % Selects mission according to Stability of Performance variation
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_H2_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end
        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_H2_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

    case 29 % case_AC = 24 - Future 1
        % Selects mission according to Stability of Performance variation
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_MILVUS2_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end
        if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
            [conditions OUTPUT_read_XLSX] = mission_changes_MILVUS2_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
        end

    case 30 % case_AC = 24 - Future 1
        [conditions OUTPUT_read_XLSX] = mission_changes_Future7_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
end