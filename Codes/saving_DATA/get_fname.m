function prefixa = get_fname(OUTPUT_read_XLSX)

case_AC = OUTPUT_read_XLSX.AC_Data_flags.case_AC;
switch case_AC
    case 1
        prefixa = strcat('Results\01_EMERGENTIA');
    case 2
        prefixa = strcat('Results\01_EMERGENTIA');
    case 3
        prefixa = strcat('Results\03_PEPINO_XXXL');
    case 4
        prefixa = strcat('Results\04_A400');
    case 5
        prefixa = strcat('Results\05_MISC');
    case 6
        prefixa = strcat('Results\06_MILVUS');
    case 7
        prefixa = strcat('Results\05_MISC');
    case 8
        prefixa = strcat('Results\08_E26TAMIZ');
    case 9
        prefixa = strcat('Results\01_EMERGENTIA_WT2');
    case 10
        prefixa = strcat('Results\01_EMERGENTIA_WT1');
    case 11
        prefixa = strcat('Results\11_EMERGENTIA_Manufactured');
    case 12
        prefixa = strcat('Results\12_ALO');
    case 13
        prefixa = strcat('Results\13_ALO_Fuel_Cell');
    case 14
        prefixa = strcat('Results\14_ONEiRE');
    case 15
        prefixa = strcat('Results\15_EMERGENTIA_v2');
    case 16
        prefixa = strcat('Results\16_TARSIS_75');
    case 17
        prefixa = strcat('Results\17_TARSIS_120');
    case 18
        prefixa = strcat('Results\18_BAT');
    case 19
        prefixa = strcat('Results\19_FALCON2000');
    case 20
        prefixa = strcat('Results\20_SOLARTII');
    case 21
        prefixa = strcat('Results\21_VANTUS');
    case 22
        prefixa = strcat('Results\22_Cessna208');
    case 23
        prefixa = strcat('Results\23_KingAir350');
    case 24
        prefixa = strcat('Results\24_MK84');
    case 25
        prefixa = strcat('Results\25_MK84B');
    case 26
        prefixa = strcat('Results\26_HA');
    case 27
        prefixa = strcat('Results\27_FR');
    case 28
        prefixa = strcat('Results\28_Future5');
    case 29
        prefixa = strcat('Results\29_Future6');
    case 30
        prefixa = strcat('Results\30_Future7');
end
