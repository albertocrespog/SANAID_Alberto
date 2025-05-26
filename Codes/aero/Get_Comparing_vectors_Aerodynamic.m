function VECTOR_XFLR5 = Get_Comparing_vectors_Aerodynamic(case_AC,OUTPUT_read_XLSX)

compare_plot_aero = OUTPUT_read_XLSX.Aerodynamic_Data_flags.compare_plot_aero;
VECTOR_XFLR5.compare = compare_plot_aero;

switch compare_plot_aero
    case 1
        % First Surface
        prompt1 = sprintf('PLOTS THE DIFFERENT AERODYNAMIC SURFACES\nIntroduce the number of files as a succession of numbers separated by a blank without any commas or brackets\nENTER THE NUMBER OF FILES FOR THE FIRST SURFACE (WING PROBABLY)');
        dlgtitle = 'Compared Files # for 1st surface';
        dims = [1 70];
        answer_SEGMENTS = inputdlg(prompt1,dlgtitle,dims);
        VECTOR_XFLR5.v1 = str2num(answer_SEGMENTS{1});
    case 2
        % First Surface
        prompt1 = sprintf('PLOTS THE DIFFERENT AERODYNAMIC SURFACES\nIntroduce the number of files as a succession of numbers separated by a blank without any commas or brackets\nENTER THE NUMBER OF FILES FOR THE FIRST SURFACE (WING PROBABLY)');
        dlgtitle = 'Compared Files # for 1st surface';
        dims = [1 70];
        answer_SEGMENTS = inputdlg(prompt1,dlgtitle,dims);
        VECTOR_XFLR5.v1 = str2num(answer_SEGMENTS{1});
        % Second Surface
        prompt2 = sprintf('ENTER THE NUMBER OF FILES FOR THE SECOND SURFACE (WING2,HTP,VEE-TAIL,CANNARD...)');
        dlgtitle = 'Compared Files # for 2nd surface';
        dims = [1 70];
        answer_SEGMENTS = inputdlg(prompt2,dlgtitle,dims);
        VECTOR_XFLR5.v2 = str2num(answer_SEGMENTS{1});
    case 3
        
        % First Surface
        prompt1 = sprintf('PLOTS THE DIFFERENT AERODYNAMIC SURFACES\nIntroduce the number of files as a succession of numbers separated by a blank without any commas or brackets\nENTER THE NUMBER OF FILES FOR THE FIRST SURFACE (WING PROBABLY)');
        dlgtitle = 'Compared Files # for 1st surface';
        dims = [1 70];
        answer_SEGMENTS = inputdlg(prompt1,dlgtitle,dims);
        VECTOR_XFLR5.v1 = str2num(answer_SEGMENTS{1});
        % Second Surface
        prompt2 = sprintf('ENTER THE NUMBER OF FILES FOR THE SECOND SURFACE (WING2,HTP,VEE-TAIL,CANNARD...)');
        dlgtitle = 'Compared Files # for 2nd surface';
        dims = [1 70];
        answer_SEGMENTS = inputdlg(prompt2,dlgtitle,dims);
        VECTOR_XFLR5.v2 = str2num(answer_SEGMENTS{1});
        % Third Surface Surface
        prompt3 = sprintf('ENTER THE NUMBER OF FILES FOR THE THIRD SURFACE (VTP)');
        dlgtitle = 'Compared Files # for 3rd surface';
        dims = [1 70];
        answer_SEGMENTS = inputdlg(prompt3,dlgtitle,dims);
        VECTOR_XFLR5.v3 = str2num(answer_SEGMENTS{1});
end


% %% Generates the file for Aerodynamic Data
% switch case_AC
%     case 1 % case_AC = 1 - EMERGENTIA 1:1
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [14];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [15,16];
%                 VECTOR_XFLR5.v2 = [20,21,22,23,24,25];
%             case 3 % compare = 3 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [15];
%                 VECTOR_XFLR5.v2 = [16];
%                 VECTOR_XFLR5.v3 = [18];
%         end
%     case 2 % case_AC = 2 - EMERGENTIA 1:2
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [14];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [15,16];
%                 VECTOR_XFLR5.v2 = [20,21,22,23,24,25];
%             case 3 % compare = 3 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [15];
%                 VECTOR_XFLR5.v2 = [16];
%                 VECTOR_XFLR5.v3 = [18];
%         end
%     case 3 % case_AC = 3 - PEPIÑO XXL
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [1,2];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1];
%                 VECTOR_XFLR5.v2 = [2];
%             case 3 % compare = 3 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [];
%                 VECTOR_XFLR5.v2 = [];
%                 VECTOR_XFLR5.v3 = [];
%         end
%     case 4 % case_AC = 4 COMERCIAL
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [1,2];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1];
%                 VECTOR_XFLR5.v2 = [2];
%             case 3 % compare = 3 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [];
%                 VECTOR_XFLR5.v2 = [];
%                 VECTOR_XFLR5.v3 = [];
%         end
%     case 5 % case_AC = 5 WIG
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [1,2];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1];
%                 VECTOR_XFLR5.v2 = [2];
%             case 3 % compare = 3 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1];
%                 VECTOR_XFLR5.v2 = [2];
%                 VECTOR_XFLR5.v3 = [3];
%         end
%     case 6 % case_AC = 1 - CERVERA
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [14];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [15,16];
%                 VECTOR_XFLR5.v2 = [20,21,22,23,24,25];
%             case 3 % compare = 3 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [15];
%                 VECTOR_XFLR5.v2 = [16];
%                 VECTOR_XFLR5.v3 = [18];
%         end
%     case 11 % case_AC = 11 - EMERGENTIA Manufacturing
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [37];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [37];
%                 VECTOR_XFLR5.v2 = [34];
%             case 3 % compare = 3 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [34];
%                 VECTOR_XFLR5.v2 = [36];
%                 VECTOR_XFLR5.v3 = [37];
%         end
%     case 12 % case_AC = 1 - ALO
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [1];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1];
%                 VECTOR_XFLR5.v2 = [2];
%             case 3 % compare = 3 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1];
%                 VECTOR_XFLR5.v2 = [2];
%                 VECTOR_XFLR5.v3 = [1];
%         end
%     case 13 % case_AC = 13 - ALO Fuel CEll
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [1];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1];
%                 VECTOR_XFLR5.v2 = [2];
%             case 3 % compare = 3 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1];
%                 VECTOR_XFLR5.v2 = [2];
%                 VECTOR_XFLR5.v3 = [1];
%         end
%     case 14 % case_AC = 13 - ALO Fuel CEll
%         switch compare_plot_aero
%             case 1 % compare = 1 elements: w1,w2, vtp... alone
%                 % W1 cases to be eplotted
%                 VECTOR_XFLR5.v1 = [1,2,3];
%             case 2 % compare = 2 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1,2,3];
%                 VECTOR_XFLR5.v2 = [7,8,9,10,11];
%             case 3 % compare = 3 elements: w1,w2, &  vtp... alone
%                 VECTOR_XFLR5.v1 = [1,2,3];
%                 VECTOR_XFLR5.v2 = [7,8,9,10,11];
%                 VECTOR_XFLR5.v3 = [14,14]; % only for VTP
%             case 4 % compare = 4 elements: w1,w2, vtp... alone
%                 VECTOR_XFLR5.v1 = [1,2,3];
%                 VECTOR_XFLR5.v2 = [7,8,9,10,11];
%                 VECTOR_XFLR5.v3 = [14, 14];
%                 VECTOR_XFLR5.v4 = [14];
%         end
% end