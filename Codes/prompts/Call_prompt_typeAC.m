function answer_typeAC = Call_prompt_typeAC

% Defines de aircraft to be analyzed
% case_AC = 1 - EMERGENTIA 1:1
% case_AC = 2 - EMERGENTIA 1:2
% case_AC = 3 - PEPIÑOXXL
% case_AC = 1;
% Function call Type AC

prompt_typeAC = sprintf('Enter the Type of Aircraft \n CASE 1: EMERGENTIA \n CASE 2: EMERGENTIA scaled \n CASE 3: PEPIÑOXXL');
dlgtitle_typeAC = 'Input Aircraft Type';
dims = [1 35];
definput_typeAC = {'1'};
opts.Interpreter = 'tex';
answer_typeAC = inputdlg(prompt_typeAC,dlgtitle_typeAC,dims,definput_typeAC);
