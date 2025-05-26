function [answer_mission1, segment] = Call_prompt_mission_t1(i,type_mission)

% case 1 TakeOff

prompt1 = sprintf('Enter the initial altitude (m):');
prompt2 = sprintf('Enter the final altitude (m):');
prompt3 = sprintf('Enter Take Off Speed (m/s) (It is Just an Estimae):');
prompt_mission1 = {prompt1,prompt2,prompt3};
dlgtitle_mission1 = 'Segment Type 1 - Take off';
dims = [1 60];
opts.Interpreter = 'tex';
% Examples of initialization
h_initial_str = '0';
h_final_str = '0';
V_TO_str = '15';
definput_mission1 = {h_initial_str,h_final_str,V_TO_str};
answer_mission1 = inputdlg(prompt_mission1,dlgtitle_mission1,dims,definput_mission1);
segment{i}.data.mision = type_mission;
segment{i}.data.h_initial = str2num(answer_mission1{1});
segment{i}.data.h_final = str2num(answer_mission1{2});
segment{i}.data.V_TO = str2num(answer_mission1{3});
[Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_TO);
segment{i}.data.Data_ATM = Data_ATM;
segment{i}.data.Performance = Performance;