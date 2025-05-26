function [answer_mission5, segment] = Call_prompt_mission_t5(i,type_mission)

prompt1 = sprintf('Enter the initial altitude (m):');
prompt2 = sprintf('Enter the final altitude (m):');
prompt3 = sprintf('Enter the cruise speed (m/s):');
prompt4 = sprintf('Enter the Endurance Time (min):');
prompt_mission5 = {prompt1,prompt2,prompt3,prompt4};
dlgtitle_mission5 = 'Segment Type 5 - Cruise Endurance';
dims = [1 60];
opts.Interpreter = 'tex';
% Examples of initialization
h_initial_str = '0';
h_final_str = '0';
V_cr_str = '5';
Endurance_str = '10';
definput_mission5 = {h_initial_str,h_final_str,V_cr_str,Endurance_str};
answer_mission5 = inputdlg(prompt_mission5,dlgtitle_mission5,dims,definput_mission5);
segment{i}.data.mision = type_mission;
segment{i}.data.h_initial = str2num(answer_mission5{1});
segment{i}.data.h_final = str2num(answer_mission5{2});
segment{i}.data.V_cr = str2num(answer_mission5{3});
segment{i}.data.Endurance = str2num(answer_mission5{4});
[Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_cr);
segment{i}.data.Data_ATM = Data_ATM;
segment{i}.data.Performance = Performance;
