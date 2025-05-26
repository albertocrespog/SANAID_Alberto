function [answer_mission6, segment] = Call_prompt_mission_t6(i,type_mission);

prompt1 = sprintf('Enter the initial altitude (m):');
prompt2 = sprintf('Enter the final altitude (m):');
prompt3 = sprintf('Enter the descent speed (m/s):');
prompt4 = sprintf('Enter the flight path angle (min):');
prompt_mission6 = {prompt1,prompt2,prompt3,prompt4};
dlgtitle_mission6 = 'Segment Type 6 - Descent';
dims = [1 60];
opts.Interpreter = 'tex';
% Examples of initialization
h_initial_str = '0';
h_final_str = '0';
V_dsc_str = '5';
gamma_dsc_str = '5';
definput_mission6 = {h_initial_str,h_final_str,V_dsc_str,gamma_dsc_str};
answer_mission6 = inputdlg(prompt_mission6,dlgtitle_mission6,dims,definput_mission6);
segment{i}.data.mision = type_mission;
segment{i}.data.h_initial = str2num(answer_mission6{1});
segment{i}.data.h_final = str2num(answer_mission6{2});
segment{i}.data.V_dsc = str2num(answer_mission6{3});
segment{i}.data.gamma_dsc = str2num(answer_mission6{4});
[Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_dsc);
segment{i}.data.V_V = segment{i}.data.V_dsc*sin(segment{i}.data.gamma_dsc);
segment{i}.data.Data_ATM = Data_ATM;
segment{i}.data.Performance = Performance;