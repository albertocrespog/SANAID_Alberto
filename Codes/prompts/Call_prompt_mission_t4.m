function [answer_mission4, segment] = Call_prompt_mission_t4(i,type_mission);

prompt1 = sprintf('Enter the initial altitude (m):');
prompt2 = sprintf('Enter the final altitude (m):');
prompt3 = sprintf('Enter the cruise speed (m/s):');
prompt4 = sprintf('Enter the Range (Km):');
prompt_mission4 = {prompt1,prompt2,prompt3,prompt4};
dlgtitle_mission4 = 'Segment Type 4 - Cruise Range';
dims = [1 60];
opts.Interpreter = 'tex';
% Examples of initialization
h_initial_str = '200';
h_final_str = '200';
V_cr_str = '26';
Range_str = '20';
definput_mission4 = {h_initial_str,h_final_str,V_cr_str,Range_str};
answer_mission4 = inputdlg(prompt_mission4,dlgtitle_mission4,dims,definput_mission4);
segment{i}.data.mision = type_mission;
segment{i}.data.h_initial = str2num(answer_mission4{1});
segment{i}.data.h_final = str2num(answer_mission4{2});
segment{i}.data.V_cr = str2num(answer_mission4{3});
segment{i}.data.Range = str2num(answer_mission4{4});
[Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_cr);
segment{i}.data.Data_ATM = Data_ATM;
segment{i}.data.Performance = Performance;