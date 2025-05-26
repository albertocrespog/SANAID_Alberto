function [answer_mission7, segment] = Call_prompt_mission_t7(i,type_mission);

prompt1 = sprintf('Enter the initial altitude (m):');
prompt2 = sprintf('Enter the final altitude (m):');
prompt3 = sprintf('Enter the descent speed (m/s) - VTOL:');
prompt_mission7 = {prompt1,prompt2,prompt3};
dlgtitle_mission7 = 'Segment Type 7 - Descent - VTOL';
dims = [1 60];
opts.Interpreter = 'tex';
% Examples of initialization
h_initial_str = '0';
h_final_str = '0';
V_VTOL_str = '-2';
definput_mission7 = {h_initial_str,h_final_str,V_VTOL_str};
answer_mission7 = inputdlg(prompt_mission7,dlgtitle_mission7,dims,definput_mission7);
segment{i}.data.mision = type_mission;
segment{i}.data.h_initial = str2num(answer_mission7{1});
segment{i}.data.h_final = str2num(answer_mission7{2});
segment{i}.data.V_VTOL = str2num(answer_mission7{3});
[Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_VTOL);
segment{i}.data.Data_ATM = Data_ATM;
segment{i}.data.Performance = Performance;
