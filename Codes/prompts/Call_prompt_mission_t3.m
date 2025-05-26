function [answer_mission3, segment] = Call_prompt_mission_t3(i,type_mission)

prompt1 = sprintf('Enter the initial altitude (m) - VTOL Flight:');
prompt2 = sprintf('Enter the final altitude (m) - VTOL Flight:');
prompt3 = sprintf('Enter the climb speed (m/s) - VTOL Flight:');
prompt_mission3 = {prompt1,prompt2,prompt3};
dlgtitle_mission3 = 'Segment Type 3 - VTOL Climb';
dims = [1 60];
opts.Interpreter = 'tex';
% Examples of initialization
h_initial_str = '0';
h_final_str = '200';
V_VTOL_str = '5';
definput_mission3 = {h_initial_str,h_final_str,V_VTOL_str};
answer_mission3 = inputdlg(prompt_mission3,dlgtitle_mission3,dims,definput_mission3);
segment{i}.data.mision = type_mission;
segment{i}.data.h_initial = str2num(answer_mission3{1});
segment{i}.data.h_final = str2num(answer_mission3{2});
segment{i}.data.V_VTOL = str2num(answer_mission3{3});
[Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_VTOL);
segment{i}.data.Data_ATM = Data_ATM;
segment{i}.data.Performance = Performance;