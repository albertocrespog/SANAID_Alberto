function [answer_mission2, segment] = Call_prompt_mission_t2(i,type_mission)

prompt1 = sprintf('Enter the initial altitude (m):');
prompt2 = sprintf('Enter the final altitude (m):');
prompt3 = sprintf('Enter the climb speed (m/s):');
prompt4 = sprintf('Enter the flight ascent path angle (deg):');
prompt_mission2 = {prompt1,prompt2,prompt3,prompt4};
dlgtitle_mission2 = 'Segment Type 2 - Climb';
dims = [1 60];
opts.Interpreter = 'tex';
% Examples of initialization
h_initial_str = '0';
h_final_str = '0';
V_cl_str = '22';
gamma_cl_str = '3';
definput_mission2 = {h_initial_str,h_final_str,V_cl_str,gamma_cl_str};
answer_mission2 = inputdlg(prompt_mission2,dlgtitle_mission2,dims,definput_mission2);
segment{i}.data.mision = type_mission;
segment{i}.data.h_initial = str2num(answer_mission2{1});
segment{i}.data.h_final = str2num(answer_mission2{2});
segment{i}.data.V_cl = str2num(answer_mission2{3});
segment{i}.data.gamma_cl = str2num(answer_mission2{4});
[Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_cl);
segment{i}.data.V_V = segment{i}.data.V_cl*sin(segment{i}.data.gamma_cl);
segment{i}.data.Data_ATM = Data_ATM;
segment{i}.data.Performance = Performance;