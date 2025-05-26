function answer_mission = Call_prompt_mission
% Asks user for the input data        
% Enter type of mission segments betwee brackets being
% - 1 TakeOff
% - 2 = Climb
% - 3 = VTOL Climb
% - 4 = Cruise (Range)
% - 5 = Cruise (Endurance)
% - 6 = Descent
% - 7 = Descent (VTOL)

prompt1 = sprintf('Enter the number of mission segments:');
prompt2 = sprintf('Enter type of mission segments betwee brackets being \n TakeOff = 1 \n Climb = 2 \n VTOL Climb = 3 \n Cruise (Range) = 4 \n Cruise (Endurance) = 5 \n Descent = 6 \n Descent (VTOL) = 7');
prompt_mission = {prompt1,prompt2};
dlgtitle = 'Input # Missions';
dims = [1 50];
opts.Interpreter = 'tex';
% Examples of initialization
num_missions = 2; 
num_missions_str = num2str(num_missions);
missions_str = '[3 4]';
definput_mission = {num_missions_str,missions_str};
answer_mission = inputdlg(prompt_mission,dlgtitle,dims,definput_mission);

% Asks user for the input data        
% Introduce data for the different segments
N_missions = str2num(answer_mission{1});
type_missions = str2num(answer_mission{2});
for i=1:N_missions
    type_mission = type_missions(i);
    switch type_mission
        case 1 % case 1 TakeOff
            mission_tex = 'TakeOff';
        case 2 % case 2 Climb
            mission_tex = 'Climb';
        case 3 % case 3 VTOL Climb
            mission_tex = 'VTOL Climb';
        case 4 % case 4 Cruise (Range)
            mission_tex = 'Cruise (Range)';
        case 5 % case 5 Cruise (Endurance)
            mission_tex = 'Cruise (Endurance)';
        case 6 % case 6 Descent
            mission_tex = 'Descent';
        case 7 % case 7 Descent (VTOL)
            mission_tex = 'Descent (VTOL)';
    end
    
    uiwait(msgbox({'Enter Properties for Mission ',num2str(type_mission),mission_tex},'input Mission Properties','modal'));
    
    switch type_mission
        case 1 % case 1 TakeOff
%             prompt1 = sprintf('Enter the initial altitude (m):');
%             prompt2 = sprintf('Enter the final altitude (m):');
%             prompt3 = sprintf('Enter Take Off Speed (m/s) (It is Just an Estimae):');
%             prompt_mission1 = {prompt1,prompt2,prompt3};
%             dlgtitle_mission1 = 'Segment Type 1 - Take off';
%             dims = [1 60];
%             opts.Interpreter = 'tex';
%             % Examples of initialization
%             h_initial_str = '0';
%             h_final_str = '0';
%             V_TO_str = '15';
%             definput_mission1 = {h_initial_str,h_final_str,V_TO_str};
%             answer_mission1 = inputdlg(prompt_mission1,dlgtitle_mission1,dims,definput_mission1);
%             segment{i}.data.mision = type_mission; 
%             segment{i}.data.h_initial = str2num(answer_mission1{1});
%             segment{i}.data.h_final = str2num(answer_mission1{2});
%             segment{i}.data.V_TO = str2num(answer_mission1{3});
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_TO);
%             segment{i}.data.Data_ATM = Data_ATM;
%             segment{i}.data.Performance = Performance;
            [answer_mission1, segment] = Call_prompt_mission_t1(i,type_mission);
        case 2 % case 2 Climb
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
            [answer_mission1, segment] = Call_prompt_mission_t1(i,type_mission);
        case 3 % case 3 VTOL Climb
%             prompt1 = sprintf('Enter the initial altitude (m) - VTOL Flight:');
%             prompt2 = sprintf('Enter the final altitude (m) - VTOL Flight:');
%             prompt3 = sprintf('Enter the climb speed (m/s) - VTOL Flight:');
%             prompt_mission3 = {prompt1,prompt2,prompt3};
%             dlgtitle_mission3 = 'Segment Type 3 - VTOL Climb';
%             dims = [1 60];
%             opts.Interpreter = 'tex';
%             % Examples of initialization
%             h_initial_str = '0';
%             h_final_str = '200';
%             V_VTOL_str = '5';
%             definput_mission3 = {h_initial_str,h_final_str,V_VTOL_str};
%             answer_mission3 = inputdlg(prompt_mission3,dlgtitle_mission3,dims,definput_mission3);
%             segment{i}.data.mision = type_mission; 
%             segment{i}.data.h_initial = str2num(answer_mission3{1});
%             segment{i}.data.h_final = str2num(answer_mission3{2});
%             segment{i}.data.V_VTOL = str2num(answer_mission3{3});
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_VTOL);
%             segment{i}.data.Data_ATM = Data_ATM;
%             segment{i}.data.Performance = Performance;
            [answer_mission3, segment] = Call_prompt_mission_t3(i,type_mission);
        case 4 % case 4 Cruise(Range)
%             prompt1 = sprintf('Enter the initial altitude (m):');
%             prompt2 = sprintf('Enter the final altitude (m):');
%             prompt3 = sprintf('Enter the cruise speed (m/s):');
%             prompt4 = sprintf('Enter the Range (Km):');
%             prompt_mission4 = {prompt1,prompt2,prompt3,prompt4};
%             dlgtitle_mission4 = 'Segment Type 4 - Cruise Range';
%             dims = [1 60];
%             opts.Interpreter = 'tex';
%             % Examples of initialization
%             h_initial_str = '200';
%             h_final_str = '200';
%             V_cr_str = '26';
%             Range_str = '20';
%             definput_mission4 = {h_initial_str,h_final_str,V_cr_str,Range_str};
%             answer_mission4 = inputdlg(prompt_mission4,dlgtitle_mission4,dims,definput_mission4);
%             segment{i}.data.mision = type_mission; 
%             segment{i}.data.h_initial = str2num(answer_mission4{1});
%             segment{i}.data.h_final = str2num(answer_mission4{2});
%             segment{i}.data.V_cr = str2num(answer_mission4{3});
%             segment{i}.data.Range = str2num(answer_mission4{4});
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_cr);
%             segment{i}.data.Data_ATM = Data_ATM;
%             segment{i}.data.Performance = Performance;
            [answer_mission4, segment] = Call_prompt_mission_t4(i,type_mission);
        case 5 % case 5 Cruise (Endurance)
%             prompt1 = sprintf('Enter the initial altitude (m):');
%             prompt2 = sprintf('Enter the final altitude (m):');
%             prompt3 = sprintf('Enter the cruise speed (m/s):');
%             prompt4 = sprintf('Enter the Endurance Time (min):');
%             prompt_mission5 = {prompt1,prompt2,prompt3,prompt4};
%             dlgtitle_mission5 = 'Segment Type 5 - Cruise Endurance';
%             dims = [1 60];
%             opts.Interpreter = 'tex';
%             % Examples of initialization
%             h_initial_str = '0';
%             h_final_str = '0';
%             V_cr_str = '5';
%             Endurance_str = '10';
%             definput_mission5 = {h_initial_str,h_final_str,V_cr_str,Endurance_str};
%             answer_mission5 = inputdlg(prompt_mission5,dlgtitle_mission5,dims,definput_mission5);
%             segment{i}.data.mision = type_mission; 
%             segment{i}.data.h_initial = str2num(answer_mission5{1});
%             segment{i}.data.h_final = str2num(answer_mission5{2});
%             segment{i}.data.V_cr = str2num(answer_mission5{3});
%             segment{i}.data.Endurance = str2num(answer_mission5{4});
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_cr);
%             segment{i}.data.Data_ATM = Data_ATM;
%             segment{i}.data.Performance = Performance;
            [answer_mission5, segment] = Call_prompt_mission_t5(i,type_mission);
        case 6 % case 6 Descent
%             prompt1 = sprintf('Enter the initial altitude (m):');
%             prompt2 = sprintf('Enter the final altitude (m):');
%             prompt3 = sprintf('Enter the descent speed (m/s):');
%             prompt4 = sprintf('Enter the flight path angle (min):');
%             prompt_mission6 = {prompt1,prompt2,prompt3,prompt4};
%             dlgtitle_mission6 = 'Segment Type 6 - Descent';
%             dims = [1 60];
%             opts.Interpreter = 'tex';
%             % Examples of initialization
%             h_initial_str = '0';
%             h_final_str = '0';
%             V_dsc_str = '5';
%             gamma_dsc_str = '5';
%             definput_mission6 = {h_initial_str,h_final_str,V_dsc_str,gamma_dsc_str};
%             answer_mission6 = inputdlg(prompt_mission6,dlgtitle_mission6,dims,definput_mission6);
%             segment{i}.data.mision = type_mission; 
%             segment{i}.data.h_initial = str2num(answer_mission6{1});
%             segment{i}.data.h_final = str2num(answer_mission6{2});
%             segment{i}.data.V_dsc = str2num(answer_mission6{3});
%             segment{i}.data.gamma_dsc = str2num(answer_mission6{4});
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_dsc);
%             segment{i}.data.V_V = segment{i}.data.V_dsc*sin(segment{i}.data.gamma_dsc);
%             segment{i}.data.Data_ATM = Data_ATM;
%             segment{i}.data.Performance = Performance;
            [answer_mission6, segment] = Call_prompt_mission_t6(i,type_mission);
        case 7 % case 7 Descent (VTOL)
%             prompt1 = sprintf('Enter the initial altitude (m):');
%             prompt2 = sprintf('Enter the final altitude (m):');
%             prompt3 = sprintf('Enter the descent speed (m/s) - VTOL:');
%             prompt_mission7 = {prompt1,prompt2,prompt3};
%             dlgtitle_mission7 = 'Segment Type 7 - Descent - VTOL';
%             dims = [1 60];
%             opts.Interpreter = 'tex';
%             % Examples of initialization
%             h_initial_str = '0';
%             h_final_str = '0';
%             V_VTOL_str = '-2';
%             definput_mission7 = {h_initial_str,h_final_str,V_VTOL_str};
%             answer_mission7 = inputdlg(prompt_mission7,dlgtitle_mission7,dims,definput_mission7);
%             segment{i}.data.mision = type_mission; 
%             segment{i}.data.h_initial = str2num(answer_mission5{1});
%             segment{i}.data.h_final = str2num(answer_mission5{2});
%             segment{i}.data.V_VTOL = str2num(answer_mission5{3});
%             [Data_ATM Performance] = Flight_Conditions_2020_v1(segment{i}.data.h_initial,segment{i}.data.V_VTOL);
%             segment{i}.data.Data_ATM = Data_ATM;
%             segment{i}.data.Performance = Performance;
            [answer_mission7, segment] = Call_prompt_mission_t7(i,type_mission);
    end
end

