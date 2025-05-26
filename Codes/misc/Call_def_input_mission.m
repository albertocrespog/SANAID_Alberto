function def_input_mission = Call_definput_mission
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
def_input_mission = {num_missions_str,missions_str};
