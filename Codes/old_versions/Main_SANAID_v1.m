%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main SANAID - Stability ANAlysis Interactive Design Tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

% Defines Generic Path so that files can be organized in folders
% change for the different versions
addpath(genpath('../../MATLAB/src'))
addpath(genpath('../../MATLAB/XFLR5'))

%% Lets user select if wants to conduct study for the propellers
answer = questdlg('Would you like execute with message dialog boxes?', ...
	'Message Dialog Boxes', ...
	'Yes','No','No');
% Handle response
switch answer
    case 'Yes'
        SANAID_prompt_user;
    case 'No'
        SANAID_command_line;
end