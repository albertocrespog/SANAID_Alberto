%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANAID - Stability ANAlysis Interactive Design Tool
% Function that adds paths for the file to be used
% Date October 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_add_path

% Defines Generic Path so that files can be organized in folders
% change for the different versions
addpath(genpath('src2')) % source codes for stability analysis
addpath(genpath('AERO')) % aerodynamic codes from XFLR5
addpath(genpath('CAD')) % CAD files
addpath(genpath('Prop_data')) % Propdata
addpath(genpath('Fuselage')) % Fuselage models
addpath(genpath('AIRCRAFT')) % Selection of input data
addpath(genpath('igesToolBox')) % Selection of input data
addpath(genpath('Results')) % Selection of input data
% Organizing files
addpath(genpath('Codes')) % Aircraft Models
addpath(genpath('Codes/6DOF_model')) % Aircraft Models
addpath(genpath('Codes/ac_models')) % Aircraft Models
addpath(genpath('Codes/aero')) % Aerodynamics
addpath(genpath('Codes/atmospheric')) % Atmospheric Models
addpath(genpath('Codes/data')) % SAving mat structures
addpath(genpath('Codes/geometry')) % generation of geometry
addpath(genpath('Codes/missions')) % generation of geometry
addpath(genpath('Codes/performance')) % Performance models
addpath(genpath('Codes/plots')) % Plots models
addpath(genpath('Codes/prompts')) % Prompts users for information
addpath(genpath('Codes/propulsion')) % Propulsion Models
addpath(genpath('Codes/read_DATA')) % Reads Data from files
addpath(genpath('Codes/results')) % Reads Data from files
addpath(genpath('Codes/saving_DATA')) % Saving data files
addpath(genpath('Codes/stability')) % Stability analysis functions
addpath(genpath('Codes/tools')) % Coding Tools and smapling codes
addpath(genpath('Codes/weights')) % Models
% addpath(genpath('misc')) % Functions that are thought that can be eliminated