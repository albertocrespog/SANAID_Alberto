%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANAID - Stability ANAlysis Interactive Design Tool
% Function that adds paths for the file to be used
% Date October 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_add_path

% Defines Generic Path so that files can be organized in folders
% change for the different versions
addpath(genpath('../src2')) % source codes for stability analysis
addpath(genpath('../AERO')) % aerodynamic codes from XFLR5
addpath(genpath('../CAD')) % CAD files
addpath(genpath('../Prop_data')) % Propdata
addpath(genpath('../Fuselage')) % Fuselage models
addpath(genpath('../AIRCRAFT')) % Selection of input data
addpath(genpath('../igesToolBox')) % Selection of input data
addpath(genpath('../Results')) % Selection of input data
% Organizing files
addpath(genpath('6DOF_model')) % Aircraft Models
addpath(genpath('ac_models')) % Aircraft Models
addpath(genpath('aero')) % Aerodynamics
addpath(genpath('atmospheric')) % Atmospheric Models
addpath(genpath('data')) % SAving mat structures
addpath(genpath('geometry')) % generation of geometry
addpath(genpath('missions')) % generation of geometry
addpath(genpath('performance')) % Performance models
addpath(genpath('plots')) % Plots models
addpath(genpath('prompts')) % Prompts users for information
addpath(genpath('propulsion')) % Propulsion Models
addpath(genpath('read_DATA')) % Reads Data from files
addpath(genpath('results')) % Reads Data from files
addpath(genpath('saving_DATA')) % Saving data files
addpath(genpath('stability')) % Stability analysis functions
addpath(genpath('tools')) % Coding Tools and smapling codes
addpath(genpath('weights')) % Models
% addpath(genpath('misc')) % Functions that are thought that can be eliminated