function [XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_CERVERA
%% Loads fuselage geometry
%--------------------------- FUSELAJE ---------------------------------
nSections = 9;% number of sectins (x coordinate) (eliminated 1 for convergence)
nPoints = 50; % number of elements per section (y and z coordinates)
lecture_Geo = 5; % lineas a partir desde donde empieza a leer
STL_PLOT = 1; % Generates Fuselage mesh from CAD STL is flag STL_PLOT = 1;
XFLR5_DATA.nSections = nSections;
XFLR5_DATA.nPoints = nPoints;
XFLR5_DATA.lecture_Geo = lecture_Geo;
XFLR5_DATA.STL_PLOT = STL_PLOT;
%--------------------------- CALCULATES GEOMETRY DATA ---------------------------------
% Defines flies to be used for the Fuselage geometry
XFLR5_file = strcat('Fuselage_DATA_scaled.txt');
STL_file = strcat('Fus_ProVant_v4.0.stl');
STL_file = strcat('ProVant.stl');