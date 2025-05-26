function [XFLR5_DATA,XFLR5_file,STL_files] = Extraction_fuselage_data(OUTPUT_read_XLSX)
%% Loads fuselage geometry
%--------------------------- FUSELAJE ---------------------------------
nSections = OUTPUT_read_XLSX.Fuselage_flags.nSections;% number of sectins (x coordinate) (eliminated 1 for convergence)
nPoints = OUTPUT_read_XLSX.Fuselage_flags.nPoints; % number of elements per section (y and z coordinates)
lecture_Geo = OUTPUT_read_XLSX.Fuselage_flags.lecture_Geo; % lineas a partir desde donde empieza a leer
STL_PLOT = OUTPUT_read_XLSX.Fuselage_flags.STL_PLOT; % Generates Fuselage mesh from CAD STL is flag STL_PLOT = 1;
XFLR5_DATA.nSections = nSections;
XFLR5_DATA.nPoints = nPoints;
XFLR5_DATA.lecture_Geo = lecture_Geo;
XFLR5_DATA.STL_PLOT = STL_PLOT;
%--------------------------- CALCULATES GEOMETRY DATA ---------------------------------
% Defines flies to be used for the Fuselage geometry
% XFLR5_file = strcat('Fuselage_DATA_scaled.txt');
% STL_file = strcat('Fus_ProVant_v4.0.stl');
% STL_file = strcat('ProVant.stl');
XFLR5_file = strcat(OUTPUT_read_XLSX.Fuselage_flags.XFLR5_file);
STL_files.STL_ac =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_ac);
STL_files.STL_fus =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_fus);
STL_files.STL_wing =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_wing);
STL_files.STL_canard =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_canard);
STL_files.STL_HTP =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_HTP);
STL_files.STL_VTP =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_VTP);
STL_files.STL_Vee =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_Vee);
STL_files.STL_engine =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_engine);
STL_files.STL_nacelle =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_nacelle);
