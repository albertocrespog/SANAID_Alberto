function [casos prefix mark_legend] = read_prop_files_May2020

% Estudios
casE = 1;

casos{casE} = 'PER3_22x12WE_RPM7000.dat'; 
st{casE} = strcat('APC 22x12WE, RPM 7000');
casE = casE + 1;

casos{casE} = 'PER3_22x12WE_RPM8000.dat'; 
st{casE} = strcat('APC 22x12WE, RPM 8000');
casE = casE + 1;

prefix = strcat('EMERGENTIA_');
mark_legend = st;


