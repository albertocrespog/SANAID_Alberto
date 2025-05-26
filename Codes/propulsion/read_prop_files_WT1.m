function [casos_prop_WT1 prefix mark_legend] = read_prop_files_WT1

% 1 - APC 20x8
% 2 - APC 22x10 
% 3 - APC 22x12 
% 4 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
% 5 - APC 21x13W 
% 6 - APC 21x14  
% Estudios
casE = 1;

casos_prop_WT1{casE} = 'WT1 - APC 20x8.dat'; 
st{casE} = strcat('WT1 - APC 20x8');
casE = casE + 1;

casos_prop_WT1{casE} = 'WT1 - APC 22x10.dat'; 
st{casE} = strcat('WT1 - APC 22x10');
casE = casE + 1;

casos_prop_WT1{casE} = 'WT1 - APC 22x12'; 
st{casE} = strcat('WT1 - APC 22x12');
casE = casE + 1;

casos_prop_WT1{casE} = 'PER3_22x12WE_RPM7000.dat'; 
st{casE} = strcat('WT1 - APC 22x12W');
casE = casE + 1;

casos_prop_WT1{casE} = 'WT1 - APC 21x13'; 
st{casE} = strcat('WT1 - APC 21x13');
casE = casE + 1;

casos_prop_WT1{casE} = 'WT1 - APC 21x14'; 
st{casE} = strcat('WT1 - APC 21x14');
casE = casE + 1;

prefix = strcat('EMERGENTIA_');
mark_legend = st;