function Geometry_changes = ONEiRE_geometry_v0

excel_file = 'matriz_geomHTP_vble.xlsx';

% Sheet for geometry
Sheet_case{1} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','K17:K21'));
Sheet_case{2} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','K23:K27'));
Sheet_case{3} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','K11:K15'));
 
Sheet_case{4} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','I17:I21'));
Sheet_case{5} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','I23:I27'));
Sheet_case{6} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','I11:I15'));

Sheet_case{7} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','M17:M21'));
Sheet_case{8} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','M23:M27'));
Sheet_case{9} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','M11:M15'));

Sheet_case{10} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','B2:B2'));

Geometry_changes.Sheet_case = Sheet_case;