function Geometry_changes = EMERGENTIA_vee_tails_geometry_v0

excel_file = 'VTAIL_GEOM.xlsx';

% Sheet for geometry

%case1
Sheet_case{1} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','C18:C22'));
Sheet_case{2} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','D18:D22'));
Sheet_case{3} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','E18:E22'));
Sheet_case{4} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','F18:F22'));
Sheet_case{5} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','G18:G22'));
Sheet_case{6} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','H18:H22'));


%case2a
Sheet_case{7} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','C27:C31'));
Sheet_case{8} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','D27:D31'));
Sheet_case{9} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','E27:E31'));
Sheet_case{10} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','F27:F31'));
Sheet_case{11} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','G27:G31'));
Sheet_case{12} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','H27:H31'));


%case2b
Sheet_case{13} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','C36:C40'));
Sheet_case{14} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','D36:D40'));
Sheet_case{15} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','E36:E40'));
Sheet_case{16} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','F36:F40'));
Sheet_case{17} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','G36:G40'));
Sheet_case{18} = cell2mat(readcell(excel_file,'Sheet','Sheet1','Range','H36:H40'));


Geometry_changes.Sheet_case = Sheet_case;