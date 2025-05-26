function OUTPUT_read_XLSX_AAA =  read_DATA_XLSX_AAA
%% Pre-load variables stored in Excel
% Function that reads all the data for the different aircraft

% get_add_path
% filename = '../AIRCRAFT/AAA/Cirrus_Visio_SF50.xlsx';
% excel_file = 'Cirrus_Visio_SF50.xlsx';
% 
% %% MATLAB Flags
% Sheet1_COL1 = readcell(excel_file,'Sheet','Sheet1','Range','A1:A78');
% Sheet1_COL2 = readcell(excel_file,'Sheet','Sheet1','Range','B1:B78');
% 
% OUTPUT_read_XLSX_AAA.Sheet1_COL1 = Sheet1_COL1;
% OUTPUT_read_XLSX_AAA.Sheet1_COL2 = Sheet1_COL2;
% 
% DER.Sheet1_COL1{1} = 1;
% 
% function generateVariableOnFly
%    % lets tic/toc to compare the use of eval and assignin
%    tic
%    eval ( 'a = zeros(10,10);' )
%    toc
%    % an alternate method is to use a 
%    % sub function which assigns vars in the caller function:
%    tic
%    variableCreator ( 'b', zeros(10,10) )
%    toc
%    % validate that a and b both exist and are the same:
%    isequal ( a, b )  
%  end
% 
%  % Use a sub function which assigns variable in the caller function:
%  function variableCreator ( newVar, variable )
%    assignin ( 'caller', newVar, variable );
%  end