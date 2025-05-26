%% Reads the data from the txt files gnerated from XFLR5
function [Data_P] = Read_data_prop_APC(casos_prop_APC)

for i =1:length(casos_prop_APC)
    % Extracción de datos APC
    file_id = fopen(casos_prop_APC{1});
    aux = textscan(file_id,'%n%n%n%n%n%n%n%n','Headerlines',0);
    % New results 13 entries
    Data_P(i).raw_data= aux;
    Data_P(i).label = casos_prop_APC{i};
    Data_P(i).V = aux{1}*0.44704;      % m/s
    Data_P(i).J = aux{2};
    Data_P(i).Pe    = aux{3};
    Data_P(i).Ct   = aux{4};
    Data_P(i).Cp   = aux{5};
    Data_P(i).Hp    = aux{6}*745.6999; % Watts
    Data_P(i).Q    = aux{7};
    Data_P(i).T    = aux{8}*4.448222;  % Newtons
    fclose(file_id);
end
