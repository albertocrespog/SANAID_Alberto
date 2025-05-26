function Storing_AERO_DATA_1 = Saving_data_Aero(VECTOR_XFLR5,Design_criteria,DATA_Ae,casos,prefix,mark_legend,X_OC,Aero_TH,Aero,DATA_PL,Performance,OUTPUT_read_XLSX,filenameS)

% prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
prefixa = get_fname(OUTPUT_read_XLSX);
%st0 = strcat('data\');
%st1 = strcat(st0,prefixa);
%st1 = strcat(prefixa);
st1 = filenameS.filename_DATA;
st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st3 = strcat(st1,st2);
% status = mkdir(st3);
st2A = strcat('\Aero.mat');
st2B = strcat('\Aero_TH.mat');
st2C = strcat('\Design_criteria.mat');

nameA   = strcat(st3,st2A);
nameB   = strcat(st3,st2B);
nameC   = strcat(st3,st2C);

% Verificar si la carpeta para `name` existe y crearla si no
folderA = fileparts(nameA); % Obtiene solo la ruta de la carpeta
if ~exist(folderA, 'dir') % Verifica si la carpeta no existe
    mkdir(folderA); % Crea la carpeta
end

folderB = fileparts(nameB); % Obtiene solo la ruta de la carpeta
if ~exist(folderB, 'dir') % Verifica si la carpeta no existe
    mkdir(folderB); % Crea la carpeta
end

folderC = fileparts(nameC); % Obtiene solo la ruta de la carpeta
if ~exist(folderC, 'dir') % Verifica si la carpeta no existe
    mkdir(folderC); % Crea la carpeta
end

save(nameA, 'Aero')
save(nameB, 'Aero_TH')
save(nameC, 'Design_criteria')

Storing_AERO_DATA_1.VECTOR_XFLR5 = VECTOR_XFLR5;
Storing_AERO_DATA_1.Design_criteria = Design_criteria;
Storing_AERO_DATA_1.DATA_Ae = DATA_Ae;
Storing_AERO_DATA_1.casos = casos;
Storing_AERO_DATA_1.prefix = prefix;
Storing_AERO_DATA_1.mark_legend = mark_legend;
Storing_AERO_DATA_1.X_OC = X_OC;
Storing_AERO_DATA_1.Aero_TH = Aero_TH;
Storing_AERO_DATA_1.Aero = Aero;
Storing_AERO_DATA_1.DATA_PL = DATA_PL;
Storing_AERO_DATA_1.Performance = Performance;


