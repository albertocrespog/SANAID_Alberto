function Storing_PERFORMANCE_DATA_1 = Saving_data_Performance_v0(Segments,handles,Weights_AP,Total_Datos,datos,Storing_PROPULSION_DATA,seg,filenameS,OUTPUT_read_XLSX)

Storing_PERFORMANCE_DATA_1.Segments = Segments;
Storing_PERFORMANCE_DATA_1.handles = handles;
Storing_PERFORMANCE_DATA_1.Weights_AP = Weights_AP;
Storing_PERFORMANCE_DATA_1.Total_Datos = Total_Datos;
Storing_PERFORMANCE_DATA_1.datos = datos;
Storing_PERFORMANCE_DATA_1.Storing_PROPULSION_DATA = Storing_PROPULSION_DATA;
Storing_PERFORMANCE_DATA_1.seg = seg;
% 
% Storing_PERFORMANCE_DATA_1.Weights_AP_var = Weights_AP_var;
% Storing_PERFORMANCE_DATA_1.fuel_total_var = fuel_total_var;
% Storing_PERFORMANCE_DATA_1.tiempo_total_var = tiempo_total_var;
% Storing_PERFORMANCE_DATA_1.distancia_total_var = distancia_total_var;
% Storing_PERFORMANCE_DATA_1.W_var = W_var;
% Storing_PERFORMANCE_DATA_1.datos_var = datos_var;
% Storing_PERFORMANCE_DATA_1.Plots_performance = Plots_performance;
% Storing_PERFORMANCE_DATA_1.Post_processing_PERFORMANCE = Post_processing_PERFORMANCE;

% %% Save Mat file
% %prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
% prefixa = get_fname(OUTPUT_read_XLSX);
% %st0 = strcat('data\');
% %st1 = strcat(st0,prefixa);
% %st1 = strcat(prefixa);
% st1 = filenameS.filename_DATA;
% st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
% st3 = strcat(st1,st2);
% st3A = strcat('\Study_var_m_varV.mat');
% name = strcat(st3,st3A);
% 
% % Verificar si la carpeta para `name` existe y crearla si no
% folder = fileparts(name); % Obtiene solo la ruta de la carpeta
% if ~exist(folder, 'dir') % Verifica si la carpeta no existe
%     mkdir(folder); % Crea la carpeta
% % end
% 
% save(name,'Segments','handles','seg','Weights_AP','Total_Datos','datos','Storing_PROPULSION_DATA','seg')