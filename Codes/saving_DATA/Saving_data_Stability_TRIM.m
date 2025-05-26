function Storing_STABILITY_DATA_1 = Saving_data_Stability_TRIM(TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts,Stab_Dyn_Long,Stab_Dyn_LatDir,Propulsion_values,OUTPUT_read_XLSX,Effects,filenameS)

% Save Mat file
%prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
prefixa = get_fname(OUTPUT_read_XLSX);
%st0 = strcat('data\');
%st1 = strcat(st0,prefixa);
%st1 = strcat(prefixa);
st1 = filenameS.filename_DATA;

st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st3 = strcat(st1,st2);
st3A = strcat('\Stability_Derivatives.mat');
name = strcat(st3,st3A);

% Verificar si la carpeta para `name` existe y crearla si no
folder = fileparts(name); % Obtiene solo la ruta de la carpeta
if ~exist(folder, 'dir') % Verifica si la carpeta no existe
    mkdir(folder); % Crea la carpeta
end

save(name,'TRIM_RESULTS','Trim_ITER','Stab_Der','Stab_Der_parts','Stab_Dyn_Long','Stab_Dyn_LatDir','Propulsion_values','Effects')

Storing_STABILITY_DATA_1.TRIM_RESULTS = TRIM_RESULTS;
Storing_STABILITY_DATA_1.Trim_ITER = Trim_ITER;
Storing_STABILITY_DATA_1.Stab_Der = Stab_Der;
Storing_STABILITY_DATA_1.Stab_Der_parts = Stab_Der_parts;
Storing_STABILITY_DATA_1.Stab_Dyn_Long = Stab_Dyn_Long;
Storing_STABILITY_DATA_1.Stab_Dyn_LatDir = Stab_Dyn_LatDir;
Storing_STABILITY_DATA_1.Propulsion_values = Propulsion_values;