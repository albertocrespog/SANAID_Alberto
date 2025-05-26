function Storing_STABILITY_DATA_2 = Saving_data_Stability_varV_vargamma(TRIM_RESULTS_var_V_gamma,Trim_ITER_var_V_gamma,Stab_Der_var_V_gamma,Stab_Der_parts_V_gamma,...
    Stab_Dyn_Long_var_V_gamma,Stab_Dyn_LatDir_var_V_gamma,case_AC,Propulsion_var_V_gamma,Restrictions_var_V_gamma,Plot_Options,OUTPUT_read_XLSX,x_XCG_VAR,filenameS)

Storing_STABILITY_DATA_2.TRIM_RESULTS_var_V_gamma = TRIM_RESULTS_var_V_gamma;
Storing_STABILITY_DATA_2.Trim_ITER_var_V_gamma = Trim_ITER_var_V_gamma;
Storing_STABILITY_DATA_2.Stab_Der_var_V_gamma = Stab_Der_var_V_gamma;
Storing_STABILITY_DATA_2.Stab_Der_parts_V_gamma = Stab_Der_parts_V_gamma;
Storing_STABILITY_DATA_2.Stab_Dyn_Long_var_V_gamma = Stab_Dyn_Long_var_V_gamma;
Storing_STABILITY_DATA_2.Stab_Dyn_LatDir_var_V_gamma = Stab_Dyn_LatDir_var_V_gamma;
Storing_STABILITY_DATA_2.Propulsion_var_V_gamma = Propulsion_var_V_gamma;
Storing_STABILITY_DATA_2.Restrictions_var_V_gamma = Restrictions_var_V_gamma;
Storing_STABILITY_DATA_2.Plot_Options = Plot_Options;
Storing_STABILITY_DATA_2.x_XCG_VAR = x_XCG_VAR;

% Save Mat file

%prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
prefixa = get_fname(OUTPUT_read_XLSX);
%st0 = strcat('data\');
%st1 = strcat(st0,prefixa);
%st1 = strcat(prefixa);
st1 = filenameS.filename_DATA;
st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st3 = strcat(st1,st2);
st3A = strcat('\Study_var_gamma_varV.mat');
name = strcat(st3,st3A);

% Verificar si la carpeta para `name` existe y crearla si no
folder = fileparts(name); % Obtiene solo la ruta de la carpeta
if ~exist(folder, 'dir') % Verifica si la carpeta no existe
    mkdir(folder); % Crea la carpeta
end

save(name,'TRIM_RESULTS_var_V_gamma','Trim_ITER_var_V_gamma','Stab_Der_var_V_gamma','Stab_Der_parts_V_gamma','Stab_Dyn_Long_var_V_gamma','Stab_Dyn_LatDir_var_V_gamma','Propulsion_var_V_gamma','Restrictions_var_V_gamma','Plot_Options','x_XCG_VAR')