function Storing_PERFORMANCE_DATA_2 = Saving_data_Performance_var_h_Max_Glide(Storing_PERFORMANCE_5,Plot_Options,Restrictions_var_h_MaxGlide,OUTPUT_read_XLSX,Storing_PROPULSION_DATA);
         
Storing_PERFORMANCE_DATA_2.Storing_PERFORMANCE_5 = Storing_PERFORMANCE_5;
Storing_PERFORMANCE_DATA_2.Plot_Options = Plot_Options;
Storing_PERFORMANCE_DATA_2.Restrictions_var_h_MaxGlide = Restrictions_var_h_MaxGlide;
Storing_PERFORMANCE_DATA_2.Storing_PROPULSION_DATA = Storing_PROPULSION_DATA;

% Save Mat file
prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
st0 = strcat('data\');
st1 = strcat(st0,prefixa);
st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st3 = strcat(st1,st2);
st3A = strcat('\Study_Performance_var_h_MaxGlide.mat');
name = strcat(st3,st3A);
save(name,'Storing_PERFORMANCE_5','Plot_Options','Restrictions_var_h_MaxGlide','Storing_PROPULSION_DATA');