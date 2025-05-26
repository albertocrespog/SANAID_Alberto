function Storing_STABILITY_DATA_5 = Conduct_Stability_Study5_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
            OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA,Storing_WEIGHT_DATA,...
            Storing_AERO_DATA,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS)

XCG_FF = OUTPUT_read_XLSX.Stability_flags.XCG_FF;
Dxcg = Modifications.Dxcg;
 
Geo_tier = Storing_GEO_DATA.Geo_tier;
Performance = Storing_AERO_DATA.Performance;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;

% %%%%%%%%%%%%%%%%%%%% Turning conditioms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant phi
phi = OUTPUT_read_XLSX.Stability_flags.phi;
% variable phi
phi_i = OUTPUT_read_XLSX.Stability_flags.phi_i;
phi_f = OUTPUT_read_XLSX.Stability_flags.phi_f;
N_Delta_phi = OUTPUT_read_XLSX.Stability_flags.N_Delta_phi;
n_viraje = OUTPUT_read_XLSX.Stability_flags.n_viraje;
phi_vec = linspace(phi_i,phi_f,N_Delta_phi);
% Storing study conditions
conditions_TRIM_turning.phi = phi;
conditions_TRIM_turning.phi_vec = phi_vec;
conditions_TRIM_turning.n_viraje = n_viraje;
Stab_Der = Storing_STABILITY_DATA_1.Stab_Der;

% rho = 1.225; % At sea level
% rho = 0.96296; % At 8000 ft
% rho = 0.74628; % At 16000 ft
conditions_TRIM_turning.rho = Performance.rho;

[Trim_ITER_LAT_Viraje] = Calculo_Trim_ITER_LAT_Viraje_v3(conv_UNITS,conditions_TRIM_turning,...
    Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,Storing_STABILITY_DATA_1);
%% Saves the data to a mat file
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_STABILITY_DATA_5 = Saving_Trim_Lat_Turn(Trim_ITER_LAT_Viraje,conditions_TRIM_turning,OUTPUT_read_XLSX,filenameS);
% prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
%     st0 = strcat('data\');
%     st1 = strcat(st0,prefixa);
%     st2 = strcat('\Study_Trim_ITER_LAT_Turn.mat');
%     name   = strcat(st1,st2);
%     save(name, 'Trim_ITER_LAT_Viraje','conditions_TRIM_turning')
%     Storing_STABILITY_DATA_5.Trim_ITER_LAT_Viraje = Trim_ITER_LAT_Viraje;
%     Storing_STABILITY_DATA_5.conditions_TRIM_turning = conditions_TRIM_turning;
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_5 = dummy;
end