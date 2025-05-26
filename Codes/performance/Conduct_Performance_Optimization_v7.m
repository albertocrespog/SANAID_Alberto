function Storing_PERFORMANCE_DATA_5 = Conduct_Performance_Optimization_v7(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
    OUTPUT_read_XLSX,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS)

Geo_tier = Storing_GEO_DATA.Geo_tier;
% Body_Geo = Storing_GEO_DATA.Body_Geo;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;
Aero = Storing_AERO_DATA.Aero;
Aero_TH = Storing_AERO_DATA.Aero_TH;
Design_criteria = Storing_AERO_DATA.Design_criteria;
Performance = Storing_AERO_DATA.Performance;
% DATA_Ae = Storing_AERO_DATA.DATA_Ae;     

%% Updates the Velocity Vector according to the minimum speeds
CL_max_ac = Aero.CL_max_ac;

% reduces the calculations since only for trim conditions
%% Study of variation of Altitude  

% Stores Variables altitude for sensitivity studies
h_low = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_low;
h_high = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_high;
h_terminal = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_terminal;

% Number of Vector elements for both altitude and velocity sensitivity study
N_h_VAR_perfo = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.N_h_VAR_perfo; %
        
h_VAR = linspace(h_low,h_high,N_h_VAR_perfo);

% Plot options
Plot_Options.h_VAR = h_VAR;
Plot_Options.N_h_VAR_perfo = N_h_VAR_perfo;


% loop
for j=1:length(h_VAR)

    h_inicial_d = h_VAR(j);
    [Temp_init,rho_init,p_init,a_init]=atmos_inter_mio(h_inicial_d);
    V_stall = sqrt(2*Weight_tier.m_TOW*conv_UNITS.g/(Geo_tier.S_w1*rho_init*CL_max_ac));
    V_TO = 1.2*V_stall;
    V_min_ope = 1.3*V_stall;

    Restrictions.V_stall = V_stall;
    Restrictions.V_TO = V_TO;
    Restrictions.V_min_ope = V_min_ope;

    OUTPUT_read_XLSX.IPP_flags.h_inicial_d    = h_inicial_d;
    OUTPUT_read_XLSX.IPP_flags.h_final_d      = h_terminal;

    % Calculations
    Storing_PERFORMANCE_5{j} = Conduct_Performance_Optimization_v2(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
        OUTPUT_read_XLSX,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA);

    % Storing
    Restrictions_var_h_V{j} = Restrictions;

end


%% Saves the data to a mat file
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_PERFORMANCE_DATA_5 = Saving_data_Performance_var_h_V_Glide(Storing_PERFORMANCE_5,Plot_Options,Restrictions_var_h_V,OUTPUT_read_XLSX,Storing_PROPULSION_DATA,filenameS);

else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_PERFORMANCE_DATA_5 = dummy;
end