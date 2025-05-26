function Storing_PERFORMANCE_DATA_4C = Conduct_Performance_Optimization_v4C(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
    OUTPUT_read_XLSX,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS)

Geo_tier = Storing_GEO_DATA.Geo_tier;
% Body_Geo = Storing_GEO_DATA.Body_Geo;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;
Aero = Storing_AERO_DATA.Aero;
Aero_TH = Storing_AERO_DATA.Aero_TH;
Design_criteria = Storing_AERO_DATA.Design_criteria;
Performance = Storing_AERO_DATA.Performance;
% DATA_Ae = Storing_AERO_DATA.DATA_Ae;

% Stores Variables altitude for sensitivity studies
h_low = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_low;
h_high = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_high;

V_low = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.V_low;
V_high = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.V_high;

% Number of Vector elements for both altitude and velocity sensitivity study
N_h_VAR_perfo = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.N_h_VAR_perfo; %
N_V_VAR_perfo = OUTPUT_read_XLSX.IPP_flags.Variable_Study5.N_V_VAR_perfo; %
        
%% Updates the Velocity Vector according to the minimum speeds
CL_max_ac = Aero.CL_max_ac;
V_stall = sqrt(2*Weight_tier.m_TOW*conv_UNITS.g/(Geo_tier.S_w1*Performance.rho*CL_max_ac));

% Variation of speeds but restricting speeds to flight velocities
if V_low < V_stall
    V_low = V_stall;
end
V_VAR = linspace(V_low,V_high,N_V_VAR_perfo);
h_VAR = linspace(h_low,h_high,N_h_VAR_perfo);

% Plot options
Plot_Options.h_VAR = h_VAR;
Plot_Options.V_VAR = V_VAR;
Plot_Options.N_h_VAR_perfo = N_h_VAR_perfo;
Plot_Options.N_V_VAR_perfo = N_V_VAR_perfo;


% reduces the calculations since only for trim conditions
%% Study of variation of Velocity and mass with XCG selected for trim conditions at desired Static Margin
for i=1:length(V_VAR)
    
    % loop
    for j=1:length(h_VAR)
        % actual values into de variable study
        V_cr = V_VAR(i);
        h_inicial_cr = h_VAR(j);

        OUTPUT_read_XLSX.IPP_flags.V_cr           = V_cr;
        OUTPUT_read_XLSX.IPP_flags.h_inicial_cr   = h_inicial_cr;
        
        % OUTPUT_read_XLSX.IPP_flags.dist_final_cr = dist_final_cr;


        [Temp_init,rho_init,p_init,a_init]=atmos_inter_mio(h_inicial_cr);

        V_stall(i) = sqrt(2*Weight_tier.m_TOW*conv_UNITS.g/(Geo_tier.S_w1*rho_init*CL_max_ac));
        V_TO(i) = 1.2*V_stall(i);
        V_min_ope(i) = 1.3*V_stall(i);

        Restrictions.V_stall = V_stall(i);
        Restrictions.V_TO = V_TO(i);
        Restrictions.V_min_ope = V_min_ope(i);
              
        % Calculations
        Storing_PERFORMANCE_4C{i,j} = Conduct_Performance_Optimization_v2(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
            OUTPUT_read_XLSX,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA);

        % Storing
            Restrictions_var_h_V{i,j} = Restrictions;

    end
end

%% Saves the data to a mat file
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_PERFORMANCE_DATA_4C = Saving_data_Performance_var_h_V(Storing_PERFORMANCE_4C,Plot_Options,Restrictions_var_h_V,OUTPUT_read_XLSX,Storing_PROPULSION_DATA,filenameS);

else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_PERFORMANCE_DATA_4C = dummy;
end