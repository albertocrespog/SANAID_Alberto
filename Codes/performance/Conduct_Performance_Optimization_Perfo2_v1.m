function Storing_PERFORMANCE_DATA_22 = Conduct_Performance_Optimization_Perfo2_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
    OUTPUT_read_XLSX,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS)

Geo_tier = Storing_GEO_DATA.Geo_tier;
% Body_Geo = Storing_GEO_DATA.Body_Geo;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;
Aero = Storing_AERO_DATA.Aero;
Aero_TH = Storing_AERO_DATA.Aero_TH;
Design_criteria = Storing_AERO_DATA.Design_criteria;
Performance = Storing_AERO_DATA.Performance;
% DATA_Ae = Storing_AERO_DATA.DATA_Ae;

N_V_VAR_perfo = OUTPUT_read_XLSX.Performance_pre_flags.N_V_VAR_perfo; %

ME_true = OUTPUT_read_XLSX.Weights_flags.ME_true; % Truen ME - Empty Weight
MF_true = OUTPUT_read_XLSX.Weights_flags.MF_true; % Truen MTOW - Maximum Fuel Weight
MPL_true = OUTPUT_read_XLSX.Weights_flags.MPL_true; % Truen MTOW - Maximum Paylod Weight

m_empty = ME_true; % Peso en vacío
m_payload = MPL_true; % Carga de pago

% Stores Variables mass for sensitivity studies
m_fuel_low = MF_true*0.1;
m_fuel_high = MF_true;
m_fuel_VAR = linspace(m_fuel_low,m_fuel_high,N_m_VAR_perfo);

m_VAR_min = m_empty + m_payload + m_fuel_low; % total mass
m_VAR_max = m_empty + m_payload + m_fuel_high; % total mass

%% Updates the Velocity Vector according to the minimum speeds
CL_max_ac = Aero.CL_max_ac;
V_stall_min_low = sqrt(2*m_VAR_min*conv_UNITS.g/(Geo_tier.S_w1*Performance.rho*CL_max_ac));
V_stall_min_high = sqrt(2*m_VAR_max*conv_UNITS.g/(Geo_tier.S_w1*Performance.rho*CL_max_ac));
V_min_ope_low = 1.2*V_stall_min_low;
V_min_ope_high = 1.2*V_stall_min_high;

% Variation of speeds
V_low = V_stall_min_low;
V_high = 1.25*V_min_ope_high;
V_VAR = linspace(V_low,V_high,N_V_VAR_perfo);

% Plot options
% Plot options
Plot_Options.m_fuel_VAR = m_fuel_VAR;
Plot_Options.V_VAR = V_VAR;
Plot_Options.N_m_VAR_perfo = N_m_VAR_perfo;
Plot_Options.N_V_VAR_perfo = N_V_VAR_perfo;

% reduces the calculations since only for trim conditions
%% Study of variation of Velocity and mass with XCG selected for trim conditions at desired Static Margin
for i=1:length(m_fuel_VAR)
    
    fuel_cr = m_fuel_VAR(i);
    % Saves variable for fuel maximization
    OUTPUT_read_XLSX.IPP_flags.fuel_cr = fuel_cr;

    m_VAR(i) = m_empty + m_payload + m_fuel_VAR(i); % total mass

    % % Stores values into de restrictions so that can be plotted
    % V_stall(i) = sqrt(2*m_VAR(i)*conv_UNITS.g/(Geo_tier.S_w1*Performance.rho*CL_max_ac));
    % V_TO(i) = 1.2*V_stall(i);
    % V_min_ope(i) = 1.3*V_stall(i);
    % 
    % Restrictions.V_stall = V_stall(i);
    % Restrictions.V_TO = V_TO(i);
    % Restrictions.V_min_ope = V_min_ope(i);
    % Restrictions.D_prop = OUTPUT_read_XLSX.Propulsive_flags.D_prop;
    % 
    % % loop
    % for j=1:length(V_VAR)
    %     % actual values into de variable study
    %     OUTPUT_read_XLSX.IPP_flags.V_cr = V_VAR(j);

        % Calculations
        Storing_PERFORMANCE_22{i} = Conduct_Performance_Optimization_Perfo0_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
            OUTPUT_read_XLSX,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS);

        % Storing
        Restrictions_var_V_m{i} = Restrictions;
    % end
end

%% Saves the data to a mat file
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_PERFORMANCE_DATA_22 = Saving_data_Performance_22(Storing_PERFORMANCE_22,Plot_Options,Restrictions_var_V_m,OUTPUT_read_XLSX,Storing_PROPULSION_DATA,filenameS);

else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_PERFORMANCE_DATA_22 = dummy;
end

   
