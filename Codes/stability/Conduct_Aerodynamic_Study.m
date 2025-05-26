function [Storing_AERO_DATA_1] = Conduct_Aerodynamic_Study(OUTPUT_read_XLSX,AC_CONFIGURATION,conv_UNITS,Performance_preliminar,case_AC,...
    Storing_WEIGHT_DATA,Storing_GEO_DATA,conditions,filenameS)

Weight_tier = Storing_WEIGHT_DATA.Weight_tier;
Geo_tier = Storing_GEO_DATA.Geo_tier;
Body_Geo = Storing_GEO_DATA.Body_Geo;
% meshData = Storing_GEO_DATA.meshData;

% Data for analysis of XFLR5 DATA
% Generates the file for Aerodynamic Data
if OUTPUT_read_XLSX.PLOT_flags.plot_aero1 == 1 || OUTPUT_read_XLSX.PLOT_flags.plot_aero2 == 1 % Decided to plot Aero and Aero Polar
    VECTOR_XFLR5 = Get_Comparing_vectors_Aerodynamic(case_AC,OUTPUT_read_XLSX);
else
    DUMMY = 1;
    VECTOR_XFLR5 = DUMMY;
end

% Selection of areodynamic properties (incidence angles)
Design_criteria = Generate_aerodynamic_data(OUTPUT_read_XLSX,AC_CONFIGURATION,conv_UNITS);

% Read Aerodynamic Data
[DATA_Ae,casos,prefix,mark_legend,X_OC] = Read_data_Aero_v2(Performance_preliminar,case_AC,OUTPUT_read_XLSX);

% Premiliminary Aerodynamic Design
% Runs twice, the first with LLT for CLmax and polar, the second for the
% rest of aerdynamic properties

% Aerodynamic Resuts from Codes (CFLR5, FLOW5)
% [Aero,DATA_PL,Performance] = Generate_Aero_Data_2021_v1(DATA_Ae,Design_criteria,Performance_preliminar,Geo_tier...
%     ,Weight_tier,conv_UNITS,AC_CONFIGURATION, Body_Geo,OUTPUT_read_XLSX);
[Aero,DATA_PL,Performance] = Generate_Aero_Data_2021_v3(DATA_Ae,Design_criteria,Performance_preliminar,Geo_tier...
    ,Weight_tier,conv_UNITS,AC_CONFIGURATION, Body_Geo,OUTPUT_read_XLSX);

% Fusion Aerodynamic Properties using theoretica results and theoretic
% [Aero_TH,Aero] = Polar_Literatura_v1(Performance_preliminar,Geo_tier,conv_UNITS,Body_Geo,AC_CONFIGURATION,OUTPUT_read_XLSX,Aero,conditions);
[Aero_TH,Aero] = Polar_Literatura_v2(Performance_preliminar,Geo_tier,conv_UNITS,Body_Geo,AC_CONFIGURATION,OUTPUT_read_XLSX,Aero,conditions);

%   Xac_cre_W_B = get_xac_cre_wing(Geo_tier.lambda_w1, Geo_tier.AR_w1,Geo_tier.Lambda_LE_w1,Performance.Mach);
%   Geo_tier.xbar_w1 = Xac_cre_W_B*cR_w1;  %Consultar!!!
%   Geo_tier.xbar_w1_e = Xac_cre_W_B*cR_w1;  %Consultar!!!

% Selects the Polar estimates that will be used including configuration
% of Cruise, TakeOff or Climb
Polar = get_Aero_Polar(Aero,Aero_TH,OUTPUT_read_XLSX,Performance,Geo_tier);
Aero.Polar = Polar;

%% Saves the data to a mat file
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_AERO_DATA_1 = Saving_data_Aero(VECTOR_XFLR5,Design_criteria,DATA_Ae,casos,prefix,mark_legend,X_OC,Aero_TH,Aero,DATA_PL,Performance,OUTPUT_read_XLSX,filenameS);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_AERO_DATA_1 = dummy;
end
