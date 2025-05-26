function [Aero,DATA_PL,Fig] = Aerodynamic_Design_Aug2018(Geo_tier,...
    Weight_tier,conv_UNITS,Design_criteria,Performance,CASOS,degrees_XAC,fig,XFLR5_DATA,...
    MAC_Estimation,DATA_Ae,X_OC,Plot_Options,VECTOR_XFLR5)

% Generates Aerodynamic Data
[Aero,DATA_PL] = Generate_Aero_Data(DATA_Ae,Design_criteria,Performance,Geo_tier,Weight_tier,conv_UNITS);

if MAC_Estimation == 1
    % Routine that determines the Aerodynamic Center  
    [XAC,Fig] = Generate_XAC(DATA_Ae,Aero,DATA_PL,CASOS,degrees_XAC,fig,...
    XFLR5_DATA,Performance,X_OC,Plot_Options,VECTOR_XFLR5)
else
    Fig = fig;
end
