function CM_alpha = CM_f_alpha(x,EMERGENTIA,Data_ref,WT_Aero_Model)

% Geo_tier = EMERGENTIA.Data_Geo.Geo_tier;
TRIM_RESULTS = EMERGENTIA.Data_Der.TRIM_RESULTS;
Stab_Der = EMERGENTIA.Data_Der.Stab_Der;
% Data_Weight = EMERGENTIA.Data_Weight.Weight_tier;
% Aero = EMERGENTIA.Data_Aero.Aero;
% Aero_TH = EMERGENTIA.Data_Aero_TH.Aero;
% Data_prop = EMERGENTIA.Data_prop;

aero_model_CM = Data_ref.aero_model_CM;

%--------------------------State Variables--------------------------------
% Velocity (wind-axis)          - % V       = x(1); alpha    = x(2); beta     = x(3);
% V       = x(1);
alpha    = x(2);
% beta     = x(3);

% Recalls the values of the derivatives - Trim
CL0        = TRIM_RESULTS.CL0_ac;
CLalpha    = Stab_Der.CL_alpha_ac;

% Aerodynamic model CL, CD, CM, CY, CLbar CN
% aero_model = 1; % Linear model
% aero_model = 2; % Non-linear
% aero_model = 3; % Wind tunnel
switch aero_model_CM
    case 1 % Linear Model
        CM_alpha = CM0 + CMalpha*alpha; 
    case 2 % Non-linear model
        CM_alpha = CM0 + CMalpha*alpha;
    case 3 % Wind tunel model
        alpha_vec = WT_Aero_Model.Final_MODEL_LON2.alpha;
        CM_alpha_vec = WT_Aero_Model.Final_MODEL_LON2.CM_ac;
        CM_alpha = interp1(alpha_vec,CM_alpha_vec,alpha,'spline');
end