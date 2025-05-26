function CL_bar_beta = CL_bar_f_beta(x,EMERGENTIA,Data_ref,WT_Aero_Model)

% Geo_tier = EMERGENTIA.Data_Geo.Geo_tier;
% TRIM_RESULTS = EMERGENTIA.Data_Der.TRIM_RESULTS;
Stab_Der = EMERGENTIA.Data_Der.Stab_Der;
% Data_Weight = EMERGENTIA.Data_Weight.Weight_tier;
% Aero = EMERGENTIA.Data_Aero.Aero;
% Aero_TH = EMERGENTIA.Data_Aero_TH.Aero;
% Data_prop = EMERGENTIA.Data_prop;

aero_model_CL_bar = Data_ref.aero_model_CL_bar;

%--------------------------State Variables--------------------------------
% Velocity (wind-axis)          - % V       = x(1); alpha    = x(2); beta     = x(3);
% V       = x(1);
% alpha    = x(2);
beta     = x(3);

% Recalls the values of the derivatives lateral-directional
Cyb        = Stab_Der.Cyb;
Clb        = Stab_Der.Clb;
Cnb        = Stab_Der.Cnb;

% Aerodynamic model CL, CD, CM, CY, CLbar CN
% aero_model = 1; % Linear model
% aero_model = 2; % Non-linear
% aero_model = 3; % Wind tunnel
switch aero_model_CY
    case 1 % Linear Model
        CY_beta = Clb*beta; 
    case 2 % Non-linear model
        CY_beta = Clb*beta; 
    case 3 % Wind tunel model
        % Lateral Contribution
        CY_beta = Clb*beta;
        % No wind tunnel data
%         beta_vec = WT_Aero_Model.Final_MODEL_LAT2.beta;
%         CY_beta_vec = WT_Aero_Model.Final_MODEL_LAT2.CY_ac;
%         CY_beta = interp1(beta_vec,CY_beta_vec,beta,'spline');
end