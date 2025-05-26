function CM = CM_model(x,Ctrl_Sig,EMERGENTIA,Data_ref,WT_Aero_Model,DATA_int,Data_ATM,conv_UNITS)
                       
% Geo_tier = EMERGENTIA.Data_Geo.Geo_tier;
% TRIM_RESULTS = EMERGENTIA.Data_Der.TRIM_RESULTS;
Stab_Der = EMERGENTIA.Data_Der.Stab_Der;
% Data_Weight = EMERGENTIA.Data_Weight.Weight_tier;
% Aero = EMERGENTIA.Data_Aero.Aero;
% Aero_TH = EMERGENTIA.Data_Aero_TH.Aero;
% Data_prop = EMERGENTIA.Data_prop;

%--------------------------State Variables--------------------------------
% Velocity (wind-axis)          - % V       = x(1); alpha    = x(2); beta     = x(3);
V       = x(1);
alpha    = x(2);
beta     = x(3);

% Velocity (body axis)          - % u_p      = x(4); v_p      = x(5); w_p      = x(6);
u_p      = x(4);
v_p      = x(5);
w_p      = x(6);
% Orientation (body-axis)       - % phi      = x(7); theta    = x(8); psi      = x(9);
phi      = x(7);
theta    = x(8);
psi      = x(9);
% Angular Velocity (body-axis)  - % p        = x(10); q        = x(11); r        = x(12);
p        = x(10);
q        = x(11);
r        = x(12);
% Orientation (wind-axis)       - % mu_angle = x(13); gamma    = x(14); xi       = x(15);
mu_angle = x(13);
gamma    = x(14);
chi       = x(15);
% Position (body-axis)          - % x_p      = x(16); y_p      = x(17); z_p      = x(18);
x_p      = x(16);
y_p      = x(17);
z_p      = x(18);
% Position (wind-axis)          - % x_p      = x(19); y_p      = x(20); z_p      = x(21);
x_w      = x(19);
y_w      = x(20);
z_w      = x(21);
h = z_w;

delta_e = Ctrl_Sig.delta_e;
% delta_T = Ctrl_Sig.delta_T;
% delta_a = Ctrl_Sig.delta_a;
% delta_r = Ctrl_Sig.delta_r;

% Recalls the values of the derivatives - Trim
% CL0        = TRIM_RESULTS.CL0_ac;
% CM0        = TRIM_RESULTS.CM0_ac;
% CD0        = Aero_TH.CD0;
% CD1        = Aero_TH.CD1;
% CD2        = Aero_TH.CD2;
% CLalpha    = TRIM_RESULTS.CL_alpha_ac;
CMalpha    = TRIM_RESULTS.CM_alpha_ac;

% CDalpha        = Stab_Der.CDalpha;
% CLalpha        = Stab_Der.CL_alpha_ac;
CMalpha        = Stab_Der.CM_alpha_ac;
% CXalfa        = Stab_Der.CXalfa;
% CZalfa        = Stab_Der.CZalfa;

% CDalphapunto        = Stab_Der.CDalfapunto;
% CLalphapunto        = Stab_Der.CLalphapunto;
CMalphapunto        = Stab_Der.CMalphapunto;
% CXalfapunto        = Stab_Der.CXalfapunto;
% CZalfapunto        = Stab_Der.CZalfapunto;

% CDq        = Stab_Der.CDq;
% CLq        = Stab_Der.CLq;
CMq        = Stab_Der.CMq;
    
% CXteta        = Stab_Der.CXteta;
% CZteta        = Stab_Der.CZteta;
% CLteta        = Stab_Der.CLteta;
% CDteta        = Stab_Der.CDteta;
CMteta        = Stab_Der.CMteta;

% CDu        = Stab_Der.CDu;
% CLu        = Stab_Der.CLu;
CMu        = Stab_Der.CMu;

% Longitudinal control
% CLdelta_e        = Stab_Der.CL_delta_e;
% CDdelta_e        = -Stab_Der.CXdeltae;
CMdelta_e        = Stab_Der.CMdeltae;

% Recalls the values of the derivatives lateral-directional
% Cyb        = Stab_Der.Cyb;
% Clb        = Stab_Der.Clb;
% Cnb        = Stab_Der.Cnb;

% Cybpunto        = Stab_Der.Cybpunto;
% Clbpunto        = Stab_Der.Clbpunto;
% Cnbpunto        = Stab_Der.Cnbpunto;

% Cyr        = Stab_Der.Cyr;
% Clr        = Stab_Der.Clr;
% Cnr        = Stab_Der.Cnr;

% Cyp        = Stab_Der.Cyp;
% Clp        = Stab_Der.Clp;
% Cnp        = Stab_Der.Cnp;

% Cydeltaa   = Stab_Der.Cydeltaa;
% Cldeltaa   = Stab_Der.Cldeltaa;
% Cndeltaa   = Stab_Der.Cndeltaa;

% Cydeltar   = Stab_Der.Cydeltar;
% Cldeltar   = Stab_Der.Cldeltar;
% Cndeltar   = Stab_Der.Cndeltar;

% Cydeltarv   = Stab_Der.Cydeltarv;
% Cldeltarv   = Stab_Der.Cldeltarv;
% Cndeltarv   = Stab_Der.Cndeltarv;

%--------------------------Constants--------------------------------
Data_const = get_const_EMERGENTIA(DATA_int,Data_ATM,h,V,1/cos(phi),EMERGENTIA,conv_UNITS);
C_der_ad_lon = Data_const.C_der_ad_lon;
C_der_ad_lat = Data_const.C_der_ad_lat;

% Adimensionalize states longitudinal and lateral-directional
KL1 = C_der_ad_lon;
KL2 = C_der_ad_lat;

% Defines Aerodynamic model
approx_Coeff = Data_ref.approx_Coeff;
switch approx_Coeff
    case 0
        % Aerodynamic Force Coeffcients in Wind Axis
        CM_alpha = CM_f_alpha(x,EMERGENTIA,Data_ref,WT_Aero_Model);
        CM = CM_alpha  + CMq*KL1*q + CMdelta_e*delta_e;
    case 1
        % Aerodynamic Force Coeffcients in Wind Axis
        CM_alpha = CM_f_alpha(x,EMERGENTIA,Data_ref,WT_Aero_Model)
        CM = CM_alpha + CMu*u + CMteta*theta + CMq*KL1*q + CMdelta_e*delta_e;
    case 2
        % Aerodynamic Force Coeffcients in Wind Axis
        CM_alpha = CM_f_alpha(x,EMERGENTIA,Data_ref,WT_Aero_Model)
        CM = CM_alpha + CMu*u + CMteta*theta + CMq*KL1*q + CMalphapunto*KL1*alphadot + CMdelta_e*delta_e;
end