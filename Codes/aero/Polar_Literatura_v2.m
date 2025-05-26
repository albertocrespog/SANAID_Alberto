function  [Aero_TH,Aero] = Polar_Literatura_v2(Performance_preliminar,Geo_tier,conv_UNITS,Body_Geo,AC_CONFIGURATION,OUTPUT_read_XLSX,Aero,conditions)

h = Performance_preliminar.h;
M = Performance_preliminar.Mach;
rho = Performance_preliminar.rho;
V = Performance_preliminar.V;

m2ft = conv_UNITS.m2ft;


% Flag that determines if fuse between CBM and numerical methods for the
% drag coefficients
fuse_aero_FLOW_and_CBM = OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM; %

Num_eng = AC_CONFIGURATION.n_eng;
l_nc = AC_CONFIGURATION.l_nc;
d_nc = AC_CONFIGURATION.d_nc;

Conf = OUTPUT_read_XLSX.Aerodynamic_Data_flags.Conf;

% Modélos teóricos de Aircraft Design - D.P. Raymer
mu = 0.000016938;
S_ref = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
w_Area_b_max = Body_Geo.w_Area_b_max;

% Skink el skin roughness coefficient (ft)
k_sfr1 = 1.015E-5; % camouflage paint on aluminium
k_sfr2 = 0.634E-5; % smooth paint
k_sfr3 = 0.405E-5; % production sheet metal
k_sfr4 = 0.152E-5; % polished sheet metal
k_sfr5 = 0.052E-5; % smooth molded composite

% correction Factor skin friction coefficient
K_corr = 1.0;

% Initialize parasite drag coefficient of elements
CD0_w1     = 0;
CD0_HTP     = 0;
CD0_v      = 0;
CD0_v2     = 0;
CD0_vee  = 0;
CD0_vee2  = 0;
CD0_can    = 0;
CD0_tb    = 0;
CD0_nac    = 0;
CD0_lnd    = 0;
CD0_missile = 0;
CD0_pod    = 0;
CD0_flap    = 0;


% Initialize parasite drag coefficient of elements for the FUSe between CMB
% and numerical
CD0_w1_fuse = 0;
CD1_w1_fuse = 0;
CD2_w1_fuse = 0;
% CD0_HTP_fuse = 0;
% CD1_HTP_fuse = 0;
% CD2_HTP_fuse = 0;
CD0_HTP_fuse = 0;
CD1_HTP_fuse = 0;
CD2_HTP_fuse = 0;
CD0_can_fuse = 0;
CD1_can_fuse = 0;
CD2_can_fuse = 0;
CD0_vee_fuse = 0;
CD1_vee_fuse = 0;
CD2_vee_fuse = 0;
CD0_vee2_fuse = 0;
CD1_vee2_fuse = 0;
CD2_vee2_fuse = 0;

% initailizes induced drag coefficient 
K_w1       = 0;
K_HTP       = 0;
K_can       = 0;
K_vee    = 0;
K_vee2    = 0;


% Select according the airplane configuration
if Conf.w1 == 1
    S_w1 = Geo_tier.S_w1;
    AR_w1 = Geo_tier.AR_w1;
    AR_w1_e = Geo_tier.AR_w1_e;
    S_w1_e = Geo_tier.S_w1_e;
    cR_w1 = Geo_tier.cR_w1;
    x_w1_LE = Geo_tier.x_w1_LE;
    cmac_w1 = Geo_tier.cmac_w1;
    wingspan2bodydiam = b_w1/w_Area_b_max;

    %%%%%%%%%%%%%%%
    % Wing 1
    % Dimensions
    % Percentage of laminar and turbulent
    c_w = cmac_w1;
    % Reynolds
    Re_w_lm = rho*V*c_w/(mu);
    % Cutoff Reynolds
    k_sfr_w = K_corr*k_sfr5;
    R_cutoff_w = 38.21*(c_w*m2ft/k_sfr_w)^1.053;
    
    % selects the smaller of both Re
    Re_w_tb = min(Re_w_lm,R_cutoff_w);
    
    % Percentage of laminar and turbulent
    lam_w = 0.25;
    tur_w = 1 - lam_w;
    
    % laminar flow
    Cf_lam_w = 1.328/sqrt(Re_w_lm);
    % turbulent flow
    Cf_turb_w = 0.455/(((log10(Re_w_tb))^2.58)*(1+0.144*M^2));
    
    x_c_w = 0.2917; % Ubicación del máximo grosor del perfil con respecto a la cuerda
    t_c_w = 0.1236; %Extraer del perfil
    
    Lambda_max_t_w = 0 ; % Flecha en el máximo espesor
    FF_w = (1+ (0.6/(x_c_w))*t_c_w + 100*(t_c_w^4))*((1.34*M^0.18)*(cos(Lambda_max_t_w))^0.28);
    % Factor de interferencia
    Qc_w = 1.05;
    S_wet_w = 2*S_w1_e;
    CD0_w1 = (lam_w*Cf_lam_w + tur_w*Cf_turb_w)*FF_w*Qc_w*S_wet_w/S_ref;
    Aero_TH.CD0_w1    = CD0_w1;
    
    % Oswalds coefficient
    e_w1 = 1.78*(1-0.045*AR_w1^(0.68)) - 0.64;
    K_w1 = (S_w1_e/S_ref)*(1/(pi*e_w1*AR_w1));

    % Fuse aerodynamic models obtained from CMB and numerical methods
    if fuse_aero_FLOW_and_CBM == 1 || 2 
        if wingspan2bodydiam <= 2
            CD0_w1_fuse = (S_w1_e/S_ref)*Aero.CD0_w1;
            CD1_w1_fuse = (S_w1_e/S_ref)*Aero.CD1_w1;
            CD2_w1_fuse = (S_w1_e/S_ref)*Aero.CD2_w1;
        else
            CD0_w1_fuse = (S_w1/S_ref)*Aero.CD0_w1;
            CD1_w1_fuse = (S_w1/S_ref)*Aero.CD1_w1;
            CD2_w1_fuse = (S_w1/S_ref)*Aero.CD2_w1;
        end
        
    else
        CD0_w1_fuse = 0;
        CD1_w1_fuse = 0;
        CD2_w1_fuse = 0;        
    end
else
    Aero_TH.CD0_w1    = 0;
end      

% Select according the airplane configuration
if Conf.can == 1
    S_can = Geo_tier.S_can;
    b_can = Geo_tier.b_can;
    AR_can = Geo_tier.AR_can;
    AR_can_e = Geo_tier.AR_can_e;
    S_can_e = Geo_tier.S_can_e;
    cR_can = Geo_tier.cR_can;
    x_can_LE = Geo_tier.x_can_LE;
    cmac_can = Geo_tier.cmac_can;
    wingspan2bodydiam = b_can/w_Area_b_max;

    %%%%%%%%%%%%%%%
    % Wing 1
    % Dimensions
    % Percentage of laminar and turbulent
    c_w = cmac_can;
    % Reynolds
    Re_w_lm = rho*V*c_w/(mu);
    % Cutoff Reynolds
    k_sfr_w = K_corr*k_sfr5;
    R_cutoff_w = 38.21*(c_w*m2ft/k_sfr_w)^1.053;
    
    % selects the smaller of both Re
    Re_w_tb = min(Re_w_lm,R_cutoff_w);
    
    % Percentage of laminar and turbulent
    lam_w = 0.25;
    tur_w = 1 - lam_w;
    
    % laminar flow
    Cf_lam_w = 1.328/sqrt(Re_w_lm);
    % turbulent flow
    Cf_turb_w = 0.455/(((log10(Re_w_tb))^2.58)*(1+0.144*M^2));
    
    x_c_w = 0.2917; % Ubicación del máximo grosor del perfil con respecto a la cuerda
    t_c_w = 0.1236;
    Lambda_max_t_w = 0 ; % Flecha en el máximo espesor
    FF_w = (1+ (0.6/(x_c_w))*t_c_w + 100*(t_c_w^4))*((1.34*M^0.18)*(cos(Lambda_max_t_w))^0.28);
    % Factor de interferencia
    Qc_w = 1.05;
    S_wet_w = 2*S_can_e;
    CD0_can = (lam_w*Cf_lam_w + tur_w*Cf_turb_w)*FF_w*Qc_w*S_wet_w/S_ref;
    Aero_TH.CD0_can    = CD0_can;
    
    % Oswalds coefficient
    e_can = 1.78*(1-0.045*AR_can^(0.68)) - 0.64;
    K_can = (S_can_e/S_ref)*(1/(pi*e_can*AR_can));

    % Fuse aerodynamic models obtained from CMB and numerical methods
    if fuse_aero_FLOW_and_CBM == 1
        if wingspan2bodydiam <= 2
            CD0_can_fuse = (S_can_e/S_ref)*Aero.CD0_can;
            CD1_can_fuse = (S_can_e/S_ref)*Aero.CD1_can;
            CD2_can_fuse = (S_can_e/S_ref)*Aero.CD2_can;
        else
            CD0_can_fuse = (S_can/S_ref)*Aero.CD0_can;
            CD1_can_fuse = (S_can/S_ref)*Aero.CD1_can;
            CD2_can_fuse = (S_can/S_ref)*Aero.CD2_can;
        end
        
    else
        CD0_can_fuse = 0;
        CD1_can_fuse = 0;
        CD2_can_fuse = 0;        
    end
else
    Aero_TH.CD0_can    = 0;
end   

if Conf.h == 1
    S_HTP = Geo_tier.S_HTP;
    b_HTP  = Geo_tier.b_HTP;
    AR_HTP = Geo_tier.AR_HTP;
    AR_HTP_e = Geo_tier.AR_HTP_e;
    S_HTP_e = Geo_tier.S_HTP_e;
    S_HTP_s = Geo_tier.S_HTP_s;
    cR_HTP = Geo_tier.cR_HTP;
    x_HTP_LE = Geo_tier.x_HTP_LE;
    cmac_HTP = Geo_tier.cmac_HTP;
    wingspan2bodydiam = b_HTP/w_Area_b_max;
    %%%%%%%%%%%%%%%
    % Wing 2
    % Dimensions
    % Percentage of laminar and turbulent
    c_w = cmac_HTP;
    % Reynolds
    Re_w_lm = rho*V*c_w/(mu);
    % Cutoff Reynolds
    k_sfr_w = K_corr*k_sfr5;
    R_cutoff_w = 38.21*(c_w*m2ft/k_sfr_w)^1.053;
    
    % selects the smaller of both Re
    Re_w_tb = min(Re_w_lm,R_cutoff_w);
    
    % Percentage of laminar and turbulent
    lam_w = 0.25;
    tur_w = 1 - lam_w;
    
    % laminar flow
    Cf_lam_w = 1.328/sqrt(Re_w_lm);
    % turbulent flow
    Cf_turb_w = 0.455/(((log10(Re_w_tb))^2.58)*(1+0.144*M^2));
    
    x_c_w = 0.25; % Ubicación del máximo grosor del perfil con respecto a la cuerda
    t_c_w = 0.12;
    Lambda_max_t_w = 0 ; % Flecha en el máximo espesor
    FF_w = (1+ (0.6/(x_c_w))*t_c_w + 100*(t_c_w^4))*((1.34*M^0.18)*(cos(Lambda_max_t_w))^0.28);
    Qc_w = 1.05;
    S_wet_w = 2*S_HTP_s;
    CD0_HTP = (lam_w*Cf_lam_w + tur_w*Cf_turb_w)*FF_w*Qc_w*S_wet_w/S_ref;
    Aero_TH.CD0_HTP    = CD0_HTP;
    
    e_HTP = 1.78*(1-0.045*AR_HTP^(0.68)) - 0.64;
    K_HTP = (S_HTP_e/S_ref)*(1/(pi*e_HTP*AR_HTP));

    % Fuse aerodynamic models obtained from CMB and numerical methods
    if fuse_aero_FLOW_and_CBM == 1
        if wingspan2bodydiam <= 2
            CD0_HTP_fuse = (S_HTP_e/S_ref)*Aero.CD0_HTP;
            CD1_HTP_fuse = (S_HTP_e/S_ref)*Aero.CD1_HTP;
            CD2_HTP_fuse = (S_HTP_e/S_ref)*Aero.CD2_HTP;
        else
            CD0_HTP_fuse = (S_HTP/S_ref)*Aero.CD0_HTP;
            CD1_HTP_fuse = (S_HTP/S_ref)*Aero.CD1_HTP;
            CD2_HTP_fuse = (S_HTP/S_ref)*Aero.CD2_HTP;
        end
        
    else
        CD0_HTP_fuse = 0;
        CD1_HTP_fuse = 0;
        CD2_HTP_fuse = 0;
    end
else
    Aero_TH.CD0_HTP    = 0;
end

if Conf.v == 1
    S_VTP = Geo_tier.S_VTP;
    AR_VTP = Geo_tier.AR_VTP;
    AR_VTP_e = Geo_tier.AR_VTP_e;
    S_VTP_e = Geo_tier.S_VTP_e;
    S_VTP_s = Geo_tier.S_VTP_s;
    cR_VTP = Geo_tier.cR_VTP;
    x_VTP_LE = Geo_tier.x_VTP_LE;
    cmac_VTP = Geo_tier.cmac_VTP;
    %%%%%%%%%%%%%%%
    % Vertical
    % Dimensions
    % Percentage of laminar and turbulent
    c_v = cmac_VTP;
    
    % Reynolds
    Re_v_lm = rho*V*c_v/(mu);
    
    % Cutoff Reynolds
    k_sfr_v = K_corr*k_sfr5;
    R_cutoff_v = 38.21*(c_v*m2ft/k_sfr_v)^1.053;
    
    % selects the smaller of both Re
    Re_v_tb = min(Re_v_lm,R_cutoff_v);
    
    % Percentage of laminar and turbulent
    lam_v = 0.25;
    tur_v = 1 - lam_v;
    
    % laminar flow
    Cf_lam_v = 1.328/sqrt(Re_v_lm);
    % turbulent flow
    Cf_turb_v = 0.455/(((log10(Re_v_tb))^2.58)*(1+0.144*M^2));
    
    x_c_v = 0.40; % Ubicación del máximo grosor del perfil con respecto a la cuerda
    t_c_v = 0.12;
    Lambda_max_t_v = 0 ; % Flecha en el máximo espesor
    FF_v = (1+ (0.6/(x_c_v))*t_c_v + 100*(t_c_v^4))*((1.34*M^0.18)*(cos(Lambda_max_t_v))^0.28);
    Qc_v = 1.08;
    S_wet_v = 2*S_VTP;
    CD0_v = (lam_v*Cf_lam_v + tur_v*Cf_turb_v)*FF_v*Qc_v*S_wet_v/S_ref;
    Aero_TH.CD0_v = CD0_v;
else
    Aero_TH.CD0_VTP    = 0;
end 

if Conf.v2 == 1
    S_VTP = Geo_tier.S_VTP;
    AR_VTP = Geo_tier.AR_VTP;
    AR_VTP_e = Geo_tier.AR_VTP_e;
    S_VTP_e = Geo_tier.S_VTP_e;
    S_VTP_s = Geo_tier.S_VTP_s;
    cR_VTP = Geo_tier.cR_VTP;
    x_VTP_LE = Geo_tier.x_VTP_LE;
    cmac_VTP = Geo_tier.cmac_VTP;
 
    %%%%%%%%%%%%%%%
    % Vertical
    % Dimensions
    % Percentage of laminar and turbulent
    c_v = cmac_VTP;

    % Reynolds
    Re_v_lm = rho*V*c_v/(mu);
    
    % Cutoff Reynolds
    k_sfr_v = K_corr*k_sfr5;
    R_cutoff_v = 38.21*(c_v*m2ft/k_sfr_v)^1.053;
    
    % selects the smaller of both Re
    Re_v_tb = min(Re_v_lm,R_cutoff_v);
    
    % Percentage of laminar and turbulent
    lam_v = 0.25;
    tur_v = 1 - lam_v;
    
    % laminar flow
    Cf_lam_v = 1.328/sqrt(Re_v_lm);
    % turbulent flow
    Cf_turb_v = 0.455/(((log10(Re_v_tb))^2.58)*(1+0.144*M^2));
    
    x_c_v = 0.40; % Ubicación del máximo grosor del perfil con respecto a la cuerda
    t_c_v = 0.12;
    Lambda_max_t_v = 0 ; % Flecha en el máximo espesor
    FF_v = (1+ (0.6/(x_c_v))*t_c_v + 100*(t_c_v^4))*((1.34*M^0.18)*(cos(Lambda_max_t_v))^0.28);
    Qc_v = 1.08;
    S_wet_v = 4*S_VTP; % Twin
    CD0_v2 = (lam_v*Cf_lam_v + tur_v*Cf_turb_v)*FF_v*Qc_v*S_wet_v/S_ref;
    Aero_TH.CD0_v2 = CD0_v2;
else
    Aero_TH.CD0_2VTP    = 0;
end

if Conf.vtail == 1
    S_vee = Geo_tier.S_vee;
    b_vee  = Geo_tier.b_vee;
    AR_vee = Geo_tier.AR_vee;
    AR_vee_e = Geo_tier.AR_vee_e;
    S_vee_e = Geo_tier.S_vee_e;
    S_vee_s = Geo_tier.S_vee_s;
    cR_vee = Geo_tier.cR_vee;
    x_vee_LE = Geo_tier.x_vee_LE;
    cmac_vee = Geo_tier.cmac_vee;
    wingspan2bodydiam = b_vee/w_Area_b_max;

    %%%%%%%%%%%%%%%
    % V-Tail
    % Dimensions
    % Percentage of laminar and turbulent
    c_vee = cmac_vee;
    
    % Reynolds
    Re_vee_lm = rho*V*c_vee/(mu);
    
    % Cutoff Reynolds
    k_sfr_vee = K_corr*k_sfr5;
    R_cutoff_vee = 38.21*(c_vee*m2ft/k_sfr_vee)^1.053;
    
    % selects the smaller of both Re
    Re_vee_tb = min(Re_vee_lm,R_cutoff_vee);
    
    % Percentage of laminar and turbulent
    lam_vee = 0.25;
    tur_vee = 1 - lam_vee;
    
    % laminar flow
    Cf_lam_vee = 1.328/sqrt(Re_vee_lm);
    % turbulent flow
    Cf_turb_vee = 0.455/(((log10(Re_vee_tb))^2.58)*(1+0.144*M^2));
    
    x_c_vee = 0.25; % Ubicación del máximo grosor del perfil con respecto a la cuerda
    t_c_vee = 0.12;
    Lambda_max_t_vee = 0 ; % Flecha en el máximo espesor
    FF_vee = (1+ (0.6/(x_c_vee))*t_c_vee + 100*(t_c_vee^4))*((1.34*M^0.18)*(cos(Lambda_max_t_vee))^0.28);
    Qc_vee = 1.03;
    S_wet_vee = 2*S_vee_s;

    CD0_vee = (lam_vee*Cf_lam_vee + tur_vee*Cf_turb_vee)*FF_vee*Qc_vee*S_wet_vee/S_ref;
    Aero_TH.CD0_vee = CD0_vee;
    
    e_vee = 1.78*(1-0.045*AR_vee_e^(0.68)) - 0.64;
    K_vee = (S_vee_e/S_ref)*(1/(pi*e_vee*AR_vee_e));

    % Fuse aerodynamic models obtained from CMB and numerical methods
    if fuse_aero_FLOW_and_CBM == 1
        if wingspan2bodydiam <= 2
            CD0_vee_fuse = (S_vee_e/S_ref)*Aero.CD0_vee;
            CD1_vee_fuse = (S_vee_e/S_ref)*Aero.CD1_vee;
            CD2_vee_fuse = (S_vee_e/S_ref)*Aero.CD2_vee;
        else
            CD0_vee_fuse = (S_vee/S_ref)*Aero.CD0_vee;
            CD1_vee_fuse = (S_vee/S_ref)*Aero.CD1_vee;
            CD2_vee_fuse = (S_vee/S_ref)*Aero.CD2_vee;
        end
        
    else
        CD0_vee_fuse = 0;
        CD1_vee_fuse = 0;
        CD2_vee_fuse = 0;
    end
else
    Aero_TH.CD0_vee    = 0;
end 

if Conf.vtail2 == 1
    S_vee2 = Geo_tier.S_vee2;
    b_vee2  = Geo_tier.b_vee2;
    AR_vee2 = Geo_tier.AR_vee2;
    AR_vee2_e = Geo_tier.AR_vee2_e;
    S_vee2_e = Geo_tier.S_vee2_e;
    S_vee2_s = Geo_tier.S_vee2_s;
    cR_vee2 = Geo_tier.cR_vee2;
    x_vee2_LE = Geo_tier.x_vee2_LE;
    cmac_vee2 = Geo_tier.cmac_vee2;
    wingspan2bodydiam = b_vee2/w_Area_b_max;

    %%%%%%%%%%%%%%%
    % V-Tail
    % Dimensions
    % Percentage of laminar and turbulent
    c_vee2 = cmac_vee2;
    
    % Reynolds
    Re_vee2_lm = rho*V*c_vee2/(mu);
    
    % Cutoff Reynolds
    k_sfr_vee2 = K_corr*k_sfr5;
    R_cutoff_vee2 = 38.21*(c_vee2*m2ft/k_sfr_vee2)^1.053;
    
    % selects the smaller of both Re
    Re_vee2_tb = min(Re_vee2_lm,R_cutoff_vee2);
    
    % Percentage of laminar and turbulent
    lam_vee2 = 0.25;
    tur_vee2 = 1 - lam_vee2;
    
    % laminar flow
    Cf_lam_vee2 = 1.328/sqrt(Re_vee2_lm);
    % turbulent flow
    Cf_turb_vee2 = 0.455/(((log10(Re_vee2_tb))^2.58)*(1+0.144*M^2));
    
    x_c_vee2 = 0.25; % Ubicación del máximo grosor del perfil con respecto a la cuerda
    t_c_vee2 = 0.12;
    Lambda_max_t_vee2 = 0 ; % Flecha en el máximo espesor
    FF_vee2 = (1+ (0.6/(x_c_vee2))*t_c_vee2 + 100*(t_c_vee2^4))*((1.34*M^0.18)*(cos(Lambda_max_t_vee2))^0.28);
    Qc_vee2 = 1.03;
    S_wet_vee2 = 2*S_vee2_s;

    CD0_vee2 = (lam_vee2*Cf_lam_vee2 + tur_vee2*Cf_turb_vee2)*FF_vee2*Qc_vee2*S_wet_vee2/S_ref;
    Aero_TH.CD0_vee2 = CD0_vee2;
    
    e_vee2 = 1.78*(1-0.045*AR_vee2_e^(0.68)) - 0.64;
    K_vee2 = (S_vee2_e/S_ref)*(1/(pi*e_vee2*AR_vee2_e));

    % Fuse aerodynamic models obtained from CMB and numerical methods
    if fuse_aero_FLOW_and_CBM == 1
        if wingspan2bodydiam <= 2
            CD0_vee2_fuse = (S_vee2_e/S_ref)*Aero.CD0_vee2;
            CD1_vee2_fuse = (S_vee2_e/S_ref)*Aero.CD1_vee2;
            CD2_vee2_fuse = (S_vee2_e/S_ref)*Aero.CD2_vee2;
        else
            CD0_vee2_fuse = (S_vee2/S_ref)*Aero.CD0_vee2;
            CD1_vee2_fuse = (S_vee2/S_ref)*Aero.CD1_vee2;
            CD2_vee2_fuse = (S_vee2/S_ref)*Aero.CD2_vee2;
        end
        
    else
        CD0_vee2_fuse = 0;
        CD1_vee2_fuse = 0;
        CD2_vee2_fuse = 0;
    end
else
    Aero_TH.CD0_vee2    = 0;
end 

if Conf.fus == 1
    % Read data from the geometry file
    % [Body_Geo] = Calc_Body_Geometry;
    Surf_TOT = Body_Geo.Surf_TOT;
    
    %%%%%%%%%%%%%%%
    % Fuselaje
    % Dimensions
    w_fus = Geo_tier.w_fus;
    h_fus = Geo_tier.h_fus;
    l_fus = Geo_tier.l_fus;
    d_fus = Geo_tier.d_fus;
    
    % Reynolds
    Re_fus_lm = rho*V*l_fus/(mu);
    
    % Cutoff Reynolds
    k_sfr_fus = K_corr*k_sfr5; % smooth molded composite
    R_cutoff_fus = 38.21*(l_fus*m2ft/k_sfr_fus)^1.053;
    
    % selects the smaller of both Re
    Re_fus_tb = min(Re_fus_lm,R_cutoff_fus);
    
    % Percentage of laminar and turbulent
    lam_fus = 0.25;
    tur_fus = 1 - lam_fus;
    
    % laminar flow
    Cf_lam_fus = 1.328/sqrt(Re_fus_lm);
    % turbulent flow
    Cf_turb_fus = 0.455/(((log10(Re_fus_tb))^2.58)*(1+0.144*M^2));
    
    f_fus = l_fus/d_fus;
    FF_fus = 1+ (60/(f_fus^3)) + (f_fus/400);
    Qc_fus = 1;
    S_wet_fus = Surf_TOT;
    CD0_fus = (lam_fus*Cf_lam_fus + tur_fus*Cf_turb_fus)*FF_fus*Qc_fus*S_wet_fus/S_ref;

    if Conf.m_fus == 1
        n_m_fus = Conf.n_m_fus;% number of multiple
        CD0_fus = CD0_fus*n_m_fus;
    end
    Aero_TH.CD0_fus    = CD0_fus;
else
    Aero_TH.CD0_fus    = 0;
end 

if Conf.nac == 1
    
    %%%%%%%%%%%%%%%
    % Engine Nacelles
    % Dimensions
    % Percentage of laminar and turbulent
    c_nac = AC_CONFIGURATION.l_nc;
    
    % Reynolds
    Re_nac_lm = rho*V*c_nac/(mu);
    
    % Cutoff Reynolds
    k_sfr_v = K_corr*k_sfr5;
    R_cutoff_nac = 38.21*(c_nac*m2ft/k_sfr_v)^1.053;
    
    % selects the smaller of both Re
    Re_nac_tb = min(Re_nac_lm,R_cutoff_nac);
    
    % Percentage of laminar and turbulent
    lam_nac = 0.05;
    tur_nac = 1 - lam_nac;
    
    % laminar flow
    Cf_lam_nac = 1.328/sqrt(Re_nac_lm);
    % turbulent flow
    Cf_turb_nac = 0.455/(((log10(Re_nac_tb))^2.58)*(1+0.144*M^2));
    
    f_nac = l_nc/d_nc;
    FF_nac = 1+ (0.35/(f_nac));
    Qc_nac = 1.30;
    S_wet_nac = 2*pi*l_nc*(d_nc/2) + 2*pi*((d_nc/2)^2);
    CD0_nac = Num_eng*(lam_nac*Cf_lam_nac + tur_nac*Cf_turb_nac)*FF_nac*Qc_nac*S_wet_nac/S_ref;
    Aero_TH.CD0_nac    = CD0_nac;
else
    Aero_TH.CD0_nac    = 0;
end


% Tailboom geometry
        
if Conf.tailboom == 1
    l_tailboom = OUTPUT_read_XLSX.InputGeometry_Data_flags.l_tailboom;
    w_tailboom = OUTPUT_read_XLSX.InputGeometry_Data_flags.w_tailboom;
    h_tailboom = OUTPUT_read_XLSX.InputGeometry_Data_flags.h_tailboom;
        
    Surf_TOT = 2*l_tailboom*w_tailboom +  2*l_tailboom*h_tailboom + 2*h_tailboom*h_tailboom;
    
    % Reynolds
    Re_tb_lm = rho*V*l_tailboom/(mu);
    
    % Cutoff Reynolds
    k_sfr_tb = K_corr*k_sfr5; % smooth molded composite
    R_cutoff_tb = 38.21*(l_tailboom*m2ft/k_sfr_tb)^1.053;
    
    % selects the smaller of both Re
    Re_tb_tb = min(Re_tb_lm,R_cutoff_tb);
    
    % Percentage of laminar and turbulent
    lam_tb = 0.0;
    tur_tb = 1 - lam_tb;
    
    % laminar flow
    Cf_lam_tb = 1.328/sqrt(Re_tb_lm);
    % turbulent flow
    Cf_turb_tb = 0.455/(((log10(Re_tb_tb))^2.58)*(1+0.144*M^2));
    
    f_tb = l_tailboom/((w_tailboom+h_tailboom)/2);
    FF_tb = 1+ (60/(f_tb^3)) + (f_tb/400);
    Qc_tb = 1;
    S_wet_tb = Surf_TOT;
    CD0_tb = (lam_tb*Cf_lam_tb + tur_tb*Cf_turb_tb)*FF_tb*Qc_tb*S_wet_tb/S_ref;

    if Conf.m_tailboom == 1
        n_m_tb = Conf.n_m_tailboom;% number of multiple
        CD0_tb = CD0_tb*n_m_tb;
    end
    Aero_TH.CD0_tailboom    = CD0_tb;
else
    Aero_TH.CD0_tailboom    = 0;
end


if Conf.pod == 1
    % Read data from the geometry file
    
    %%%%%%%%%%%%%%%
    % Fuselaje
    % Dimensions
    w_fus = Geo_tier.w_fus;
    h_fus = Geo_tier.h_fus;
    l_fus = Geo_tier.l_fus;
    d_fus = Geo_tier.d_fus;

    l_pod = OUTPUT_read_XLSX.InputGeometry_Data_flags.l_pod; % Length of pod
    w_pod = OUTPUT_read_XLSX.InputGeometry_Data_flags.w_pod; % Width of pod
    h_pod = OUTPUT_read_XLSX.InputGeometry_Data_flags.w_pod;% Height of pod

    Surf_TOT = 2*l_pod*w_pod +  2*l_pod*h_pod + 2*h_pod*w_pod;

    % Reynolds
    Re_pod_lm = rho*V*l_pod/(mu);
    
    % Cutoff Reynolds
    k_sfr_pod = K_corr*k_sfr5; % smooth molded composite
    R_cutoff_pod = 38.21*(l_pod*m2ft/k_sfr_pod)^1.053;
    
    % selects the smaller of both Re
    Re_pod_tb = min(Re_pod_lm,R_cutoff_pod);
    
    % Percentage of laminar and turbulent
    lam_pod = 0.0;
    tur_pod = 1 - lam_pod;
    
    % laminar flow
    Cf_lam_pod = 1.328/sqrt(Re_pod_lm);
    % turbulent flow
    Cf_turb_pod = 0.455/(((log10(Re_pod_tb))^2.58)*(1+0.144*M^2));
    
    d_pod = (w_pod + h_pod)/2;
    
    f_pod = l_pod/d_pod;
    FF_pod = 1+ (60/(f_pod^3)) + (f_pod/400);
    Qc_pod = 1;
    S_wet_pod = Surf_TOT;
    CD0_pod = (lam_pod*Cf_lam_pod + tur_pod*Cf_turb_pod)*FF_pod*Qc_pod*S_wet_pod/S_ref;

    n_m_pod = 2;% number of multiple pods assume to be one on each wing
    CD0_pod = CD0_pod*n_m_pod;
else
    Aero_TH.CD0_pod    = 0;
end 

if Conf.landgear == 1
    CD0_lnd = 0.0011;
    % type AC: GA 1, medium/small 2, large transport 3
    type_AC_lnd = 2;
    if  type_AC_lnd == 1
        CD0_lnd = 0.01;
    elseif type_AC_lnd == 2
        CD0_lnd = 0.015;
    elseif type_AC_lnd == 3
        CD0_lnd = 0.025;
    end
else
    CD0_lnd = 0;
end

% Missile Configuration
if Conf.missile == 1
    if conditions.T120 == 1 % Tarsis 120
        S_ref_missile = 0.0021;
        CD0_missile_ref = 0.4594; % Value from AERTEC
        CD0_missile = Conf.n_missile*CD0_missile_ref*S_ref_missile/S_ref;
    elseif conditions.T120 == 2 || conditions.T120 == 3 % Tarsis 120 with KSA
        S_ref_KSA = 0.0095; % FROM Defensa
        Delta0 = 28.23; % Delta in Figure to determine increment in CD0
        Delta0CD0 = 0.2; % Value of CD for the associated Delta0
        if M<=0.1
            Delta1 = 20.93; % MEasured from figure provided by Defensa
            CD0_missile_ref = Delta0CD0 + Delta0CD0*(Delta1/Delta0); % Value from Defensa Graph
        elseif M<=0.2
            Delta1 = 16.69; % MEasured from figure provided by Defensa
            CD0_missile_ref = Delta0CD0 + Delta0CD0*(Delta1/Delta0); % Value from Defensa Graph
        elseif M<=0.3
            Delta1 = 14.85;
            CD0_missile_ref = Delta0CD0 + Delta0CD0*(Delta1/Delta0); % Value from Defensa Graph
        elseif M<=0.4
            Delta1 = 14.30; % MEasured from figure provided by Defensa
            CD0_missile_ref = Delta0CD0 + Delta0CD0*(Delta1/Delta0); % Value from Defensa Graph
        elseif M>=0.5
            Delta1 = 13.75; % MEasured from figure provided by Defensa
            CD0_missile_ref = Delta0CD0 + Delta0CD0*(Delta1/Delta0); % Value from Defensa Graph
        end
        CD0_missile = Conf.n_missile*CD0_missile_ref*S_ref_KSA/S_ref;
    end
else
    CD0_missile = 0;
end

% Landing Gear
if Conf.flap == 1
    dflap_TO = 20;
    dflap_LD = 40;
    F_flap = 0.0074;

    cf_flap = Geo_tier.cf_flap;
    K_y1_flap_w1 = Geo_tier.K_y1_flap_w1;
    K_y2_flap_w1 = Geo_tier.K_y2_flap_w1;
    K_y1y2_flap_w1 = K_y2_flap_w1 - K_y1_flap_w1;

    CD_flap_TO = F_flap*cf_flap*K_y1y2_flap_w1*(dflap_TO-10); 
    CD_flap_LN = F_flap*cf_flap*K_y1y2_flap_w1*(dflap_LD-10);
else
    CD_flap_TO = 0; 
    CD_flap_LN = 0;
end

% Correction CDi in ground effect
CDi_IGE = (16*(h/b_w1)^2)/(1 + 15*(h/b_w1)^2);
CDi_IGE = 1; % eliminates the ground effect correction

%% Results only from CBM
% CD0_fus=0;

CD0_par = CD0_w1 + CD0_can + CD0_HTP + CD0_v + CD0_v2 + CD0_vee + CD0_vee2 + CD0_fus + CD0_tb + CD0_nac + CD0_missile + CD0_pod;
CD0_misc = CD0_lnd;%[-] Coef de resistencia en configuración de crucero - miscelaneos
CD0_lp = (CD0_par + CD0_misc)*0.10; %[-] Coef de resistencia en configuración de crucero miscelanea

CD0_CBM = CD0_par + CD0_misc + CD0_lp;
CD1_CBM = 0;
CD2_CBM = CDi_IGE*(K_w1 + K_can + K_HTP + K_vee);

% Actualización de los resultados de FLOW con los componentes no analizados
Aero.CD0_ac = Aero.CD0_ac + CD0_v + CD0_v2 + CD0_vee+ CD0_vee2 + CD0_fus + CD0_tb + CD0_nac + CD0_lnd + CD0_missile + CD0_misc + CD0_pod + CD0_lp;

% Polar for the rest of elements - una vez obtenidos los resultados de
% Fuse aerodynamic models obtained from CMB and numerical methods
if fuse_aero_FLOW_and_CBM == 1

    % CRUISE CONDITIONS
    CD0_par = CD0_w1_fuse + CD0_can_fuse + CD0_HTP_fuse + CD0_v + CD0_v2 + CD0_vee_fuse+ CD0_vee2_fuse + CD0_fus + CD0_tb + CD0_nac + CD0_pod +  CD0_lnd + CD0_missile;
    CD0_misc = 0;%[-] Coef de resistencia en configuración de crucero - miscelaneos
    CD0_lp = (CD0_par + CD0_misc)*0.10; %[-] Coef de resistencia en configuración de crucero miscelanea

    % Total Drag coefficients 
    CD0 = CD0_par + CD0_misc + CD0_lp;    
    CD1 = CD1_w1_fuse + CD1_can_fuse + CD1_HTP_fuse + CD1_vee_fuse + CD1_vee2_fuse;
    CD2 = CDi_IGE*(CD2_w1_fuse + CD2_can_fuse + CD2_HTP_fuse + CD2_vee_fuse + CD2_vee2_fuse);

    % TAKE OFF
    CD0_par_TO = CD0_w1_fuse + CD0_can_fuse + CD0_HTP_fuse + CD0_v + CD0_v2 + CD0_vee_fuse+ CD0_vee2_fuse + CD0_fus + CD0_tb + CD0_nac + CD0_pod +  CD0_lnd + CD0_missile;
    CD0_misc_TO = CD_flap_TO;%[-] Coef de resistencia en configuración de crucero - miscelaneos
    CD0_lp_TO = (CD0_par_TO + CD0_misc_TO)*0.10; %[-] Coef de resistencia en configuración de crucero miscelanea

    CD0_TO = CD0_par_TO + CD0_misc_TO + CD0_lp_TO;    
    CD1_TO = CD1_w1_fuse + CD1_can_fuse + CD1_HTP_fuse + CD1_vee_fuse + CD1_vee2_fuse;
    CD2_TO = CDi_IGE*(CD2_w1_fuse + CD2_can_fuse + CD2_HTP_fuse + CD2_vee_fuse + CD2_vee2_fuse);

    % LANDING
    CD0_par_LN = CD0_w1_fuse + CD0_can_fuse + CD0_HTP_fuse + CD0_v + CD0_v2 + CD0_vee_fuse+ CD0_vee2_fuse + CD0_fus + CD0_tb + CD0_nac + CD0_pod +  CD0_lnd + CD0_missile;
    CD0_misc_LN = CD_flap_LN;%[-] Coef de resistencia en configuración de crucero - miscelaneos
    CD0_lp_LN = (CD0_par_LN + CD0_misc_LN)*0.10; %[-] Coef de resistencia en configuración de crucero miscelanea

    CD0_LN = CD0_par_LN + CD0_misc_LN + CD0_lp_LN;    
    CD1_LN = CD1_w1_fuse + CD1_can_fuse + CD1_HTP_fuse + CD1_vee_fuse + CD1_vee2_fuse;
    CD2_LN = CDi_IGE*(CD2_w1_fuse + CD2_can_fuse + CD2_HTP_fuse + CD2_vee_fuse + CD2_vee2_fuse);

elseif fuse_aero_FLOW_and_CBM == 2 % If only CBM results is used

    CD0 = CD0_CBM;
    CD1 = CD1_w1_fuse + CD1_can_fuse + CD1_HTP_fuse + CD1_vee_fuse + CD1_vee2_fuse;
    CD2 = CDi_IGE*(CD2_w1_fuse + CD2_can_fuse + CD2_HTP_fuse + CD2_vee_fuse + CD2_vee2_fuse);

    CD0_TO = CD0_CBM;
    CD1_TO = CD1_w1_fuse + CD1_can_fuse + CD1_HTP_fuse + CD1_vee_fuse + CD1_vee2_fuse;
    CD2_TO = CDi_IGE*(CD2_w1_fuse + CD2_can_fuse + CD2_HTP_fuse + CD2_vee_fuse + CD2_vee2_fuse);

    CD0_LN = CD0_CBM;
    CD1_LN = CD1_w1_fuse + CD1_can_fuse + CD1_HTP_fuse + CD1_vee_fuse + CD1_vee2_fuse;
    CD2_LN = CDi_IGE*(CD2_w1_fuse + CD2_can_fuse + CD2_HTP_fuse + CD2_vee_fuse + CD2_vee2_fuse);
   
else

    CD0 = CD0_CBM;
    CD1 = CD1_CBM;
    CD2 = CD2_CBM;
    
    % TAKE OFF
    CD0_par_TO = CD0;
    CD0_misc_TO = CD_flap_TO; %[-] Coef de resistencia en configuración de crucero - miscelaneos
    CD0_lp_TO = (CD0_par_TO + CD0_misc_TO)*0.10; %[-] Coef de resistencia en configuración de crucero miscelanea

    CD0_TO = CD0_par_TO + CD0_misc_TO + CD0_lp_TO;    
    CD1_TO = CD1_w1_fuse + CD1_can_fuse + CD1_HTP_fuse + CD1_vee_fuse + CD1_vee2_fuse;
    CD2_TO = CDi_IGE*(CD2_w1_fuse + CD2_can_fuse + CD2_HTP_fuse + CD2_vee_fuse + CD2_vee2_fuse);

    % LANDING
    CD0_par_LN = CD0;
    CD0_misc_LN = CD_flap_LN; %[-] Coef de resistencia en configuración de crucero - miscelaneos
    CD0_lp_LN = (CD0_par_LN + CD0_misc_LN)*0.10; %[-] Coef de resistencia en configuración de crucero miscelanea

    CD0_LN = CD0_par_LN + CD0_misc_LN + CD0_lp_LN;    
    CD1_LN = CD1_w1_fuse + CD1_can_fuse + CD1_HTP_fuse + CD1_vee_fuse + CD1_vee2_fuse;
    CD2_LN = CDi_IGE*(CD2_w1_fuse + CD2_can_fuse + CD2_HTP_fuse + CD2_vee_fuse + CD2_vee2_fuse);
end

% % Oswalds coefficient
% e_w1 = 1.78*(1-0.045*AR_w1^(0.68)) - 0.64;
% e_HTP = 1.78*(1-0.045*AR_HTP^(0.68)) - 0.64;
% K_w1 = 1/(pi*e_w1*AR_w1);
% K_HTP = 1/(pi*e_HTP*AR_HTP);

% CD2 =
% Almacenado de valores
Aero_TH.CD0_w1    = CD0_w1;
Aero_TH.CD0_can    = CD0_can;
Aero_TH.CD0_HTP    = CD0_HTP;
Aero_TH.CD0_v    = CD0_v;
Aero_TH.CD0_v2  = CD0_v2;
Aero_TH.CD0_vee  = CD0_vee;
Aero_TH.CD0_vee2  = CD0_vee2;
Aero_TH.CD0_fus  = CD0_fus;
Aero_TH.CD0_tb  = CD0_tb;
Aero_TH.CD0_nac  = CD0_nac;
Aero_TH.CD0_missile = CD0_missile;
Aero_TH.CD0_pod = CD0_pod;

Aero_TH.CD0_par  = CD0_par; %Si fuse_aero_FLOW_and_CBM == 1, está calculado con CD012_element_fuse, y falta CD0,1,2_HTP_fuse pq es 0.
Aero_TH.CD0_misc  = CD0_misc;
Aero_TH.CD0_lp  = CD0_lp;
Aero_TH.CD0  = CD0;
Aero_TH.CD1  = CD1;
Aero_TH.CD2  = CD2;

Aero_TH.CD0_CBM  = CD0_CBM;
Aero_TH.CD1_CBM  = CD1_CBM;
Aero_TH.CD2_CBM  = CD2_CBM;

% Cruise conditions
Aero_TH.CR.CD0 = CD0;
Aero_TH.CR.CD1 = CD1;
Aero_TH.CR.CD2 = CD2;

% Cruise conditions
Aero_TH.TO.CD0 = CD0_TO;
Aero_TH.TO.CD1 = CD1_TO;
Aero_TH.TO.CD2 = CD2_TO;

% Cruise conditions
Aero_TH.LND.CD0 = CD0_LN;
Aero_TH.LND.CD1 = CD1_LN;
Aero_TH.LND.CD2 = CD2_LN;

% Cruise conditions
Aero_TH.CR.K_w1 = K_w1;
Aero_TH.CR.K_can = K_can;
Aero_TH.CR.K_HTP = K_HTP;
Aero_TH.CR.K_vee = K_vee;
Aero_TH.CR.K_vee2 = K_vee2;