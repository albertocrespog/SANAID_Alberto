function  [Aero] = Polar_Literatura(Performance_preliminar,Geo_tier,conv_UNITS,Body_Geo,AC_CONFIGURATION,Conf)

h = Performance_preliminar.h;
M = Performance_preliminar.Mach;
rho = Performance_preliminar.rho;
V = Performance_preliminar.V;

m2ft = conv_UNITS.m2ft;

Num_eng = AC_CONFIGURATION.n_eng;
l_nc = AC_CONFIGURATION.l_nc;
d_nc = AC_CONFIGURATION.d_nc;

% Modélos teóricos de Aircraft Design - D.P. Raymer
mu = 0.000016938;
S_ref = Geo_tier.S_w1;

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
CD0_w2     = 0;
CD0_v      = 0;
CD0_v2     = 0;
CD0_vtail  = 0;
CD0_fus    = 0;
CD0_nac    = 0;
CD0_lnd    = 0;

% initailizes induced drag coefficient 
K_w1       = 0;
K_w2       = 0;
K_vtail    = 0;

% Select according the airplane configuration
if Conf.w1 == 1
    S_w1 = Geo_tier.S_w1;
    AR_w1 = Geo_tier.AR_w1;
    AR_w1_e = Geo_tier.AR_w1_e;
    S_w1_e = Geo_tier.S_w1_e;
    cR_w1 = Geo_tier.cR_w1;
    cR_w2 = Geo_tier.cR_w2;
    x_w1_LE = Geo_tier.x_w1_LE;
    cmac_w1 = Geo_tier.cmac_w1;
    
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
    t_c_w = 0.1236;
    Lambda_max_t_w = 0 ; % Flecha en el máximo espesor
    FF_w = (1+ (0.6/(x_c_w))*t_c_w + 100*(t_c_w^4))*((1.34*M^0.18)*(cos(Lambda_max_t_w))^0.28);
    % Factor de interferencia
    Qc_w = 1.05;
    S_wet_w = 2*S_w1_e;
    CD0_w1 = (lam_w*Cf_lam_w + tur_w*Cf_turb_w)*FF_w*Qc_w*S_wet_w/S_ref;
    Aero.CD0_w1    = CD0_w1;
    
    % Oswalds coefficient
    e_w1 = 1.78*(1-0.045*AR_w1^(0.68)) - 0.64;
    K_w1 = (S_w1_e/S_ref)*(1/(pi*e_w1*AR_w1));
end      
if Conf.h == 1
    S_w2 = Geo_tier.S_w2;
    AR_w2 = Geo_tier.AR_w2;
    AR_w2_e = Geo_tier.AR_w2_e;
    S_w2_e = Geo_tier.S_w2_e;
    S_w2_s = Geo_tier.S_w2_s;
    cR_w2 = Geo_tier.cR_w2;
    x_w2_LE = Geo_tier.x_w2_LE;
    cmac_w2 = Geo_tier.cmac_w2;
    
    %%%%%%%%%%%%%%%
    % Wing 2
    % Dimensions
    % Percentage of laminar and turbulent
    c_w = cmac_w2;
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
    S_wet_w = 2*S_w2_s;
    CD0_w2 = (lam_w*Cf_lam_w + tur_w*Cf_turb_w)*FF_w*Qc_w*S_wet_w/S_ref;
    Aero.CD0_w2    = CD0_w2;
    
    e_w2 = 1.78*(1-0.045*AR_w2^(0.68)) - 0.64;
    K_w2 = (S_w2_e/S_ref)*(1/(pi*e_w2*AR_w2));
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
    Aero.CD0_v = CD0_v;
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
    Aero.CD0_v2 = CD0_v2;
end
if Conf.vtail == 1
    S_w2 = Geo_tier.S_w2;
    AR_w2 = Geo_tier.AR_w2;
    AR_w2_e = Geo_tier.AR_w2_e;
    S_w2_e = Geo_tier.S_w2_e;
    S_w2_s = Geo_tier.S_w2_s;
    cR_w2 = Geo_tier.cR_w2;
    x_w2_LE = Geo_tier.x_w2_LE;
    cmac_w2 = Geo_tier.cmac_w2;
    
    %%%%%%%%%%%%%%%
    % V-Tail
    % Dimensions
    % Percentage of laminar and turbulent
    c_vtail = cmac_w2;
    
    % Reynolds
    Re_vtail_lm = rho*V*c_vtail/(mu);
    
    % Cutoff Reynolds
    k_sfr_vtail = K_corr*k_sfr5;
    R_cutoff_vtail = 38.21*(c_vtail*m2ft/k_sfr_vtail)^1.053;
    
    % selects the smaller of both Re
    Re_vtail_tb = min(Re_vtail_lm,R_cutoff_vtail);
    
    % Percentage of laminar and turbulent
    lam_vtail = 0.25;
    tur_vtail = 1 - lam_vtail;
    
    % laminar flow
    Cf_lam_vtail = 1.328/sqrt(Re_vtail_lm);
    % turbulent flow
    Cf_turb_vtail = 0.455/(((log10(Re_vtail_tb))^2.58)*(1+0.144*M^2));
    
    x_c_vtail = 0.25; % Ubicación del máximo grosor del perfil con respecto a la cuerda
    t_c_vtail = 0.12;
    Lambda_max_t_vtail = 0 ; % Flecha en el máximo espesor
    FF_vtail = (1+ (0.6/(x_c_vtail))*t_c_vtail + 100*(t_c_vtail^4))*((1.34*M^0.18)*(cos(Lambda_max_t_vtail))^0.28);
    Qc_vtail = 1.03;
    S_wet_vtail = 2*S_w2_s;

    CD0_vtail = (lam_vtail*Cf_lam_vtail + tur_vtail*Cf_turb_vtail)*FF_vtail*Qc_vtail*S_wet_vtail/S_ref;
    Aero.CD0_vtail = CD0_vtail;
    
    e_vtail = 1.78*(1-0.045*AR_w2_e^(0.68)) - 0.64;
    K_vtail = (S_w2_e/S_ref)*(1/(pi*e_vtail*AR_w2_e));
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
    lam_fus = 0.20;
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
    Aero.CD0_fus    = CD0_fus;
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
    lam_nac = 0.20;
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
    Aero.CD0_nac    = CD0_nac;
end 
if Conf.landgear == 1
    
end

% Polar for the rest of elements - una vez obtenidos los resultados de
CD0_par = CD0_w1 + CD0_w2 + CD0_v + CD0_v2 + CD0_vtail + CD0_fus + CD0_nac + CD0_lnd;
CD0_misc = 0;%[-] Coef de resistencia en configuración de crucero - miscelaneos
CD0_lp = (CD0_par + CD0_misc)*0.10; %[-] Coef de resistencia en configuración de crucero miscelanea

CD0 = CD0_par + CD0_misc + CD0_lp;
CD1 = 0;
CD2 = K_w1 + K_w2 + K_vtail;

% % Oswalds coefficient
% e_w1 = 1.78*(1-0.045*AR_w1^(0.68)) - 0.64;
% e_w2 = 1.78*(1-0.045*AR_w2^(0.68)) - 0.64;
% K_w1 = 1/(pi*e_w1*AR_w1);
% K_w2 = 1/(pi*e_w2*AR_w2);

% CD2 =
% Almacenado de valores
Aero.CD0_w1    = CD0_w1;
Aero.CD0_w2    = CD0_w2;
Aero.CD0_v    = CD0_v;
Aero.CD0_v2  = CD0_v2;
Aero.CD0_vtail  = CD0_vtail;
Aero.CD0_fus  = CD0_fus;
Aero.CD0_nac  = CD0_nac;


Aero.CD0_par  = CD0_par;
Aero.CD0_misc  = CD0_misc;
Aero.CD0_lp  = CD0_lp;
Aero.CD0  = CD0;
Aero.CD1  = CD1;
Aero.CD2  = CD2;