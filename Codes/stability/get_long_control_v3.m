function [Stab_Der,Trim_ITER] = get_long_control_v3(AC_CONFIGURATION,Geo_tier,Stab_Der_parts,Aero_TH,conditions,afe,Performance,Stab_Der,Trim_ITER,OUTPUT_read_XLSX)
W1  = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

prop_wash_effect = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;

twin_VTP = AC_CONFIGURATION.twin_VTP;
d_ail    = AC_CONFIGURATION.d_ail;
d_ele    = AC_CONFIGURATION.d_ele;
d_elevon = AC_CONFIGURATION.d_elevon;
d_flap   = AC_CONFIGURATION.d_flap;
d_rudder = AC_CONFIGURATION.d_rudder;
d_rudvtr = AC_CONFIGURATION.d_rudvtr;
d_rudvtr2 = AC_CONFIGURATION.d_rudvtr2;
d_can    = AC_CONFIGURATION.d_can;
AC_type  = AC_CONFIGURATION.AC_type;

% Available control surfaces
d_ail    = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_ail; %Definition of available control surface - aileron
d_ele    = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_ele; %Definition of available control surface - elevator
d_elevon = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_elevon; %Definition of available control surface - elevon
d_flap   = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_flap; %Definition of available control surface - flap
d_rudder = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudder; %Definition of available control surface - rudder
d_rudvtr = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudvtr; %Definition of available control surface - ruddervator
d_rudvtr2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudvtr2; %Definition of available control surface - ruddervator
d_can    = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_can; %Definition of available control surface - canard

x_w1_LE      = Geo_tier.x_w1_LE;
cmac_w1      = Geo_tier.cmac_w1;
CL0_ac       = Stab_Der_parts.CL0_ac;
C_D2         = Aero_TH.CD2;
x_XCG        = conditions.x_XCG;
cR_w1        = Geo_tier.cR_w1;
Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
lambda_w1    = Geo_tier.lambda_w1;
AR_w1_e      = Geo_tier.AR_w1_e;

if prop_wash_effect == 1
    if W1 == 1
        eta_w1_afe_S_w1_afe_S_ref       = afe.eta_w1_afe_S_w1_afe_S_ref;
        eta_w1_no_afe_S_w1_no_afe_S_ref = afe.eta_w1_no_afe_S_w1_no_afe_S_ref;
    end
    if HTP == 1 
        eta_HTP_afe_S_HTP_afe_S_ref       = afe.eta_HTP_afe_S_HTP_afe_S_ref;
        eta_HTP_no_afe_S_HTP_no_afe_S_ref = afe.eta_HTP_no_afe_S_HTP_no_afe_S_ref;
    end
    if Vee ==1
        eta_vee_afe_S_vee_afe_S_ref       = afe.eta_vee_afe_S_vee_afe_S_ref;
        eta_vee_no_afe_S_vee_no_afe_S_ref = afe.eta_vee_no_afe_S_vee_no_afe_S_ref;
    end
    if Vee2 ==1
        eta_vee2_afe_S_vee2_afe_S_ref       = afe.eta_vee2_afe_S_vee2_afe_S_ref;
        eta_vee2_no_afe_S_vee2_no_afe_S_ref = afe.eta_vee2_no_afe_S_vee2_no_afe_S_ref;
    end
end

CL_alpha_w1 = Stab_Der_parts.CL_alpha_w1;
V           = conditions.V;
a           = Performance.a;
Mach        = V/a;
% Xac_cre_W_B = get_xac_cre_wing(lambda_w1, AR_w1_e,Lambda_LE_w1,Mach);
% % Xac_cre_W_B = (0.25*cmac_w1)/cR_w1; %eq 3.38 PAMADI % assumes the influence of the body on the location of the wing aerodynamic center is small; 25% de la cuerda hipótesis    
% xac_w_bar = Xac_cre_W_B*cR_w1/cmac_w1;
xac_w_bar = Geo_tier.xbar_w1/cmac_w1;
x_xbar_w1 = xac_w_bar + x_w1_LE/cmac_w1;

if HTP == 1 
    AR_HTP_e   = Geo_tier.AR_HTP_e;
    x_xbar_HTP = Geo_tier.x_xbar_HTP;
    lambda_HTP = Geo_tier.lambda_HTP;
    diedro_HTP = Geo_tier.dihedral_HTP;
end

if Vee ==1
    AR_vee_e   = Geo_tier.AR_vee_e;
    x_xbar_vee = Geo_tier.x_xbar_vee;
    lambda_vee = Geo_tier.lambda_vee;
    diedro_vee = Geo_tier.dihedral_vee;
end

if Vee2 ==1
    AR_vee2_e   = Geo_tier.AR_vee2_e;
    x_xbar_vee2 = Geo_tier.x_xbar_vee2;
    lambda_vee2 = Geo_tier.lambda_vee2;
    diedro_vee2 = Geo_tier.dihedral_vee2;
end

if Can == 1 
    AR_can_e   = Geo_tier.AR_can_e;
    x_xbar_can = Geo_tier.x_xbar_can;
    lambda_can = Geo_tier.lambda_can;
    diedro_can = Geo_tier.dihedral_can;
end

% ASpro CLdeltae, CDdeltae, CMdeltae
%% Estimación con proyecciones en el plano horizontal
if d_rudvtr == 1
    cf_rudvtr = Geo_tier.cf_rudvtr;
    cf_d_rudvtr = cf_rudvtr;
    t_c_rudvtr = Geo_tier.t_c_rudvtr;
    t_c_d_rudvtr = t_c_rudvtr;
    K_y1_rudvtr_vee = Geo_tier.K_y1_rudvtr_vee;
    K_y2_rudvtr_vee = Geo_tier.K_y2_rudvtr_vee;
    K_y1 = K_y1_rudvtr_vee;
    K_y2 = K_y2_rudvtr_vee;
end

if d_rudvtr2 == 1
    cf_rudvtr2 = Geo_tier.cf_rudvtr2;
    cf_d_rudvtr2 = cf_rudvtr2;
    t_c_rudvtr2 = Geo_tier.t_c_rudvtr2;
    t_c_d_rudvtr2 = t_c_rudvtr2;
    K_y1_rudvtr_vee2 = Geo_tier.K_y1_rudvtr_vee2;
    K_y2_rudvtr_vee2 = Geo_tier.K_y2_rudvtr_vee2;
    K_y1 = K_y1_rudvtr_vee2;
    K_y2 = K_y2_rudvtr_vee2;
end

if d_ele == 1
    cf_ele = Geo_tier.cf_ele;
    cf_d_ele = cf_ele;
    t_c_ele = Geo_tier.t_c_ele;
    t_c_d_ele = t_c_ele;
    K_y1_ele_HTP = Geo_tier.K_y1_ele_HTP;
    K_y2_ele_HTP = Geo_tier.K_y2_ele_HTP;
    K_y1 = K_y1_ele_HTP;
    K_y2 = K_y2_ele_HTP;
%     K_y1 = 0.02;
%     K_y2 = 0.9541;
end

if d_elevon == 1
    cf_elevon = Geo_tier.cf_elevon;
    cf_d_elevon = cf_elevon;
    t_c_elevon = Geo_tier.t_c_elevon;
    t_c_d_elevon = t_c_elevon;
    K_y1_elevon_w1 = Geo_tier.K_y1_elevon_w1;
    K_y2_elevon_w1 = Geo_tier.K_y2_elevon_w1;
    K_y1 = K_y1_elevon_w1;
    K_y2 = K_y2_elevon_w1;
%     K_y1 = 0.02;
%     K_y2 = 0.9541;
end

if d_rudder == 1
    cf_rudder = Geo_tier.cf_rudder;
    cf_d_rudder = cf_rudder;
    t_c_rudder = Geo_tier.t_c_rudder;
    t_c_d_rudder = t_c_rudder;
    K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
    K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
    K_y1 = K_y1_rudder_VTP;
    K_y2 = K_y2_rudder_VTP;
end

if d_can == 1
    cf_can = Geo_tier.cf_canard;
    cf_d_can = cf_can;
    t_c_can = Geo_tier.t_c_canard;
    t_c_d_can = t_c_can;
    K_y1_canard_can = Geo_tier.K_y1_canard_can;
    K_y2_canard_can = Geo_tier.K_y2_canard_can;
    K_y1 = K_y1_canard_can;
    K_y2 = K_y2_canard_can;
end

if d_ail == 1
    cf_ail = Geo_tier.cf_ail;
    cf_d_ail = cf_ail;
    t_c_ail = Geo_tier.t_c_ail;
    t_c_d_ail = t_c_ail;
    K_y1_ail_w1 = Geo_tier.K_y1_ail_w1;
    K_y2_ail_w1 = Geo_tier.K_y2_ail_w1;
    K_y1 = K_y1_ail_w1;
    K_y2 = K_y2_ail_w1;
end

if d_flap == 1
    cf_flap = Geo_tier.cf_flap;
    cf_d_flap = cf_flap;
    t_c_flap = Geo_tier.t_c_flap;
    t_c_d_flap = t_c_flap;
    K_y1_flap_w1 = Geo_tier.K_y1_flap_w1;
    K_y2_flap_w1 = Geo_tier.K_y2_flap_w1;
    K_y1 = K_y1_flap_w1;
    K_y2 = K_y2_flap_w1;
end

% % theoretical V-Tail sectional lift curve slope at Mach equal to zero
% k_prima              = 1;
% Cla_M0               = 2*pi + 5.0525*t_c_d;
% Kb_HTP                = Kb_calc(K_y1, K_y2, lambda_HTP);              % Fig 3.36 Pamadi
% Clde_Cldetheo_HTP     = Cldelta_Cldeltatheory_calc(cf_d, Cla_M0);% Fig 3.37b Pamadi
% Cldetheo_HTP          = Cldeltatheory_calc(cf_d, t_c_d);               % Fig 3.37a Pamadi
% alphaCL_alphaCl_HTP   = alphaCL_alphaCl_calc(cf_d,AR_HTP_e);              % Fig 3.35 Pamadi
% alpha_delta_e       = Kb_HTP*Clde_Cldetheo_HTP*Cldetheo_HTP/Cla_M0*alphaCL_alphaCl_HTP;

% Revision to determine CLalpha of the horizontal control surface
if d_rudvtr == 1
    %     CL_alpha_wb_long = CL_alpha_wb_vee;
    %Debería ser:
    CL_alpha_wb_vee = Stab_Der_parts.CL_alpha_vee;
%     CL_alpha_wb_long = Stab_Der_parts.CL_alpha_vee;
end

% Revision to determine CLalpha of the horizontal control surface
if d_rudvtr2 == 1
    %     CL_alpha_wb_long = CL_alpha_wb_vee;
    %Debería ser:
    CL_alpha_wb_vee2 = Stab_Der_parts.CL_alpha_vee2;
%     CL_alpha_wb_long = Stab_Der_parts.CL_alpha_vee;
end

if d_ele == 1
    %     CL_alpha_wb_long = CL_alpha_wb_HTP;
    %Debería ser:
    CL_alpha_wb_HTP = Stab_Der_parts.CL_alpha_HTP;
%     CL_alpha_wb_long = Stab_Der_parts.CL_alpha_HTP;
end

if d_elevon == 1
    %     CL_alpha_wb_long = CL_alpha_wb_HTP;
    %Debería ser:
    CL_alpha_wb_elevon = Stab_Der_parts.CL_alpha_w1;
%     CL_alpha_wb_long = Stab_Der_parts.CL_alpha_HTP;
end

if d_can == 1
    %     CL_alpha_wb_long = CL_alpha_wb_HTP;
    %Debería ser:
    CL_alpha_wb_can = Stab_Der_parts.CL_alpha_can;    
%     CL_alpha_wb_long = Stab_Der_parts.CL_alpha_can;
end
%
% CL_delta_e          = alpha_delta_e*CL_alpha_wb_long*(S_HTP_e/S_ref);
% CD_delta_e          = 2*CL0_ac*C_D2*(eta_HTP_afe_S_HTP_afe_S_ref*CL_delta_e + eta_HTP_no_afe_S_HTP_no_afe_S_ref*CL_delta_e);
% CM_delta_e            = - CL_delta_e*((x_xbar_HTP - x_XCG)/cmac_w1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPLET FOR DIFFERENT CONTROL SURFACES!
%% CL_delta_rv
if d_rudvtr == 1
    % theoretical V-Tail sectional lift curve slope at Mach equal to zero
    k_prima              = 1;
    %% Perfil revisado para el ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_vee11);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c_d_rudvtr); %Pamadi 3.17
    clalpha_vee_M0 = 1.05*cla_clatheo*clatheo; %Pamadi 3.17
    clalpha_vee = clalpha_vee_M0/sqrt(1-Mach^2); %Pamadi 3.17
    Cla_theo_vee = 2*pi + 5.0525*t_c_d_rudvtr;
    Cla_M0 = clalpha_vee_M0/Cla_theo_vee;
    Kb_vee                = Kb_calc(K_y1_rudvtr_vee, K_y2_rudvtr_vee, lambda_vee);              % Fig 3.36 Pamadi
    Clde_Cldetheo_vee     = Cldelta_Cldeltatheory_calc(cf_d_rudvtr, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_vee          = Cldeltatheory_calc(cf_d_rudvtr, t_c_d_rudvtr);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_vee   = alphaCL_alphaCl_calc(cf_d_rudvtr,AR_vee_e);              % Fig 3.35 Pamadi
    alpha_delta_e       = Kb_vee*Clde_Cldetheo_vee*Cldetheo_vee/clalpha_vee*alphaCL_alphaCl_vee;
    %% CL_delta_rv
    Balance_rv = cf_rudvtr;
    %% CL_delta_rv
    Balance_rv = cf_rudvtr;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_rv = 0.83 + 0.30714*Balance_rv; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    %     CL_i_vee = (eta_vee_afe_S_vee_afe_S_ref*CL_alpha_wb_vee + eta_vee_no_afe_S_vee_no_afe_S_ref*CL_alpha_wb_vee);
    CL_i_vee = CL_alpha_wb_vee; %ya multiplicado por eta*S/S_ref en propwash_influence;
    CL_delta_rv_0_darc = f_val_rv*CL_i_vee*alpha_delta_e;
    CL_delta_rv_0_naca = f_val_rv*(CL_i_vee/cos(diedro_vee))*alpha_delta_e; %Preferencia por NACA a expensas de comparar.
    CL_delta_rv_0 = CL_delta_rv_0_naca;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    
    CL_delta_rv = CL_delta_rv_0*K_prima2;
    CL_delta_e = CL_delta_rv;
    
    
    %% CD_delta_rv %
    CD_i_vee          = 2*CL0_ac*C_D2*(CL_alpha_wb_vee); %C_D2 = 1/(e*AR*pi) %ya multiplicado por eta*S/S_ref en propwash_influence;
    CD_delta_rv       = alpha_delta_e*CD_i_vee;
    CD_delta_e       = CD_delta_rv;
    
    %% CMdelta_rv
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_vee = -CL_i_vee*((x_xbar_vee - x_XCG)/cmac_w1);
    CM_delta_rv_0 = f_val_rv*CM_i_vee*alpha_delta_e;
    CM_delta_rv = CM_delta_rv_0*K_prima2;
    CM_delta_e = CM_delta_rv;
    
    % Store real component of control effort
    CL_delta_rv = CL_delta_e;
    CD_delta_rv = CD_delta_e;
    CM_delta_rv = CM_delta_e;
else
    % assigns 0 if no present
    CL_delta_rv = 0;
    CD_delta_rv = 0;
    CM_delta_rv = 0;
end

if d_rudvtr2 == 1
    % theoretical V-Tail sectional lift curve slope at Mach equal to zero
    k_prima              = 1;
    %% Perfil revisado para el ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_vee21);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c_d_rudvtr); %Pamadi 3.17
    clalpha_vee2_M0 = 1.05*cla_clatheo*clatheo; %Pamadi 3.17
    clalpha_vee2 = clalpha_vee2_M0/sqrt(1-Mach^2); %Pamadi 3.17
    Cla_theo_vee2 = 2*pi + 5.0525*t_c_d_rudvtr;
    Cla_M0 = clalpha_vee2_M0/Cla_theo_vee2;
    Kb_vee2                = Kb_calc(K_y1_rudvtr_vee2, K_y2_rudvtr_vee2, lambda_vee2);              % Fig 3.36 Pamadi
    Clde_Cldetheo_vee2     = Cldelta_Cldeltatheory_calc(cf_d_rudvtr, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_vee2          = Cldeltatheory_calc(cf_d_rudvtr, t_c_d_rudvtr);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_vee2   = alphaCL_alphaCl_calc(cf_d_rudvtr,AR_vee2_e);              % Fig 3.35 Pamadi
    alpha_delta_e       = Kb_vee2*Clde_Cldetheo_vee2*Cldetheo_vee2/clalpha_vee2*alphaCL_alphaCl_vee2;
    %% CL_delta_rv
    Balance_rv = cf_rudvtr;
    %% CL_delta_rv
    Balance_rv = cf_rudvtr;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_rv = 0.83 + 0.30714*Balance_rv; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    %     CL_i_vee2 = (eta_vee2_afe_S_vee2_afe_S_ref*CL_alpha_wb_vee2 + eta_vee2_no_afe_S_vee2_no_afe_S_ref*CL_alpha_wb_vee2);
    CL_i_vee2 = CL_alpha_wb_vee2; %ya multiplicado por eta*S/S_ref en propwash_influence;
    CL_delta_rv_0_darc = f_val_rv*CL_i_vee2*alpha_delta_e;
    CL_delta_rv_0_naca = f_val_rv*(CL_i_vee2/cos(diedro_vee2))*alpha_delta_e; %Preferencia por NACA a expensas de comparar.
    CL_delta_rv_0 = CL_delta_rv_0_naca;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    
    CL_delta_rv = CL_delta_rv_0*K_prima2;
    CL_delta_e_tmp = CL_delta_rv;
    
    
    %% CD_delta_rv %
    CD_i_vee2          = 2*CL0_ac*C_D2*(CL_alpha_wb_vee2); %C_D2 = 1/(e*AR*pi) %ya multiplicado por eta*S/S_ref en propwash_influence;
    CD_delta_rv       = alpha_delta_e*CD_i_vee2;
    CD_delta_e_tmp       = CD_delta_rv;
    
    %% CMdelta_rv
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_vee2 = -CL_i_vee2*((x_xbar_vee2 - x_XCG)/cmac_w1);
    CM_delta_rv_0 = f_val_rv*CM_i_vee2*alpha_delta_e;
    CM_delta_rv = CM_delta_rv_0*K_prima2;
    CM_delta_e_tmp = CM_delta_rv;
    
    % Store real component of control effort
    CL_delta_rv2 = CL_delta_e_tmp;
    CD_delta_rv2 = CD_delta_e_tmp;
    CM_delta_rv2 = CM_delta_e_tmp;
else
    % assigns 0 if no present
    CL_delta_rv2 = 0;
    CD_delta_rv2 = 0;
    CM_delta_rv2 = 0;
end

if d_ele == 1
    % theoretical V-Tail sectional lift curve slope at Mach equal to zero
    k_prima = 1;
    %% Perfil revisado para el ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_HTP1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c_ele);
    clalpha_HTP_M0 = 1.05*cla_clatheo*clatheo;
    clalpha_HTP = clalpha_HTP_M0/sqrt(1-Mach^2);
    Cla_theo_HTP = 2*pi + 5.0525*t_c_d_ele;
    Cla_M0 = clalpha_HTP_M0/Cla_theo_HTP;
    Kb_HTP                = Kb_calc(K_y1_ele_HTP, K_y2_ele_HTP, lambda_HTP);      % Fig 3.36 Pamadi
    Clde_Cldetheo_HTP     = Cldelta_Cldeltatheory_calc(cf_d_ele, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_HTP          = Cldeltatheory_calc(cf_d_ele, t_c_d_ele);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_HTP   = alphaCL_alphaCl_calc(cf_d_ele,AR_HTP_e);              % Fig 3.35 Pamadi
    alpha_delta_e       = Kb_HTP*Clde_Cldetheo_HTP*Cldetheo_HTP/clalpha_HTP*alphaCL_alphaCl_HTP;

    Balance_e = cf_ele;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_e = 0.83 + 0.30714*Balance_e; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    CL_i_HTP = CL_alpha_wb_HTP; %ya multiplicado por eta*S/S_ref en propwash_influence;
    
    CL_delta_e_0 = f_val_e*CL_i_HTP*alpha_delta_e;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    CL_delta_e = CL_delta_e_0*K_prima2;
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    %% CD_delta
    CD_i_h          = 2*CL0_ac*C_D2*(CL_alpha_wb_HTP); %ya multiplicado por eta*S/S_ref en propwash_influence;
    CD_delta_e       = alpha_delta_e*CD_i_h;
    
    %% CMdelta_r
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_HTP = -CL_i_HTP*((x_xbar_HTP - x_XCG)/cmac_w1);
    CM_delta_e_0 = f_val_e*CM_i_HTP*alpha_delta_e;
    CM_delta_e = CM_delta_e_0*K_prima2;
    
else
    % assigns 0 if no present
    CL_delta_e = 0;
    CD_delta_e = 0;
    CM_delta_e = 0;
end

if d_elevon == 1
    % theoretical sectional lift curve slope at Mach equal to zero
    k_prima = 1;
    %% Perfil revisado para el ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_w1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c_elevon);
    clalpha_w1_M0 = 1.05*cla_clatheo*clatheo;
    clalpha_w1 = clalpha_w1_M0/sqrt(1-Mach^2);
    Cla_theo_w1 = 2*pi + 5.0525*t_c_d_elevon;
    Cla_M0 = clalpha_w1_M0/Cla_theo_w1;
    Kb_w1                = Kb_calc(K_y1_elevon_w1, K_y2_elevon_w1, lambda_w1);      % Fig 3.36 Pamadi
    Clde_Cldetheo_w1     = Cldelta_Cldeltatheory_calc(cf_d_elevon, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_w1          = Cldeltatheory_calc(cf_d_elevon, t_c_d_elevon);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_w1   = alphaCL_alphaCl_calc(cf_d_elevon,AR_w1_e);              % Fig 3.35 Pamadi
    alpha_delta_e       = Kb_w1*Clde_Cldetheo_w1*Cldetheo_w1/clalpha_w1*alphaCL_alphaCl_w1;
    
    Balance_e = cf_elevon;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_e = 0.83 + 0.30714*Balance_e; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    CL_i_w1 = CL_alpha_wb_elevon; %ya multiplicado por eta*S/S_ref en propwash_influence;
    
    CL_delta_e_0 = f_val_e*CL_i_w1*alpha_delta_e;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    CL_delta_e_tmp = CL_delta_e_0*K_prima2;
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    %% CD_delta
    CD_i_h          = 2*CL0_ac*C_D2*(CL_alpha_wb_elevon); %ya multiplicado por eta*S/S_ref en propwash_influence;
    CD_delta_e_tmp       = alpha_delta_e*CD_i_h;
    
    %% CMdelta_r
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_w1 = -CL_i_w1*((x_xbar_w1 - x_XCG)/cmac_w1);
    CM_delta_e_0 = f_val_e*CM_i_w1*alpha_delta_e;
    CM_delta_e_tmp = CM_delta_e_0*K_prima2;    

    % Save as elevons
    CL_delta_elevon = CL_delta_e_tmp;
    CD_delta_elevon = CD_delta_e_tmp;
    CM_delta_elevon = CM_delta_e_tmp;
else
    % assigns 0 if no present
    CL_delta_elevon = 0;
    CD_delta_elevon = 0;
    CM_delta_elevon = 0;
end

if d_can == 1
    % theoretical Canard sectional lift curve slope at Mach equal to zero
    k_prima              = 1;
    %% Perfil revisado para el ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_c1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c_can);
    clalpha_can_M0 = 1.05*cla_clatheo*clatheo;
    clalpha_can = clalpha_can_M0/sqrt(1-Mach^2);
    Cla_theo_can = 2*pi + 5.0525*t_c_d_can;
    Cla_M0 = clalpha_can_M0/Cla_theo_can;
    Kb_can                = Kb_calc(K_y1_canard_can, K_y2_canard_can, lambda_can);      % Fig 3.36 Pamadi
    Clde_Cldetheo_can     = Cldelta_Cldeltatheory_calc(cf_d_can, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_can          = Cldeltatheory_calc(cf_d_can, t_c_d_can);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_can   = alphaCL_alphaCl_calc(cf_d_can,AR_can_e);              % Fig 3.35 Pamadi
    alpha_delta_can       = Kb_can*Clde_Cldetheo_can*Cldetheo_can/clalpha_can*alphaCL_alphaCl_can;
    
    Balance_can = cf_can;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_canard = 0.83 + 0.30714*Balance_can; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    CL_i_can = CL_alpha_wb_can; %ya multiplicado por eta*S/S_ref en propwash_influence;
    
    CL_delta_can_0 = f_val_canard*CL_i_can*alpha_delta_can;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    CL_delta_e2 = CL_delta_can_0*K_prima2;
    CL_delta_can = CL_delta_e2;
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    %% CD_delta
    CD_i_h          = 2*CL0_ac*C_D2*(CL_alpha_wb_can); %ya multiplicado por eta*S/S_ref en propwash_influence;
    CD_delta_can       = alpha_delta_e*CD_i_h;
    
    %% CMdelta_r
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_can = -CL_i_can*((x_xbar_can - x_XCG)/cmac_w1);
    CM_delta_can_0 = f_val_canard*CM_i_can*alpha_delta_can;
    CM_delta_e = CM_delta_can_0*K_prima2;
    CM_delta_can = CM_delta_e;    
else
    % assigns 0 if no present
    CL_delta_can = 0;
    CD_delta_can = 0;
    CM_delta_can = 0;
end

if d_ail ==1
    % theoretical wing sectional lift curve slope at Mach equal to zero
    k_prima              = 1;
    %% Perfil revisado para el ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_w1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c_d_ail);
    clalpha_w1_M0 = 1.05*cla_clatheo*clatheo;
    clalpha_w1 = clalpha_w1_M0/sqrt(1-Mach^2);
    Cla_theo_w1 = 2*pi + 5.0525*t_c_d_ail;
    Cla_M0 = clalpha_w1_M0/Cla_theo_w1;
    Kb_w1                = Kb_calc(K_y1_ail_w1, K_y2_ail_w1, lambda_w1);              % Fig 3.36 Pamadi
    Clde_Cldetheo_w1     = Cldelta_Cldeltatheory_calc(cf_d_ail, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_w1          = Cldeltatheory_calc(cf_d_ail, t_c_d_ail);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_w1   = alphaCL_alphaCl_calc(cf_d_ail,AR_w1_e);              % Fig 3.35 Pamadi
    alpha_delta_ail       = Kb_w1*Clde_Cldetheo_w1*Cldetheo_w1/clalpha_w1*alphaCL_alphaCl_w1;
    
    Balance_ail = cf_ail;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_ail = 0.83 + 0.30714*Balance_ail; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    CL_i_w1 = CL_alpha_w1;
    
    CL_delta_ail_0 = f_val_ail*CL_i_w1*alpha_delta_ail;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    CL_delta_ail = CL_delta_ail_0*K_prima2;
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    
    %% CD_delta
    CD_i_h          = 2*CL0_ac*C_D2*(CL_alpha_w1); %ya multiplicado por eta*S/S_ref en propwash_influence;
    CD_delta_ail       = alpha_delta_ail*CD_i_h;
    
    %% CMdelta_rv %% por qué rv?!?!?!
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_w1 = -CL_i_w1*((x_xbar_w1 - x_XCG)/cmac_w1); %OJO CON ESTO, NO DEBERÍA SER EN EL ALA, Y NO EN LA COLA? ROLLO_ x_xbar_w1!!!
    CM_delta_ail_0 = f_val_ail*CM_i_w1*alpha_delta_ail;
    CM_delta_ail = CM_delta_ail_0*K_prima2;
else
    % assigns 0 if no present
    CL_delta_ail = 0;
    CD_delta_ail = 0;
    CM_delta_ail = 0;
end

if d_flap ==1
    % theoretical wing sectional lift curve slope at Mach equal to zero
    k_prima              = 1;
    %% Perfil revisado para el ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_w1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c_d_flap);
    clalpha_w1_M0 = 1.05*cla_clatheo*clatheo;
    clalpha_w1 = clalpha_w1_M0/sqrt(1-Mach^2);
    Cla_theo_w1 = 2*pi + 5.0525*t_c_d_flap;
    Cla_M0 = clalpha_w1_M0/Cla_theo_w1;
    Kb_w1                = Kb_calc(K_y1_flap_w1, K_y2_flap_w1, lambda_w1);              % Fig 3.36 Pamadi
    Clde_Cldetheo_w1     = Cldelta_Cldeltatheory_calc(cf_d_flap, Cla_M0);% Fig 3.37b Pamadi
    Cldetheo_w1          = Cldeltatheory_calc(cf_d_flap, t_c_d_flap);               % Fig 3.37a Pamadi
    alphaCL_alphaCl_w1   = alphaCL_alphaCl_calc(cf_d_flap,AR_w1_e);              % Fig 3.35 Pamadi
    alpha_delta_flap       = Kb_w1*Clde_Cldetheo_w1*Cldetheo_w1/clalpha_w1*alphaCL_alphaCl_w1;
    
    Balance_flap = cf_flap;
    % f_val_rv = 0.83 + 0.19857*Balance_rv; % Elliptic nose
    f_val_flap = 0.83 + 0.30714*Balance_flap; % Round nose
    % The airplane lift-coefficient-due-to-V-Tail-incidence derivative
    CL_i_w1 = CL_alpha_w1;
    
    CL_delta_flap_0 = f_val_flap*CL_i_w1*alpha_delta_flap;
    % The airplane lift-coefficient-due-to-ruddervator-deflection
    K_prima2 = 1; % APPROXIMATED > the slope at the the given ruddervator deflection.:
    CL_delta_flap = CL_delta_flap_0*K_prima2;
    
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI and is a function of
    % ruddervator deflection angle and the average ruddervator chord to V-Tail chord ratio aft of hinge line.
    % K_prima = K_prima_calc(delta_f,cf_rudvtr)
    
    %% CD_delta
    CD_i_h          = 2*CL0_ac*C_D2*(CL_alpha_w1);
    CD_delta_flap       = alpha_delta_flap*CD_i_h;
    
    %% CMdelta_rv
    % correction for nonlinear pitching moment behavior of plain flaps is found from Figure 8.13 in Airplane Design Part VI
    % and is a function of elevator deflection angle and the average elevator chord to horizontal tail chord ratio aft of hinge line.
    CM_i_w1 = -CL_i_w1*((x_xbar_w1 - x_XCG)/cmac_w1);
    CM_delta_flap_0 = f_val_ail*CM_i_w1*alpha_delta_flap;
    CM_delta_flap = CM_delta_flap_0*K_prima2;
else
    % assigns 0 if no present
    CL_delta_flap = 0;
    CD_delta_flap = 0;
    CM_delta_flap = 0;
end

Trim_ITER.CL_delta_e = CL_delta_e;
Trim_ITER.CD_delta_e = CD_delta_e;
Trim_ITER.CM_delta_e = CM_delta_e;

Trim_ITER.CL_delta_elevon = CL_delta_elevon;
Trim_ITER.CD_delta_elevon = CD_delta_elevon;
Trim_ITER.CM_delta_elevon = CM_delta_elevon;

Trim_ITER.CL_delta_rv = CL_delta_rv;
Trim_ITER.CD_delta_rv = CD_delta_rv;
Trim_ITER.CM_delta_rv = CM_delta_rv;

Trim_ITER.CL_delta_rv2 = CL_delta_rv2;
Trim_ITER.CD_delta_rv2 = CD_delta_rv2;
Trim_ITER.CM_delta_rv2 = CM_delta_rv2;

Trim_ITER.CL_delta_can = CL_delta_can;
Trim_ITER.CD_delta_can = CD_delta_can;
Trim_ITER.CM_delta_can = CM_delta_can;

%  Elevator control effectiveness in longitudinal
Stab_Der.CL_delta_e = CL_delta_e + CL_delta_elevon;
Stab_Der.CD_delta_e = CD_delta_e + CD_delta_elevon;
Stab_Der.CM_delta_e = CM_delta_e + CM_delta_elevon;

%  Elevator control effectiveness in longitudinal
Stab_Der.CL_delta_elevon = CL_delta_elevon;
Stab_Der.CD_delta_elevon = CD_delta_elevon;
Stab_Der.CM_delta_elevon = CM_delta_elevon;

%  ruddervator - Vtail control effectiveness
Stab_Der.CL_delta_rv = CL_delta_rv;
Stab_Der.CD_delta_rv = CD_delta_rv;
Stab_Der.CM_delta_rv = CM_delta_rv;

Stab_Der.CL_delta_rv2 = CL_delta_rv2;
Stab_Der.CD_delta_rv2 = CD_delta_rv2;
Stab_Der.CM_delta_rv2 = CM_delta_rv2;

%  Canard control effectiveness
Stab_Der.CL_delta_can = CL_delta_can;
Stab_Der.CD_delta_can = CD_delta_can;
Stab_Der.CM_delta_can = CM_delta_can;
%  Aileron control effectiveness
Stab_Der.CL_delta_ail = CL_delta_ail;
Stab_Der.CD_delta_ail = CD_delta_ail;
Stab_Der.CM_delta_ail = CM_delta_ail;
%  Flap control effectiveness
Stab_Der.CL_delta_flap = CL_delta_flap;
Stab_Der.CD_delta_flap = CD_delta_flap;
Stab_Der.CM_delta_flap = CM_delta_flap;

% Control derivatives
CX_delta_e = -CD_delta_e;
CZ_delta_e = -CL_delta_e;

Stab_Der.CX_delta_e = CX_delta_e;
Stab_Der.CZ_delta_e = CZ_delta_e;
end