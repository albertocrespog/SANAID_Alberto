function [Stab_Der_parts Stab_Der] = getCyp(AC_CONFIGURATION,modelo,alpha,Stab_Der,Stab_Der_parts,Trim_ITER, Geo_tier, OUTPUT_read_XLSX)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

% CL          = modelo.ala.CLa; %% CAMBIAR
CL_w1       = Trim_ITER.CL_w1;
CLa_we      = Stab_Der_parts.CLalpha_w1_e_pw;
AR_w        = modelo.ala.AR;
AR_we       = modelo.ala.ARwe;
eOswald_w   = modelo.ala.oswald;
LAMc4_w     = modelo.ala.LAMc4;
TR_w        = modelo.ala.TR;
% Cla_w       = modelo.ala.Cla;
diedro_w    = modelo.ala.diedro;
Zca_w       = -Geo_tier.z_w1_LE; % positive if the wing is below thefuselage centerline.

b_w         = modelo.ala.b;
Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.ala.Xca, 'pchip');


% Vee & VTP
if VTP == 1 
    Zca_v       = modelo.vertical.Zca;
    Xca_v       = modelo.vertical.Xca;
    S_v         = modelo.vertical.S;
    b_v         = modelo.vertical.b;
    eta_v       = modelo.vertical.eta;
end

if Vee == 1
    Zca_v       = modelo.vertical.Zca;
    Xca_v       = modelo.vertical.Xca;
    S_v         = modelo.vertical.S;
    b_v         = modelo.vertical.b;
    eta_v       = modelo.vertical.eta;
end

% Vee 2
if Vee2 == 1
    Zca_v2       = modelo.vertical2.Zca;
    Xca_v2       = modelo.vertical2.Xca;
    S_v2         = modelo.vertical2.S;
    b_v2         = modelo.vertical2.b;
    eta_v2       = modelo.vertical2.eta;
end

Sref        = modelo.general.Sref;
Minf        = modelo.general.Minf;
Xcg         = modelo.general.Xcg;
Zcg         = modelo.general.Zcg;
Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca, 'pchip');


beta = sqrt(1 - Minf^2);
B_prandtl = sqrt(1 - (Minf^2)*(cos(LAMc4_w))^2);

% Vee & VTP
if VTP == 1 
    Cy_beta_VTP = Stab_Der_parts.Cy_beta_VTP;
    
    % Twin Vertical Tail configuration
    if AC_CONFIGURATION.twin_VTP == 1
        Cy_beta_VTP = 2*Cy_beta_VTP;
    end
    
    if isnan(Cy_beta_VTP)
        Cy_beta_VTP = 0;
    end
else
        Cy_beta_VTP = 0;
end

% Vee
if Vee == 1
    Cy_beta_vee = Stab_Der_parts.Cy_beta_vee;
    
    if isnan(Cy_beta_vee)
        Cy_beta_vee = 0;
    end
else
        Cy_beta_vee = 0;
end

% Vee 2
if Vee2 == 1
    Cy_beta_vee2 = Stab_Der_parts.Cy_beta_vee2;
    if isnan(Cy_beta_vee2)
        Cy_beta_vee2 = 0;
    end
else
    Cy_beta_vee2 = 0;
end

% Can
if Can == 1
    Cy_beta_can = Stab_Der_parts.Cy_beta_can;
    if isnan(Cy_beta_can)
        Cy_beta_can = 0;
    end
else
    Cy_beta_can = 0;
end

% HTP
if HTP == 1
    Cy_beta_HTP = Stab_Der_parts.Cy_beta_HTP;
    if isnan(Cy_beta_HTP)
        Cy_beta_HTP = 0;
    end
else
    Cy_beta_HTP = 0;
end

%% CY_p
%% Wing contributiuon
if W1 == 1
    % Ala
    CL_w = Trim_ITER.CL_w1;
    CLa_we      = Stab_Der_parts.CLalpha_w1_e_pw;
    AR_w    = Geo_tier.AR_w1;
    AR_we       =  Geo_tier.AR_w1_e;
    b_w     = Geo_tier.b_w1;
    LAMLE_w = Geo_tier.Lambda_LE_w1;
    LAMc2_w = Geo_tier.Lambda_c2_w1;
    LAMc4_w = Geo_tier.Lambda_c4_w1;
    TR_w    = Geo_tier.lambda_w1;
    diedro_w  = Geo_tier.dihedral_w1;
    Zca_w       = -Geo_tier.z_w1_LE; % positive if the wing is below thefuselage centerline.

    taper_e = 1;
    a1 = 0.0004;
    a2 = -0.0080;
    a3 = 0.0501;
    a4 =0.8642;
    lambda1 = AR_we*taper_e/(cos(LAMLE_w));
    R = a1*lambda1^3 + a2*lambda1^2 + a3*lambda1 + a4;
    eOswald_w = 1.1*CLa_we/(R*CLa_we + (1-R)*pi*AR_we);

    if eOswald_w > 1
        eOswald_w = 0.95;
    end

% Se han mezclado los métodos del DATCOM y los de Pamadi (pag 406)
    aw1         = CLa_we/(pi*AR_we*eOswald_w);
    K           = (1-aw1)/(1-eOswald_w*aw1);    % Pamadi Ecuacion 4.549
    %K           = ((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - (k1*CLa + 2*k2*CL*CLa))/((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - 2*CL*CLa/pi/ARexp_w);
    Cyp_CL_M0   = Cyp_CL_M0_calc(LAMc4_w*180/pi, TR_w);
    B_prandtl = sqrt(1 - (Minf^2)*(cos(LAMc4_w))^2);
    Cyp_CL_CL0  = (AR_we + 4*cos(LAMc4_w))/(AR_we*B_prandtl + 4*cos(LAMc4_w))*(AR_we*B_prandtl + cos(LAMc4_w))/(AR_we + cos(LAMc4_w))*Cyp_CL_M0;
    t_c = Geo_tier.t_c_flap;

    %% Perfil revisado para el Ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_w1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
    clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
    k_Clp       = clalpha_w2_M0/2/pi;
    betaClp_k   = C_lp_calc(TR_w, AR_w, LAMc4_w, Minf);
    Clp_die0    = betaClp_k*k_Clp/beta;
    DCy_p_die   = (3*sin(diedro_w)*(1 - 4*(Zca_w-Zcg)/b_w*sin(diedro_w)))*Clp_die0;
    Cy_p_w1      = K*Cyp_CL_CL0*CL_w + DCy_p_die; %DATCOM 7.1.2.1
    % z is the vertical distance between the c.g. and the wing root
    % quarter-chord point, positive for the c.g. above the wing root chord.
    % DATCOM 7.1.2.1 (pagina 2523)
    %DATCOM  7.4.2.1 pag.332/444

    if isnan(Cy_p_w1)
        Cy_p_w1 = 0;
    end
else
    Cy_p_w1 = 0;
end



%% Canard contributiuon
if Can == 1
    % Ala
    CL_w = Trim_ITER.CL_can;
    CLa_we      = Stab_Der_parts.CLalpha_can_e_pw;
    AR_w    = Geo_tier.AR_can;
    AR_we       =  Geo_tier.AR_can_e;
    b_w     = Geo_tier.b_can;
    LAMLE_w = Geo_tier.Lambda_LE_can;
    LAMc2_w = Geo_tier.Lambda_c2_can;
    LAMc4_w = Geo_tier.Lambda_c4_can;
    TR_w    = Geo_tier.lambda_can;
    diedro_w  = Geo_tier.dihedral_can;
    Zca_w       = -Geo_tier.z_can_LE; % positive if the wing is below thefuselage centerline.

    taper_e = 1;
    a1 = 0.0004;
    a2 = -0.0080;
    a3 = 0.0501;
    a4 =0.8642;
    lambda1 = AR_we*taper_e/(cos(LAMLE_w));
    R = a1*lambda1^3 + a2*lambda1^2 + a3*lambda1 + a4;
    eOswald_w = 1.1*CLa_we/(R*CLa_we + (1-R)*pi*AR_we);

    if eOswald_w > 1
        eOswald_w = 0.95;
    end

% Se han mezclado los métodos del DATCOM y los de Pamadi (pag 406)
    aw1         = CLa_we/(pi*AR_we*eOswald_w);
    K           = (1-aw1)/(1-eOswald_w*aw1);    % Pamadi Ecuacion 4.549
    %K           = ((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - (k1*CLa + 2*k2*CL*CLa))/((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - 2*CL*CLa/pi/ARexp_w);
    Cyp_CL_M0   = Cyp_CL_M0_calc(LAMc4_w*180/pi, TR_w);
    B_prandtl = sqrt(1 - (Minf^2)*(cos(LAMc4_w))^2);
    Cyp_CL_CL0  = (AR_we + 4*cos(LAMc4_w))/(AR_we*B_prandtl + 4*cos(LAMc4_w))*(AR_we*B_prandtl + cos(LAMc4_w))/(AR_we + cos(LAMc4_w))*Cyp_CL_M0;
    t_c = Geo_tier.t_c_flap;

    %% Perfil revisado para el Ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_can1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
    clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
    k_Clp       = clalpha_w2_M0/2/pi;
    betaClp_k   = C_lp_calc(TR_w, AR_w, LAMc4_w, Minf);
    Clp_die0    = betaClp_k*k_Clp/beta;
    DCy_p_die   = (3*sin(diedro_w)*(1 - 4*(Zca_w-Zcg)/b_w*sin(diedro_w)))*Clp_die0;
    Cy_p_can      = K*Cyp_CL_CL0*CL_w + DCy_p_die; %DATCOM 7.1.2.1
    % z is the vertical distance between the c.g. and the wing root
    % quarter-chord point, positive for the c.g. above the wing root chord.
    % DATCOM 7.1.2.1 (pagina 2523)
    %DATCOM  7.4.2.1 pag.332/444

    if isnan(Cy_p_w1)
        Cy_p_can = 0;
    end
else
    Cy_p_can = 0;
end


%% Canard contributiuon
if HTP == 1
    % Ala
    CL_w = Trim_ITER.CL_HTP;
    CLa_we      = Stab_Der_parts.CLalpha_HTP_e_pw;
    AR_w    = Geo_tier.AR_HTP;
    AR_we       =  Geo_tier.AR_HTP_e;
    b_w     = Geo_tier.b_HTP;
    LAMLE_w = Geo_tier.Lambda_LE_HTP;
    LAMc2_w = Geo_tier.Lambda_c2_HTP;
    LAMc4_w = Geo_tier.Lambda_c4_HTP;
    TR_w    = Geo_tier.lambda_HTP;
    diedro_w  = Geo_tier.dihedral_HTP;
    Zca_w       = -Geo_tier.z_HTP_LE; % positive if the wing is below thefuselage centerline.

    taper_e = 1;
    a1 = 0.0004;
    a2 = -0.0080;
    a3 = 0.0501;
    a4 =0.8642;
    lambda1 = AR_we*taper_e/(cos(LAMLE_w));
    R = a1*lambda1^3 + a2*lambda1^2 + a3*lambda1 + a4;
    eOswald_w = 1.1*CLa_we/(R*CLa_we + (1-R)*pi*AR_we);

    if eOswald_w > 1
        eOswald_w = 0.95;
    end

% Se han mezclado los métodos del DATCOM y los de Pamadi (pag 406)
    aw1         = CLa_we/(pi*AR_we*eOswald_w);
    K           = (1-aw1)/(1-eOswald_w*aw1);    % Pamadi Ecuacion 4.549
    %K           = ((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - (k1*CLa + 2*k2*CL*CLa))/((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - 2*CL*CLa/pi/ARexp_w);
    Cyp_CL_M0   = Cyp_CL_M0_calc(LAMc4_w*180/pi, TR_w);
    B_prandtl = sqrt(1 - (Minf^2)*(cos(LAMc4_w))^2);
    Cyp_CL_CL0  = (AR_we + 4*cos(LAMc4_w))/(AR_we*B_prandtl + 4*cos(LAMc4_w))*(AR_we*B_prandtl + cos(LAMc4_w))/(AR_we + cos(LAMc4_w))*Cyp_CL_M0;
    t_c = Geo_tier.t_c_flap;

    %% Perfil revisado para el Ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_HTP1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
    clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
    k_Clp       = clalpha_w2_M0/2/pi;
    betaClp_k   = C_lp_calc(TR_w, AR_w, LAMc4_w, Minf);
    Clp_die0    = betaClp_k*k_Clp/beta;
    DCy_p_die   = (3*sin(diedro_w)*(1 - 4*(Zca_w-Zcg)/b_w*sin(diedro_w)))*Clp_die0;
    Cy_p_HTP      = K*Cyp_CL_CL0*CL_w + DCy_p_die; %DATCOM 7.1.2.1
    % z is the vertical distance between the c.g. and the wing root
    % quarter-chord point, positive for the c.g. above the wing root chord.
    % DATCOM 7.1.2.1 (pagina 2523)
    %DATCOM  7.4.2.1 pag.332/444

    if isnan(Cy_p_w1)
        Cy_p_HTP = 0;
    end
else
    Cy_p_HTP = 0;
end

% Vee & VTP
if VTP == 1
    % Vertical
    Cy_p_VTP   = 2/b_w*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)*Cy_beta_VTP;

    % Twin Vertical Tail configuration
    if AC_CONFIGURATION.twin_VTP == 1
        Cy_p_VTP = 2*Cy_p_VTP;
    end

    if isnan(Cy_p_VTP)
        Cy_p_VTP = 0;
    end
else
    Cy_p_VTP = 0;
end

if Vee == 1
    % Vertical
    Cy_p_vee   = 2/b_w*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)*Cy_beta_vee;

    if isnan(Cy_p_vee)
        Cy_p_vee = 0;
    end
else
    Cy_p_vee = 0;
end


% Vee 2
if Vee2 == 1
    % Vertical
    Cy_p_vee2   = 2/b_w*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)*Cy_beta_vee2;

    if isnan(Cy_p_vee2)
        Cy_p_vee2 = 0;
    end
else
    Cy_p_vee2 = 0;
end

% Vee 2
if Nac == 1
    Cy_p_nac = 0;
else
    Cy_p_nac = 0;
end

% % Vertical
% Cy_p_vert   = 2/b_w*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)*Cy_beta_VTP;
% 
% % Twin Vertical Tail configuration
% if AC_CONFIGURATION.twin_VTP == 1
%     Cy_p_VTP = 2*Cy_p_VTP;
% end
% 
% if isnan(Cy_p_vert)
%     Cy_p_VTP = 0;
% end

%% DERIVADA TOTAL
Cy_p = Cy_p_w1 + Cy_p_VTP + Cy_p_vee + Cy_p_vee2 + + Cy_p_HTP + Cy_p_can + Cy_p_nac;
Stab_Der.Cyp_w_diedro = DCy_p_die;
Stab_Der.Cyp_w1 = Cy_p_w1;
Stab_Der.Cyp_VTP = Cy_p_VTP;
Stab_Der.Cyp_vee = Cy_p_vee;
Stab_Der.Cyp_vee2 = Cy_p_vee2;
Stab_Der.Cyp_can = Cy_p_can;
Stab_Der.Cyp_nac = Cy_p_nac;
Stab_Der.Cyp_HTP = Cy_p_HTP;
Stab_Der.Cyp = Cy_p;

% Stores Derivatives per parts
Stab_Der_parts.Cy_p = Cy_p;
Stab_Der_parts.Cy_p_w1 = Cy_p_w1;
Stab_Der_parts.Cy_p_VTP = Cy_p_VTP;
Stab_Der_parts.Cy_p_vee = Cy_p_vee;
Stab_Der_parts.Cy_p_vee2 = Cy_p_vee2;
Stab_Der_parts.Cy_p_can = Cy_p_can;
Stab_Der_parts.Cy_p_HTP = Cy_p_HTP;
Stab_Der_parts.Cy_p_nac = Cy_p_nac;
end