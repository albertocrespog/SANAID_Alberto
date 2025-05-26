function [Stab_Der_parts Stab_Der] = getCnp_v2(AC_CONFIGURATION,modelo,alpha,Stab_Der,Stab_Der_parts,Trim_ITER, Geo_tier, OUTPUT_read_XLSX)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

AR_we      = modelo.ala.ARwe;
LAMc4_w    = modelo.ala.LAMc4;
CL_w1         = Trim_ITER.CL_w1;
CLa_we      = Stab_Der_parts.CLalpha_w1_e_pw;
b_w        = modelo.ala.b;
eOswald_w  = modelo.ala.oswald;
Xca_w      = modelo.ala.Xca;
cMAC_w     = modelo.ala.MAC;

Zca_v      = modelo.vertical.Zca;
Xca_v      = modelo.vertical.Xca;

Xcg        = modelo.general.Xcg;
Minf       = modelo.general.Minf;

Cy_beta_VTP = Stab_Der_parts.Cy_beta_VTP;
Cy_beta_vee2 = Stab_Der_parts.Cy_beta_vee2;
Cy_beta_vee = Stab_Der_parts.Cy_beta_vee;


beta = sqrt(1 - Minf^2);


%% Cn_p

%% Ala
if W1 == 1

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


    Xca_w     = Geo_tier.x_xbar_w1;
    cMAC_w   = Geo_tier.cmac_w1;

    % Ala
    aw1         = CLa_we/(pi*AR_we*eOswald_w);
    K           = (1-aw1)/(1-eOswald_w*aw1);    % Pamadi Ecuacion 4.549
    xi         = (Xca_w - Xcg)/cMAC_w;

    Cnp_CL_CL0_M0   = -(AR_we + 6*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)*(xi/AR_we + tan(LAMc4_w)/12))/(6*(AR_we + 4*cos(LAMc4_w)));
    Cnp_CL_CL0      = ((AR_we + 4*cos(LAMc4_w))/((AR_we*beta + 4*cos(LAMc4_w))))*...
        ((AR_we*beta + 0.5*((AR_we*beta + cos(LAMc4_w)))*(tan(LAMc4_w)^2))/...
        (AR_we + 0.5*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)^2))*Cnp_CL_CL0_M0;

    Cl_p_w       = Stab_Der.Clp_w1;

    Cn_p_w1          = Cl_p_w*tan(alpha)*(K-1) + K*Cnp_CL_CL0*CL_w; %Pamadi 4.594
    if isnan(Cn_p_w1)
        Cn_p_w1 = 0;
    end
else
    Cn_p_w1 = 0;
end


%% CAnard
if Can == 1

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


    Xca_w     = Geo_tier.x_xbar_can;
    cMAC_w   = Geo_tier.cmac_can;

    % Ala
    aw1         = CLa_we/(pi*AR_we*eOswald_w);
    K           = (1-aw1)/(1-eOswald_w*aw1);    % Pamadi Ecuacion 4.549
    xi         = (Xca_w - Xcg)/cMAC_w;

    Cnp_CL_CL0_M0   = -(AR_we + 6*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)*(xi/AR_we + tan(LAMc4_w)/12))/(6*(AR_we + 4*cos(LAMc4_w)));
    Cnp_CL_CL0      = ((AR_we + 4*cos(LAMc4_w))/((AR_we*beta + 4*cos(LAMc4_w))))*...
        ((AR_we*beta + 0.5*((AR_we*beta + cos(LAMc4_w)))*(tan(LAMc4_w)^2))/...
        (AR_we + 0.5*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)^2))*Cnp_CL_CL0_M0;

    Cl_p_w       = Stab_Der.Clp_can;
    Cn_p_can          = Cl_p_w*tan(alpha)*(K-1) + K*Cnp_CL_CL0*CL_w; %Pamadi 4.594

    if isnan(Cn_p_can)
        Cn_p_can = 0;
    end
else
    Cn_p_can = 0;
end


%% HTP
if HTP == 1

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


    Xca_w     = Geo_tier.x_xbar_HTP;
    cMAC_w   = Geo_tier.cmac_HTP;

    % Ala
    aw1         = CLa_we/(pi*AR_we*eOswald_w);
    K           = (1-aw1)/(1-eOswald_w*aw1);    % Pamadi Ecuacion 4.549
    xi         = (Xca_w - Xcg)/cMAC_w;

    Cnp_CL_CL0_M0   = -(AR_we + 6*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)*(xi/AR_we + tan(LAMc4_w)/12))/(6*(AR_we + 4*cos(LAMc4_w)));
    Cnp_CL_CL0      = ((AR_we + 4*cos(LAMc4_w))/((AR_we*beta + 4*cos(LAMc4_w))))*...
        ((AR_we*beta + 0.5*((AR_we*beta + cos(LAMc4_w)))*(tan(LAMc4_w)^2))/...
        (AR_we + 0.5*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)^2))*Cnp_CL_CL0_M0;

    Cl_p_w       = Stab_Der.Clp_HTP;
    Cn_p_HTP          = Cl_p_w*tan(alpha)*(K-1) + K*Cnp_CL_CL0*CL_w; %Pamadi 4.594

    if isnan(Cn_p_HTP)
        Cn_p_HTP = 0;
    end
else
    Cn_p_HTP = 0;
end

if Nac == 1
    Cn_p_nac = 0;
else
    Cn_p_nac = 0;
end

if VTP == 1
    % Vertical
    Cn_p_VTP       = -2/b_w*(Zca_v*sin(alpha) + (Xca_v - Xcg)*cos(alpha))*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)/b_w*Cy_beta_VTP;

    if isnan(Cn_p_VTP)
        Cn_p_VTP = 0;
    end
    
    % Twin Vertical Tail configuration
    if AC_CONFIGURATION.twin_VTP == 1
        Cn_p_VTP = 2*Cn_p_VTP;
    end

else
    Cn_p_VTP = 0;
end

if Vee == 1
    % Vertical
    Cn_p_vee       = -2/b_w*(Zca_v*sin(alpha) + (Xca_v - Xcg)*cos(alpha))*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)/b_w*Cy_beta_vee;

    if isnan(Cn_p_vee)
        Cn_p_vee = 0;
    end
else
    Cn_p_vee = 0;
end

if Vee2 == 1
    Zca_v2      = modelo.vertical2.Zca;
    Xca_v2      = modelo.vertical2.Xca;
    % Vertical
    Cn_p_vee2       = -2/b_w*(Zca_v2*sin(alpha) + (Xca_v2 - Xcg)*cos(alpha))*((Zca_v2*cos(alpha) - (Xca_v2 - Xcg)*sin(alpha)) - Zca_v2)/b_w*Cy_beta_vee2;

    if isnan(Cn_p_vee)
        Cn_p_vee2 = 0;
    end
else
    Cn_p_vee2 = 0;
end




% DERIVADA TOTAL
Cn_p            = Cn_p_w1 + Cn_p_VTP + Cn_p_vee + Cn_p_vee2 + Cn_p_nac + Cn_p_can + Cn_p_HTP;

Stab_Der.Cnp_w1 = Cn_p_w1;
Stab_Der.Cnp_VTP = Cn_p_VTP;
Stab_Der.Cnp_vee = Cn_p_vee;
Stab_Der.Cnp_vee2 = Cn_p_vee2;
Stab_Der.Cnp_nac = Cn_p_nac;
Stab_Der.Cnp_can = Cn_p_can;
Stab_Der.Cnp_HTP = Cn_p_HTP;
Stab_Der.Cnp = Cn_p;

% Stores Derivatives per parts
Stab_Der_parts.Cn_p = Cn_p;
Stab_Der_parts.Cn_p_w1 = Cn_p_w1;
Stab_Der_parts.Cn_p_VTP = Cn_p_VTP;
Stab_Der_parts.Cn_p_vee = Cn_p_vee;
Stab_Der_parts.Cn_p_vee2 = Cn_p_vee2;
Stab_Der_parts.Cn_p_nac = Cn_p_nac;
Stab_Der_parts.Cn_p_can = Cn_p_can;
Stab_Der_parts.Cn_p_HTP = Cn_p_HTP;

end
