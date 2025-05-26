function [Stab_Der_parts Stab_Der] = getClp_v2(AC_CONFIGURATION,modelo,alpha,Stab_Der,Stab_Der_parts, Geo_tier, OUTPUT_read_XLSX)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

Zca_v     = modelo.vertical.Zca;
Xca_v     = modelo.vertical.Xca;

b_w       = modelo.ala.b;
Zca_w     = -Geo_tier.z_w1_LE;
diedro_w  = modelo.ala.diedro;
TR_w      = modelo.ala.TR;
AR_w      = modelo.ala.AR;
LAMc4_w   = modelo.ala.LAMc4;

Xcg       = modelo.general.Xcg;
Minf      = modelo.general.Minf;

% Cy_beta_vert   = Stab_Der.Cyb_v;
% Cy_p_vert      = Stab_Der.Cyp_v;

beta = sqrt(1 - Minf^2);
%% Cl_p

if VTP == 1
    Cy_beta_VTP = Stab_Der_parts.Cy_beta_VTP;
    % Pamadi Ecuacion 4.579, Pagina 412
    Cl_p_VTP           = abs(2*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w)*(((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)/b_w))*Cy_beta_VTP;

    % Twin Vertical Tail configuration
    if AC_CONFIGURATION.twin_VTP == 1
        Cl_p_VTP = 2*Cl_p_VTP;
    end

    if isnan(Cl_p_VTP)
        Cl_p_VTP = 0;
    end
else
    Cl_p_VTP = 0;
end

if Vee == 1
    Cy_beta_vee = Stab_Der_parts.Cy_beta_vee;
    % Pamadi Ecuacion 4.579, Pagina 412
    Cl_p_vee           = abs(2*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w)*(((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)/b_w))*Cy_beta_vee;

    if isnan(Cl_p_vee)
        Cl_p_vee = 0;
    end
else
    Cl_p_vee = 0;
end

if Vee2 == 1
    Zca_v     = modelo.vertical2.Zca;
    Xca_v     = modelo.vertical2.Xca;

    Cy_beta_vee2 = Stab_Der_parts.Cy_beta_vee2;
    % Pamadi Ecuacion 4.579, Pagina 412
    Cl_p_vee2           = abs(2*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w)*(((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)/b_w))*Cy_beta_vee2;

    if isnan(Cl_p_vee2)
        Cl_p_vee2 = 0;
    end
else
    Cl_p_vee2 = 0;
end

% Ala
% Pamadi Ecuacion 4.576, Pagina 412
% Effect of drag on rolling moment is ignored

if W1 == 1

    % Ala
    
    LAMc4_w = Geo_tier.Lambda_c4_w1;
    AR_w    = Geo_tier.AR_w1;
    TR_w    = Geo_tier.lambda_w1;
    b_w     = Geo_tier.b_w1;
    diedro_w  = Geo_tier.dihedral_w1;
    Zca_w       = -Geo_tier.z_w1_LE; % positive if the wing is below thefuselage centerline.
    t_c = Geo_tier.t_c_flap;
    
    %% Perfil revisado para el Ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_w1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
    clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
    k_Clp       = clalpha_w2_M0/2/pi;
    betaClp_k   = C_lp_calc(TR_w, AR_w, LAMc4_w, Minf);
    Clp_ClpDiedro0  = 1 - 4*Zca_w/b_w*sin(diedro_w) + 3*(2*Zca_w/b_w)^2*(sin(diedro_w))^2;
    Cl_p_w1          = betaClp_k*k_Clp/beta*Clp_ClpDiedro0;

    if isnan(Cl_p_w1)
        Cl_p_w1 = 0;
    end
else
    Cl_p_w1 = 0;

end

%% Canard
if Can == 1

    % Ala
    
    LAMc4_w = Geo_tier.Lambda_c4_can;
    AR_w    = Geo_tier.AR_can;
    TR_w    = Geo_tier.lambda_can;
    b_w     = Geo_tier.b_can;
    diedro_w  = Geo_tier.dihedral_can;
    Zca_w       = -Geo_tier.z_can_LE; % positive if the wing is below thefuselage centerline.
    t_c = Geo_tier.t_c_flap;
    
    %% Perfil revisado para el Ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_can1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
    clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
    k_Clp       = clalpha_w2_M0/2/pi;
    betaClp_k   = C_lp_calc(TR_w, AR_w, LAMc4_w, Minf);
    Clp_ClpDiedro0  = 1 - 4*Zca_w/b_w*sin(diedro_w) + 3*(2*Zca_w/b_w)^2*(sin(diedro_w))^2;
    Cl_p_can          = betaClp_k*k_Clp/beta*Clp_ClpDiedro0;
    
    if isnan(Cl_p_can)
        Cl_p_can = 0;
    end
else
    Cl_p_can = 0;

end

%% HTP
if HTP == 1

    % Ala
    
    LAMc4_w = Geo_tier.Lambda_c4_HTP;
    AR_w    = Geo_tier.AR_HTP;
    TR_w    = Geo_tier.lambda_HTP;
    b_w     = Geo_tier.b_HTP;
    diedro_w  = Geo_tier.dihedral_HTP;
    Zca_w       = -Geo_tier.z_HTP_LE; % positive if the wing is below thefuselage centerline.
    t_c = Geo_tier.t_c_flap;
    
    %% Perfil revisado para el Ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_HTP1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
    clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
    k_Clp       = clalpha_w2_M0/2/pi;
    betaClp_k   = C_lp_calc(TR_w, AR_w, LAMc4_w, Minf);
    Clp_ClpDiedro0  = 1 - 4*Zca_w/b_w*sin(diedro_w) + 3*(2*Zca_w/b_w)^2*(sin(diedro_w))^2;
    Cl_p_HTP          = betaClp_k*k_Clp/beta*Clp_ClpDiedro0;
    
    if isnan(Cl_p_HTP)
        Cl_p_HTP = 0;
    end
else
    Cl_p_HTP = 0;

end

% Nacelle
if Nac == 1
    Cl_p_nac          = 0;
    
    if isnan(Cl_p_HTP)
        Cl_p_nac = 0;
    end
else
    Cl_p_nac = 0;

end

% DERIVADA TOTAL

Cl_p            = Cl_p_w1 + Cl_p_VTP + Cl_p_vee + Cl_p_vee2 + Cl_p_HTP + Cl_p_can + Cl_p_nac;

Stab_Der.Clp_w1 = Cl_p_w1;

Stab_Der.Clp_w1 = Cl_p_w1;
Stab_Der.Clp_can = Cl_p_can;
Stab_Der.Clp_vee = Cl_p_vee;
Stab_Der.Clp_vee2 = Cl_p_vee2;
Stab_Der.Clp_nac = Cl_p_nac;
Stab_Der.Clp_VTP = Cl_p_VTP;
Stab_Der.Clp_HTP = Cl_p_HTP;
Stab_Der.Clp = Cl_p;

% Stores Derivatives per parts
Stab_Der_parts.Cl_p = Cl_p;
Stab_Der_parts.Cl_p_w1 = Cl_p_w1;
Stab_Der_parts.Cl_p_can = Cl_p_can;
Stab_Der_parts.Cl_p_HTP = Cl_p_HTP;
Stab_Der_parts.Cl_p_vee = Cl_p_vee;
Stab_Der_parts.Cl_p_vee2 = Cl_p_vee2;
Stab_Der_parts.Cl_p_VTP = Cl_p_VTP;
Stab_Der_parts.Cl_p_nac = Cl_p_nac;

end