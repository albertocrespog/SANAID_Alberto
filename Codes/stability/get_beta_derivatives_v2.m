function [Stab_Der_parts,Stab_Der] = get_beta_derivatives_v2(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,Geo_tier,TRIM_RESULTS,Trim_ITER,Body_Geo,Aero,conditions,OUTPUT_read_XLSX)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

%%%%%%%%%%%%%%%%%%%%%derivadas en funcion de beta%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cyb
% ASPro Cyb
%     Stab_Der = getCybeta(modelo,Stab_Der);
%     Cyb_w = Stab_Der.Cyb_w;
%     Cyb_fus = Stab_Der.Cyb_fus;
%     Cyb_v = Stab_Der.Cyb_v;
%     Cyb = Stab_Der.Cyb;

qinf    = modelo.general.qinf;
W       = modelo.general.mtow*modelo.general.w_w0*9.8065;
Sref    = modelo.general.Sref;
C_L      = modelo.general.CL;
C_Lw   = modelo.general.CL_w;
%     C_Lh   = modelo.general.CL_h;

AR_w    = modelo.ala.AR;

z_w     = modelo.ala.Zca;
LAM_w   = modelo.ala.LAMc2;
LAMc4_w = modelo.ala.LAMc4;
diedro  = modelo.ala.diedro;
dihedral_w1 = diedro;

S_v     = modelo.vertical.S;
eta_v   = modelo.vertical.eta;
CLa_v   = modelo.vertical.CLa;
b_v     = modelo.vertical.b;
Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca, 'pchip');

vol_fus = modelo.fuselaje.vol;
CLa_fus = modelo.fuselaje.CLa;


CL_w1 = C_Lw;
Lambda_c2_w1 = LAM_w;
AR_w1 = modelo.ala.AR;
x_Area_body = modelo.fuselaje.x;
Area_body = modelo.fuselaje.D_x;
xbar_w1 = modelo.ala.xca;
x_xbar_w1 = Geo_tier.x_xbar_w1;
%     z_zbar_w1 = Geo_tier.z_zbar_w1;
z_w1_LE1 = Geo_tier.z_zbar_w1;
% Positrive bellow fuselage center line
z_w1_LE1 = -Geo_tier.z_w1_LE;

Vol_TOT = modelo.fuselaje.vol;
CLalpha_fus = Stab_Der_parts.CL_alpha_fus;
S_ref = Geo_tier.S_ref;
trim_alpha = TRIM_RESULTS.trim_alpha;
b_w1 = Geo_tier.b_w1;
lambda_w1 = Geo_tier.lambda_w1;
x_XCG = conditions.x_XCG;
alpha = TRIM_RESULTS.trim_alpha;
Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
length_fus = Body_Geo.l_fus;
height_x_position = Body_Geo.height_x_position;
width_x_position = Body_Geo.width_x_position;
length_x_position = Body_Geo.length_x_position;
%% Cy_beta = Cy_beta_fus + Cy_beta_wing + Cy_beta_VTP

if W1 ==1
    if OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL ==1
        Cy_beta_w1 = Aero.CYbeta_w1;
        Cl_beta_w1 = Aero.Clbeta_w1;
        Cn_beta_w1 = Aero.Cnbeta_w1;

        Cy_beta_w1_CL = 0;
        Cy_beta_w1_dihedral = 0;
    else
        %% APORTE ALA
        % - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
        Cy_beta_w1_dihedral       = -0.00573*abs(dihedral_w1*180/pi);

        % - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
        C_L = W/qinf/S_ref;
        CL_w1 = Trim_ITER.CL_w1;
        Cy_beta_w1_CL       = (CL_w1^2)*(6*tan(Lambda_c2_w1)*sin(Lambda_c2_w1))/(pi*AR_w1*(AR_w1 + 4*cos(Lambda_c2_w1)));

        % total contribution wing
        Cy_beta_w1 = Cy_beta_w1_CL + Cy_beta_w1_dihedral;
    end
else
    Cy_beta_w1_CL = 0;
    Cy_beta_w1_dihedral = 0;
    Cy_beta_w1 = 0;
end
% W1

Stab_Der_parts.Cy_beta_w1 = Cy_beta_w1;
Stab_Der_parts.Cy_beta_w1_CL = Cy_beta_w1_CL;
Stab_Der_parts.Cy_beta_w1_dihedral = Cy_beta_w1_dihedral;


if Can ==1
    if OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL ==1
        Cy_beta_can = Aero.CYbeta_can;
        Cl_beta_can = Aero.Clbeta_can;
        Cn_beta_can = Aero.Cnbeta_can;

        Cy_beta_can_CL = 0;
        Cy_beta_can_dihedral = 0;
    else
        % APORTE Canard
        Lambda_c2_can = Geo_tier.Lambda_c2_can;
        Lambda_c4_can = Geo_tier.Lambda_c4_can;
        dihedral_can = Geo_tier.dihedral_can;
        AR_can = Geo_tier.AR_can;
        CL_can = Trim_ITER.CL_can;

        % - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
        Cy_beta_can_dihedral       = -0.00573*abs(dihedral_can*180/pi);

        % - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
        C_L = W/qinf/S_ref;
        Cy_beta_can_CL       = (CL_can^2)*(6*tan(Lambda_c2_can)*sin(Lambda_c2_can))/(pi*AR_can*(AR_can + 4*cos(Lambda_c2_can)));

        % total contribution
        Cy_beta_can = Cy_beta_can_CL + Cy_beta_can_dihedral;
    end
else
    Cy_beta_can_CL = 0;
    Cy_beta_can_dihedral = 0;
    Cy_beta_can = 0;
end

% Can
Stab_Der_parts.Cy_beta_can = Cy_beta_can;
Stab_Der_parts.Cy_beta_can_CL = Cy_beta_can_CL;
Stab_Der_parts.Cy_beta_can_dihedral = Cy_beta_can_dihedral;

if HTP ==1
    if OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL ==1
        Cy_beta_HTP = Aero.CYbeta_HTP;
        Cl_beta_HTP = Aero.Clbeta_HTP;
        Cn_beta_HTP = Aero.Cnbeta_HTP;

        Cy_beta_HTP_CL = 0;
        Cy_beta_HTP_dihedral = 0;

    else
        %% APORTE Canard
        Lambda_c2_HTP = Geo_tier.Lambda_c2_HTP;
        Lambda_c4_HTP = Geo_tier.Lambda_c4_HTP;
        dihedral_HTP = Geo_tier.dihedral_HTP;
        AR_HTP = Geo_tier.AR_HTP;
        CL_HTP = Trim_ITER.CL_HTP;

        % - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
        Cy_beta_HTP_dihedral       = -0.00573*abs(dihedral_HTP*180/pi);

        % - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
        C_L = W/qinf/S_ref;
        Cy_beta_HTP_CL       = (CL_HTP^2)*(6*tan(Lambda_c2_HTP)*sin(Lambda_c2_HTP))/(pi*AR_HTP*(AR_HTP + 4*cos(Lambda_c2_HTP)));
        % total contribution
        Cy_beta_HTP = Cy_beta_HTP_CL + Cy_beta_HTP_dihedral;
    end

else
    Cy_beta_HTP_CL = 0;
    Cy_beta_HTP_dihedral = 0;
    Cy_beta_HTP = 0;
end
% HTP
Stab_Der_parts.Cy_beta_HTP = Cy_beta_HTP;
Stab_Der_parts.Cy_beta_HTP_CL = Cy_beta_HTP_CL;
Stab_Der_parts.Cy_beta_HTP_dihedral = Cy_beta_HTP_dihedral;

%% APORTE FUSELAJE
Dfus_w1  = interp1(length_x_position,height_x_position,x_xbar_w1, 'pchip');
% vertical distance between the fuselage line and wing root (positive if wing below center line)
%     Ki  = Ki_calc(-z_zbar_w1, Dfus_w1);
Ki  = Ki_calc(-z_w1_LE1, Dfus_w1);

% - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
%Cy_beta_fus     = -2*Ki*S0/Sref;

% - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
A_refFus        = Vol_TOT^(2/3);
Cy_beta_fus     = -Ki*CLalpha_fus*A_refFus/S_ref;

if VTP == 1
    if OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL ==1
        Cy_beta_VTP = Aero.CYbeta_VTP;
        Cl_beta_VTP = Aero.Clbeta_VTP;
        Cn_beta_VTP = Aero.Cnbeta_VTP;

    else
        CLalpha_VTP = Aero.CL_alpha_VTP_CR;
        Cy_beta_VTP = CLalpha_VTP;
    end

    twin_VTP = AC_CONFIGURATION.twin_VTP;
    %% APORTE VERTICAL
    % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
    % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

    x_xbar_VTP = Geo_tier.x_xbar_VTP;
    if twin_VTP ==1
        S_VTP = 2*Geo_tier.S_VTP;
    else
        S_VTP = Geo_tier.S_VTP;
    end

    % CLalpha_VTP =  Aero.CL_alpha_VTP_CR;
    b_VTP_s = Geo_tier.b_VTP_s;
    b_VTP = Geo_tier.b_VTP;
    b_w2 = Geo_tier.b_w2;
    % x_bar_VTP = Geo_tier.xbar_VTP;
    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash = 0.724 + 3.06*(S_VTP/S_ref)/(1+cos(Lambda_c4_w1)) + 0.4*(z_w1_LE1)/Dfus_w1 + 0.009*AR_w1;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

    % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
    % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
    Dfus_VTP  = interp1(length_x_position,height_x_position, x_xbar_VTP, 'pchip');
    % Avoids negative relation
    if Dfus_VTP <0
        Dfus_VTP = 0;
    end

    k       = k_calc(b_VTP_s, Dfus_VTP);

    % Single Vertical tail
    Cy_beta_VTP    = -k*Cy_beta_VTP*sidewash*S_VTP/S_ref;

    if twin_VTP ==1
        % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
        % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
        if AC_CONFIGURATION.AC_type == 1 % AC_type = 1 - flying wing
            b_surf = b_w1;
        else
            b_surf = b_w2;
        end
        Cybeta_v_Cybeta_eff = Cybeta_v_Cybeta_eff_calc(b_surf,length_fus,Dfus_VTP/2,b_VTP);
        %         Cy_beta_VTP = 2*Cy_beta_VTP*Cybeta_v_Cybeta_eff*(S_VTP/S_ref);
        Cy_beta_VTP = 2*Cy_beta_VTP*Cybeta_v_Cybeta_eff;
    end
else
    Cy_beta_VTP = 0;
end
% VTP
Stab_Der_parts.Cy_beta_VTP = Cy_beta_VTP;

if Vee == 1
    if OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL ==1
        Cy_beta_vee = Aero.CYbeta_vee;
        Cl_beta_vee = Aero.Clbeta_vee;
        Cn_beta_vee = Aero.Cnbeta_vee;

    else
        CYbeta_vee = Aero.CYbeta_vee;
        Cy_beta_vee = CYbeta_vee;

    end

    S_vee_pv = Geo_tier.S_vee_pv;
    S_vee = Geo_tier.S_vee;
    S_vee_s = Geo_tier.S_vee_s;
    x_xbar_vee = Geo_tier.x_xbar_vee;
    b_vee_s = Geo_tier.b_vee_s;
    dihedral_vee = Geo_tier.dihedral_vee;
    % CYbeta_vee = Aero.CYbeta_vee;

    %% APORTE VERTICAL
    % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
    % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash_vee = 0.724 + 3.06*(S_vee_pv/S_ref)/(1 + cos(Lambda_c4_w1)) + 0.4*(z_w1_LE1)/Dfus_w1 + 0.009*AR_w1;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))
    % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
    % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
    Dfus_vee  = interp1(length_x_position,height_x_position, x_xbar_vee, 'pchip');
    Dfus_vee  = interp1(length_x_position,width_x_position, x_xbar_vee, 'pchip');
    % Avoids negative relation
    if Dfus_vee <0
        Dfus_vee = 0;
    end

    k       = k_calc(b_vee_s, Dfus_vee);
    % Single Vertical tail
    %     Cy_beta_vee2    = -k*abs(CYbeta_vee)*sidewash*S_vee/S_ref;
    Cy_beta_vee    = -k*abs(Cy_beta_vee)*sidewash_vee*S_vee_s/S_ref;

    %     Cy_beta_vee = -0.483133428328602;
else
    Cy_beta_vee = 0;
end
% Vee
Stab_Der_parts.Cy_beta_vee = Cy_beta_vee;

if Vee2 == 1
    if OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL ==1
        Cy_beta_vee2 = Aero.CYbeta_vee2;
        Cl_beta_vee2 = Aero.Clbeta_vee2;
        Cn_beta_vee2 = Aero.Cnbeta_vee2;

    else
        CYbeta_vee2 = Aero.CYbeta_vee;
        Cy_beta_vee2 = CYbeta_vee2;

    end

    S_vee2_pv = Geo_tier.S_vee2_pv;
    S_vee2 = Geo_tier.S_vee2;
    S_vee2_s = Geo_tier.S_vee2_s;
    x_xbar_vee2 = Geo_tier.x_xbar_vee2;
    b_vee2_s = Geo_tier.b_vee2_s;
    dihedral_vee2 = Geo_tier.dihedral_vee2;
    % CYbeta_vee2 = Aero.CYbeta_vee2;
    %% APORTE VERTICAL
    % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
    % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash_vee2 = 0.724 + 3.06*(S_vee2_pv/S_ref)/(1 + cos(Lambda_c4_w1)) + 0.4*(z_w1_LE1)/Dfus_w1 + 0.009*AR_w1;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))
    % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
    % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
    Dfus_vee2  = interp1(length_x_position,height_x_position, x_xbar_vee2, 'pchip');
    Dfus_vee2  = interp1(length_x_position,width_x_position, x_xbar_vee2, 'pchip');
    % Avoids negative relation
    if Dfus_vee2 <0
        Dfus_vee2 = 0;
    end

    k       = k_calc(b_vee2_s, Dfus_vee2);
    % Single Vertical tail
%     Cy_beta_vee22    = -k*abs(CYbeta_vee2)*sidewash*S_vee2/S_ref;
    Cy_beta_vee2    = -k*abs(Cy_beta_vee2)*sidewash_vee2*S_vee2_s/S_ref;

%     Cy_beta_vee2 = -0.483133428328602;
else
    Cy_beta_vee2 = 0;
end
% Vee
Stab_Der_parts.Cy_beta_vee2 = Cy_beta_vee2;


if Nac == 1
    Cy_beta_nac = 0;
else
    Cy_beta_nac = 0;
end

%% DERIVADA TOTAL
Cy_beta         = Cy_beta_w1 + Cy_beta_fus + Cy_beta_VTP + Cy_beta_vee + Cy_beta_HTP + Cy_beta_can; 
Stab_Der.Cyb_w = Cy_beta_w1;
Stab_Der.Cyb_fus = Cy_beta_fus;
Stab_Der.Cyb_VTP = Cy_beta_VTP;
Stab_Der.Cyb_vee = Cy_beta_vee;
Stab_Der.Cyb_vee2 = Cy_beta_vee2;
Stab_Der.Cyb_HTP = Cy_beta_HTP;
Stab_Der.Cyb_can = Cy_beta_can;
Stab_Der.Cyb = Cy_beta;

% Stores Derivatives per parts
Stab_Der_parts.Cy_beta = Cy_beta;

% Fus
Stab_Der_parts.Cy_beta_fus = Cy_beta_fus;
% NAC
Stab_Der_parts.Cy_beta_nac = Cy_beta_nac;

if VTP == 1
    Cy_beta_VTP = Cy_beta_VTP;
    z_zbar_VTP = Geo_tier.z_zbar_VTP;
    Stab_Der_parts.Cy_beta_VTP = Cy_beta_VTP;
elseif Vee == 1
    Cy_beta_vee = Cy_beta_vee;
    z_zbar_vee = Geo_tier.z_zbar_vee;
    Stab_Der_parts.Cy_beta_vee = Cy_beta_vee;
elseif Vee2 == 1
    Cy_beta_vee2 = Cy_beta_vee2;
    z_zbar_vee2 = Geo_tier.z_zbar_vee2;
    Stab_Der_parts.Cy_beta_vee2 = Cy_beta_vee2;
end

%% Clb
% ASPro Clb
% IMPORTANT MISSING TO INCLUDE
% The wing-fuselage, horizontal tail-fuselage and/or canard-fuselage contributions to the dihedral effect are found from:
%     Stab_Der = getClbeta(modelo,trim_alpha,Stab_Der);
%     Clb_wb = Stab_Der.Clb_wb;
%     Clb_v = Stab_Der.Clb_v;
%     Clb = Stab_Der.Clb;

Mach       = modelo.general.Minf;
% qinf    = modelo.general.qinf;
% W       = modelo.general.mtow*modelo.general.w_w0*9.8065;
Sref    = modelo.general.Sref;
% C_L      = modelo.general.CL;
C_Lw   = modelo.general.CL_w;
%     C_Lh   = modelo.general.CL_h;

D_fus   = modelo.fuselaje.D;
if modelo.vertical.Xca < modelo.fuselaje.l
    Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca);
else
    Dfus_v  = 0;
end

AR_w    = modelo.ala.AR;
b_w     = modelo.ala.b;
Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.ala.Xca, 'pchip');
z_w     = modelo.ala.Zca;
% Positive bellow fuselage center line
z_w1_LE1 = -Geo_tier.z_w1_LE;
TR_w    = modelo.ala.TR;
LAMc2_w = modelo.ala.LAMc2;
LAMc4_w = modelo.ala.LAMc4;
lt_w    = modelo.ala.Xca - modelo.ala.xca + modelo.ala.le_y(end) + modelo.ala.ct/2;
diedro  = modelo.ala.diedro;

S_v     = modelo.vertical.S;
b_v     = modelo.vertical.b;
Zca_v     = modelo.vertical.Zca;
CLa_v   = modelo.vertical.CLa;
Xca_v    = modelo.vertical.Xca;
eta_v   = modelo.vertical.eta;

x_XCG = modelo.general.Xcg;
z_XCG = modelo.general.Zcg;

%% FUSELAJE
%     Cl_beta_fus = 1.2*sqrt(AR_w1)*(-z_zbar_w1)/b_w1*2*Dfus_w1/b_w1; % 1/rad

if W1 ==1
    CL_w = Trim_ITER.CL_w1;
    %%  APORTE DEL ALA-FUSELAJE
    AR_w    = Geo_tier.AR_w1;
    b_w     = Geo_tier.b_w1;
    Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, Geo_tier.x_xbar_w1, 'pchip');
    % Positive bellow fuselage center line
    z_w = -Geo_tier.z_w1_LE;
    TR_w    = Geo_tier.lambda_w1;
    LAMc2_w = Geo_tier.Lambda_c2_w1;
    LAMc4_w = Geo_tier.Lambda_c4_w1;
    lt_w    = Geo_tier.x_xbar_w1 - Geo_tier.xbar_w1 + (Geo_tier.b_w1/2)*tan(Geo_tier.Lambda_LE_w1) + Geo_tier.cT_w1/2;
    diedro  = Geo_tier.dihedral_w1;
    lambda_w = Geo_tier.lambda_w1;

    % Estimacion de Cl_beta_wf extraida de:
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
    %   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)
    DClbeta_zw          = 1.2*sqrt(AR_w)*z_w/b_w*2*Dfus_w/b_w; % 1/rad
    DClbeta_tau         = -180/pi*0.0005*sqrt(AR_w)*(Dfus_w/b_w)^2; % 1/rad
    Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(lambda_w, AR_w, LAMc2_w); % 1/rad

    K_MLAM              = K_MLAM_calc(Mach, AR_w, Lambda_c2_w1);

    Kf                  = Kf_calc(lt_w, b_w, AR_w, LAMc2_w); % Corrected with DATCOM
    % Saturates to 1 if NaN
    if isnan(Kf)
        Kf = 1;
    end
    Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_w, lambda_w); % 1/rad
    Clbeta_die          = 180/pi*Clbeta_die_calc(lambda_w, AR_w, LAMc2_w); % 1/rad
    K_Mdie              = K_Mdie_calc(Mach, AR_w, LAMc2_w); % Corrected with DATCOM

    % C_L = W/qinf/Sref;
    %% APORTE ALA
    Cl_beta_wf_CL = C_Lw*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR);
    Cl_beta_wf_dihedral = diedro*(Clbeta_die*K_Mdie + DClbeta_tau);
    Cl_beta_wf_z = DClbeta_zw;

    Cl_beta_wf = Cl_beta_wf_CL + Cl_beta_wf_dihedral + Cl_beta_wf_z;
    % die: Angulo de diedro en radianes
else
    Cl_beta_wf = 0;
    Cl_beta_wf_CL = 0;
    Cl_beta_wf_dihedral = 0;
    Cl_beta_wf_z = 0;
end

if Can ==1
    %%  APORTE DEL ALA-FUSELAJE
    CL_w = Trim_ITER.CL_can;
    AR_w    = Geo_tier.AR_can;
    b_w     = Geo_tier.b_can;
    Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, Geo_tier.x_xbar_can, 'pchip');
    % Positive bellow fuselage center line
    z_w = -Geo_tier.z_can_LE;
    TR_w    = Geo_tier.lambda_can;
    LAMc2_w = Geo_tier.Lambda_c2_can;
    LAMc4_w = Geo_tier.Lambda_c4_can;
    lt_w    = Geo_tier.x_xbar_can - Geo_tier.xbar_can + (Geo_tier.b_can/2)*tan(Geo_tier.Lambda_LE_can) + Geo_tier.cT_can/2;
    diedro  = Geo_tier.dihedral_can;
    lambda_w = Geo_tier.lambda_can;

    % Estimacion de Cl_beta_wf extraida de:
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
    %   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)
    DClbeta_zw          = 1.2*sqrt(AR_w)*z_w/b_w*2*Dfus_w/b_w; % 1/rad
    DClbeta_tau         = -180/pi*0.0005*sqrt(AR_w)*(Dfus_w/b_w)^2; % 1/rad
    Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(lambda_w, AR_w, LAMc2_w); % 1/rad

    K_MLAM              = K_MLAM_calc(Mach, AR_w, Lambda_c2_w1);

    Kf                  = Kf_calc(lt_w, b_w, AR_w, LAMc2_w); % Corrected with DATCOM
    % Saturates to 1 if NaN
    if isnan(Kf)
        Kf = 1;
    end
    Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_w, lambda_w); % 1/rad
    Clbeta_die          = 180/pi*Clbeta_die_calc(lambda_w, AR_w, LAMc2_w); % 1/rad
    K_Mdie              = K_Mdie_calc(Mach, AR_w, LAMc2_w); % Corrected with DATCOM

    % C_L = W/qinf/Sref;
    %% Canard Contribution
    Cl_beta_can_CL = C_Lw*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR);
    Cl_beta_can_dihedral = diedro*(Clbeta_die*K_Mdie + DClbeta_tau);
    Cl_beta_can_z = DClbeta_zw;

    Cl_beta_can = Cl_beta_can_CL + Cl_beta_can_dihedral + Cl_beta_can_z;
    % die: Angulo de diedro en radianes
else
    Cl_beta_can = 0;
    Cl_beta_can_CL = 0;
    Cl_beta_can_dihedral = 0;
    Cl_beta_can_z = 0;
end

if HTP ==1
    %%  APORTE DEL ALA-FUSELAJE
    CL_w = Trim_ITER.CL_HTP;
    AR_w    = Geo_tier.AR_HTP;
    b_w     = Geo_tier.b_HTP;
    Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, Geo_tier.x_xbar_HTP, 'pchip');
    % Positive bellow fuselage center line
    z_w = -Geo_tier.z_HTP_LE;
    TR_w    = Geo_tier.lambda_HTP;
    LAMc2_w = Geo_tier.Lambda_c2_HTP;
    LAMc4_w = Geo_tier.Lambda_c4_HTP;
    lt_w    = Geo_tier.x_xbar_HTP - Geo_tier.xbar_HTP + (Geo_tier.b_HTP/2)*tan(Geo_tier.Lambda_LE_HTP) + Geo_tier.cT_HTP/2;
    diedro  = Geo_tier.dihedral_HTP;
    lambda_w = Geo_tier.lambda_HTP;

    % Estimacion de Cl_beta_wf extraida de:
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
    %   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)
    DClbeta_zw          = 1.2*sqrt(AR_w)*z_w/b_w*2*Dfus_w/b_w; % 1/rad
    DClbeta_tau         = -180/pi*0.0005*sqrt(AR_w)*(Dfus_w/b_w)^2; % 1/rad
    Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(lambda_w, AR_w, LAMc2_w); % 1/rad

    K_MLAM              = K_MLAM_calc(Mach, AR_w, Lambda_c2_w1);

    Kf                  = Kf_calc(lt_w, b_w, AR_w, LAMc2_w); % Corrected with DATCOM

    % Saturates to 1 if NaN
    if isnan(Kf)
        Kf = 1;
    end

    Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_w, lambda_w); % 1/rad
    Clbeta_die          = 180/pi*Clbeta_die_calc(lambda_w, AR_w, LAMc2_w); % 1/rad
    K_Mdie              = K_Mdie_calc(Mach, AR_w, LAMc2_w); % Corrected with DATCOM

    %% HTP Contribution
    Cl_beta_HTP_CL = C_Lw*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR);
    Cl_beta_HTP_dihedral = diedro*(Clbeta_die*K_Mdie + DClbeta_tau);
    Cl_beta_HTP_z = DClbeta_zw;
    Cl_beta_HTP = Cl_beta_HTP_CL + Cl_beta_HTP_dihedral + Cl_beta_HTP_z;
else
    Cl_beta_HTP = 0;
    Cl_beta_HTP_CL = 0;
    Cl_beta_HTP_dihedral = 0;
    Cl_beta_HTP_z = 0;
end

%% DERIVADA DE ESTABILIDAD COMPLETA
if VTP == 1

    %% APORTE VERTICAL
    % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
    % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash = 0.724 + 3.06*(S_VTP/S_ref)/(1+cos(Lambda_c4_w1)) + 0.4*(z_w1_LE1)/Dfus_w1 + 0.009*AR_w1;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

    % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
    % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
    Dfus_VTP  = interp1(length_x_position,height_x_position, x_xbar_VTP, 'pchip');
    % Avoids negative values
    if Dfus_VTP <0
        Dfus_VTP = 0;
    end
    k       = k_calc(b_VTP_s, Dfus_VTP);

    %         % Single Vertical tail
    %         Cy_beta_VTP    = -k*CLalpha_VTP*sidewash*S_VTP/S_ref;
    %         if twin_VTP ==1
    %             % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
    %             % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
    %             Dfus_VTP  = interp1(x_Area_body, Area_body, x_xbar_VTP, 'pchip');
    %             Cybeta_v_Cybeta_eff = Cybeta_v_Cybeta_eff_calc(b_w2,length_fus,Dfus_VTP/2,b_VTP);
    %             Cy_beta_VTP = 2*Cy_beta_VTP*Cybeta_v_Cybeta_eff*(S_VTP/S_ref);
    %         end
    z_zbar_VTP = Geo_tier.z_zbar_VTP;
    Cl_beta_VTP    = Cy_beta_VTP*((z_zbar_VTP - z_XCG)*cos(alpha)-(x_xbar_VTP - x_XCG)*sin(alpha))/b_w1;
else
    Cl_beta_VTP = 0;
end

if Vee == 1
    
    %% Contribution Horizontal projection
    CL_w = Trim_ITER.CL_vee;
    AR_w    = Geo_tier.AR_vee;
    b_w     = Geo_tier.b_vee;
    Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, Geo_tier.x_xbar_vee, 'pchip');
    % Positive bellow fuselage center line
    z_w = -Geo_tier.z_vee_LE;
    TR_w    = Geo_tier.lambda_vee;
    LAMc2_w = Geo_tier.Lambda_c2_vee;
    LAMc4_w = Geo_tier.Lambda_c4_vee;
    lt_w    = Geo_tier.x_xbar_vee - Geo_tier.xbar_vee + (Geo_tier.b_vee/2)*tan(Geo_tier.Lambda_LE_vee) + Geo_tier.cT_vee/2;
    diedro  = Geo_tier.dihedral_vee;
    lambda_w = Geo_tier.lambda_vee;

    % Estimacion de Cl_beta_wf extraida de:
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
    %   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)
    DClbeta_zw          = 1.2*sqrt(AR_w)*z_w/b_w*2*Dfus_w/b_w; % 1/rad
    DClbeta_tau         = -180/pi*0.0005*sqrt(AR_w)*(Dfus_w/b_w)^2; % 1/rad
    Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(lambda_w, AR_w, LAMc2_w); % 1/rad

    K_MLAM              = K_MLAM_calc(Mach, AR_w, Lambda_c2_w1);

    Kf                  = Kf_calc(lt_w, b_w, AR_w, LAMc2_w); % Corrected with DATCOM
    % Colapses Kf in case that the results is out of the scope
    if isnan(Kf)
        Kf  = 1;
    end
    
    Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_w, lambda_w); % 1/rad
    Clbeta_die          = 180/pi*Clbeta_die_calc(lambda_w, AR_w, LAMc2_w); % 1/rad
    K_Mdie              = K_Mdie_calc(Mach, AR_w, LAMc2_w); % Corrected with DATCOM

    % C_L = W/qinf/Sref;
    %% Canard Contribution
    Cl_beta_vee_CL = C_Lw*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR);
    Cl_beta_vee_dihedral = diedro*(Clbeta_die*K_Mdie + DClbeta_tau);
    Cl_beta_vee_z = DClbeta_zw;

    % Do not consider the dihedral contribution since that it is taken into
    % account wth CYbeta from FLOW or from NACA Report 823
    Cl_beta_vee_dihedral = 0;
    Cl_beta_vee_HP = Cl_beta_vee_CL + Cl_beta_vee_dihedral + Cl_beta_vee_z; % Horizontal projection
            
    z_zbar_vee = Geo_tier.z_zbar_vee;
    x_xbar_vee = Geo_tier.x_xbar_vee;
    %% APORTE VERTICAL
    % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
    % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash = 0.724 + 3.06*(S_vee_pv/S_ref)/(1 + cos(Lambda_c4_w1)) + 0.4*(z_w1_LE1)/Dfus_w1 + 0.009*AR_w1;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

    % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
    % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
    Dfus_vee  = interp1(length_x_position,width_x_position, x_xbar_vee, 'pchip');
    % Avoids negative values
    if Dfus_vee <0
        Dfus_vee = 0;
    end
    
    k       = k_calc(b_vee_s, Dfus_vee);

    % V tail
    %         Cy_beta_vee    = -k*(CLalpha_vee_e_pw*(tan(dihedral_vee)))*sidewash*S_vee/S_ref;
    Cl_beta_vee_VP = Cy_beta_vee*((z_zbar_vee - z_XCG)*cos(alpha)-(x_xbar_vee - x_XCG)*sin(alpha))/b_w1; % Vertical Projection
    Cl_beta_vee = Cl_beta_vee_VP + Cl_beta_vee_HP;
    % Single Vertical tail

else
    Cl_beta_vee_HP = 0;
    Cl_beta_vee_VP = 0;
    Cl_beta_vee_CL = 0;
    Cl_beta_vee_dihedral = 0;
    Cl_beta_vee_z = 0;
    Cl_beta_vee = 0;
end

if Vee2 == 1
     %% Contribution Horizontal projection
    CL_w = Trim_ITER.CL_vee2;
    AR_w    = Geo_tier.AR_vee2;
    b_w     = Geo_tier.b_vee2;
    Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, Geo_tier.x_xbar_vee2, 'pchip');
    % Positive bellow fuselage center line
    z_w = -Geo_tier.z_vee2_LE;
    TR_w    = Geo_tier.lambda_vee2;
    LAMc2_w = Geo_tier.Lambda_c2_vee2;
    LAMc4_w = Geo_tier.Lambda_c4_vee2;
    lt_w    = Geo_tier.x_xbar_vee2 - Geo_tier.xbar_vee2 + (Geo_tier.b_vee2/2)*tan(Geo_tier.Lambda_LE_vee2) + Geo_tier.cT_vee2/2;
    diedro  = Geo_tier.dihedral_vee2;
    lambda_w = Geo_tier.lambda_vee2;

    % Estimacion de Cl_beta_wf extraida de:
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
    %   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)
    DClbeta_zw          = 1.2*sqrt(AR_w)*z_w/b_w*2*Dfus_w/b_w; % 1/rad
    DClbeta_tau         = -180/pi*0.0005*sqrt(AR_w)*(Dfus_w/b_w)^2; % 1/rad
    Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(lambda_w, AR_w, LAMc2_w); % 1/rad

    K_MLAM              = K_MLAM_calc(Mach, AR_w, Lambda_c2_w1);

    Kf                  = Kf_calc(lt_w, b_w, AR_w, LAMc2_w); % Corrected with DATCOM
    % Colapses Kf in case that the results is out of the scope
    if isnan(Kf)
        Kf  = 1;
    end

    Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_w, lambda_w); % 1/rad
    Clbeta_die          = 180/pi*Clbeta_die_calc(lambda_w, AR_w, LAMc2_w); % 1/rad
    K_Mdie              = K_Mdie_calc(Mach, AR_w, LAMc2_w); % Corrected with DATCOM

    % C_L = W/qinf/Sref;
    %% Canard Contribution
    Cl_beta_vee2_CL = C_Lw*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR);
    Cl_beta_vee2_dihedral = diedro*(Clbeta_die*K_Mdie + DClbeta_tau);
    Cl_beta_vee2_z = DClbeta_zw;

    % Do not consider the dihedral contribution since that it is taken into
    % account wth CYbeta from FLOW or from NACA Report 823
    Cl_beta_vee2_dihedral = 0;
    Cl_beta_vee2_HP = Cl_beta_vee2_CL + Cl_beta_vee2_dihedral + Cl_beta_vee2_z; % Horizontal projection
     
    z_zbar_vee2 = Geo_tier.z_zbar_vee2;
    x_xbar_vee2 = Geo_tier.x_xbar_vee2;
    %% APORTE VERTICAL
    % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
    % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash = 0.724 + 3.06*(S_vee2_pv/S_ref)/(1 + cos(Lambda_c4_w1)) + 0.4*(z_w1_LE1)/Dfus_w1 + 0.009*AR_w1;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

    % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
    % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
    Dfus_vee2  = interp1(length_x_position,width_x_position, x_xbar_vee2, 'pchip');
    % Avoids negative values
    if Dfus_vee2 <0
        Dfus_vee2 = 0;
    end
    
    k       = k_calc(b_vee2_s, Dfus_vee2);

    % V tail
    %         Cy_beta_vee2    = -k*(CLalpha_vee2_e_pw*(tan(dihedral_vee2)))*sidewash*S_vee2/S_ref;
    Cl_beta_vee2_VP    = Cy_beta_vee2*((z_zbar_vee2 - z_XCG)*cos(alpha)-(x_xbar_vee2 - x_XCG)*sin(alpha))/b_w1; % Vertical Projection
    Cl_beta_vee2 = Cl_beta_vee2_VP + Cl_beta_vee2_HP;
    % Single Vertical tail

else
    Cl_beta_vee2_HP = 0;
    Cl_beta_vee2_VP = 0;
    Cl_beta_vee2_CL = 0;
    Cl_beta_vee2_dihedral = 0;
    Cl_beta_vee2_z = 0;
    Cl_beta_vee2 = 0;
end

if Nac == 1
    Cl_beta_nac = 0;
else
    Cl_beta_nac = 0;
end

%% DERIVADA TOTAL
Cl_beta         = Cl_beta_wf + Cl_beta_VTP + Cl_beta_vee + Cl_beta_vee2 + Cl_beta_can + Cl_beta_HTP;

Stab_Der.Clb_wf = Cl_beta_wf;
Stab_Der.Clb_VTP = Cl_beta_VTP;
Stab_Der.Clb_vee = Cl_beta_vee;
Stab_Der.Clb_vee2 = Cl_beta_vee2;
Stab_Der.Clb_can = Cl_beta_can;
Stab_Der.Clb_HTP = Cl_beta_HTP;
Stab_Der.Clb = Cl_beta;
% Stores Derivatives per parts
Stab_Der_parts.Cl_beta = Cl_beta;
% w1
Stab_Der_parts.Cl_beta_wf = Cl_beta_wf;
Stab_Der_parts.Cl_beta_wf_CL = Cl_beta_wf_CL;
Stab_Der_parts.Cl_beta_wf_dihedral = Cl_beta_wf_dihedral;
Stab_Der_parts.Cl_beta_wf_z = Cl_beta_wf_z;
% VTP
Stab_Der_parts.Cl_beta_VTP = Cl_beta_VTP;
% Vee
Stab_Der_parts.Cl_beta_vee = Cl_beta_vee;
Stab_Der_parts.Cl_beta_vee_VP = Cl_beta_vee_VP;
Stab_Der_parts.Cl_beta_vee_HP = Cl_beta_vee_HP;
Stab_Der_parts.Cl_beta_vee_CL = Cl_beta_vee_CL;
Stab_Der_parts.Cl_beta_vee_dihedral = Cl_beta_vee_dihedral;
Stab_Der_parts.Cl_beta_vee_z = Cl_beta_vee_z;
% Vee2
Stab_Der_parts.Cl_beta_vee2 = Cl_beta_vee2;
Stab_Der_parts.Cl_beta_vee2_VP = Cl_beta_vee2_VP;
Stab_Der_parts.Cl_beta_vee2_HP = Cl_beta_vee2_HP;
Stab_Der_parts.Cl_beta_vee2_CL = Cl_beta_vee2_CL;
Stab_Der_parts.Cl_beta_vee2_dihedral = Cl_beta_vee2_dihedral;
Stab_Der_parts.Cl_beta_vee2_z = Cl_beta_vee2_z;
% NAC
Stab_Der_parts.Cl_beta_nac = Cl_beta_nac;
% Can
Stab_Der_parts.Cl_beta_can = Cl_beta_can;
Stab_Der_parts.Cl_beta_can_CL = Cl_beta_can_CL;
Stab_Der_parts.Cl_beta_can_dihedral = Cl_beta_can_dihedral;
Stab_Der_parts.Cl_beta_can_z = Cl_beta_can_z;
% HTP
Stab_Der_parts.Cl_beta_HTP = Cl_beta_HTP;
Stab_Der_parts.Cl_beta_HTP_CL = Cl_beta_HTP_CL;
Stab_Der_parts.Cl_beta_HTP_dihedral = Cl_beta_HTP_dihedral;
Stab_Der_parts.Cl_beta_HTP_z = Cl_beta_HTP_z;

%% Cnb
length_x_position = Body_Geo.length_x_position;
height_x_position = Body_Geo.height_x_position;
width_x_position = Body_Geo.width_x_position;
alpha = TRIM_RESULTS.trim_alpha;

Sref    = modelo.general.Sref;
W       = modelo.general.mtow*modelo.general.w_w0*9.8065;
qinf    = modelo.general.qinf;
vinf    = modelo.general.Vinf;
rhoinf  = modelo.general.rhoinf;
Xcg     = modelo.general.Xcg;
Zcg     = modelo.general.Zcg;
C_L     = modelo.general.CL;
C_Lw   = modelo.general.CL_w;
% C_Lh   = modelo.general.CL_h;

x_w     = modelo.ala.Xca;
MAC_w   = modelo.ala.MAC;

AR_w    = modelo.ala.AR;
b_w     = modelo.ala.b;
LAMc4_w = modelo.ala.LAMc4;
% z_w     = modelo.ala.Zca;
% Positrive bellow fuselage center line
z_w1_LE1 = -Geo_tier.z_w1_LE;
x_w     = modelo.ala.Xca;
MAC_w   = modelo.ala.MAC;
Dfus_w  = interp1(length_x_position,height_x_position, x_w, 'pchip');

S_v     = modelo.vertical.S;
b_v     = modelo.vertical.b;
CLa_v   = modelo.vertical.CLa;
x_v     = modelo.vertical.Xca;
Dfus_v = interp1(length_x_position,width_x_position,x_v,'pchip');


l_fus   = modelo.fuselaje.l;
S_fusS  = modelo.fuselaje.Sside;
D_fus   = interp1(length_x_position,height_x_position, x_w, 'pchip');
W_fus   = modelo.fuselaje.W;
h1      = interp1(length_x_position,height_x_position, 0.25*l_fus, 'pchip');
h2      = interp1(length_x_position,height_x_position, 0.75*l_fus, 'pchip');

% ASPro Cnb
if W1 ==1
    %%  APORTE DEL ALA-FUSELAJE
    AR_w    = Geo_tier.AR_w1;
    b_w     = Geo_tier.b_w1;
    Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, Geo_tier.x_xbar_w1, 'pchip');
    % Positive bellow fuselage center line
    z_w = -Geo_tier.z_w1_LE;
    TR_w    = Geo_tier.lambda_w1;
    LAMc2_w = Geo_tier.Lambda_c2_w1;
    LAMc4_w = Geo_tier.Lambda_c4_w1;
    lt_w    = Geo_tier.x_xbar_w1 - Geo_tier.xbar_w1 + (Geo_tier.b_w1/2)*tan(Geo_tier.Lambda_LE_w1) + Geo_tier.cT_w1/2;
    diedro  = Geo_tier.dihedral_w1;
    lambda_w = Geo_tier.lambda_w1;

    Lambda_c2_w = Geo_tier.Lambda_c2_w1;
    Lambda_c4_w = Geo_tier.Lambda_c4_w1;
    dihedral_w = Geo_tier.dihedral_w1;
    AR_w = Geo_tier.AR_w1;
    CL_w = Trim_ITER.CL_w1;

    x_w     = Geo_tier.x_xbar_w1;
    MAC_w   = Geo_tier.cmac_w1;

    Cn_beta_w1 = CL_w^2 * (1/4/pi/AR_w - (tan(LAMc4_w)/pi/AR_w/(AR_w + 4*cos(LAMc4_w)))*(cos(LAMc4_w)...
        - AR_w/2 - AR_w^2/8/cos(LAMc4_w) + (6*(x_w - Xcg)*sin(LAMc4_w))/MAC_w/AR_w));

else
    Cn_beta_w1 = 0;
end

if Can ==1
    %%  APORTE DEL Canard
    AR_w    = Geo_tier.AR_can;
    b_w     = Geo_tier.b_can;
    Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, Geo_tier.x_xbar_can, 'pchip');
    % Positive bellow fuselage center line
    z_w = -Geo_tier.z_can_LE;
    TR_w    = Geo_tier.lambda_can;
    LAMc2_w = Geo_tier.Lambda_c2_can;
    LAMc4_w = Geo_tier.Lambda_c4_can;
    lt_w    = Geo_tier.x_xbar_can - Geo_tier.xbar_can + (Geo_tier.b_can/2)*tan(Geo_tier.Lambda_LE_can) + Geo_tier.cT_can/2;
    diedro  = Geo_tier.dihedral_can;
    lambda_w = Geo_tier.lambda_can;

    Lambda_c2_w = Geo_tier.Lambda_c2_can;
    Lambda_c4_w = Geo_tier.Lambda_c4_can;
    dihedral_w = Geo_tier.dihedral_can;
    AR_w = Geo_tier.AR_can;
    CL_w = Trim_ITER.CL_can;

    x_w     = Geo_tier.x_xbar_can;
    MAC_w   = Geo_tier.cmac_can;

    Cn_beta_can = CL_w^2 * (1/4/pi/AR_w - (tan(LAMc4_w)/pi/AR_w/(AR_w + 4*cos(LAMc4_w)))*(cos(LAMc4_w)...
        - AR_w/2 - AR_w^2/8/cos(LAMc4_w) + (6*(x_w - Xcg)*sin(LAMc4_w))/MAC_w/AR_w));
    
else
    Cn_beta_can = 0;
end

if HTP ==1
    %%  APORTE DEL Canard
    AR_w    = Geo_tier.AR_HTP;
    b_w     = Geo_tier.b_HTP;
    Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, Geo_tier.x_xbar_HTP, 'pchip');
    % Positive bellow fuselage center line
    z_w = -Geo_tier.z_HTP_LE;
    TR_w    = Geo_tier.lambda_HTP;
    LAMc2_w = Geo_tier.Lambda_c2_HTP;
    LAMc4_w = Geo_tier.Lambda_c4_HTP;
    lt_w    = Geo_tier.x_xbar_HTP - Geo_tier.xbar_HTP + (Geo_tier.b_HTP/2)*tan(Geo_tier.Lambda_LE_HTP) + Geo_tier.cT_HTP/2;
    diedro  = Geo_tier.dihedral_HTP;
    lambda_w = Geo_tier.lambda_HTP;

    Lambda_c2_w = Geo_tier.Lambda_c2_HTP;
    Lambda_c4_w = Geo_tier.Lambda_c4_HTP;
    dihedral_w = Geo_tier.dihedral_HTP;
    AR_w = Geo_tier.AR_HTP;
    CL_w = Trim_ITER.CL_HTP;

    x_w     = Geo_tier.x_xbar_HTP;
    MAC_w   = Geo_tier.cmac_HTP;

    Cn_beta_HTP = CL_w^2 * (1/4/pi/AR_w - (tan(LAMc4_w)/pi/AR_w/(AR_w + 4*cos(LAMc4_w)))*(cos(LAMc4_w)...
        - AR_w/2 - AR_w^2/8/cos(LAMc4_w) + (6*(x_w - Xcg)*sin(LAMc4_w))/MAC_w/AR_w));
    
else
    Cn_beta_HTP = 0;
end

%% APORTE FUSELAJE
% 2 FORMAS:
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 398, 2088 PDF)
%   PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu; 3.5 (pag 271)

mu_air      = 1.7e-5;  % Viscosidad dinámica del aire
Re          = rhoinf*vinf*l_fus/mu_air;

K_N     = K_N_calc(l_fus, S_fusS, Xcg, D_fus, h1, h2, W_fus);
K_Rl    = K_Rl_calc(Re);

Cn_beta_fus = -180/pi*K_N*K_Rl*S_fusS*l_fus/Sref/b_w;

%   - AIRCRAFT DESIGN: A CONCEPTUAL APPROACH 5th EDITION; Raymer, Daniel P.; 16.4.5 (pag 634)
% Cn_beta_fus = -1.3*modelo.fuselaje.vol/modelo.general.Sref/modelo.ala.b*modelo.fuselaje.D/modelo.fuselaje.W;
% D_fus: Fuselage depth
% W_fus: Fuselage width
% V_fus: Fuselage volume

%% DERIVADA DE ESTABILIDAD COMPLETA
if VTP == 1

    %% APORTE VERTICAL
    %% APORTE VERTICAL
    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w1_LE1/D_fus + 0.009*AR_w;
    %% REVISE I think the correct is Dfus_w
    % sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w/Dfus_w + 0.009*AR_w;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))


    % Calculo de la pendiente de sustentacion del vertical

    % X_le        = modelo.vertical.Xac - modelo.vertical.xac + intep1(modelo.vertical.y, modelo.vertical.le_y, modelo.horizontal.Zac,'pchip');
    % dx_vle2hac  = horizontal.Xac - X_le;
    % dz_fus2hac  = horizontal.Zac;
    %
    % ARv_eff = getARv_eff(S_v, modelo.vertical.b, modelo.vertical.TR,...
    % modelo.horizontal.S, dx_vle2hac, dz_fus2hac, Dfus_v, modelo.general.vert_conf);
    %
    % CLa_v = getCLa_geometryBased(ARv_eff, modelo.vertical.LAMc2, modelo.general.Minf, modelo.vertical.Cla);

    k       = k_calc(b_v, Dfus_v);
    Cn_beta_VTP    = - Cy_beta_VTP*((z_zbar_VTP - Zcg)*sin(alpha)+(x_v - Xcg)*cos(alpha))/b_w;

    if isnan(Cn_beta_VTP)
        Cn_beta_VTP = 0;
    end

else
    Cn_beta_VTP = 0;
end

%% DERIVADA DE ESTABILIDAD COMPLETA
if Vee == 1

    %% APORTE VERTICAL
    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w1_LE1/D_fus + 0.009*AR_w;
    %% REVISE I think the correct is Dfus_w
    % sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w/Dfus_w + 0.009*AR_w;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

    % Calculo de la pendiente de sustentacion del vertical

    % X_le        = modelo.vertical.Xac - modelo.vertical.xac + intep1(modelo.vertical.y, modelo.vertical.le_y, modelo.horizontal.Zac,'pchip');
    % dx_vle2hac  = horizontal.Xac - X_le;
    % dz_fus2hac  = horizontal.Zac;
    %
    % ARv_eff = getARv_eff(S_v, modelo.vertical.b, modelo.vertical.TR,...
    % modelo.horizontal.S, dx_vle2hac, dz_fus2hac, Dfus_v, modelo.general.vert_conf);
    %
    % CLa_v = getCLa_geometryBased(ARv_eff, modelo.vertical.LAMc2, modelo.general.Minf, modelo.vertical.Cla);

    k       = k_calc(b_v, Dfus_v);
    Cn_beta_vee    = - Cy_beta_vee*((z_zbar_vee - Zcg)*sin(alpha)+(x_v - Xcg)*cos(alpha))/b_w;

    if isnan(Cn_beta_vee)
        Cn_beta_vee = 0;
    end

else
    Cn_beta_vee = 0;
end

%% DERIVADA DE ESTABILIDAD COMPLETA
if Vee2 == 1

    %% APORTE VERTICAL
    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w1_LE1/D_fus + 0.009*AR_w;
    %% REVISE I think the correct is Dfus_w
    % sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w/Dfus_w + 0.009*AR_w;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

    % Calculo de la pendiente de sustentacion del vertical

    % X_le        = modelo.vertical.Xac - modelo.vertical.xac + intep1(modelo.vertical.y, modelo.vertical.le_y, modelo.horizontal.Zac,'pchip');
    % dx_vle2hac  = horizontal.Xac - X_le;
    % dz_fus2hac  = horizontal.Zac;
    %
    % ARv_eff = getARv_eff(S_v, modelo.vertical.b, modelo.vertical.TR,...
    % modelo.horizontal.S, dx_vle2hac, dz_fus2hac, Dfus_v, modelo.general.vert_conf);
    %
    % CLa_v = getCLa_geometryBased(ARv_eff, modelo.vertical.LAMc2, modelo.general.Minf, modelo.vertical.Cla);

    k       = k_calc(b_v, Dfus_v);
    Cn_beta_vee2    = - Cy_beta_vee2*((z_zbar_vee2 - Zcg)*sin(alpha)+(x_v - Xcg)*cos(alpha))/b_w;

    if isnan(Cn_beta_vee2)
        Cn_beta_vee2 = 0;
    end

else
    Cn_beta_vee2 = 0;
end

if Nac == 1
    Cn_beta_nac = 0;
else
    Cn_beta_nac = 0;
end
Cn_beta = Cn_beta_w1 + Cn_beta_fus + Cn_beta_VTP + Cn_beta_vee + Cn_beta_vee2 + Cn_beta_HTP + Cn_beta_can;

% Stab_Der = getCnbeta(modelo,Stab_Der,Body_Geo,TRIM_RESULTS,Cy_beta_VTP,z_zbar_VTP,Geo_tier);

Stab_Der.Cnb_w = Cn_beta_w1;
Stab_Der.Cnb_b = Cn_beta_fus;
Stab_Der.Cnb_wb = Cn_beta_w1 + Cn_beta_fus;    
Stab_Der.Cnb_VTP = Cn_beta_VTP;    
Stab_Der.Cnb_vee = Cn_beta_vee;
Stab_Der.Cnb_vee2 = Cn_beta_vee2;
Stab_Der.Cnb_can = Cn_beta_can;
Stab_Der.Cnb_HTP = Cn_beta_HTP;
Stab_Der.Cnb = Cn_beta;

% Stores Derivatives per parts
Stab_Der_parts.Cn_beta = Cn_beta;
Stab_Der_parts.Cn_beta_w1 = Cn_beta_w1;
Stab_Der_parts.Cn_beta_fus = Cn_beta_fus;
Stab_Der_parts.Cn_beta_VTP = Cn_beta_VTP;
Stab_Der_parts.Cn_beta_vee = Cn_beta_vee;
Stab_Der_parts.Cn_beta_vee2 = Cn_beta_vee2;
Stab_Der_parts.Cn_beta_can = Cn_beta_can;
Stab_Der_parts.Cn_beta_HTP = Cn_beta_HTP;
Stab_Der_parts.Cn_beta_nac = Cn_beta_nac;

% Cnb_b = Stab_Der.Cnb_b;
% Cnb_w = Stab_Der.Cnb_w;
% Cnb_wb = Stab_Der.Cnb_wb;
% Cnb_v = Stab_Der.Cnb_v;
% Cnb = Stab_Der.Cnb;

% %Cnbdiedro
% % Strip theory approximation
% Cnb_diedro = -Sigma_w*D2R*(CL-CD_alpha)/4; % Rectangular wings with constant Chord - Pamadi Eq 3.270
% % approximation that solves the problem of strip theory of ignoring the
% % induced drag effects
% Cnb_diedro = -0.075*Sigma_w*D2R*CL; % Rectangular wings with constant Chord - Pamadi Eq 3.271
% Cnb_flecha = CL^2/(4*pi*AR_w); % Appproximation Pamadi 3.299 with all terms with Sweep ignored
% % Input data to the subfigures in Fig 3.73 Pamadi
% x_m = (x_XCG + mamparo_F_front)/length_fus;
% inFig_3_73 = (length_fus^2)/Area_side;
% h_1 = z_Area_b_1_4;
% h_2 = z_Area_b_3_4;
% in_h1_h2 = sqrt(h_1/h_2);
% in_h_bf_max = h_Area_b_max/w_Area_b_max;
% Knn=0.00125;                       %de grafica de altura/bfmax
% % Figura 3.74 - Pamadi
% Re_fus = rho*V*length_fus/(mu);
% Re_fus_f_Re_fus = [1e6,7e6, 30e6, 50e6, 80e6 350e6];
% K_Ri_f_Re_fus = [1,1.4, 1.7 1.8, 1.9,2.2];
% Kri_f_Re_fus  = interp1(Re_fus_f_Re_fus,K_Ri_f_Re_fus,Re_fus,'spline');
% Kri=Kri_f_Re_fus;                         %sale del numero de reynolds
% Cnb_b = -Knn*Kri*(Area_side/S_w)*(length_fus/b_w)*R2D;
% % Cnbw_b = -1.3*(Vol_TOT/(S_w*b_w))*(h_Area_b_max/w_Area_b_max);
% Cnb_w = Cnb_diedro + Cnb_flecha;
% Cnb_wb = Cnb_diedro + Cnb_flecha + Cnb_b;
% Cnb_v = - Cyb_v*((x_v_xbar - x_XCG)*cos(trim_alpha) + (z_v_xbar - Zcg)*sin(trim_alpha))/b_w;
% Cnb =  Cnb_wb + Cnb_v;

%%%%%%%%%%%%%%%derivadas propulsivas laterales-direccionales%%%%%%%%%%%%%%%
% % Data from NACA REport 640
% w_R_30 = 0.0525*2;
% w_R_60 = 0.073*2;
% w_R_90 = 0.045*2;
% Nprop = 1;
% CNalpha_p_KN = 0.1;
% KN = 262*w_R_30 + 262*w_R_60 + 135*w_R_90;
% dCN_dalpha = CNalpha_p_KN*(1 + 0.8*(KN/80.7) - 1);
% CyTb = (pi/4)*Nprop*(D_prop^2)*dCN_dalpha/S_w1
% % Thrust lines inclination angle
% psi_T = 0;
% l_prop = (x_XCG-x_m_propeller)*cos(psi_T);
% CNTb = (pi/4)*Nprop*l_prop*(D_prop^2)*dCN_dalpha/(S_w1*b_w1)
% Stab_Der.CyTb = CyTb;
% Stab_Der.CNTb = CNTb;
% pause

%% REVIEW
%ASpro CyTb
Stab_Der = getCyTbeta(modelo,Stab_Der);
CyTb = Stab_Der.CyTb;

%% REVIEW
%Aspro CNTb
Stab_Der = getCnTbeta(modelo,Stab_Der);
CNTb = Stab_Der.CNTb;

Stab_Der_parts.sidewash = sidewash;
Stab_Der.sidewash = sidewash;
end