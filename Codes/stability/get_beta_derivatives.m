function [Stab_Der_parts,Stab_Der] = get_beta_derivatives(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,conditions)
    W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
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
%% Cy_beta = Cy_beta_fus + Cy_beta_wing + Cy_beta_vert

if W1 ==1
    %% APORTE ALA
    % - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
    Cy_beta_w1       = -0.00573*abs(dihedral_w1*180/pi);

    % - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
    C_L = W/qinf/S_ref;
    Cy_beta_w1       = (CL_w1^2)*(6*tan(Lambda_c2_w1)*sin(Lambda_c2_w1))/(pi*AR_w1*(AR_w1 + 4*cos(Lambda_c2_w1)));
else
    Cy_beta_w1 = 0;
end
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
    CLalpha_VTP =  Aero.CL_alpha_VTP_CR;
    b_VTP_s = Geo_tier.b_VTP_s;
    b_VTP = Geo_tier.b_VTP;
    b_w2 = Geo_tier.b_w2;
    xbar_VTP = Geo_tier.xbar_VTP;
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
    Cy_beta_VTP    = -k*CLalpha_VTP*sidewash*S_VTP/S_ref;

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

if Vee == 1
    S_w2_pv = Geo_tier.S_w2_pv;
    S_w2 = Geo_tier.S_w2;
    S_w2_s = Geo_tier.S_w2_s;
    x_xbar_w2 = Geo_tier.x_xbar_w2;
    b_w2_s = Geo_tier.b_w2_s;
    dihedral_w2 = Geo_tier.dihedral_w2;
    CYbeta_vee = Aero.CYbeta_vee;
    %% APORTE VERTICAL
    % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
    % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash = 0.724 + 3.06*(S_w2_pv/S_ref)/(1 + cos(Lambda_c4_w1)) + 0.4*(z_w1_LE1)/Dfus_w1 + 0.009*AR_w1;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))
    % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
    % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
    Dfus_Vee  = interp1(length_x_position,height_x_position, x_xbar_w2, 'pchip');
    Dfus_Vee  = interp1(length_x_position,width_x_position, x_xbar_w2, 'pchip');
    % Avoids negative relation
    if Dfus_Vee <0
        Dfus_Vee = 0;
    end

    k       = k_calc(b_w2_s, Dfus_Vee);
    % Single Vertical tail
%     Cy_beta_Vee2    = -k*abs(CYbeta_vee)*sidewash*S_w2/S_ref;
    Cy_beta_Vee    = -k*abs(CYbeta_vee)*sidewash*S_w2_s/S_ref;

%     Cy_beta_Vee = -0.483133428328602;
else
    Cy_beta_Vee = 0;
end



%% DERIVADA TOTAL
Cy_beta         = Cy_beta_w1 + Cy_beta_fus + Cy_beta_VTP + Cy_beta_Vee;
Stab_Der.Cyb_w = Cy_beta_w1;
Stab_Der.Cyb_fus = Cy_beta_fus;
Stab_Der.Cyb_VTP = Cy_beta_VTP;
Stab_Der.Cyb_Vee = Cy_beta_Vee;
Stab_Der.Cyb = Cy_beta;
% Stores Derivatives per parts
Stab_Der_parts.Cy_beta = Cy_beta;
Stab_Der_parts.Cy_beta_w1 = Cy_beta_w1;
Stab_Der_parts.Cy_beta_fus = Cy_beta_fus;
%     Stab_Der_parts.Cy_beta_HTP = Cy_beta_HTP;
Stab_Der_parts.Cy_beta_Vee = Cy_beta_Vee;
%     Stab_Der_parts.Cy_beta_can = Cy_beta_can;
%     Stab_Der_parts.Cy_beta_nac = Cy_beta_nac;
if VTP == 1
    Cy_beta_vert = Cy_beta_VTP;
    z_zbar_vert = Geo_tier.z_zbar_VTP;
    Stab_Der_parts.Cy_beta_vert = Cy_beta_vert;
elseif Vee == 1
    Cy_beta_vert = Cy_beta_Vee;
    z_zbar_vert = Geo_tier.z_zbar_w2;
    Stab_Der_parts.Cy_beta_vert = Cy_beta_vert;
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
% Positrive bellow fuselage center line
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
    %%  APORTE DEL ALA-FUSELAJE
    % Estimacion de Cl_beta_wf extraida de:
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
    %   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)
    DClbeta_zw          = 1.2*sqrt(AR_w)*z_w1_LE1/b_w*2*Dfus_w1/b_w; % 1/rad
    DClbeta_tau         = -180/pi*0.0005*sqrt(AR_w1)*(Dfus_w1/b_w1)^2; % 1/rad
    Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(lambda_w1, AR_w1, Lambda_c2_w1); % 1/rad
    
    K_MLAM              = K_MLAM_calc(Mach, AR_w1, Lambda_c2_w1);

    Kf                  = Kf_calc(lt_w, b_w1, AR_w1, Lambda_c2_w1); % Corrected with DATCOM
    Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_w1, lambda_w1); % 1/rad
    Clbeta_die          = 180/pi*Clbeta_die_calc(lambda_w1, AR_w1, Lambda_c2_w1); % 1/rad
    K_Mdie              = K_Mdie_calc(Mach, AR_w, Lambda_c2_w1); % Corrected with DATCOM

    % C_L = W/qinf/Sref;
    %% APORTE ALA
    Cl_beta_wf = C_Lw*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR) + diedro*(Clbeta_die*K_Mdie + DClbeta_tau) + DClbeta_zw;
    % die: Angulo de diedro en radianes
else
    Cl_beta_wf = 0;
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
    z_zbar_w2 = Geo_tier.z_zbar_w2;
    x_xbar_w2 = Geo_tier.x_xbar_w2;
    %% APORTE VERTICAL
    % AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
    % FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

    % Calculo del side-wash: deflexion de la corriente debida a la presencia
    % del ala
    sidewash = 0.724 + 3.06*(S_w2_pv/S_ref)/(1 + cos(Lambda_c4_w1)) + 0.4*(z_w1_LE1)/Dfus_w1 + 0.009*AR_w1;
    % S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
    % LAMc4: flecha del ala en c/4
    % z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

    % The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
    % and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
    Dfus_Vee  = interp1(length_x_position,width_x_position, x_xbar_w2, 'pchip');
    % Avoids negative values
    if Dfus_Vee <0
        Dfus_Vee = 0;
    end
    
    k       = k_calc(b_w2_s, Dfus_Vee);

    % V tail
    %         Cy_beta_Vee    = -k*(CLalpha_Vee_e_pw*(tan(dihedral_w2)))*sidewash*S_w2/S_ref;
    Cl_beta_Vee    = Cy_beta_Vee*((z_zbar_w2 - z_XCG)*cos(alpha)-(x_xbar_w2 - x_XCG)*sin(alpha))/b_w1;
    % Single Vertical tail

else
    Cl_beta_Vee = 0;
end

%% DERIVADA TOTAL
Cl_beta         = Cl_beta_wf  + Cl_beta_VTP + Cl_beta_Vee;

Stab_Der.Clb_wf = Cl_beta_wf;
Stab_Der.Clb_VTP = Cl_beta_VTP;
Stab_Der.Clb_Vee = Cl_beta_Vee;
Stab_Der.Clb = Cl_beta;
% Stores Derivatives per parts
Stab_Der_parts.Clb = Cl_beta;
Stab_Der_parts.Clb_wf = Cl_beta_wf;
Stab_Der_parts.Clb_VTP = Cl_beta_VTP;
Stab_Der_parts.Clb_Vee = Cl_beta_Vee;


%% Cnb
% ASPro Cnb
Stab_Der = getCnbeta(modelo,Stab_Der,Body_Geo,TRIM_RESULTS,Cy_beta_vert,z_zbar_vert,Geo_tier);
Cnb_b = Stab_Der.Cnb_b;
Cnb_w = Stab_Der.Cnb_w;
Cnb_wb = Stab_Der.Cnb_wb;
Cnb_v = Stab_Der.Cnb_v;
Cnb = Stab_Der.Cnb;

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