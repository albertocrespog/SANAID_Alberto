function [Stab_Der_parts Stab_Der] = get_tab_derivatives(modelo,OUTPUT_read_XLSX, AC_CONFIGURATION,Geo_tier,afe,alpha,Aero,Stab_Der,Stab_Der_parts)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

% Available control surfaces
d_TT_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_TT_ail; %Definition of available control surface - aileron trim tab
d_TT_ele = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_TT_ele; %Definition of available control surface - elevator trim tab
d_TT_rudder = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_TT_rudder; %Definition of available control surface - rudder  trim tab
d_TT_sp = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_TT_sp; %Definition of available control surface - Spoiler

if d_TT_ail == 1
    % Data
    Minf    = modelo.general.Minf;

    W       = modelo.general.mtow*modelo.general.w_w0*9.8065;
    qinf    = modelo.general.qinf;
    Sref    = modelo.general.Sref;

    yf_b2   = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_TT_ail_w1;
    y0_b2   = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_TT_ail_w1;
    TR      = modelo.ala.TR;
    LAMc4   = modelo.ala.LAMc4;
    % Relation between aileron and wing chord
    ca_c_tmp    = modelo.ala.cm_c;
    % relation between aileron trim tab and aileron
    ca_TT   = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_TT_ail;
    % relation between aileron trim tab and wing chord
    ca_c = ca_c_tmp*ca_TT;
    t_c     = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_TT_ail;
    AR      = modelo.ala.AR;
    CL      = W/qinf/Sref;
    %% Aileron Tab Control Derivatives
    % Sideforce Coefficient due to Aileron Tab Deflection Derivative
    % The airplane sideforce-coefficient-due-to-aileron-tab-deflection derivative is negligible for most conventional aileron arragements.

    %% CALCULO DE LA DERIVADA DE CONTROL C_y_da
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 442, 2137 PDF)
    Cy_delta_a_Tab = 0;

    %% CALCULO DE LA DERIVADA DE CONTROL C_l_da
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 442, 2137 PDF)

    beta = sqrt(1-Minf^2);

    betaCldelta_kf  = betaCldelta_k_calc(yf_b2,TR,AR,LAMc4,Minf);
    betaCldelta_k0  = betaCldelta_k_calc(y0_b2,TR,AR,LAMc4,Minf);

    %% Perfil revisado para el Ala
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_w1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
    clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
    clalpha_w2 = clalpha_w2_M0/beta;
    Cla_theo_w2 = 2*pi + 5.0525*t_c;
    Cla_M0 = clalpha_w2_M0/Cla_theo_w2;

    k = clalpha_w2_M0/(2*pi);

    Cldelta_prima   = k/beta*(betaCldelta_kf - betaCldelta_k0);

    Cldelta_int     = Cldelta_Cldeltatheory_calc(ca_c,Cla_M0)*Cldeltatheory_calc(ca_c,t_c);

    alpha_delta     = Cldelta_int/clalpha_w2;

    Cl_da = Cldelta_prima*alpha_delta;

    Cl_delta_a_Tab = Cl_da;

    %% CALCULO DE Cn_da_TT
    %Yawing Moment Coefficient due to Aileron Tab Deflection Derivative
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 448, 2138 PDF)
    Ka      = Ka_calc(yf_b2, AR, TR);

    Cn_da   = 2*Ka*CL*Cl_da;
    Cn_delta_a_Tab = Cn_da;

else
    Cy_delta_a_Tab = 0;
    Cl_delta_a_Tab = 0;
    Cn_delta_a_Tab = 0;
end

%% Elevator Tab Control Derivatives

if d_TT_ele == 1
    CD_delta_e_Tab = 0;
    CL_delta_e_Tab = 0;
    CM_delta_e_Tab = 0;
else
    CD_delta_e_Tab = 0;
    CL_delta_e_Tab = 0;
    CM_delta_e_Tab = 0;
end

%% Rudder Tab Control Derivatives
if d_TT_rudder == 1

    % Data
    y0_b2   = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_TT_rudder_VTP;
    y1_b2   = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_TT_rudder_VTP;
    AR      = modelo.vertical.AR;
    TR      = modelo.vertical.TR;
    % Relation between rudder and VTP Chord chord
    cf_c_tmp    = modelo.vertical.cm_c;
    % relation between aileron trim tab and VTP
    crudder_TT   = OUTPUT_read_XLSX.InputGeometry_Data_flags.cf_TT_rudder;
    % relation between aileron trim tab and wing chord
    cf_c = cf_c_tmp*crudder_TT;

    t_c     = OUTPUT_read_XLSX.InputGeometry_Data_flags.t_c_TT_rudder;
    CLa     = Aero.CL_alpha_VTP_CR; %el de aero, sin multiplicar por eta*S/Sref
    S_v     = modelo.vertical.S;
    Xca_v   = Geo_tier.x_xbar_VTP;
    Zca_v   = Geo_tier.z_zbar_VTP;
    b       = Geo_tier.b_w1;

    prop_wash_effect = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;

    Sref    = modelo.general.Sref;
    % %Added
    if prop_wash_effect == 1
        rend    = afe.eta_w2_no_afe;
    else
        rend=1;
    end

    x_XCG = modelo.general.Xcg;
    z_XCG = modelo.general.Zcg;

    Mach = modelo.general.Minf;

    %% CÁLCULO DE LA DERIVADA DE ESTABILIDAD
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)
    %% CÁLCULO DE LA DERIVADA DE ESTABILIDAD
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)
    k_prima                 = 1;

    %% Perfil revisado para el VTP
    airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_VTP1);
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
    clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
    clalpha_w2 = clalpha_w2_M0/sqrt(1-Mach^2);
    Cla_theo_w2 = 2*pi + 5.0525*t_c;
    Cla_M0 = clalpha_w2_M0/Cla_theo_w2;
    % The three dimensional ruddervator effectiveness parameter is determined from Figure 8.53 in Airplane Design Part VI
    % and is a function of the V-tail aspect ratio, and ruddervator chord to V-tail chord ratio:
    alphaCL_alphaCl         = alphaCL_alphaCl_calc(cf_c,AR);
    % Calculo de Cldelta_Cldeltatheory: Correction Factor for Plain Flap Lift. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.15 (pag 230, 1920 PDF))
    Cldelta_Cldeltatheory   = Cldelta_Cldeltatheory_calc(cf_c,Cla_M0);
    % Calculo de Cldeltatheory: Lift Effectiveness. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.14 (pag 228, 1918 PDF)
    Cldelta_theory          = Cldeltatheory_calc(cf_c,t_c);
    % Calculo de Kb: Effect of the taper ratio and flap span. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.4.2, Fig. 8.52 (pag 260, 1950 PDF))
    K_b                     = Kb_calc(y0_b2, y1_b2, TR);

    % The change in sideslip due to ruddervator directional deflection
    beta_delta_r = -K_b*Cldelta_Cldeltatheory*Cldelta_theory*(k_prima/clalpha_w2)*alphaCL_alphaCl;


    Cy_dr   = -rend*CLa*(S_v/Sref)*beta_delta_r;
    if AC_CONFIGURATION.twin_VTP == 1
        Cy_dr = 2*Cy_dr;
    end
    Cy_delta_r_Tab = Cy_dr;

    % Cy_dr   = K_b*CLa*S_v/Sref*rend*k_prima/Cla*alphaCL_alphaCl*Cldelta_Cldeltatheory*Cldelta_theory

    if isnan(Cy_dr)
        Cy_dr = 0;
        Cy_delta_r_Tab = Cy_dr;
    end

    %% Cálculo preliminar de la derivada de estabilidad C_l_dr.
    Cl_dr   = ((Zca_v - z_XCG)*cos(alpha)-(Xca_v - x_XCG)*sin(alpha))/b*Cy_dr;
    Cl_delta_r_Tab = Cl_dr;

    %% CÁLCULO DE LA DERIVADA DE CONTROL
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 462, 2152 PDF)

    Cn_dr   = -((Zca_v - z_XCG)*sin(alpha)+ (Xca_v - x_XCG)*cos(alpha))/b*Cy_dr;
    Cn_delta_r_Tab = Cn_dr;

else
    Cy_delta_r_Tab = 0;
    Cl_delta_r_Tab = 0;
    Cn_delta_r_Tab = 0;
end

%% Spoiler Tab Control Derivatives
if d_TT_sp == 1
    CD_delta_sp = 0;
    CL_delta_sp = 0;
    CM_delta_sp = 0;
else
    CD_delta_sp = 0;
    CL_delta_sp = 0;
    CM_delta_sp = 0;
end

%% Storing
% Elevator
Stab_Der.CL_delta_e_Tab = CL_delta_e_Tab;
Stab_Der.CD_delta_e_Tab = CD_delta_e_Tab;
Stab_Der.CM_delta_e_Tab = CM_delta_e_Tab;

% Aileron
% Stab_Der.Ch_deltaa_Tab = Ch_deltaa_Tab;
Stab_Der.Cy_delta_a_Tab = Cy_delta_a_Tab;
Stab_Der.Cl_delta_a_Tab = Cl_delta_a_Tab;
Stab_Der.Cn_delta_a_Tab = Cn_delta_a_Tab;

% Rudder
% Stab_Der.Ch_deltaa_Tab = Ch_deltaa_Tab;
Stab_Der.Cy_delta_r_Tab = Cy_delta_r_Tab;
Stab_Der.Cl_delta_r_Tab = Cl_delta_r_Tab;
Stab_Der.Cn_delta_r_Tab = Cn_delta_r_Tab;

% Spoiler
Stab_Der.CL_delta_sp = CL_delta_sp;
Stab_Der.CD_delta_sp = CD_delta_sp;
Stab_Der.CM_delta_sp = CM_delta_sp;

% Elevator
Stab_Der_parts.CL_delta_e_Tab = CL_delta_e_Tab;
Stab_Der_parts.CD_delta_e_Tab = CD_delta_e_Tab;
Stab_Der_parts.CM_delta_e_Tab = CM_delta_e_Tab;

% Aileron
% Stab_Der.Ch_deltaa_Tab = Ch_deltaa_Tab;
Stab_Der_parts.Cy_delta_a_Tab = Cy_delta_a_Tab;
Stab_Der_parts.Cl_delta_a_Tab = Cl_delta_a_Tab;
Stab_Der_parts.Cn_delta_a_Tab = Cn_delta_a_Tab;

% Rudder
% Stab_Der.Ch_deltaa_Tab = Ch_deltaa_Tab;
Stab_Der_parts.Cy_delta_r_Tab = Cy_delta_r_Tab;
Stab_Der_parts.Cl_delta_r_Tab = Cl_delta_r_Tab;
Stab_Der_parts.Cn_delta_r_Tab = Cn_delta_r_Tab;

% Spoiler
Stab_Der_parts.CL_delta_sp = CL_delta_sp;
Stab_Der_parts.CD_delta_sp = CD_delta_sp;
Stab_Der_parts.CM_delta_sp = CM_delta_sp;
