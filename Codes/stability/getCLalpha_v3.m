function [Stab_Der_parts,  Trim_ITER] = getCLalpha_v3(AC_CONFIGURATION,Body_Geo,Geo_tier,body_interference,Stab_Der_parts, Effects,Design_criteria, OUTPUT_read_XLSX)

%% Input
W1      = AC_CONFIGURATION.W1;
HTP     = AC_CONFIGURATION.HTP;
VTP     = AC_CONFIGURATION.VTP;
Can     = AC_CONFIGURATION.Can;
Vee     = AC_CONFIGURATION.Vee;
Vee2     = AC_CONFIGURATION.Vee2;
Nac     = AC_CONFIGURATION.Nac;
AC_type = AC_CONFIGURATION.AC_type;

length_fus         = Body_Geo.l_fus;
w_Area_b_max       = Body_Geo.w_Area_b_max;
x_Area_b_max       = Body_Geo.x_Area_b_max;
Area_b_max         = Body_Geo.Area_b_max;
S_ref              = Geo_tier.S_ref;

CLalpha_w1_e_pw    = Stab_Der_parts.CLalpha_w1_e_pw;
CL0_w1_e_corrected = Stab_Der_parts.CL0_w1_e_corrected;

S_w1   = Geo_tier.S_w1;
S_w1_e = Geo_tier.S_w1_e;
b_w1   = Geo_tier.b_w1;
i_w1   = Design_criteria.i_w1;

wingspan2bodydiam = b_w1/w_Area_b_max;
flagwingspan2bodydiam = OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam; 
Munk_fuselage_constribution = OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution;

%% Contribution of Fuselage
if   wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1 || Munk_fuselage_constribution == 1
    CLa_fus            = Body_Geo.CLa_fus;
else
    CLa_fus            = 0;
end

if Can == 1
    upwash              = Effects.upwash;
    CL0_can_e_corrected = Stab_Der_parts.CL0_can_e_corrected;
    CLalpha_can_e_pw    = Stab_Der_parts.CLalpha_can_e_pw;
    S_can               = Geo_tier.S_can;
    S_can_e             = Geo_tier.S_can_e; % effective wing area: proyected in the y-plane
    b_can               = Geo_tier.b_can;
    i_can               = Design_criteria.i_can;
    eps_can             = Effects.eps_can;
    % When CAnard  in configuration, calculates the downwash in wing due to
    % canard
    if W1 == 1
        downwash_w1          = Effects.downwash_w1;
        CL0_w1_e_corrected = Stab_Der_parts.CL0_w1_e_corrected;
        CLalpha_w1_e_pw    = Stab_Der_parts.CLalpha_w1_e_pw;
        S_w1                = Geo_tier.S_w1;
        S_w1_s              = Geo_tier.S_w1_s; % real wing area: proyected along the surface of wing
        b_w1                = Geo_tier.b_w1;
        i_w1                = Design_criteria.i_w1;
        eps_w1              = Effects.eps_w1;
    end
end

if HTP == 1
    downwash_HTP          = Effects.downwash_HTP;
    CL0_HTP_e_corrected = Stab_Der_parts.CL0_HTP_e_corrected;
    CLalpha_HTP_e_pw    = Stab_Der_parts.CLalpha_HTP_e_pw;
    S_HTP                = Geo_tier.S_HTP;
    S_HTP_s              = Geo_tier.S_HTP_s; % real wing area: proyected along the surface of wing
    b_HTP                = Geo_tier.b_HTP;
    i_HTP                = Design_criteria.i_HTP;
    eps_HTP              = Effects.eps_HTP;
end

if Vee == 1
    downwash_vee            = Effects.downwash_vee;
    CL0_vee_e_corrected = Stab_Der_parts.CL0_vee_e_corrected;
    CLalpha_vee_e_pw    = Stab_Der_parts.CLalpha_vee_e_pw;
    S_vee                = Geo_tier.S_vee;
    S_vee_s              = Geo_tier.S_vee_s; % real wing area: proyected along the surface of wing
    b_vee                = Geo_tier.b_vee;
    i_vee                = Design_criteria.i_vee;
    eps_vee              = Effects.eps_vee;
end

if Vee2 == 1
    downwash_vee2            = Effects.downwash_vee2;
    CL0_vee2_e_corrected = Stab_Der_parts.CL0_vee2_e_corrected;
    CLalpha_vee2_e_pw    = Stab_Der_parts.CLalpha_vee2_e_pw;
    S_vee2                = Geo_tier.S_vee2;
    S_vee2_s              = Geo_tier.S_vee2_s; % real wing area: proyected along the surface of wing
    b_vee2                = Geo_tier.b_vee2;
    i_vee2                = Design_criteria.i_vee2;
    eps_vee2              = Effects.eps_vee2;
end


%% Correction of Lift curve slope of w1 associated to fuselage interaction
% Pendiente de sustentación del morro. No se usa, porque CLa_fus se coge
% directamente de XFLR5.

Fineness_Ratio = length_fus/w_Area_b_max;
%digitaliazacion figura PAMADI CAP3 Fig 3.6
x_f_k2_k1 = [4.,5.,6.,8.,10.,12.,14.,16.,18.,20.];
y_f_k2_k1 = [.77,.825,.865,.91,.94,.955,.965,.97,.973,.975];
f_k2_k1  = interp1(x_f_k2_k1,y_f_k2_k1,Fineness_Ratio,'spline');
x1 = x_Area_b_max;
Stab_Der_parts.f_k2_k1 = f_k2_k1;

if W1 == 1
    %% W1 body interference
    aN_w1 = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    %     aN_w1 = CLa_fus; %  corrigiendo con la integración del fuselaje a partir de modelo XFLR5

    %     f_k2_k1         = k2_k1_calc(Body_Geo.l_fus, Geo_tier.d_fus);      % Fuselage apparent mass coefficient. Pamadi, Figure 3.6
    %     Body_Geo.CLa_fus    = 2*º*(Body_Geo.S_front/Sref);         % Pendiente de sustentación del morro aislado


    KN_w1 = (aN_w1/CLalpha_w1_e_pw)*(S_w1/S_w1_e);  %eq 3.25 PAMADI
    KWB_w1 = 0.1714*(w_Area_b_max/b_w1)^2 + 0.8326*(w_Area_b_max/b_w1) + 0.9974; %eq 3.27 PAMADI
    KBW_w1 = 0.7810*(w_Area_b_max/b_w1)^2 + 1.1976*(w_Area_b_max/b_w1) + 0.0088; %eq 3.28 PAMADI

    % %     C*********       KW(B) VS D/B  DATA FIGURE 4.3.1.2-10-A
    %
    %        X10A=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %        Y10A=[1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0];
    %        KWB_w1 = interp1(X10A, Y10A, w_Area_b_max/b_w1);
    %
    % % C ****         KB(W) VS D/B  DATA  FIGURE 4.3.1.2-10-B
    %
    %       X10B=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %       Y10B=[0.0,0.13,0.29,0.45,0.62,0.80,1.0,1.22,1.45,1.70,2.0];
    %         KBW_w1 = interp1(X10B, Y10B, w_Area_b_max/b_w1);

    % Queda preguntar si se añaden el resto de gráficas para casos supersónicos
    % haciendo diferenciación entre formas triangulares y no triangulares:
    % DATCOM Cap4-PartB pag 103/592.


    a_bw_w1 = KBW_w1*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.34 PAMADI
    a_wb_w1 = KWB_w1*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.33 PAMADI
    a_WB_w1 = (KN_w1 + KWB_w1 + KBW_w1)*CLalpha_w1_e_pw*(S_w1_e/S_w1); %eq 3.24 PAMADI
    CLalpha_WB_w1 = a_WB_w1; % contribution of w1 to fuselage and fuselage to w1
    CL_alpha_wb_w1 = CLalpha_WB_w1;
    CL_alpha_w1_corrected = CL_alpha_wb_w1;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_w1 = CL_alpha_wb_w1;
    Stab_Der_parts.a_WB_w1 = a_WB_w1;
    Stab_Der_parts.a_wb_w1 = a_wb_w1;
    Stab_Der_parts.KWB_w1 = KWB_w1;
    Stab_Der_parts.KBW_w1 = KBW_w1;

    CLalpha_fus = 0;
end

%% Canard body interference
if Can == 1
    aN_can = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    aN_can = CLa_fus; %  corrigiendo con la integrtación del fuselaje a partir de modelo XFLR5

    KN_can = (aN_can/CLalpha_can_e_pw)*(S_can/S_can_e);  %eq 3.25 PAMADI
    CLalpha_fus = aN_can;

    KWB_can = 0.1714*(w_Area_b_max/b_can)^2 + 0.8326*(w_Area_b_max/b_can) + 0.9974; %eq 3.27 PAMADI
    KBW_can = 0.7810*(w_Area_b_max/b_can)^2 + 1.1976*(w_Area_b_max/b_can) + 0.0088; %eq 3.28 PAMADI

    % %     C*********       KW(B) VS D/B  DATA FIGURE 4.3.1.2-10-A
    %
    %        X10A=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %        Y10A=[1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0];
    %        KWB_can = interp1(X10A, Y10A, w_Area_b_max/b_can);
    %
    % % C ****         KB(W) VS D/B  DATA  FIGURE 4.3.1.2-10-B
    %
    %       X10B=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %       Y10B=[0.0,0.13,0.29,0.45,0.62,0.80,1.0,1.22,1.45,1.70,2.0];
    %         KBW_can = interp1(X10B, Y10B, w_Area_b_max/b_can);

    % Queda preguntar si se añaden el resto de gráficas para casos supersónicos
    % haciendo diferenciación entre formas triangulares y no triangulares:
    % DATCOM Cap4-PartB pag 103/592.


    %     a_bw_can = KBW_can*CLalpha_can_e_pw*(S_can_e/S_can); %eq 3.34 PAMADI
    %     a_wb_can = KWB_can*CLalpha_can_e_pw*(S_can_e/S_can); %eq 3.33 PAMADI
    a_WB_can = (KN_can + KWB_can + KBW_can)*CLalpha_can_e_pw*(S_can_e/S_can); %eq 3.24 PAMADI
    CLalpha_WB_can = a_WB_can; % contribution of ww1 to fuselage and fuselage to w1
    CL_alpha_wb_can = CLalpha_WB_can;
    %% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
    CL_alpha_wb_can = CLalpha_can_e_pw;
    CL_alpha_can_corrected = CL_alpha_wb_can;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_can = CL_alpha_wb_can;
else
    Stab_Der_parts.CL_alpha_wb_can = 0;
end

%% W2 body interference (HTP or Vee tail)
if HTP == 1
    aN_HTP = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    aN_HTP = CLa_fus; %  corrigiendo con la integrtación del fuselaje a partir de modelo XFLR5
    KN_HTP = (aN_HTP/CLalpha_HTP_e_pw)*(S_HTP/S_HTP_s);  %eq 3.25 PAMADI
    CLalpha_fus = aN_HTP;

    KWB_HTP = 0.1714*(w_Area_b_max/b_HTP)^2 + 0.8326*(w_Area_b_max/b_HTP) + 0.9974; %eq 3.27 PAMADI
    KBW_HTP = 0.7810*(w_Area_b_max/b_HTP)^2 + 1.1976*(w_Area_b_max/b_HTP) + 0.0088; %eq 3.28 PAMADI
    % %     C*********       KW(B) VS D/B  DATA FIGURE 4.3.1.2-10-A
    %
    %        X10A=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %        Y10A=[1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0];
    %        KWB_HTP = interp1(X10A, Y10A, w_Area_b_max/b_HTP);
    %
    % % C ****         KB(W) VS D/B  DATA  FIGURE 4.3.1.2-10-B
    %
    %       X10B=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %       Y10B=[0.0,0.13,0.29,0.45,0.62,0.80,1.0,1.22,1.45,1.70,2.0];
    %         KBW_HTP = interp1(X10B, Y10B, w_Area_b_max/b_HTP);

    % Queda preguntar si se añaden el resto de gráficas para casos supersónicos
    % haciendo diferenciación entre formas triangulares y no triangulares:
    % DATCOM Cap4-PartB pag 103/592.
    %     a_bw_HTP = KBW_w1*CLalpha_HTP_e_pw*(S_HTP_s/S_HTP); %eq 3.34 PAMADI
    %     a_wb_HTP = KWB_w1*CLalpha_HTP_e_pw*(S_HTP_s/S_HTP); %eq 3.33 PAMADI
    a_WB_HTP = (KN_HTP + KWB_HTP + KBW_HTP)*CLalpha_HTP_e_pw*(S_HTP_s/S_HTP); %eq 3.24 PAMADI
    CLalpha_WB_HTP = a_WB_HTP; % contribution of HTP to fuselage and fuselage to HTP
    CLalpha_wb_HTP = CLalpha_WB_HTP;
    %% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
    CL_alpha_wb_HTP = CLalpha_HTP_e_pw;
    CL_alpha_HTP_corrected = CLalpha_wb_HTP;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_HTP = CL_alpha_wb_HTP;
else
        Stab_Der_parts.CL_alpha_wb_HTP = 0;
end

%% Vee body interference (HTP or Vee tail)
if Vee==1
    %     aN_vee = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    aN_vee = CLa_fus; %  corrigiendo con la integrtación del fuselaje a partir de modelo XFLR5
    KN_vee = (aN_vee/CLalpha_vee_e_pw)*(S_vee/S_vee_s);  %eq 3.25 PAMADI
    CLalpha_fus = aN_vee;

    KWB_vee = 0.1714*(w_Area_b_max/b_vee)^2 + 0.8326*(w_Area_b_max/b_vee) + 0.9974; %eq 3.27 PAMADI
    KBW_vee = 0.7810*(w_Area_b_max/b_vee)^2 + 1.1976*(w_Area_b_max/b_vee) + 0.0088; %eq 3.28 PAMADI
    % %     C*********       KW(B) VS D/B  DATA FIGURE 4.3.1.2-10-A
    %
    %        X10A=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %        Y10A=[1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0];
    %        KWB_vee = interp1(X10A, Y10A, w_Area_b_max/b_vee);
    %
    % % C ****         KB(W) VS D/B  DATA  FIGURE 4.3.1.2-10-B
    %
    %       X10B=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %       Y10B=[0.0,0.13,0.29,0.45,0.62,0.80,1.0,1.22,1.45,1.70,2.0];
    %         KBW_vee = interp1(X10B, Y10B, w_Area_b_max/b_vee);

    % Queda preguntar si se añaden el resto de gráficas para casos supersónicos
    % haciendo diferenciación entre formas triangulares y no triangulares:
    % DATCOM Cap4-PartB pag 103/592.
    %     a_bw_vee = KBW_w1*CLalpha_vee_e_pw*(S_vee_s/S_vee); %eq 3.34 PAMADI
    %     a_wb_vee = KWB_w1*CLalpha_vee_e_pw*(S_vee_s/S_vee); %eq 3.33 PAMADI
    a_WB_vee = (KN_vee + KWB_vee + KBW_vee)*CLalpha_vee_e_pw*(S_vee_s/S_vee); %eq 3.24 PAMADI
    CLalpha_WB_vee = a_WB_vee; % contribution of vee to fuselage and fuselage to vee
    CLalpha_wb_vee = CLalpha_WB_vee;
    %% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
    CL_alpha_wb_vee = CLalpha_vee_e_pw;
    CL_alpha_vee_corrected = CLalpha_wb_vee;
    % Correction of CLalphais for wing with no dihedral
    %     CL_alpha_wb_vee = CLalpha_vee_e_pw*(cos(dihedral_vee_e))^2;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_vee = CL_alpha_wb_vee;

else
    Stab_Der_parts.CL_alpha_wb_vee = 0;
end

%% vee2 body interference (HTP or vee2 tail)
if Vee2==1
    %     aN_vee2 = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    aN_vee2 = CLa_fus; %  corrigiendo con la integrtación del fuselaje a partir de modelo XFLR5
    KN_vee2 = (aN_vee2/CLalpha_vee2_e_pw)*(S_vee2/S_vee2_s);  %eq 3.25 PAMADI
    CLalpha_fus = aN_vee2;

    KWB_vee2 = 0.1714*(w_Area_b_max/b_vee2)^2 + 0.8326*(w_Area_b_max/b_vee2) + 0.9974; %eq 3.27 PAMADI
    KBW_vee2 = 0.7810*(w_Area_b_max/b_vee2)^2 + 1.1976*(w_Area_b_max/b_vee2) + 0.0088; %eq 3.28 PAMADI
    % %     C*********       KW(B) VS D/B  DATA FIGURE 4.3.1.2-10-A
    %
    %        X10A=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %        Y10A=[1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0];
    %        KWB_vee2 = interp1(X10A, Y10A, w_Area_b_max/b_vee2);
    %
    % % C ****         KB(W) VS D/B  DATA  FIGURE 4.3.1.2-10-B
    %
    %       X10B=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %       Y10B=[0.0,0.13,0.29,0.45,0.62,0.80,1.0,1.22,1.45,1.70,2.0];
    %         KBW_vee2 = interp1(X10B, Y10B, w_Area_b_max/b_vee2);

    % Queda preguntar si se añaden el resto de gráficas para casos supersónicos
    % haciendo diferenciación entre formas triangulares y no triangulares:
    % DATCOM Cap4-PartB pag 103/592.
    %     a_bw_vee2 = KBW_w1*CLalpha_vee2_e_pw*(S_vee2_s/S_vee2); %eq 3.34 PAMADI
    %     a_wb_vee2 = KWB_w1*CLalpha_vee2_e_pw*(S_vee2_s/S_vee2); %eq 3.33 PAMADI
    a_WB_vee2 = (KN_vee2 + KWB_vee2 + KBW_vee2)*CLalpha_vee2_e_pw*(S_vee2_s/S_vee2); %eq 3.24 PAMADI
    CLalpha_WB_vee2 = a_WB_vee2; % contribution of vee2 to fuselage and fuselage to vee2
    CLalpha_wb_vee2 = CLalpha_WB_vee2;
    %% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
    CL_alpha_wb_vee2 = CLalpha_vee2_e_pw;
    CL_alpha_vee2_corrected = CLalpha_wb_vee2;
    % Correction of CLalphais for wing with no dihedral
    %     CL_alpha_wb_vee2 = CLalpha_vee2_e_pw*(cos(dihedral_vee2_e))^2;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_vee2 = CL_alpha_wb_vee2;

else
        Stab_Der_parts.CL_alpha_wb_vee2 = 0;
end

if Nac == 1
        Stab_Der_parts.CL_alpha_wb_nac = 0;
else
        Stab_Der_parts.CL_alpha_wb_nac = 0;
end

CL_alpha_wb_nac = 0;
% Storing Values
Stab_Der_parts.CL_alpha_fus = CLalpha_fus;
Stab_Der_parts.CL_alpha_nac = CL_alpha_wb_nac;

% Identifies which CLalpha is used, the prop-wash correction or the body
% influence
body_interference = 0;

if W1 == 1
    switch body_interference
        case 0
            CL_alpha_w1 = CLalpha_w1_e_pw;
        case 1
            CL_alpha_w1 = CL_alpha_wb_w1;
    end
else
    CL_alpha_w1 = 0;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
end

if Can == 1
    switch body_interference
        case 0
            CL_alpha_can = CLalpha_can_e_pw;
        case 1
            CL_alpha_can = CL_alpha_wb_can;
    end
else
    CL_alpha_can = 0;
    Stab_Der_parts.CL_alpha_can = CL_alpha_can;
end

if HTP == 1
    switch body_interference
        case 0
            CL_alpha_HTP = CLalpha_HTP_e_pw;
        case 1
            CL_alpha_HTP = CL_alpha_wb_HTP;
    end
else
    CL_alpha_HTP = 0;
    Stab_Der_parts.CL_alpha_HTP = CL_alpha_HTP;
end

if Vee == 1
    switch body_interference
        case 0
            CL_alpha_vee = CLalpha_vee_e_pw;
        case 1
            CL_alpha_vee = CL_alpha_wb_vee;
    end
else
    CL_alpha_vee = 0;
    Stab_Der_parts.CL_alpha_vee = CL_alpha_vee;
end

if Vee2 == 1
    switch body_interference
        case 0
            CL_alpha_vee2 = CLalpha_vee2_e_pw;
        case 1
            CL_alpha_vee2 = CL_alpha_wb_vee2;
    end
else
    CL_alpha_vee2 = 0;
    Stab_Der_parts.CL_alpha_vee2 = CL_alpha_vee2;
end

if VTP == 1
    CLalpha_VTP_e_pw = Stab_Der_parts.CLalpha_VTP_e_pw ; %modificar en caso de que body_interference = 1;

    switch body_interference
        case 0
            CLalpha_VTP = CLalpha_VTP_e_pw;
        case 1
            CLalpha_VTP = CL_alpha_wb_VTP;
    end
    Stab_Der_parts.CLalpha_VTP = CLalpha_VTP;
else
    CLalpha_VTP = 0;
    Stab_Der_parts.CLalpha_VTP = CLalpha_VTP;
end


if Nac == 1
    CL_alpha_nac = 0;
    Stab_Der_parts.CL_alpha_nac = CL_alpha_nac;
else
    CL_alpha_nac = 0;
    Stab_Der_parts.CL_alpha_nac = CL_alpha_nac;
end


if AC_type == 1
    CL_alpha_ac = CL_alpha_w1;
    %     X_NP = (CL_alpha_w1*XAC_WB)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1;
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
elseif AC_type == 2
   if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
       CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_HTP*(downwash_HTP);
   else 
       CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_HTP*(downwash_HTP);
   end
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_HTP*(downwash_HTP)*x_xbar_HTP)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_HTP_e_corrected + CL_alpha_HTP*(i_HTP - eps_HTP);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_HTP_e_corrected = CL0_HTP_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_HTP = CL0_HTP_e_corrected;
    Stab_Der_parts.CL_alpha_HTP = CL_alpha_HTP;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
elseif AC_type == 3
    if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
        % CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_can*(upwash) + CL_alpha_HTP*(downwash_HTP);
        CL_alpha_ac = CL_alpha_wb_w1*(downwash_w1) + CL_alpha_can*(upwash) + CL_alpha_HTP*(downwash_HTP);
    else
        CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_can*(upwash) + CL_alpha_HTP*(downwash_HTP);
        CL_alpha_ac = CL_alpha_w1*(downwash_w1) + CLa_fus + CL_alpha_can*(upwash) + CL_alpha_HTP*(downwash_HTP);
    end
    
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can + CL_alpha_HTP*(downwash_HTP)*x_xbar_HTP)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_HTP_e_corrected + CL_alpha_HTP*(i_HTP - eps_HTP) + ...
        CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*(i_w1-eps_w1) + CL0_HTP_e_corrected + CL_alpha_HTP*(i_HTP - eps_HTP) + ...
        CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);


    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_HTP_e_corrected = CL0_HTP_e_corrected;
    Trim_ITER.CL0_can_e_corrected = CL0_can_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_HTP = CL0_HTP_e_corrected;
    Stab_Der_parts.CL_alpha_HTP = CL_alpha_HTP;
    Stab_Der_parts.CL0_can = CL0_can_e_corrected;
    Stab_Der_parts.CL_alpha_can = CL_alpha_can;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;

elseif AC_type == 4
    if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
        CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_vee*(downwash_vee);
    else
        CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_vee*(downwash_vee);
    end
    
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_vee*(downwash_vee)*x_xbar_HTP)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_vee_e_corrected + CL_alpha_vee*(i_vee - eps_vee);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_vee_e_corrected = CL0_vee_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_vee = CL0_vee_e_corrected;
    Stab_Der_parts.CL_alpha_vee = CL_alpha_vee;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;

elseif AC_type == 5
    if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
        CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_can*(upwash) + CL_alpha_vee*(downwash_vee);
        CL_alpha_ac = CL_alpha_wb_w1*(downwash_w1) + CL_alpha_can*(upwash) + CL_alpha_vee*(downwash_vee);
    else
       CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_can*(upwash) + CL_alpha_vee*(downwash_vee);
       CL_alpha_ac = CL_alpha_w1*(downwash_w1) + CLa_fus + CL_alpha_can*(upwash) + CL_alpha_vee*(downwash_vee);
    end
    
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can + CL_alpha_vee*(downwash_vee)*x_xbar_HTP)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_vee_e_corrected + CL_alpha_vee*(i_vee - eps_vee) + ...
        CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*(i_w1-eps_w1) + CL0_vee_e_corrected + CL_alpha_vee*(i_vee - eps_vee) + ...
        CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);

    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_vee_e_corrected = CL0_vee_e_corrected;
    Trim_ITER.CL0_can_e_corrected = CL0_can_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_vee = CL0_vee_e_corrected;
    Stab_Der_parts.CL_alpha_vee = CL_alpha_vee;
    Stab_Der_parts.CL0_can = CL0_can_e_corrected;
    Stab_Der_parts.CL_alpha_can = CL_alpha_can;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;

elseif AC_type == 6
    if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
        CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_can*(upwash);
        CL_alpha_ac = CL_alpha_wb_w1*(downwash_w1) + CL_alpha_can*(upwash);
    else
        CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_can*(upwash);
        CL_alpha_ac = CL_alpha_w1*(downwash_w1) + CLa_fus + CL_alpha_can*(upwash);
    end
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*(i_w1-eps_w1) + CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);

    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_can_e_corrected = CL0_can_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_can = CL0_can_e_corrected;
    Stab_Der_parts.CL_alpha_can = CL_alpha_can;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;

elseif AC_type == 7
    if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
        CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_vee*(downwash_vee);
    else
        CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_vee*(downwash_vee);
    end
    
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_vee*(downwash_vee)*x_xbar_HTP)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_vee_e_corrected + CL_alpha_vee*(i_vee - eps_vee);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_vee_e_corrected = CL0_vee_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_vee = CL0_vee_e_corrected;
    Stab_Der_parts.CL_alpha_vee = CL_alpha_vee;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;

elseif AC_type == 8
    if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
        CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_vee*(downwash_vee) + CL_alpha_vee2*(downwash_vee2);
    else
        CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_vee*(downwash_vee) + CL_alpha_vee2*(downwash_vee2);
    end
    
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_vee*(downwash_vee)*x_xbar_HTP)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_vee_e_corrected + CL_alpha_vee*(i_vee - eps_vee) + CL0_vee2_e_corrected + CL_alpha_vee2*(i_vee2 - eps_vee2);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_vee_e_corrected = CL0_vee_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_vee = CL0_vee_e_corrected;
    Stab_Der_parts.CL0_vee2 = CL0_vee2_e_corrected;
    Stab_Der_parts.CL_alpha_vee = CL_alpha_vee;
    Stab_Der_parts.CL_alpha_vee2 = CL_alpha_vee2;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
end

end
