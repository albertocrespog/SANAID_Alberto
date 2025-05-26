function [Stab_Der_parts,  Trim_ITER] = getCLalpha(AC_CONFIGURATION,Body_Geo,Geo_tier,body_interference,Stab_Der_parts, Effects,Design_criteria, OUTPUT_read_XLSX)

%% Input
W1      = AC_CONFIGURATION.W1;
HTP     = AC_CONFIGURATION.HTP;
VTP     = AC_CONFIGURATION.VTP;
Can     = AC_CONFIGURATION.Can;
Vee     = AC_CONFIGURATION.Vee;
Nac     = AC_CONFIGURATION.Nac;
AC_type = AC_CONFIGURATION.AC_type;

length_fus         = Body_Geo.l_fus;
w_Area_b_max       = Body_Geo.w_Area_b_max;
x_Area_b_max       = Body_Geo.x_Area_b_max;
Area_b_max         = Body_Geo.Area_b_max;
S_ref              = Geo_tier.S_ref;
CLa_fus            = Body_Geo.CLa_fus;
CLalpha_w1_e_pw    = Stab_Der_parts.CLalpha_w1_e_pw;
CL0_w1_e_corrected = Stab_Der_parts.CL0_w1_e_corrected;

S_w1   = Geo_tier.S_w1;
S_w1_e = Geo_tier.S_w1_e;
b_w1   = Geo_tier.b_w1;
i_w1   = Design_criteria.i_w1;

wingspan2bodydiam = b_w1/w_Area_b_max;

flagwingspan2bodydiam = OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam; 

if Can == 1
    upwash              = Effects.upwash;
    CL0_can_e_corrected = Stab_Der_parts.CL0_can_e_corrected;
    CLalpha_can_e_pw    = Stab_Der_parts.CLalpha_can_e_pw;
    S_can               = Geo_tier.S_can;
    S_can_e             = Geo_tier.S_can_e; % effective wing area: proyected in the y-plane
    b_can               = Geo_tier.b_can;
    i_can               = Design_criteria.i_w3;
    eps_can             = Effects.eps_can;
end
if HTP == 1
    downwash            = Effects.downwash;
    CL0_HTP_e_corrected = Stab_Der_parts.CL0_HTP_e_corrected;
    CLalpha_HTP_e_pw    = Stab_Der_parts.CLalpha_HTP_e_pw;
    S_w2                = Geo_tier.S_w2;
    S_w2_s              = Geo_tier.S_w2_s; % real wing area: proyected along the surface of wing
    b_w2                = Geo_tier.b_w2;
    i_w2                = Design_criteria.i_w2;
    eps_w2              = Effects.eps_w2;
end
if Vee == 1
    downwash            = Effects.downwash;
    CL0_Vee_e_corrected = Stab_Der_parts.CL0_Vee_e_corrected;
    CLalpha_Vee_e_pw    = Stab_Der_parts.CLalpha_Vee_e_pw;
    S_w2                = Geo_tier.S_w2;
    S_w2_s              = Geo_tier.S_w2_s; % real wing area: proyected along the surface of wing
    b_w2                = Geo_tier.b_w2;
    i_w2                = Design_criteria.i_w2;
    eps_w2              = Effects.eps_w2;
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
end

%% W2 body interference (HTP or Vee tail)
if HTP == 1
    aN_w2 = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    aN_w2 = CLa_fus; %  corrigiendo con la integrtación del fuselaje a partir de modelo XFLR5
    KN_w2 = (aN_w2/CLalpha_HTP_e_pw)*(S_w2/S_w2_s);  %eq 3.25 PAMADI
    CLalpha_fus = aN_w2;

    KWB_w2 = 0.1714*(w_Area_b_max/b_w2)^2 + 0.8326*(w_Area_b_max/b_w2) + 0.9974; %eq 3.27 PAMADI
    KBW_w2 = 0.7810*(w_Area_b_max/b_w2)^2 + 1.1976*(w_Area_b_max/b_w2) + 0.0088; %eq 3.28 PAMADI
    % %     C*********       KW(B) VS D/B  DATA FIGURE 4.3.1.2-10-A
    %
    %        X10A=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %        Y10A=[1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0];
    %        KWB_w2 = interp1(X10A, Y10A, w_Area_b_max/b_w2);
    %
    % % C ****         KB(W) VS D/B  DATA  FIGURE 4.3.1.2-10-B
    %
    %       X10B=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %       Y10B=[0.0,0.13,0.29,0.45,0.62,0.80,1.0,1.22,1.45,1.70,2.0];
    %         KBW_w2 = interp1(X10B, Y10B, w_Area_b_max/b_w2);

    % Queda preguntar si se añaden el resto de gráficas para casos supersónicos
    % haciendo diferenciación entre formas triangulares y no triangulares:
    % DATCOM Cap4-PartB pag 103/592.
    %     a_bw_w2 = KBW_w1*CLalpha_HTP_e_pw*(S_w2_s/S_w2); %eq 3.34 PAMADI
    %     a_wb_w2 = KWB_w1*CLalpha_HTP_e_pw*(S_w2_s/S_w2); %eq 3.33 PAMADI
    a_WB_w2 = (KN_w2 + KWB_w2 + KBW_w1)*CLalpha_HTP_e_pw*(S_w2_s/S_w2); %eq 3.24 PAMADI
    CLalpha_WB_w2 = a_WB_w2; % contribution of w2 to fuselage and fuselage to w2
    CLalpha_wb_w2 = CLalpha_WB_w2;
    %% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
    CL_alpha_wb_HTP = CLalpha_HTP_e_pw;
    CL_alpha_HTP_corrected = CLalpha_wb_w2;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_HTP = CL_alpha_wb_HTP;
end

%% W2 body interference (HTP or Vee tail)
if Vee==1
    %     aN_w2 = 2*f_k2_k1*(Area_b_max/S_ref);   %eq 3.26 PAMADI
    aN_w2 = CLa_fus; %  corrigiendo con la integrtación del fuselaje a partir de modelo XFLR5
    KN_w2 = (aN_w2/CLalpha_Vee_e_pw)*(S_w2/S_w2_s);  %eq 3.25 PAMADI
    CLalpha_fus = aN_w2;

    KWB_w2 = 0.1714*(w_Area_b_max/b_w2)^2 + 0.8326*(w_Area_b_max/b_w2) + 0.9974; %eq 3.27 PAMADI
    KBW_w2 = 0.7810*(w_Area_b_max/b_w2)^2 + 1.1976*(w_Area_b_max/b_w2) + 0.0088; %eq 3.28 PAMADI
    % %     C*********       KW(B) VS D/B  DATA FIGURE 4.3.1.2-10-A
    %
    %        X10A=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %        Y10A=[1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0];
    %        KWB_w2 = interp1(X10A, Y10A, w_Area_b_max/b_w2);
    %
    % % C ****         KB(W) VS D/B  DATA  FIGURE 4.3.1.2-10-B
    %
    %       X10B=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    %       Y10B=[0.0,0.13,0.29,0.45,0.62,0.80,1.0,1.22,1.45,1.70,2.0];
    %         KBW_w2 = interp1(X10B, Y10B, w_Area_b_max/b_w2);

    % Queda preguntar si se añaden el resto de gráficas para casos supersónicos
    % haciendo diferenciación entre formas triangulares y no triangulares:
    % DATCOM Cap4-PartB pag 103/592.
    %     a_bw_w2 = KBW_w1*CLalpha_Vee_e_pw*(S_w2_s/S_w2); %eq 3.34 PAMADI
    %     a_wb_w2 = KWB_w1*CLalpha_Vee_e_pw*(S_w2_s/S_w2); %eq 3.33 PAMADI
    a_WB_w2 = (KN_w2 + KWB_w2 + KBW_w1)*CLalpha_Vee_e_pw*(S_w2_s/S_w2); %eq 3.24 PAMADI
    CLalpha_WB_w2 = a_WB_w2; % contribution of w2 to fuselage and fuselage to w2
    CLalpha_wb_w2 = CLalpha_WB_w2;
    %% NOTE: maintains the correction with dynamic pressure but not the corrections with body interference
    CL_alpha_wb_Vee = CLalpha_Vee_e_pw;
    CL_alpha_Vee_corrected = CLalpha_wb_w2;
    % Correction of CLalphais for wing with no dihedral
    %     CL_alpha_wb_Vee = CLalpha_Vee_e_pw*(cos(dihedral_w2_e))^2;
    % Storing Values
    Stab_Der_parts.CL_alpha_wb_Vee = CL_alpha_wb_Vee;

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
end

if Can == 1
    switch body_interference
        case 0
            CL_alpha_can = CLalpha_can_e_pw;
        case 1
            CL_alpha_can = CL_alpha_wb_can;
    end
end

if HTP == 1
    switch body_interference
        case 0
            CL_alpha_HTP = CLalpha_HTP_e_pw;
        case 1
            CL_alpha_HTP = CL_alpha_wb_HTP;
    end
end

if Vee == 1
    switch body_interference
        case 0
            CL_alpha_Vee = CLalpha_Vee_e_pw;
        case 1
            CL_alpha_Vee = CL_alpha_wb_Vee;
    end
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
       CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_HTP*(downwash);
   else 
       CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_HTP*(downwash);
   end
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_HTP*(downwash)*x_xbar_w2)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_HTP_e_corrected + CL_alpha_HTP*(i_w2 - eps_w2);
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
        CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_can*(upwash) + CL_alpha_HTP*(downwash);
    else
        CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_can*(upwash) + CL_alpha_HTP*(downwash);
    end
    
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can + CL_alpha_HTP*(downwash)*x_xbar_w2)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_HTP_e_corrected + CL_alpha_HTP*(i_w2 - eps_w2) + ...
        CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_HTP_e_corrected = CL0_HTP_e_corrected;
    Trim_ITER.CL0_can_e_corrected = CL0_can_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_HTP = CL0_HTP_e_corrected;
    Stab_Der_parts.CL_alpha_HTP = CL_alpha_htp;
    Stab_Der_parts.CL0_can = CL0_can_e_corrected;
    Stab_Der_parts.CL_alpha_can = CL_alpha_can;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
elseif AC_type == 4
    if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
        CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_Vee*(downwash);
    else
        CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_Vee*(downwash);
    end
    
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_Vee*(downwash)*x_xbar_w2)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_Vee_e_corrected + CL_alpha_Vee*(i_w2 - eps_w2);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_Vee_e_corrected = CL0_Vee_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_Vee = CL0_Vee_e_corrected;
    Stab_Der_parts.CL_alpha_Vee = CL_alpha_Vee;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
elseif AC_type == 5
    if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
        CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_can*(upwash) + CL_alpha_Vee*(downwash);
    else
       CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_can*(upwash) + CL_alpha_Vee*(downwash);
    end
    
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can + CL_alpha_Vee*(downwash)*x_xbar_w2)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_Vee_e_corrected + CL_alpha_Vee*(i_w2 - eps_w2) + ...
        CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);
    Trim_ITER.CL0_w1_e_corrected = CL0_w1_e_corrected;
    Trim_ITER.CL0_Vee_e_corrected = CL0_Vee_e_corrected;
    Trim_ITER.CL0_can_e_corrected = CL0_can_e_corrected;
    Trim_ITER.CL0_ac = CL0_ac;
    % Store values
    Stab_Der_parts.CL0_w1 = CL0_w1_e_corrected;
    Stab_Der_parts.CL_alpha_w1 = CL_alpha_w1;
    Stab_Der_parts.CL0_Vee = CL0_Vee_e_corrected;
    Stab_Der_parts.CL_alpha_Vee = CL_alpha_Vee;
    Stab_Der_parts.CL0_can = CL0_can_e_corrected;
    Stab_Der_parts.CL_alpha_can = CL_alpha_can;
    Stab_Der_parts.CL0_ac = CL0_ac;
    Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
elseif AC_type == 6
    if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
        CL_alpha_ac = CL_alpha_wb_w1 + CL_alpha_can*(upwash);
        
    else
        CL_alpha_ac = CL_alpha_w1 + CLa_fus + CL_alpha_can*(upwash);
    end
    %     X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can)/CL_alpha_ac;
    CL0_ac = CL0_w1_e_corrected + CL_alpha_w1*i_w1 + CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can);
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
end


end
