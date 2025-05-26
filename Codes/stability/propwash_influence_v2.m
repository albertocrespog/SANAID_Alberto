function [Stab_Der_parts, afe] = propwash_influence(AC_CONFIGURATION,Geo_tier,Posicion_Palanca,conditions,Propulsion,Performance,Aero,Effects,OUTPUT_read_XLSX,Body_Geo)

W1  = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Nac = AC_CONFIGURATION.Nac;

AC_type  = AC_CONFIGURATION.AC_type;
twin_VTP = AC_CONFIGURATION.twin_VTP;



% Upwash influencing in prop
prop_wash_effect      = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;
flagwingspan2bodydiam = OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam;

% Location of the prop disk (source of thrust)
x_prop_cF = Geo_tier.x_prop_cF; % location of prop
AR_w1     = Geo_tier.AR_w1;
x_xbar_w1 = Geo_tier.x_xbar_w1;
x_prop_cF = Geo_tier.x_prop_cF; % location of prop
cmac_w1   = Geo_tier.cmac_w1;
S_w1_pw   = Geo_tier.S_w1_pw;
S_w1      = Geo_tier.S_w1;
S_ref     = Geo_tier.S_ref;
S_w1_e    = Geo_tier.S_w1_e;
b_w1      = Geo_tier.b_w1;

V   = conditions.V;
rho = Performance.rho;

CLalpha_w1   = Aero.CL_alpha_w1_CR;
CL0_w1       = Aero.CL_0_w1_CR;


depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);
q_inf = 0.5*rho*V^2;
afe.depsu_dalpha = depsu_dalpha;
w_Area_b_max = Body_Geo.w_Area_b_max;

wingspan2bodydiam = b_w1/w_Area_b_max;

if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
    CLalpha_w1_e = CLalpha_w1*S_w1_e/S_w1;
    CL0_w1_e     = CL0_w1*S_w1_e/S_w1;
else
    CLalpha_w1_e = CLalpha_w1;
    CL0_w1_e     = CL0_w1;
end

%% Selects if prop wash effect is calculated
%% !!!!!!!!!!!!!!!! BE CAREFUL!!!! IT SEEMS THAT THE AREAS ARE CALCULATED WRONG AS IT PROVIDES NEGATIVE VALUES!!!!
%% NEED TO REVIEW SO USE prop_wash_effect = 0 IN THE eXCEL
if prop_wash_effect == 1
    v_i = Propulsion.v_i; % Induced velocity at prop disk
    
    if W1 == 1
        S_w1_afe = S_w1_pw;
        S_w1_no_afe = S_w1 - S_w1_afe;
        q_w1_no_afe = q_inf;
        if Posicion_Palanca == 0
            V_w1 = V;
            q_w1_afe    = q_inf;
        else
            V_w1 = V + v_i;
            q_w1_afe       = 0.5*rho*V_w1^2;
        end
        eta_w1_afe = (q_w1_afe/q_inf);
        eta_w1_no_afe = (q_w1_no_afe/q_inf);
        eta_w1_afe_S_w1_afe_S_ref = eta_w1_afe*(S_w1_afe/S_ref);
        eta_w1_no_afe_S_w1_no_afe_S_ref = eta_w1_no_afe*(S_w1_no_afe/S_ref);
        afe.eta_w1_afe_S_w1_afe_S_ref = eta_w1_afe_S_w1_afe_S_ref;
        afe.eta_w1_no_afe_S_w1_no_afe_S_ref = eta_w1_no_afe_S_w1_no_afe_S_ref;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Correction of Lift curve slope of w1 associated to prop wash
        % eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
        % eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
        CLalpha_w1_e_pw = eta_w1_afe_S_w1_afe_S_ref*CLalpha_w1_e + eta_w1_no_afe_S_w1_no_afe_S_ref*CLalpha_w1_e;
        CL0_w1_e_corrected = eta_w1_afe_S_w1_afe_S_ref*CL0_w1_e + eta_w1_no_afe_S_w1_no_afe_S_ref*CL0_w1_e;
        Stab_Der_parts.CL0_w1_e_corrected = CL0_w1_e_corrected;
        Stab_Der_parts.CLalpha_w1_e_pw = CLalpha_w1_e_pw;
        
    end
    
    if Can == 1
        S_can_pw = Geo_tier.S_can_pw;
        S_can    = Geo_tier.S_can;
        b_can    = Geo_tier.b_can;
        S_can_e  = Geo_tier.S_can_e;
        
        w_can = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_can);
        wingspan2bodydiam_can = b_can/w_can;
        
        CLalpha_can   = Aero.CL_alpha_can_CR;
        CL0_can       = Aero.CL_0_can_CR;
        
        if wingspan2bodydiam_can <= 2 || flagwingspan2bodydiam==1
            CLalpha_can_e = CLalpha_can*S_can_e/Geo_tier.S_can;
            CL0_can_e     = CL0_can*S_can_e/Geo_tier.S_can;
        else
            CLalpha_can_e = CLalpha_can;
            CL0_can_e     = CL0_can;
        end
        
        S_can_afe = S_can_pw;
        S_can_no_afe = S_can - S_can_afe;
        q_can_no_afe = q_inf;
        
        if Posicion_Palanca == 0
            V_can = V;
            q_can_afe    = q_inf;
        else
            V_can = V + v_i;
            q_can_afe       = 0.5*rho*V_can^2;
        end
        
        eta_can_afe = (q_can_afe/q_inf);
        eta_can_no_afe = (q_can_no_afe/q_inf);
        afe.eta_can_afe = eta_can_afe;
        afe.eta_can_no_afe = eta_can_no_afe;
        eta_can_afe_S_can_afe_S_ref = eta_can_afe*(S_can_afe/S_ref);
        eta_can_no_afe_S_can_no_afe_S_ref = eta_can_no_afe*(S_can_no_afe/S_ref);
        afe.eta_can_afe_S_can_afe_S_ref = eta_can_afe_S_can_afe_S_ref;
        afe.eta_can_no_afe_S_can_no_afe_S_ref = eta_can_no_afe_S_can_no_afe_S_ref;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Correction of Lift curve slope of can associated to prop wash
        % eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
        % eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
        CLalpha_can_e_pw = eta_can_afe_S_can_afe_S_ref*CLalpha_can_e + eta_can_no_afe_S_can_no_afe_S_ref*CLalpha_can_e;
        CL0_can_e_corrected = eta_can_afe_S_can_afe_S_ref*CL0_can_e + eta_can_no_afe_S_can_no_afe_S_ref*CL0_can_e;
        Stab_Der_parts.CL0_can_e_corrected = CL0_can_e_corrected;
        Stab_Der_parts.CLalpha_can_e_pw = CLalpha_can_e_pw;
    end
    
    if HTP == 1
        
        S_w2_pw = Geo_tier.S_w2_pw;
        S_w2_s  = Geo_tier.S_w2_s;
        S_w2_e  = Geo_tier.S_w2_e;
        b_w2  = Geo_tier.b_w2;
        
        w_w2 = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_w2);
        wingspan2bodydiam_w2 = b_w2/w_w2;
        
        eta_w2 = Effects.eta_w2;
        
        CL0_w2      = Aero.CL_0_w2_CR;
        CL_alpha_w2 = Aero.CL_alpha_w2_CR;
        
        
        if wingspan2bodydiam_w2 <= 2 || flagwingspan2bodydiam==1
            CL0_w2_e     = CL0_w2*S_w2_e/Geo_tier.S_w2;
            CLalpha_w2_e = CL_alpha_w2*S_w2_e/Geo_tier.S_w2;
        else
            CL0_w2_e         = CL0_w2;
            CLalpha_w2_e     = CL_alpha_w2;
        end
        
        S_w2_afe = S_w2_pw;
        S_w2_no_afe = S_w2_s - S_w2_afe;
        q_w2_no_afe    = eta_w2*q_inf;
        %presion dinamica no afectada
        
        if Posicion_Palanca == 0
            V_w2 = V;
            q_w2_afe    = q_inf;
        else
            V_w2 = sqrt(eta_w2)*V + 2*v_i;
            q_w2_afe       = 0.5*rho*V_w2^2;
        end
        
        %     if Posicion_Palanca == 0
        %         V_w1 = V;
        %         q_w1_afe    = q_inf;
        %     else
        %         V_w1 = V + v_i;
        %         q_w1_afe       = 0.5*rho*V_w1^2;
        %     end
        %
        eta_w2_afe = (q_w2_afe/q_inf);
        eta_w2_no_afe = (q_w2_no_afe/q_inf);
        afe.eta_w2_afe = eta_w2_afe;
        afe.eta_w2_no_afe = eta_w2_no_afe;
        eta_w2_afe_S_w2_afe_S_ref = eta_w2_afe*(S_w2_afe/S_ref);
        eta_w2_no_afe_S_w2_no_afe_S_ref = eta_w2_no_afe*(S_w2_no_afe/S_ref);
        afe.eta_w2_afe_S_w2_afe_S_ref = eta_w2_afe_S_w2_afe_S_ref;
        afe.eta_w2_no_afe_S_w2_no_afe_S_ref = eta_w2_no_afe_S_w2_no_afe_S_ref;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Correction of Lift curve slope of HTP associated to prop wash
        % eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
        % eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
        CLalpha_HTP_e_pw = eta_w2_afe_S_w2_afe_S_ref*CLalpha_w2_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CLalpha_w2_e;
        CL0_HTP_e_corrected = eta_w2_afe_S_w2_afe_S_ref*CL0_w2_e + eta_w2_no_afe_S_w2_no_afe_S_ref*CL0_w2_e;
        Stab_Der_parts.CL0_HTP_e_corrected = CL0_HTP_e_corrected;
        Stab_Der_parts.CLalpha_HTP_e_pw = CLalpha_HTP_e_pw;
        
    end
    
    if Vee == 1
        S_vee_pw = Geo_tier.S_vee_pw;
        S_vee_s  = Geo_tier.S_vee_s;
        S_vee_e  =Geo_tier.S_vee_e;
        b_vee  = Geo_tier.b_vee;
        
        w_vee = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee);
        wingspan2bodydiam_vee = b_vee/w_vee;
        
        eta_vee = Effects.eta_vee;
        
        CL0_vee          = Aero.CL_0_vee_CR;
        CL_alpha_wb_Vee = Aero.CLalpha_vee;
        
        if wingspan2bodydiam_vee <= 2 || flagwingspan2bodydiam==1
            CL0_vee_e = CL0_vee*S_vee_e/Geo_tier.S_vee;
            CLalpha_vee_e = CL_alpha_wb_Vee*S_vee_e/Geo_tier.S_vee;
        else
            CL0_vee_e         = CL0_vee;
            CLalpha_vee_e     = CL_alpha_wb_Vee;
        end
        
        
        S_vee_afe = S_vee_pw;
        S_vee_no_afe = S_vee_s - S_vee_afe;
        q_vee_no_afe    = eta_vee*q_inf;
        %presion dinamica no afectada
        
        if Posicion_Palanca == 0
            V_vee = V;
            q_vee_afe    = q_inf;
        else
            V_vee = sqrt(eta_vee)*V + 2*v_i;
            q_vee_afe       = 0.5*rho*V_vee^2;
        end
        
        eta_vee_afe = (q_vee_afe/q_inf);
        eta_vee_no_afe = (q_vee_no_afe/q_inf);
        afe.eta_vee_afe = eta_vee_afe;
        afe.eta_vee_no_afe = eta_vee_no_afe;
        eta_vee_afe_S_vee_afe_S_ref = eta_vee_afe*(S_vee_afe/S_ref);
        eta_vee_no_afe_S_vee_no_afe_S_ref = eta_vee_no_afe*(S_vee_no_afe/S_ref);
        afe.eta_vee_afe_S_vee_afe_S_ref = eta_vee_afe_S_vee_afe_S_ref;
        afe.eta_vee_no_afe_S_vee_no_afe_S_ref = eta_vee_no_afe_S_vee_no_afe_S_ref;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Correction of Lift curve slope of Vee associated to prop wash
        % eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
        % eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
        CLalpha_Vee_e_pw = eta_vee_afe_S_vee_afe_S_ref*CLalpha_vee_e + eta_vee_no_afe_S_vee_no_afe_S_ref*CLalpha_vee_e;
        CL0_Vee_e_corrected = eta_vee_afe_S_vee_afe_S_ref*CL0_vee_e + eta_vee_no_afe_S_vee_no_afe_S_ref*CL0_vee_e;
        Stab_Der_parts.CL0_Vee_e_corrected = CL0_Vee_e_corrected;
        Stab_Der_parts.CLalpha_Vee_e_pw = CLalpha_Vee_e_pw;
    end

    if Vee2 == 1
        S_vee2_pw = Geo_tier.S_vee2_pw;
        S_vee2_s  = Geo_tier.S_vee2_s;
        S_vee2_e  =Geo_tier.S_vee2_e;
        b_vee2  = Geo_tier.b_vee2;
        
        w_vee2 = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee2);
        wingspan2bodydiam_vee2 = b_vee2/w_vee2;
        
        eta_vee2 = Effects.eta_vee2;
        
        CL0_vee2          = Aero.CL_0_vee2_CR;
        CL_alpha_wb_vee2 = Aero.CLalpha_vee2;
        
        if wingspan2bodydiam_vee2 <= 2 || flagwingspan2bodydiam==1
            CL0_vee2_e = CL0_vee2*S_vee2_e/Geo_tier.S_vee2;
            CLalpha_vee2_e = CL_alpha_wb_vee2*S_vee2_e/Geo_tier.S_vee2;
        else
            CL0_vee2_e         = CL0_vee2;
            CLalpha_vee2_e     = CL_alpha_wb_vee2;
        end
        
        
        S_vee2_afe = S_vee2_pw;
        S_vee2_no_afe = S_vee2_s - S_vee2_afe;
        q_vee2_no_afe    = eta_vee2*q_inf;
        %presion dinamica no afectada
        
        if Posicion_Palanca == 0
            V_vee2 = V;
            q_vee2_afe    = q_inf;
        else
            V_vee2 = sqrt(eta_vee2)*V + 2*v_i;
            q_vee2_afe       = 0.5*rho*V_vee2^2;
        end
        
        eta_vee2_afe = (q_vee2_afe/q_inf);
        eta_vee2_no_afe = (q_vee2_no_afe/q_inf);
        afe.eta_vee2_afe = eta_vee2_afe;
        afe.eta_vee2_no_afe = eta_vee2_no_afe;
        eta_vee2_afe_S_vee2_afe_S_ref = eta_vee2_afe*(S_vee2_afe/S_ref);
        eta_vee2_no_afe_S_vee2_no_afe_S_ref = eta_vee2_no_afe*(S_vee2_no_afe/S_ref);
        afe.eta_vee2_afe_S_vee2_afe_S_ref = eta_vee2_afe_S_vee2_afe_S_ref;
        afe.eta_vee2_no_afe_S_vee2_no_afe_S_ref = eta_vee2_no_afe_S_vee2_no_afe_S_ref;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Correction of Lift curve slope of vee2 associated to prop wash
        % eta_w1_afe = (q_w1_afe/q_inf)*(S_w1_afe/S_ref)
        % eta_w1_no_afe = (q_w1_no_afe/q_inf)*(S_w1_no_afe/S_ref)
        CLalpha_vee2_e_pw = eta_vee2_afe_S_vee2_afe_S_ref*CLalpha_vee2_e + eta_vee2_no_afe_S_vee2_no_afe_S_ref*CLalpha_vee2_e;
        CL0_vee2_e_corrected = eta_vee2_afe_S_vee2_afe_S_ref*CL0_vee2_e + eta_vee2_no_afe_S_vee2_no_afe_S_ref*CL0_vee2_e;
        Stab_Der_parts.CL0_vee2_e_corrected = CL0_vee2_e_corrected;
        Stab_Der_parts.CLalpha_vee2_e_pw = CLalpha_vee2_e_pw;
    end

    if VTP == 1
        
        % %         warndlg('En caso de aplicar esta corrección debido al propwash effect, se ha de revisar y cambiar la formulación de esta sección del código y adimensionalización en función del AC_type')
        %         S_VTP_pw = Geo_tier.S_VTP_pw;
        %         S_VTP_s = Geo_tier.S_VTP_s;
        %         eta_VTP = Effects.eta_VTP;
        %         CLalpha_VTP = Aero.CL_alpha_VTP_CR;
        %         CLalpha_VTP_e = CLalpha_VTP;
        %
        %
        %         if twin_VTP == 1
        %             S_VTP_afe = 2*S_VTP_pw;
        %             S_VTP_no_afe = S_VTP_s - S_VTP_afe;
        %             q_VTP_no_afe    = eta_VTP*q_inf;                    %presion dinamica no afectada
        %             if Posicion_Palanca == 0
        %                 V_VTP = V;
        %                 q_VTP_afe    = q_inf;
        %             else
        %                 V_VTP = sqrt(eta_VTP)*V + 2*v_i;
        %                 q_VTP_afe       = 0.5*rho*V_VTP^2;
        %             end
        %             eta_VTP_afe = (q_VTP_afe/q_inf);
        %             eta_VTP_no_afe = (q_VTP_no_afe/q_inf);
        %             afe.eta_VTP_afe = eta_VTP_afe;
        %             afe.eta_VTP_no_afe = eta_VTP_no_afe;
        %             eta_VTP_afe_S_VTP_afe_S_ref = eta_VTP_afe*(S_VTP_afe/S_ref);
        %             eta_VTP_no_afe_S_VTP_no_afe_S_ref = eta_VTP_no_afe*(S_VTP_no_afe/S_ref);
        %             afe.eta_VTP_afe_S_VTP_afe_S_ref = eta_VTP_afe_S_VTP_afe_S_ref;
        %             afe.eta_VTP_no_afe_S_VTP_no_afe_S_ref = eta_VTP_no_afe_S_VTP_no_afe_S_ref;
        %         else
        %             S_VTP_afe = S_VTP_pw;
        %             S_VTP_no_afe = S_VTP_s - S_VTP_afe;
        %             q_VTP_no_afe    = eta_VTP*q_inf;                    %presion dinamica no afectada
        %             if Posicion_Palanca == 0
        %                 V_VTP = V;
        %                 q_VTP_afe    = q_inf;
        %             else
        %                 V_VTP = sqrt(eta_VTP)*V + 2*v_i;
        %                 q_VTP_afe       = 0.5*rho*V_VTP^2;
        %             end
        %             eta_VTP_afe = (q_VTP_afe/q_inf);
        %             eta_VTP_no_afe = (q_VTP_no_afe/q_inf);
        %             eta_VTP_afe_S_VTP_afe_S_ref = eta_VTP_afe*(S_VTP_afe/S_ref);
        %             eta_VTP_no_afe_S_VTP_no_afe_S_ref = eta_VTP_no_afe*(S_VTP_no_afe/S_ref);
        %             afe.eta_VTP_afe_S_VTP_afe_S_ref = eta_VTP_afe_S_VTP_afe_S_ref;
        %             afe.eta_VTP_no_afe_S_VTP_no_afe_S_ref = eta_VTP_no_afe_S_VTP_no_afe_S_ref;
        %         end
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %% Correction of Lift curve slope of VTP1 associated to prop wash
        %         CLalpha_VTP_e_pw = eta_VTP_afe_S_VTP_afe_S_ref*CLalpha_VTP_e + eta_VTP_no_afe_S_VTP_no_afe_S_ref*CLalpha_VTP_e;
        %         CL0_VTP_e_corrected = eta_VTP_afe_S_VTP_afe_S_ref*CL0_w2_e + eta_VTP_no_afe_S_VTP_no_afe_S_ref*CL0_w2_e;
        %         Stab_Der_parts.CL0_VTP_e_corrected = CL0_VTP_e_corrected;
        %         Stab_Der_parts.CLalpha_VTP_e_pw = CLalpha_VTP_e_pw;
        
        if VTP == 1
            
            S_w1        = Geo_tier.S_w1;
            S_VTP       = Geo_tier.S_VTP;
            S_VTP_e     = Geo_tier.S_VTP_e;
            
            CLalpha_VTP   = Aero.CL_alpha_VTP_CR;
            CY_0_VTP_CR   = Aero.CY_0_VTP_CR;
            
            switch AC_type
                case 1
                    % Aerodynamic results conducted with the wing
                    %                 conv_VTP = S_w1_e/S_w1;
                    conv_VTP = 1; %ya adimensionalizados en flow5 con la Sref
                case 2
                    % If there is HTP the Aerodynamic studies in FLOW have been conducted
                    % HTP + VTP hence the area used for adimentionalizing is the one from the HTP
                    if HTP == 1
                        S_HTP_e = Geo_tier.S_w2_e;
                        S_HTP = Geo_tier.S_w2;
                        %                     conv_VTP = S_HTP_e/S_w1;
                        conv_VTP = S_HTP/S_w1; %adimensionalizados con la S_w2_total en flow5, no con la expuesta;
                    end
                case 3
                    % If there is HTP the Aerodynamic studies in FLOW have been conducted
                    % HTP + VTP hence the area used for adimentionalizing is the one from the HTP
                    if HTP == 1
                        S_HTP_e = Geo_tier.S_w2_e;
                        S_HTP = Geo_tier.S_w2;
                        %                     conv_VTP = S_HTP_e/S_w1;
                        conv_VTP = S_HTP/S_w1; %adimensionalizados con la S_w2_total en flow5, no con la expuesta;
                    end
                case 4
                    % No VTP
                case 5
                    % No VTP
                case 6
                    % Aerodynamic results conducted with the wing
                    %                 conv_VTP = S_w1_e/S_w1;
                    conv_VTP = 1; %ya adimensionalizados en flow5 con la Sref
                case 7
                    conv_VTP = 1; %ya adimensionalizados en flow5 con la Sref
                case 8
                    % No VTP
            end
            
            
            
            
            CLalpha_VTP_e_pw = CLalpha_VTP_e*conv_VTP;
            CL0_VTP_e_corrected = CY_0_VTP_CR_e*conv_VTP;
            
            Stab_Der_parts.CL0_VTP_e_corrected = CL0_VTP_e_corrected;
            Stab_Der_parts.CLalpha_VTP_e_pw = CLalpha_VTP_e_pw;
        end
    end
    
else
    v_i = 0;
    
    % NO PROP WASH CONTRIBUTION
    if W1 == 1
        CLalpha_w1_e_pw    = CLalpha_w1_e;
        CL0_w1_e_corrected = CL0_w1_e;
        
        Stab_Der_parts.CL0_w1_e_corrected = CL0_w1_e_corrected;
        Stab_Der_parts.CLalpha_w1_e_pw = CLalpha_w1_e_pw;
    end
    
    if Can == 1
        S_can = Geo_tier.S_can;
        S_can_e = Geo_tier.S_can_e;
        b_can = Geo_tier.b_can;

        CLalpha_can = Aero.CL_alpha_can_CR;
        CL0_can     = Aero.CL_0_can_CR;
        
        w_can = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_can);
        wingspan2bodydiam_can = b_can/w_can;
        
        if wingspan2bodydiam_can <= 2 || flagwingspan2bodydiam==1
            CLalpha_can_e = CLalpha_can*S_can_e/Geo_tier.S_can;
            CL0_can_e     = CL0_can*S_can_e/Geo_tier.S_can;
        else
            CLalpha_can_e = CLalpha_can;
            CL0_can_e     = CL0_can;
        end
        
        conv_can = S_can/S_w1;%adimensionalizados con la S_can total en flow5, no con la expuesta;
        
        CLalpha_can_e_pw = CLalpha_can_e*conv_can;
        CL0_can_e_corrected = CL0_can_e*conv_can;
        
        
        Stab_Der_parts.CL0_can_e_corrected = CL0_can_e_corrected;
        Stab_Der_parts.CLalpha_can_e_pw = CLalpha_can_e_pw;
    end
    
    if HTP == 1
        S_w2 = Geo_tier.S_w2;
        S_w2_e = Geo_tier.S_w2_e;
        b_w2 = Geo_tier.b_w2;
        CL0_w2 = Aero.CL_0_w2_CR;
        CL_alpha_w2 = Aero.CL_alpha_w2_CR;
        w_w2 = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_w2);
        wingspan2bodydiam_w2 = b_w2/w_w2;
        
        if wingspan2bodydiam_w2 <= 2 || flagwingspan2bodydiam==1
            CL0_w2_e     = CL0_w2*S_w2_e/Geo_tier.S_w2;
            CLalpha_w2_e = CL_alpha_w2*S_w2_e/Geo_tier.S_w2;
        else
            CL0_w2_e         = CL0_w2;
            CLalpha_w2_e     = CL_alpha_w2;
        end
        
        conv_w2 = S_w2/S_w1; %adimensionalizados con la S_w2_total en flow5, no con la expuesta;
        CLalpha_HTP_e_pw = CLalpha_w2_e*conv_w2;
        CL0_HTP_e_corrected = CL0_w2_e*conv_w2;
        
        Stab_Der_parts.CL0_HTP_e_corrected = CL0_HTP_e_corrected;
        Stab_Der_parts.CLalpha_HTP_e_pw = CLalpha_HTP_e_pw;
    end
    
    if Vee == 1
        S_vee = Geo_tier.S_vee;
        S_vee_e = Geo_tier.S_vee_e;
        b_vee = Geo_tier.b_vee;
        CL0_vee = Aero.CL_0_vee_CR;
        CL_alpha_wb_Vee = Aero.CLalpha_vee;
        w_vee = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee);
        wingspan2bodydiam_vee = b_vee/w_vee;
        
        if wingspan2bodydiam_vee <= 2 || flagwingspan2bodydiam==1
            CL0_vee_e = CL0_vee*S_vee_e/Geo_tier.S_vee;
            CLalpha_vee_e = CL_alpha_wb_Vee*S_vee_e/Geo_tier.S_vee;
        else
            CL0_vee_e         = CL0_vee;
            CLalpha_vee_e     = CL_alpha_wb_Vee;
        end
        
        conv_vee = S_vee/S_w1; %adimensionalizados con la S_vee_total en flow5, no con la expuesta;
        
        CLalpha_Vee_e_pw = CLalpha_vee_e*conv_vee;
        CL0_Vee_e_corrected = CL0_vee_e*conv_vee;
        
        Stab_Der_parts.CL0_Vee_e_corrected = CL0_Vee_e_corrected;
        Stab_Der_parts.CLalpha_Vee_e_pw = CLalpha_Vee_e_pw;
    end

    if Vee2 == 1
        S_vee2 = Geo_tier.S_vee2;
        S_vee2_e = Geo_tier.S_vee2_e;
        b_vee2 = Geo_tier.b_vee2;
        CL0_vee2 = Aero.CL_0_vee2_CR;
        CL_alpha_wb_vee2 = Aero.CLalpha_vee2;
        w_vee2 = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee2);
        wingspan2bodydiam_vee2 = b_vee2/w_vee2;

        if wingspan2bodydiam_vee2 <= 2 || flagwingspan2bodydiam==1
            CL0_vee2_e = CL0_vee2*S_vee2_e/Geo_tier.S_vee2;
            CLalpha_vee2_e = CL_alpha_wb_vee2*S_vee2_e/Geo_tier.S_vee2;
        else
            CL0_vee2_e         = CL0_vee2;
            CLalpha_vee2_e     = CL_alpha_wb_vee2;
        end

        conv_vee2 = S_vee2/S_w1; %adimensionalizados con la S_vee2_total en flow5, no con la expuesta;

        CLalpha_vee2_e_pw = CLalpha_vee2_e*conv_vee2;
        CL0_vee2_e_corrected = CL0_vee2_e*conv_vee2;

        Stab_Der_parts.CL0_vee2_e_corrected = CL0_vee2_e_corrected;
        Stab_Der_parts.CLalpha_vee2_e_pw = CLalpha_vee2_e_pw;
    end

    if VTP == 1
        
        S_VTP       = Geo_tier.S_VTP;
        S_VTP_e     = Geo_tier.S_VTP_e;
        
        CLalpha_VTP = Aero.CL_alpha_VTP_CR;
        CY_0_VTP_CR = Aero.CY_0_VTP_CR;
        CLalpha_VTP_e = CLalpha_VTP;
        CY_0_VTP_CR_e = CY_0_VTP_CR;
        switch AC_type
            case 1
                % Aerodynamic results conducted with the wing
                %                 conv_VTP = S_w1_e/S_w1;
                conv_VTP = 1; %ya adimensionalizados en flow5 con la Sref
            case 2
                % If there is HTP the Aerodynamic studies in FLOW have been conducted
                % HTP + VTP hence the area used for adimentionalizing is the one from the HTP
                if HTP == 1
                    S_HTP_e = Geo_tier.S_w2_e;
                    S_HTP = Geo_tier.S_w2;
                    %                     conv_VTP = S_HTP_e/S_w1;
                    conv_VTP = S_HTP/S_w1; %adimensionalizados con la S_w2_total en flow5, no con la expuesta;
                end
            case 3
                % If there is HTP the Aerodynamic studies in FLOW have been conducted
                % HTP + VTP hence the area used for adimentionalizing is the one from the HTP
                if HTP == 1
                    S_HTP_e = Geo_tier.S_w2_e;
                    S_HTP = Geo_tier.S_w2;
                    %                     conv_VTP = S_HTP_e/S_w1;
                    conv_VTP = S_HTP/S_w1; %adimensionalizados con la S_w2_total en flow5, no con la expuesta;
                end
            case 4
                % No VTP
            case 5
                % No VTP
            case 6
                % Aerodynamic results conducted with the wing
                %                 conv_VTP = S_w1_e/S_w1;
                conv_VTP = 1; %ya adimensionalizados en flow5 con la Sref
                           % No VTP
            case 7
                conv_VTP = 1; %ya adimensionalizados en flow5 con la Sref
            case 8
                % No VTP
         end
        
        
        CLalpha_VTP_e_pw = CLalpha_VTP_e*conv_VTP;
        CL0_VTP_e_corrected = CY_0_VTP_CR_e*conv_VTP;
        
        Stab_Der_parts.CL0_VTP_e_corrected = CL0_VTP_e_corrected;
        Stab_Der_parts.CLalpha_VTP_e_pw = CLalpha_VTP_e_pw;
    end
end


end