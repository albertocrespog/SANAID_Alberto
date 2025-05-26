function Effects = effects_v3(AC_CONFIGURATION,Geo_tier,Design_criteria,conditions,Aero_TH,Aero,conv_UNITS)

W1  = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;


z_cR_w1_TE    = Geo_tier.z_cR_w1_TE;
z_w1_LE       = Geo_tier.z_w1_LE;
x_cR_w1_LE    = Geo_tier.x_cR_w1_LE;
x_cR_w1_TE    = Geo_tier.x_cR_w1_TE;
cR_w1         = Geo_tier.cR_w1;
cmac_w1       = Geo_tier.cmac_w1;
AR_w1_e       = Geo_tier.AR_w1_e;
AR_w1         = Geo_tier.AR_w1;
b_w1          = Geo_tier.b_w1;
lambda_w1     = Geo_tier.lambda_w1;
Lambda_c4_w1  = Geo_tier.Lambda_c4_w1;
i_w1          = Design_criteria.i_w1;
alpha_f       = conditions.alpha_f;
CD0_w1        = Aero_TH.CD0_w1;
alpha_CL_0_w1 = Aero.alpha_CL_0_w1_CR;
D2R           = conv_UNITS.D2R;
CL0_w1        = Aero.CL_0_w1_CR;


%% COLA EN V
if Vee == 1

    z_vee_LE = Geo_tier.z_vee_LE;
    l_xac_w1vee = Geo_tier.l_xac_w1vee;
    z_zbar_vee = Geo_tier.z_zbar_vee;
    x_xbar_vee = Geo_tier.x_xbar_vee;

    %%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kh = (1-(z_vee_LE - z_w1_LE)/b_w1)/(2*l_xac_w1vee/b_w1)^(1/3);          %eq 3.43 PAMADI
    kl = (10-3*lambda_w1)/7;                                     %eq 3.43 PAMADI
    ka = 1/AR_w1_e -1/(1 + AR_w1_e^(1.7));                        %eq 3.43 PAMADI
    deps_dalpha_clean = 4.44*(kh*kl*ka*sqrt(cos(Lambda_c4_w1)))^1.19;  %eq 3.43 PAMADI
    deps_dalpha_poff = 0; % Assumes not affecting
    deps_dalpha_vee = deps_dalpha_clean + deps_dalpha_poff;

    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash = 1 - deps_dalpha_vee;
    % eps_vee_0 = deps_dalpha*(alpha_CL_0_w1*D2R - i_w1);
    %     eps_vee_0 = deps_dalpha_vee*(alpha_CL_0_w1*D2R - i_w1);
    eps_vee_0 = 2*CL0_w1/(pi*AR_w1_e);
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_vee = eps_vee_0 + deps_dalpha_vee*alpha_f;
    Effects.downwash_vee = downwash;
    Effects.eps_vee = eps_vee;
    Effects.eps_vee_0 = eps_vee_0;
    Effects.deps_dalpha_vee = deps_dalpha_vee;
    %% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection

    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_vee - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to vee MAC
    l = x_xbar_vee - x_cR_w1_TE; % x distance from the trailing edge of the wing to the vee MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_w1 + atan(z_TE/l);
    gamma_vee = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_w1 - eps_vee);
    x_vee_wake = a*cos(gamma_vee - alpha_f - i_w1 - eps_vee);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_vee_power = 0;
    % NOTE- need to take into account power effects
    % eta_vee = eta_vee_power_off + eta_vee_power
    eta_vee = eta_vee_power_off + eta_vee_power;
    Effects.eta_vee = eta_vee;
else
    deps_dalpha_vee = 0;
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash_vee = 1 - deps_dalpha_vee;
    eps_vee_0 = 0;
    eps_vee = eps_vee_0 + deps_dalpha_vee*alpha_f;
    eta_vee = 1;
    Effects.downwash_vee = downwash_vee;
    Effects.eps_vee = eps_vee;
    Effects.eps_vee_0 = eps_vee_0;
    Effects.deps_dalpha_vee = deps_dalpha_vee;
    Effects.eta_vee = eta_vee;
end

%% COLA EN V2
if Vee2 == 1

    z_vee2_LE = Geo_tier.z_vee2_LE;
    l_xac_w1vee2 = Geo_tier.l_xac_w1vee2;
    z_zbar_vee2 = Geo_tier.z_zbar_vee2;
    x_xbar_vee2 = Geo_tier.x_xbar_vee2;

    %%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kh = (1-(z_vee2_LE - z_w1_LE)/b_w1)/(2*l_xac_w1vee2/b_w1)^(1/3);          %eq 3.43 PAMADI
    kl = (10-3*lambda_w1)/7;                                     %eq 3.43 PAMADI
    ka = 1/AR_w1_e -1/(1 + AR_w1_e^(1.7));                        %eq 3.43 PAMADI
    deps_dalpha_clean = 4.44*(kh*kl*ka*sqrt(cos(Lambda_c4_w1)))^1.19;  %eq 3.43 PAMADI
    deps_dalpha_poff = 0; % Assumes not affecting
    deps_dalpha_vee2 = deps_dalpha_clean + deps_dalpha_poff;

    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash_vee2 = 1 - deps_dalpha_vee2;
    % eps_vee2_0 = deps_dalpha*(alpha_CL_0_w1*D2R - i_w1);
    %     eps_vee2_0 = deps_dalpha_vee2*(alpha_CL_0_w1*D2R - i_w1);
    eps_vee2_0 = 2*CL0_w1/(pi*AR_w1_e);
    % eps_vee2 = eps_vee2_0 + deps_dalpha*alpha_f;
    eps_vee2 = eps_vee2_0 + deps_dalpha_vee2*alpha_f;
    Effects.downwash_vee2 = downwash_vee2;
    Effects.eps_vee2 = eps_vee2;
    Effects.eps_vee2_0 = eps_vee2_0;
    Effects.deps_dalpha_vee2 = deps_dalpha_vee2;
    %% NOTE eps_vee2_0 missing The V-Tail downwash angle increment due to flap deflection

    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio vee2 tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_vee2 - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to vee2 MAC
    l = x_xbar_vee2 - x_cR_w1_TE; % x distance from the trailing edge of the wing to the vee2 MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee2 = i_w1 + atan(z_TE/l);
    gamma_vee2 = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee2_wake = acos(gamma_vee2 - alpha_f - i_w1 - eps_vee2);
    x_vee2_wake = a*cos(gamma_vee2 - alpha_f - i_w1 - eps_vee2);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee2_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee2_wake/cmac_w1 + 0.15));
    Deltaz_vee2_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee2_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee2_wake = a*sin(gamma_vee2 - alpha_f - i_w1 + eps_vee2);
    z_vee2_wake = a*sin(gamma_vee2 - alpha_f - i_w1 + eps_vee2);
    % eta_vee2_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee2_wake)/(2*Deltaz_vee2_wake)))^2)/((x_vee2_wake/cmac_w1)+0.3);
    eta_vee2_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee2_wake)/(2*Deltaz_vee2_wake)))^2)/((x_vee2_wake/cmac_w1)+0.3);
    % eta_vee2_power = 0;
    eta_vee2_power = 0;
    % NOTE- need to take into account power effects
    % eta_vee2 = eta_vee2_power_off + eta_vee2_power
    eta_vee2 = eta_vee2_power_off + eta_vee2_power;
    Effects.eta_vee2 = eta_vee2;
else
    deps_dalpha_vee2 = 0;
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash_vee2 = 1 - deps_dalpha_vee2;
    eps_vee2_0 = 0;
    eps_vee2 = eps_vee2_0 + deps_dalpha_vee2*alpha_f;
    eta_vee2 = 1;
    Effects.downwash_vee2 = downwash_vee2;
    Effects.eps_vee2 = eps_vee2;
    Effects.eps_vee2_0 = eps_vee2_0;
    Effects.deps_dalpha_vee2 = deps_dalpha_vee2;
    Effects.eta_vee2 = eta_vee2;
end

%% HTP

if HTP == 1

    z_HTP_LE    = Geo_tier.z_HTP_LE;
    l_xac_w1HTP = Geo_tier.l_xac_w1HTP;
    z_zbar_HTP  = Geo_tier.z_zbar_HTP;
    x_xbar_HTP  = Geo_tier.x_xbar_HTP;

    %%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kh = (1-(z_HTP_LE - z_w1_LE)/b_w1)/(2*l_xac_w1HTP/b_w1)^(1/3);          %eq 3.43 PAMADI
    kl = (10-3*lambda_w1)/7;                                     %eq 3.43 PAMADI
    ka = 1/AR_w1_e -1/(1 + AR_w1_e^(1.7));                        %eq 3.43 PAMADI

    deps_dalpha_clean = 4.44*(kh*kl*ka*sqrt(cos(Lambda_c4_w1)))^1.19;  %eq 3.43 PAMADI
    deps_dalpha_poff  = 0; % Assumes not affecting
    deps_dalpha_HTP     = deps_dalpha_clean + deps_dalpha_poff;

    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash_HTP = 1 - deps_dalpha_HTP;
    % eps_vee_0 = deps_dalpha*(alpha_CL_0_w1*D2R - i_w1);
    %     eps_HTP_0 = deps_dalpha_HTP*(alpha_CL_0_w1*D2R - i_w1);
    eps_HTP_0 = 2*CL0_w1/(pi*AR_w1_e);

    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_HTP = eps_HTP_0 + deps_dalpha_HTP*alpha_f;
    Effects.downwash_HTP = downwash_HTP;
    Effects.eps_HTP = eps_HTP;
    Effects.eps_HTP_0 = eps_HTP_0;
    Effects.deps_dalpha_HTP = deps_dalpha_HTP;

    %% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection

    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_HTP - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to HTP MAC
    l = x_xbar_HTP - x_cR_w1_TE; % x distance from the trailing edge of the wing to the HTP MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_w1 + atan(z_TE/l);
    gamma_HTP = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_w1 - eps_vee);
    x_HTP_wake = a*cos(gamma_HTP - alpha_f - i_w1 - eps_HTP);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    Deltaz_HTP_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_HTP_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    z_HTP_wake = a*sin(gamma_HTP - alpha_f - i_w1 + eps_HTP);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_HTP_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_HTP_wake)/(2*Deltaz_HTP_wake)))^2)/((x_HTP_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_HTP_power = 0;

    % NOTE- need to take into account power effects
    % eta_vee = eta_vee_power_off + eta_vee_power
    eta_HTP = eta_HTP_power_off + eta_HTP_power;
    Effects.eta_HTP = eta_HTP;
% elseif Vee == 1
%     deps_dalpha_HTP = 0;
%     % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
%     downwash = 1 - deps_dalpha_HTP;
%     eps_vee_0 = 0;
%     % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
%     eps_vee = eps_vee_0 + deps_dalpha_HTP*alpha_f;
%     Effects.eta_vee = eta_vee;
%     Effects.downwash = downwash;
%     Effects.eps_vee = eps_vee;
%     Effects.deps_dalpha_HTP = deps_dalpha_HTP;
% elseif Vee2 == 1
%     deps_dalpha_HTP = 0;
%     % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
%     downwash = 1 - deps_dalpha_HTP;
%     eps_vee2_0 = 0;
%     % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
%     eps_vee2 = eps_vee2_0 + deps_dalpha_HTP*alpha_f;
%     Effects.eta_vee2 = eta_vee2;
%     Effects.downwash = downwash;
%     Effects.eps_vee2 = eps_vee2;
%     Effects.deps_dalpha_HTP = deps_dalpha_HTP;
else
    deps_dalpha_HTP = 0;
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash_HTP = 1 - deps_dalpha_HTP;
    eps_HTP_0 = 0;
    eps_HTP = eps_HTP_0 + deps_dalpha_HTP*alpha_f;
    eta_HTP = 1;
    Effects.downwash_HTP = downwash_HTP;
    Effects.eps_HTP = eps_HTP;
    Effects.eps_HTP_0 = eps_HTP_0;
    Effects.deps_dalpha_HTP = deps_dalpha_HTP;
    Effects.eta_HTP = eta_HTP;


end

%% Canard
if Can == 1
    x_xbar_can = Geo_tier.x_xbar_can;
    z_zbar_can = Geo_tier.z_zbar_can;
    x_1R_y1_can = Geo_tier.x_1R_y1_can;
    cR_can = Geo_tier.cR_can;
    
    X_w = (x_cR_w1_LE + cR_w1/4);
    X_c = (x_1R_y1_can + cR_can/4);
    %%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL UPWASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deps_dalpha_clean = upwash_calc(AR_w1_e, X_w, X_c, cR_w1);
    deps_dalpha_poff = 0; % Assumes not affecting
    deps_dalpha_can = deps_dalpha_clean + deps_dalpha_poff;
    
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    upwash = 1 + deps_dalpha_can;
%     eps_can_0 = deps_dalpha_can*(alpha_CL_0_w1*D2R - i_w1);
    eps_can_0 = 2*CL0_w1/(pi*AR_w1_e);
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_can = eps_can_0 + deps_dalpha_can*alpha_f;
    Effects.upwash = upwash;
    Effects.eps_can = eps_can;
    Effects.deps_dalpha_can = deps_dalpha_can;

    %% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection
    
    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_can - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to HTP MAC
    l = x_xbar_can - x_cR_w1_TE; % x distance from the trailing edge of the wing to the HTP MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the Canard aerodynamic center
    gamma_can = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    x_can_wake = a*cos(gamma_can - alpha_f - i_w1 - eps_can);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    Deltaz_can_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_can_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    z_can_wake = a*sin(gamma_can - alpha_f - i_w1 + eps_can);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_can_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_can_wake)/(2*Deltaz_can_wake)))^2)/((x_can_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_can_power = 0;
    
    % NOTE- need to take into account power effects
    eta_can = eta_can_power_off + eta_can_power;
    Effects.eta_can = eta_can;
else
    deps_dalpha_can = 0;
    % NOTE deps_dalpha missing The downwash gradient at the can including flap effects
    upwash_can = 1 + deps_dalpha_can;
    eps_can_0 = 0;
    eps_can = eps_can_0 + deps_dalpha_can*alpha_f;
    eta_can = 1;
    Effects.upwash_can = upwash_can;
    Effects.eps_can = eps_can;
    Effects.eps_can_0 = eps_can_0;
    Effects.deps_dalpha_can = deps_dalpha_can;
    Effects.eta_can = eta_can;

end


if Can == 1 && W1 == 1
    % For the wing

    % Data from the Can to calculate downwash on wing
    z_cR_can_TE    = Geo_tier.z_cR_can_TE;
    z_can_LE       = Geo_tier.z_can_LE;
    x_cR_can_LE    = Geo_tier.x_cR_can_LE;
    x_cR_can_TE    = Geo_tier.x_cR_can_TE;
    cR_can         = Geo_tier.cR_can;
    cmac_can       = Geo_tier.cmac_can;
    AR_can_e       = Geo_tier.AR_can_e;
    AR_can         = Geo_tier.AR_can;
    b_can          = Geo_tier.b_can;
    lambda_can     = Geo_tier.lambda_can;
    Lambda_c4_can  = Geo_tier.Lambda_c4_can;
    i_can          = Design_criteria.i_can;
    alpha_f       = conditions.alpha_f;
    CD0_can        = Aero_TH.CD0_can;
    alpha_CL_0_can = Aero.alpha_CL_0_can_CR;
    CL0_can        = Aero.CL_0_can_CR;

    z_w1_LE    = Geo_tier.z_w1_LE;
    l_xac_canw1 = Geo_tier.l_xac_canw1;
    z_zbar_w1  = Geo_tier.z_zbar_w1;
    x_xbar_w1  = Geo_tier.x_xbar_w1;

    %%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kh = (1-(z_w1_LE - z_can_LE)/b_can)/(2*l_xac_canw1/b_can)^(1/3);          %eq 3.43 PAMADI
    kl = (10-3*lambda_can)/7;                                     %eq 3.43 PAMADI
    ka = 1/AR_can_e -1/(1 + AR_can_e^(1.7));                        %eq 3.43 PAMADI

    deps_dalpha_clean = 4.44*(kh*kl*ka*sqrt(cos(Lambda_c4_can)))^1.19;  %eq 3.43 PAMADI
    deps_dalpha_poff  = 0; % Assumes not affecting
    deps_dalpha_w1     = deps_dalpha_clean + deps_dalpha_poff;

    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash_w1 = 1 - deps_dalpha_w1;
    % eps_vee_0 = deps_dalpha*(alpha_CL_0_can*D2R - i_can);
    %     eps_w1_0 = deps_dalpha_w1*(alpha_CL_0_can*D2R - i_can);
    eps_w1_0 = 2*CL0_can/(pi*AR_can_e);

    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_w1 = eps_w1_0 + deps_dalpha_w1*alpha_f;
    Effects.downwash_w1 = downwash_w1;
    Effects.eps_w1 = eps_w1;
    Effects.eps_w1_0 = eps_w1_0;
    Effects.deps_dalpha_w1 = deps_dalpha_w1;

    %% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection

    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_w1 - z_cR_can_TE + 0.75*cR_can*tan(i_can); %z-distance from the TE can root chord to w1 MAC
    l = x_xbar_w1 - x_cR_can_TE; % x distance from the trailing edge of the wing to the w1 MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_can + atan(z_TE/l);
    gamma_w1 = i_can + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_can - eps_vee);
    x_w1_wake = a*cos(gamma_w1 - alpha_f - i_can - eps_w1);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_can*sqrt(CD0_can*(x_vee_wake/cmac_can + 0.15));
    Deltaz_w1_wake = 0.68*cmac_can*sqrt(CD0_can*(x_w1_wake/cmac_can + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_can + eps_vee);
    z_w1_wake = a*sin(gamma_w1 - alpha_f - i_can + eps_w1);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_can)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_can)+0.3);
    eta_w1_power_off = 1 - (2.42*sqrt(CD0_can)*(cos((pi*z_w1_wake)/(2*Deltaz_w1_wake)))^2)/((x_w1_wake/cmac_can)+0.3);
    % eta_vee_power = 0;
    eta_w1_power = 0;

    % NOTE- need to take into account power effects
    % eta_vee = eta_vee_power_off + eta_vee_power
    eta_w1 = eta_w1_power_off + eta_w1_power;
    Effects.eta_w1 = eta_w1;
    Effects.deps_dalpha_w1 = deps_dalpha_w1;
else
    deps_dalpha_w1 = 0;
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash_w1 = 1 - deps_dalpha_w1;
    eps_w1_0 = 0;
    eps_w1 = eps_w1_0 + deps_dalpha_w1*alpha_f;
    eta_w1 = 1;
    Effects.downwash_w1 = downwash_w1;
    Effects.eps_w1 = eps_w1;
    Effects.eps_w1_0 = eps_w1_0;
    Effects.deps_dalpha_w1 = deps_dalpha_w1;
    Effects.eta_w1 = eta_w1;
end


%% VTP
if VTP == 1
    x_xbar_VTP = Geo_tier.x_xbar_VTP;
    z_zbar_VTP = Geo_tier.z_zbar_VTP;
    
    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_VTP - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to HTP MAC
    l = x_xbar_VTP - x_cR_w1_TE; % x distance from the trailing edge of the wing to the HTP MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_w1 + atan(z_TE/l);
    gamma_VTP = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_w1 - eps_vee);
    x_VTP_wake = acos(gamma_VTP - alpha_f - i_w1 - eps_HTP);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    Deltaz_VTP_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_VTP_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    z_VTP_wake = a*sin(gamma_VTP - alpha_f - i_w1 + eps_HTP);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_VTP_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_VTP_wake)/(2*Deltaz_VTP_wake)))^2)/((x_VTP_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_VTP_power = 0;
    
    % NOTE- need to take into account power effects
    eta_VTP = eta_VTP_power_off + eta_VTP_power;
    Effects.eta_VTP = eta_VTP;
else 
    deps_dalpha_VTP = 0;
    % NOTE deps_dalpha missing The downwash gradient at the VTP including flap effects
    downwash_VTP = 1 - deps_dalpha_VTP;
    eps_VTP_0 = 0;
    eps_VTP = eps_VTP_0 + deps_dalpha_VTP*alpha_f;
    eta_VTP = 1;
    Effects.downwash_VTP = downwash_VTP;
    Effects.eps_VTP = eps_VTP;
    Effects.eps_VTP_0 = eps_VTP_0;
    Effects.deps_dalpha_VTP = deps_dalpha_VTP;
    Effects.eta_VTP = eta_VTP;

end

end
