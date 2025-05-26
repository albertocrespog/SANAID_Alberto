function Effects = effects(AC_CONFIGURATION,Geo_tier,Design_criteria,conditions,Aero_TH,Aero,conv_UNITS)

W1  = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
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
    
z_w2_LE = Geo_tier.z_w2_LE;
l_xac_w1w2 = Geo_tier.l_xac_w1w2;
z_zbar_w2 = Geo_tier.z_zbar_w2;
x_xbar_w2 = Geo_tier.x_xbar_w2;


    %%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kh = (1-(z_w2_LE - z_w1_LE)/b_w1)/(2*l_xac_w1w2/b_w1)^(1/3);          %eq 3.43 PAMADI
    kl = (10-3*lambda_w1)/7;                                     %eq 3.43 PAMADI
    ka = 1/AR_w1_e -1/(1 + AR_w1_e^(1.7));                        %eq 3.43 PAMADI
    deps_dalpha_clean = 4.44*(kh*kl*ka*sqrt(cos(Lambda_c4_w1)))^1.19;  %eq 3.43 PAMADI
    deps_dalpha_poff = 0; % Assumes not affecting
    deps_dalpha_Vee = deps_dalpha_clean + deps_dalpha_poff;
    
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash = 1 - deps_dalpha_Vee;
    % eps_vee_0 = deps_dalpha*(alpha_CL_0_w1*D2R - i_w1);
%     eps_w2_0 = deps_dalpha_Vee*(alpha_CL_0_w1*D2R - i_w1);
    eps_w2_0 = 2*CL0_w1/(pi*AR_w1_e);
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_w2 = eps_w2_0 + deps_dalpha_Vee*alpha_f;
    Effects.downwash = downwash;
    Effects.eps_w2 = eps_w2;
    Effects.eps_w2_0 = eps_w2_0;
    Effects.deps_dalpha_Vee = deps_dalpha_Vee;
    %% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection
    
    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_w2 - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to w2 MAC
    l = x_xbar_w2 - x_cR_w1_TE; % x distance from the trailing edge of the wing to the w2 MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_w1 + atan(z_TE/l);
    gamma_w2 = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_w1 - eps_vee);
    x_w2_wake = a*cos(gamma_w2 - alpha_f - i_w1 - eps_w2);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    Deltaz_w2_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_w2_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    z_w2_wake = a*sin(gamma_w2 - alpha_f - i_w1 + eps_w2);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_w2_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_w2_wake)/(2*Deltaz_w2_wake)))^2)/((x_w2_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_w2_power = 0;
    % NOTE- need to take into account power effects
    % eta_vee = eta_vee_power_off + eta_vee_power
    eta_w2 = eta_w2_power_off + eta_w2_power;
    Effects.eta_w2 = eta_w2;
else
    deps_dalpha_Vee = 0;
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
%     downwash = 1 - deps_dalpha_Vee;
    eps_w2_0 = 0;
%     eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_w2 = eps_w2_0 + deps_dalpha_Vee*alpha_f;
%     Effects.eta_w2 = eta_w2;
%     Effects.downwash = downwash;
%     Effects.eps_w2 = eps_w2;
    Effects.deps_dalpha_Vee = deps_dalpha_Vee;
end

%% HTP

if HTP == 1

z_w2_LE    = Geo_tier.z_w2_LE;
l_xac_w1w2 = Geo_tier.l_xac_w1w2;
z_zbar_w2  = Geo_tier.z_zbar_w2;
x_xbar_w2  = Geo_tier.x_xbar_w2;

    %%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kh = (1-(z_w2_LE - z_w1_LE)/b_w1)/(2*l_xac_w1w2/b_w1)^(1/3);          %eq 3.43 PAMADI
    kl = (10-3*lambda_w1)/7;                                     %eq 3.43 PAMADI
    ka = 1/AR_w1_e -1/(1 + AR_w1_e^(1.7));                        %eq 3.43 PAMADI

    deps_dalpha_clean = 4.44*(kh*kl*ka*sqrt(cos(Lambda_c4_w1)))^1.19;  %eq 3.43 PAMADI
    deps_dalpha_poff  = 0; % Assumes not affecting
    deps_dalpha_h     = deps_dalpha_clean + deps_dalpha_poff;
    
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash = 1 - deps_dalpha_h;
    % eps_vee_0 = deps_dalpha*(alpha_CL_0_w1*D2R - i_w1);
%     eps_w2_0 = deps_dalpha_h*(alpha_CL_0_w1*D2R - i_w1);
    eps_w2_0 = 2*CL0_w1/(pi*AR_w1_e);
    
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_w2 = eps_w2_0 + deps_dalpha_h*alpha_f;
    Effects.downwash = downwash;
    Effects.eps_w2 = eps_w2;
    Effects.eps_w2_0 = eps_w2_0;
    Effects.deps_dalpha_h = deps_dalpha_h;
    
    %% NOTE eps_vee_0 missing The V-Tail downwash angle increment due to flap deflection
    
    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_w2 - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to w2 MAC
    l = x_xbar_w2 - x_cR_w1_TE; % x distance from the trailing edge of the wing to the w2 MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_w1 + atan(z_TE/l);
    gamma_w2 = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_w1 - eps_vee);
    x_w2_wake = a*cos(gamma_w2 - alpha_f - i_w1 - eps_w2);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    Deltaz_w2_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_w2_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    z_w2_wake = a*sin(gamma_w2 - alpha_f - i_w1 + eps_w2);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_w2_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_w2_wake)/(2*Deltaz_w2_wake)))^2)/((x_w2_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_w2_power = 0;
    
    % NOTE- need to take into account power effects
    % eta_vee = eta_vee_power_off + eta_vee_power
    eta_w2 = eta_w2_power_off + eta_w2_power;
    Effects.eta_w2 = eta_w2;
elseif Vee == 1
    deps_dalpha_h = 0;
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    downwash = 1 - deps_dalpha_h;
    eps_w2_0 = 0;
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_w2 = eps_w2_0 + deps_dalpha_h*alpha_f;
    Effects.eta_w2 = eta_w2;
    Effects.downwash = downwash;
    Effects.eps_w2 = eps_w2;
    Effects.deps_dalpha_h = deps_dalpha_h;
else
    deps_dalpha_h = 0;
    Effects.deps_dalpha_h = deps_dalpha_h;
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
    z_TE = z_zbar_can - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to w2 MAC
    l = x_xbar_can - x_cR_w1_TE; % x distance from the trailing edge of the wing to the w2 MAC
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
    % NOTE deps_dalpha missing The downwash gradient at the V-Tail including flap effects
    upwash = 1 + deps_dalpha_can;
    eps_can_0 = 0;
    % eps_vee = eps_vee_0 + deps_dalpha*alpha_f;
    eps_can = eps_can_0 + deps_dalpha_can*alpha_f;
%     Effects.eta_can = eta_can;
    Effects.upwash = upwash;
    Effects.eps_can = eps_can;
    Effects.deps_dalpha_can = deps_dalpha_can;

end

%% VTP

if VTP == 1
    x_xbar_VTP = Geo_tier.x_xbar_VTP;
    z_zbar_VTP = Geo_tier.z_zbar_VTP;
    
    %% Determination of pressure dynamic ratios
    % Determination of dynamic pressure ratio Vee tail
    % The distance between the wing trailing edge and the V-Tail aerodynamic center
    z_TE = z_zbar_VTP - z_cR_w1_TE + 0.75*cR_w1*tan(i_w1); %z-distance from the TE w1 root chord to w2 MAC
    l = x_xbar_VTP - x_cR_w1_TE; % x distance from the trailing edge of the wing to the w2 MAC
    a = sqrt(z_TE^2 + l^2);
    % The angle between the wing chord plane and the line connecting the wing root chord trailing edge and the V-Tail aerodynamic center
    % gamma_vee = i_w1 + atan(z_TE/l);
    gamma_VTP = i_w1 + atan(z_TE/l);
    % The distance between the wing trailing edge and the V-Tail aerodynamic center along the centerline of the wake
    % x_vee_wake = acos(gamma_vee - alpha_f - i_w1 - eps_vee);
    x_VTP_wake = acos(gamma_VTP - alpha_f - i_w1 - eps_w2);
    % The half-width of the wing wake perpendicular to the centerline of the wake
    % Deltaz_vee_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_vee_wake/cmac_w1 + 0.15));
    Deltaz_VTP_wake = 0.68*cmac_w1*sqrt(CD0_w1*(x_VTP_wake/cmac_w1 + 0.15));
    % The perpendicular distance from the centerline of the wake to the V-Tail aerodynamic center
    % z_vee_wake = a*sin(gamma_vee - alpha_f - i_w1 + eps_vee);
    z_VTP_wake = a*sin(gamma_VTP - alpha_f - i_w1 + eps_w2);
    % eta_vee_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_vee_wake)/(2*Deltaz_vee_wake)))^2)/((x_vee_wake/cmac_w1)+0.3);
    eta_VTP_power_off = 1 - (2.42*sqrt(CD0_w1)*(cos((pi*z_VTP_wake)/(2*Deltaz_VTP_wake)))^2)/((x_VTP_wake/cmac_w1)+0.3);
    % eta_vee_power = 0;
    eta_VTP_power = 0;
    
    % NOTE- need to take into account power effects
    eta_VTP = eta_VTP_power_off + eta_VTP_power;
    Effects.eta_VTP = eta_VTP;
end

end
