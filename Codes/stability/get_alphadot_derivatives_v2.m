function [Stab_Der,Stab_Der_parts] = get_alphadot_derivatives_v2(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,Effects,...
    Aero_TH,OUTPUT_read_XLSX)


%% Input
W1  = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

prop_wash_effect = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;

f_k2_k1 = Stab_Der_parts.f_k2_k1;
Area_b_max = Body_Geo.Area_b_max;
Vol_TOT = Body_Geo.Vol_TOT;
length_fus = Body_Geo.l_fus;
x_Area_body = Body_Geo.x_Area_body;
Area_body = Body_Geo.Area_body;
w_Area_b_max = Body_Geo.w_Area_b_max;
x_XCG = conditions.x_XCG;
S_w1 = Geo_tier.S_w1;
AR_w1 = Geo_tier.AR_w1;
S_w1_e = Geo_tier.S_w1_e;
S_ref = Geo_tier.S_ref;
cmac_w1 = Geo_tier.cmac_w1;
lambda_w1 = Geo_tier.lambda_w1;
cmac_w1_e = Geo_tier.cmac_w1_e;
AR_w1_e = Geo_tier.AR_w1_e;
b_w1 = Geo_tier.b_w1;
xbar_w1 = Geo_tier.xbar_w1;
x_xbar_w1_e = Geo_tier.x_xbar_w1_e;
xbar_w1_e = Geo_tier.xbar_w1_e;
V = conditions.V;
a = Performance.a;
Mach = V/a;
CLalpha_w1_e = Stab_Der_parts.CL_alpha_w1;
cR_w1 = Geo_tier.cR_w1;
Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Alpha punto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% CL_alphapunto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuselage contribution
CLalpha_B = 2*f_k2_k1*(Area_b_max/Vol_TOT^(2/3));      %eq 4.532 PAMADI
CLalpha_B_p = CLalpha_B*(Vol_TOT^(2/3)/(Area_b_max));    %eq 4.531 PAMADI
CL_alphapunto_B = 2*CLalpha_B_p*(Vol_TOT/(length_fus*Area_b_max));  %eq 4.530 PAMADI
CL_alphapunto_fus = CL_alphapunto_B*((length_fus*Area_b_max)/(S_ref*cmac_w1));      %eq 4.527 PAMADI

if W1 == 1
    beta = sqrt(1 - Mach^2);
    tau = beta*AR_w1_e;

    % Only valid for tau < 4
    if tau< 4
        CL_g = (-pi*AR_w1_e/(2*beta^2))*(0.0013*tau^4 - 0.0122*tau^3 + ...
            0.0317*tau^2 + 0.0186*tau - 0.0004);                %eq 4.529 PAMADI
    else
        CL_g = 0;
    end
    %% OJO! using he distance from apex of wing LE root to Xac
    CL_alphapunto_e1 = 1.5*(xbar_w1_e/cR_w1)*CLalpha_w1_e; %eq 4.528 PAMADI
    CL_alphapunto_e2 = 3*CL_g;                              %eq 4.528 PAMADI
    % CL_alphapunto_e2 = 0; %correccion
    CL_alphapunto_e = CL_alphapunto_e1 + CL_alphapunto_e2;  %eq 4.528 PAMADI
    % W1 and w2 with no dynamic pressure correction
    % Ojo revisar por la contribución de la presión dinámica
    KWB_w1 = 0.1714*(w_Area_b_max/b_w1)^2 + 0.8326*(w_Area_b_max/b_w1) + 0.9974; %eq 3.27 PAMADI
    KBW_w1 = 0.7810*(w_Area_b_max/b_w1)^2 + 1.1976*(w_Area_b_max/b_w1) + 0.0088; %eq 3.28 PAMADI
    CL_alphapunto_w1 = (KWB_w1 + KBW_w1)*(cmac_w1_e/cmac_w1)*CL_alphapunto_e;    %eq 4.527 PAMADI
    % Dyanmic pressure correction
    if prop_wash_effect == 1
        eta_w1_afe_S_w1_afe_S_ref = afe.eta_w1_afe_S_w1_afe_S_ref;
        eta_w1_no_afe_S_w1_no_afe_S_ref = afe.eta_w1_no_afe_S_w1_no_afe_S_ref;
        CL_alphapunto_w1_pw = eta_w1_afe_S_w1_afe_S_ref*CL_alphapunto_w1 +...
            eta_w1_no_afe_S_w1_no_afe_S_ref*CL_alphapunto_w1;
    else
        CL_alphapunto_w1_pw = CL_alphapunto_w1;
    end    
else
    CL_alphapunto_w1_pw = 0;
end

if HTP == 1
    deps_dalpha_h = Effects.deps_dalpha_h;
    % Reference wing
    S_w2_e = Geo_tier.S_w2_e;
    S_w2 = Geo_tier.S_w2;
    b_w2 = Geo_tier.b_w2;
    conv_w2 = S_w2/S_w1; %adimensionalizados con la S_w2_total en flow5, no con la expuesta;
    w_w2 = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_w2);
    wingspan2bodydiam_w2 = b_w2/w_w2;
    x_xbar_w2 = Geo_tier.x_xbar_w2;
    %       Vbar_h = ((x_xbar_w2 - x_XCG)/cmac_w1)*(S_w2_e/S_ref);
    Vbar_h = ((x_xbar_w2 - x_XCG))/(cmac_w1); %para no mult dos veces *(S_w2_e/S_ref)
    CL_alpha_w2 = Aero.CL_alpha_w2_CR;
    
    if wingspan2bodydiam_w2 <= 2
        CLalpha_w2_e = CL_alpha_w2*S_w2_e/Geo_tier.S_w2;
    else
        CLalpha_w2_e = CL_alpha_w2;
    end
    
    % Modified version Pamadi to account for different dynamic pressure
    % W1 and w2 with no dynamic pressure correction
    CL_alphapunto_w2 = 2*CLalpha_w2_e*Vbar_h*deps_dalpha_h; %eq 4.525 PAMADI
    % Revisar OJO - contribuición de presiópn dinámica en ala w2
    
    % Dyanmic pressure correction
    if prop_wash_effect == 1
        eta_w2_afe_S_w2_afe_S_ref = afe.eta_w2_afe_S_w2_afe_S_ref;
        eta_w2_no_afe_S_w2_no_afe_S_ref = afe.eta_w2_no_afe_S_w2_no_afe_S_ref;
        CL_alphapunto_HTP_pw = eta_w2_afe_S_w2_afe_S_ref*CL_alphapunto_w2 + ...
            eta_w2_no_afe_S_w2_no_afe_S_ref*CL_alphapunto_w2;
    else
        CL_alphapunto_HTP_pw = CL_alphapunto_w2*conv_w2;
    end
else
    CL_alphapunto_HTP_pw = 0;
end

if Vee == 1
    % Reference wing
    S_vee_e = Geo_tier.S_vee_e;
    S_vee = Geo_tier.S_vee;
    b_vee = Geo_tier.b_vee;
    w_vee = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee);
    wingspan2bodydiam_vee = b_vee/w_vee;

    conv_vee = S_vee/S_w1; %adimensionalizados con la S_vee_total en flow5, no con la expuesta;
    
    x_xbar_vee = Geo_tier.x_xbar_vee;
    %         Vbar_vee = ((x_xbar_vee - x_XCG)/cmac_w1)*(S_vee_s/S_ref);
    Vbar_vee = ((x_xbar_vee - x_XCG)/cmac_w1);%para no mult dos veces *(S_vee_e/S_ref)
    CL_alpha_wb_Vee = Aero.CL_alpha_vee_CR;
    
    if wingspan2bodydiam_vee <= 2
        CLalpha_vee_e = CL_alpha_wb_Vee*S_vee_e/Geo_tier.S_vee;
    else
        CLalpha_vee_e = CL_alpha_wb_Vee;
    end    
    
    deps_dalpha_Vee = Effects.deps_dalpha_Vee;
    % W1 and vee with no dynamic pressure correction
    CL_alphapunto_vee = 2*CLalpha_vee_e*Vbar_vee*deps_dalpha_Vee; %eq 4.525 PAMADI
    
    % Dyanmic pressure correction
    if prop_wash_effect == 1
        eta_vee_afe_S_vee_afe_S_ref = afe.eta_vee_afe_S_vee_afe_S_ref;
        eta_vee_no_afe_S_vee_no_afe_S_ref = afe.eta_vee_no_afe_S_vee_no_afe_S_ref;
        CL_alphapunto_Vee_pw = eta_vee_afe_S_vee_afe_S_ref*CL_alphapunto_vee + ...
            eta_vee_no_afe_S_vee_no_afe_S_ref*CL_alphapunto_vee;
    else
        CL_alphapunto_Vee_pw = CL_alphapunto_vee*conv_vee;
    end
    % Revisar OJO - contribuición de presiópn dinámica en ala vee
else
    CL_alphapunto_Vee_pw = 0;
end

if Vee2 == 1
    % Reference wing
    S_vee2_e = Geo_tier.S_vee2_e;
    S_vee2 = Geo_tier.S_vee2;
    b_vee2 = Geo_tier.b_vee2;
    w_vee2 = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee2);
    wingspan2bodydiam_vee2 = b_vee2/w_vee2;

    conv_vee2 = S_vee2/S_w1; %adimensionalizados con la S_vee2_total en flow5, no con la expuesta;
    
    x_xbar_vee2 = Geo_tier.x_xbar_vee2;
    %         Vbar_vee2 = ((x_xbar_vee2 - x_XCG)/cmac_w1)*(S_vee2_s/S_ref);
    Vbar_vee2 = ((x_xbar_vee2 - x_XCG)/cmac_w1);%para no mult dos veces *(S_vee2_e/S_ref)
    CL_alpha_wb_vee2 = Aero.CL_alpha_vee2_CR;
    
    if wingspan2bodydiam_vee2 <= 2
        CLalpha_vee2_e = CL_alpha_wb_vee2*S_vee2_e/Geo_tier.S_vee2;
    else
        CLalpha_vee2_e = CL_alpha_wb_vee2;
    end    
    
    deps_dalpha_vee2 = Effects.deps_dalpha_vee2;
    % W1 and vee2 with no dynamic pressure correction
    CL_alphapunto_vee2 = 2*CLalpha_vee2_e*Vbar_vee2*deps_dalpha_vee2; %eq 4.525 PAMADI
    
    % Dyanmic pressure correction
    if prop_wash_effect == 1
        eta_vee2_afe_S_vee2_afe_S_ref = afe.eta_vee2_afe_S_vee2_afe_S_ref;
        eta_vee2_no_afe_S_vee2_no_afe_S_ref = afe.eta_vee2_no_afe_S_vee2_no_afe_S_ref;
        CL_alphapunto_vee2_pw = eta_vee2_afe_S_vee2_afe_S_ref*CL_alphapunto_vee2 + ...
            eta_vee2_no_afe_S_vee2_no_afe_S_ref*CL_alphapunto_vee2;
    else
        CL_alphapunto_vee2_pw = CL_alphapunto_vee2*conv_vee2;
    end
    % Revisar OJO - contribuición de presiópn dinámica en ala vee2
else
    CL_alphapunto_vee2_pw = 0;
end

if Can == 1
    deps_dalpha_can = Effects.deps_dalpha_can;
    x_xbar_can = Geo_tier.x_xbar_can;
    
    % Reference wing
    S_can   = Geo_tier.S_can;
    b_can   = Geo_tier.b_can;
    S_can_e = Geo_tier.S_can_e;
    w_can = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_can);
    wingspan2bodydiam_can = b_can/w_can;   
    conv_can    = S_can/S_w1; %adimensionalizados con la S_can total en flow5, no con la expuesta;
    CLalpha_can = Aero.CL_alpha_can_CR;
    
    if wingspan2bodydiam_can <= 2
        CLalpha_can_e = CLalpha_can*S_can_e/S_can;
    else
        CLalpha_can_e  = CLalpha_can;
    end
    
    Vbar_can = ((-x_xbar_can + x_XCG)/cmac_w1); %para no mult dos veces *(S_can_e/S_ref)
    CL_alphapunto_can = 2*CLalpha_can_e*Vbar_can*deps_dalpha_can; %eq 4.525 PAMADI
    % Revisar OJO - contribuición de presiópn dinámica en ala can
    % Dyanmic pressure correction
    if prop_wash_effect == 1
        eta_can_afe_S_can_afe_S_ref = afe.eta_can_afe_S_can_afe_S_ref;
        eta_can_no_afe_S_can_no_afe_S_ref = afe.eta_can_no_afe_S_can_no_afe_S_ref;
        CL_alphapunto_can_pw = eta_can_afe_S_can_afe_S_ref*CL_alphapunto_can + ...
            eta_can_no_afe_S_can_no_afe_S_ref*CL_alphapunto_can;
    else
        CL_alphapunto_can_pw = CL_alphapunto_can*conv_can;
    end
else
    CL_alphapunto_can_pw = 0;
end

% NAcelle contribution
if Nac== 1
    CL_alphapunto_nac_pw = 0;
else
    CL_alphapunto_nac_pw = 0;
end

% Total Derivative CLalphapunto
CL_alphapuntoWB = CL_alphapunto_w1_pw + CL_alphapunto_fus;
CL_alphapunto = CL_alphapunto_w1_pw + CL_alphapunto_fus + CL_alphapunto_HTP_pw + ...
    CL_alphapunto_Vee_pw + CL_alphapunto_Vee2_pw + CL_alphapunto_can_pw + CL_alphapunto_nac_pw;

% Storing DATA
Stab_Der.CL_alphapunto = CL_alphapunto;
%     Stab_Der.CL_alphapunto_w1_pw = CL_alphapunto_w1_pw;
%     Stab_Der.CL_alphapunto_w2_pw = CL_alphapunto_w2_pw;
%     Stab_Der.CL_alphapunto_can_pw = CL_alphapunto_can_pw;
%     Stab_Der.CL_alphapunto_fus = CL_alphapunto_fus;

% Stores Derivatives per parts
Stab_Der_parts.CL_alphapunto = CL_alphapunto;
Stab_Der_parts.CL_alphapunto_w1 = CL_alphapunto_w1_pw;
Stab_Der_parts.CL_alphapunto_fus = CL_alphapunto_fus;
Stab_Der_parts.CL_alphapunto_HTP = CL_alphapunto_HTP_pw;
Stab_Der_parts.CL_alphapunto_Vee = CL_alphapunto_Vee_pw;
Stab_Der_parts.CL_alphapunto_Vee2 = CL_alphapunto_Vee2_pw;
Stab_Der_parts.CL_alphapunto_can = CL_alphapunto_can_pw;
Stab_Der_parts.CL_alphapunto_nac = CL_alphapunto_nac_pw;
Stab_Der_parts.CL_alphapunto_WB = CL_alphapuntoWB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CM_alphapunto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_Sb1 = @(x) interp1(x_Area_body(2:end),Area_body(2:end),x).*(x);               %eq 4.514 PAMADI
x_c_fus = (1/Vol_TOT)*quad(@(x) int_Sb1(x),0,length_fus);           %eq 4.514 PAMADI
x_m = x_XCG;
xm_1 = x_m/length_fus;                      %eq 4.513 PAMADI
xc_1 = x_c_fus/length_fus;                  %eq 4.513 PAMADI
VB_1 = Vol_TOT/(Area_b_max*length_fus);     %eq 4.513 PAMADI
% Cálculo del centro aerodinamico del fuselaje
% Derivada del área:
% dSdX = diff(Area_body)./diff(x_Area_body);
% %  Calculates where fuselage ceases to be potential
% i=1;
% while dSdX(i) > 0
%     x_0_v = i;
%     i=i+1;
% end
dSdX = Body_Geo.dSdX;
x_0_v = 2; %  does not include the firs elements
while dSdX(x_0_v) > 0
    x_0_v=x_0_v+1;
end
x_0 = x_Area_body(x_0_v+1);
int_Sb_CMq = @(x) interp1(x_Area_body(2:end),dSdX(2:end),x).*(x_m - x);    %eq 4.516 PAMADI
CMalpha_B = (2*f_k2_k1/Vol_TOT)*quad(@(x) int_Sb_CMq(x),0,x_0);       %eq 4.516 PAMADI
CMalpha_p = CMalpha_B*VB_1;
CM_alphapunto_B = 2*CMalpha_p*((xc_1 - xm_1)/(1 - xm_1 - VB_1))*(Vol_TOT/(length_fus*Area_b_max)); %eq 4.543 PAMADI
CM_alphapunto_fus = CM_alphapunto_B*((Area_b_max*length_fus^2)/(S_w1*cmac_w1^2))  ;    %eq 4.539 PAMADI

if W1 == 1
    % Only valid for tau < 4
    if tau<4
        CMo_g = (pi*AR_w1_e/(2*beta^2))*(0.0008*tau^4 - 0.0075*tau^3 +...
            0.0185*tau^2 + 0.0128*tau - 0.0003); %eq 4.542 PAMADI
    else
        CMo_g = 0;
    end
    
    %% OJO! using he distance from apex of wing LE root to Xac
    CM_alphapunto_e_pp_1 = -(81/32)*((xbar_w1/cmac_w1_e)^2)*CLalpha_w1_e;        %eq 4.541 PAMADI
    CM_alphapunto_e_pp_2 = (9/2)*CMo_g;                                         %eq 4.541 PAMADI
    CM_alphapunto_e_pp =  CM_alphapunto_e_pp_1 + CM_alphapunto_e_pp_2;          %eq 4.541 PAMADI
    CM_alphapunto_e = CM_alphapunto_e_pp + ((x_xbar_w1_e - x_XCG)/cmac_w1_e)*CL_alphapunto_e; %eq 4.540 PAMADI
    
    % W1 and w2 with no dynamic pressure correction %eq 4.539 PAMADI
    % Ojo revisar por la contribución de la presión dinámica
    CM_alphapunto_w1 = (KWB_w1 + KBW_w1)*(cmac_w1_e^2/cmac_w1^2)*CM_alphapunto_e;    %eq 4.539 PAMADI
    % Dyanmic pressure correction
    if prop_wash_effect == 1
        eta_w1_afe_S_w1_afe_S_ref = afe.eta_w1_afe_S_w1_afe_S_ref;
        eta_w1_no_afe_S_w1_no_afe_S_ref = afe.eta_w1_no_afe_S_w1_no_afe_S_ref;
        CM_alphapunto_w1_pw = eta_w1_afe_S_w1_afe_S_ref*CM_alphapunto_w1 +...
            eta_w1_no_afe_S_w1_no_afe_S_ref*CM_alphapunto_w1;
    else
        CM_alphapunto_w1_pw = CM_alphapunto_w1;
    end
    
else
    CM_alphapunto_w1_pw = 0;
end

if HTP == 1
    CM_alphapunto_w2 = - CL_alphapunto_w2*Vbar_h;      %eq 4.537 PAMADI
    % W1 and w2 with no dynamic pressure correction
    %% Revisar OJO - contribuición de presión dinámica en ala w1
    % Dynamic pressure correction
    if prop_wash_effect == 1
        eta_w2_afe_S_w2_afe_S_ref = afe.eta_w2_afe_S_w2_afe_S_ref;
        eta_w2_no_afe_S_w2_no_afe_S_ref = afe.eta_w2_no_afe_S_w2_no_afe_S_ref;
        CM_alphapunto_HTP_pw = eta_w2_afe_S_w2_afe_S_ref*CM_alphapunto_w2 + ...
            eta_w2_no_afe_S_w2_no_afe_S_ref*CM_alphapunto_w2;
    else
        CM_alphapunto_HTP_pw = CM_alphapunto_w2*conv_w2;
    end
else
    CM_alphapunto_HTP_pw = 0;
end

if Vee == 1
    CM_alphapunto_vee = - CL_alphapunto_vee*((x_xbar_vee - x_XCG)/cmac_w1);      %eq 4.537 PAMADI

    % W1 and vee with no dynamic pressure correction
    %% Revisar OJO - contribuición de presión dinámica en ala w1
    % Dynamic pressure correction
    if prop_wash_effect == 1
        eta_vee_afe_S_vee_afe_S_ref = afe.eta_vee_afe_S_vee_afe_S_ref;
        eta_vee_no_afe_S_vee_no_afe_S_ref = afe.eta_vee_no_afe_S_vee_no_afe_S_ref;
        CM_alphapunto_Vee_pw = eta_vee_afe_S_vee_afe_S_ref*CM_alphapunto_vee + ...
            eta_vee_no_afe_S_vee_no_afe_S_ref*CM_alphapunto_vee;
    else
        CM_alphapunto_Vee_pw = CM_alphapunto_vee*conv_vee;
    end
else
    CM_alphapunto_Vee_pw = 0;
end

if Vee2 == 1
    CM_alphapunto_vee2 = - CL_alphapunto_vee2*((x_xbar_vee2 - x_XCG)/cmac_w1);      %eq 4.537 PAMADI

    % W1 and vee2 with no dynamic pressure correction
    %% Revisar OJO - contribuición de presión dinámica en ala w1
    % Dynamic pressure correction
    if prop_wash_effect == 1
        eta_vee2_afe_S_vee2_afe_S_ref = afe.eta_vee2_afe_S_vee2_afe_S_ref;
        eta_vee2_no_afe_S_vee2_no_afe_S_ref = afe.eta_vee2_no_afe_S_vee2_no_afe_S_ref;
        CM_alphapunto_vee2_pw = eta_vee2_afe_S_vee2_afe_S_ref*CM_alphapunto_vee2 + ...
            eta_vee2_no_afe_S_vee2_no_afe_S_ref*CM_alphapunto_vee2;
    else
        CM_alphapunto_vee2_pw = CM_alphapunto_vee2*conv_vee2;
    end
else
    CM_alphapunto_vee2_pw = 0;
end

if Can == 1
    CM_alphapunto_can = - CL_alphapunto_can*((x_xbar_can - x_XCG)/cmac_w1);      %eq 4.537 PAMADI
    % W1 and w2 with no dynamic pressure correction
    %% Revisar OJO - contribuición de presión dinámica en ala w1
    % Dynamic pressure correction
    if prop_wash_effect == 1
        eta_can_afe_S_can_afe_S_ref = afe.eta_can_afe_S_can_afe_S_ref;
        eta_can_no_afe_S_can_no_afe_S_ref = afe.eta_can_no_afe_S_can_no_afe_S_ref;
        CM_alphapunto_can_pw = eta_can_afe_S_can_afe_S_ref*CM_alphapunto_can + ...
            eta_can_no_afe_S_can_no_afe_S_ref*CM_alphapunto_can;
    else
        CM_alphapunto_can_pw = CM_alphapunto_can*conv_can;
    end
else
    CM_alphapunto_can_pw = 0;
end

% NAcelle contribution
if Nac== 1
    CM_alphapunto_nac_pw = 0;
else
    CM_alphapunto_nac_pw = 0;
end

% Total Derivative CMalphapunto
CM_alphapunto = CM_alphapunto_w1_pw + CM_alphapunto_fus + CM_alphapunto_HTP_pw + ...
    CM_alphapunto_Vee_pw + CM_alphapunto_Vee2_pw + CM_alphapunto_can_pw + CM_alphapunto_nac_pw;

% Storing DATA
Stab_Der.CM_alphapunto = CM_alphapunto;
%     Stab_Der.CM_alphapunto_w1_pw = CM_alphapunto_w1_pw;
%     Stab_Der.CM_alphapunto_w2_pw = CM_alphapunto_w2_pw;
%     Stab_Der.CM_alphapunto_can_pw = CM_alphapunto_can_pw;
%     Stab_Der.CM_alphapunto_fus = CM_alphapunto_fus;

% Stores Derivatives per parts
Stab_Der_parts.CM_alphapunto = CM_alphapunto;
Stab_Der_parts.CM_alphapunto_w1 = CM_alphapunto_w1_pw;
Stab_Der_parts.CM_alphapunto_fus = CM_alphapunto_fus;
Stab_Der_parts.CM_alphapunto_HTP = CM_alphapunto_HTP_pw;
Stab_Der_parts.CM_alphapunto_Vee = CM_alphapunto_Vee_pw;
Stab_Der_parts.CM_alphapunto_Vee2 = CM_alphapunto_Vee2_pw;
Stab_Der_parts.CM_alphapunto_can = CM_alphapunto_can_pw;
Stab_Der_parts.CM_alphapunto_nac = CM_alphapunto_nac_pw;

%%%%%%%%%%%%%%%%%%%%%%%%% Alpha punto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CD_alphapunto = 0;              % vale 0 para vuelo subsónico
CX_alphapunto = -CD_alphapunto;              % vale 0 para vuelo subsónico
CZ_alphapunto = -CL_alphapunto; % Roskam 3.14

Stab_Der.CD_alphapunto = CD_alphapunto;
Stab_Der.CX_alphapunto = CX_alphapunto;
Stab_Der.CZ_alphapunto = CZ_alphapunto;
Stab_Der.CL_alphapunto = CL_alphapunto;
Stab_Der.CM_alphapunto = CM_alphapunto;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CDalfa%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Selects the source of the polar model
% C_D0 = Aero.Polar.C_D0;
% C_D1 = Aero.Polar.C_D1;
% C_D2 = Aero.Polar.C_D2;
% 
% CL_alpha_ac = TRIM_RESULTS.CL_alpha_ac;
% trim_alpha = TRIM_RESULTS.trim_alpha;
% CM_alpha_ac = TRIM_RESULTS.CM_alpha_ac_des;
% CD = Stab_Der.CD;
% CD_alpha = (C_D1 + 2*C_D2*CL_needed)*CL_alpha_ac;
% 
% CX_alfa = - CD_alpha + CL_alpha_ac*trim_alpha + CL_needed; % Roskam 3.128
% CX_alfa = - CD_alpha  + CL; % Steady State Flight condition Roskam 3.128
% CX_alpha = - CD_alpha + CL; % Pamadi 4.447
% 
% % CD_trim_delta = abs(CD_delta*trim_delta_e);
% %     CD = C_D0 + C_D1*CL + C_D2*CL^2 + CD_trim_delta;  %polar aeronave
% CZ_alpha = - CL_alpha_ac -CD_alpha*trim_alpha - CD; % Roskam 3.131
% CZ_alpha = - CL_alpha_ac - CD; % Steady State Flight condition Roskam 3.131
% CM_alpha = CM_alpha_ac;
% 
% Stab_Der.CX_alpha = CX_alpha;
% Stab_Der.CZ_alpha = CZ_alpha;
% Stab_Der.CM_alpha = CM_alpha;

end
