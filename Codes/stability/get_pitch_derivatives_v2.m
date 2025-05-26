    function [Stab_Der, Stab_Der_parts] = get_pitch_derivatives_v2(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,Body_Geo,Aero,OUTPUT_read_XLSX)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%Pitch Derivatives%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W1      = AC_CONFIGURATION.W1;
HTP     = AC_CONFIGURATION.HTP;
VTP     = AC_CONFIGURATION.VTP;
Can     = AC_CONFIGURATION.Can;
Vee     = AC_CONFIGURATION.Vee;
Vee2     = AC_CONFIGURATION.Vee2;
Nac     = AC_CONFIGURATION.Nac;
AC_type = AC_CONFIGURATION.AC_type;

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
S_w1_e = Geo_tier.S_w1_e;
S_ref = Geo_tier.S_ref;
cmac_w1 = Geo_tier.cmac_w1;
cmac_w1_e = Geo_tier.cmac_w1_e;
AR_w1 = Geo_tier.AR_w1;
AR_w1_e = Geo_tier.AR_w1_e;
Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
b_w1 = Geo_tier.b_w1;
V = conditions.V;
a = Performance.a;
Mach = V/a;

CLalpha_w1_e_pw = Stab_Der_parts.CLalpha_w1_e_pw;
% CLalpha_w1 = Stab_Der_parts.CL_alpha_w1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%CLq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLalpha_fus = 2*f_k2_k1*Area_b_max/(Vol_TOT^(2/3));  %eq 4.495 PAMADI
CLalpha_Bp = CLalpha_fus*(Vol_TOT^(2/3))/Area_b_max; %eq 4.494 PAMADI
CLq_B = 2*CLalpha_Bp*(1 - x_XCG/length_fus);           %eq 4.493 PAMADI
CLq_fus = CLq_B*((length_fus*Area_b_max)/(S_w1*cmac_w1)); %eq 4.488 PAMADI
KWB_w1  = Stab_Der_parts.KWB_w1;
KBW_w1  = Stab_Der_parts.KBW_w1;

if W1 == 1
    x_xbar_w1_e = Geo_tier.x_xbar_w1_e;
    B_prandalt = sqrt(1 - (Mach^2)*(cos(Lambda_c4_w1))^2);     %eq 4.510 PAMADI
    chi = (x_xbar_w1_e - x_XCG)/cmac_w1_e; %eq 4.490 PAMADI
    % CLalpha_wb_w1 already includes (KWB_w1 + KBW_w1)
    CLq_e = CLalpha_w1_e_pw*(0.5 + 2*chi);   %eq 4.489 PAMADI
    CLq_w1 = ((S_w1_e*cmac_w1_e)/(S_ref*cmac_w1))*CLq_e;
    CLq_w  =  (KWB_w1 + KBW_w1)*CLq_w1;
%     % ROSKAM AAA
%     % Curve lift slope with no flaps
%     CL_alpha_clean = CLalpha_w1_e_pw;
%     % increment of wing-fuselage lift curve slope due to power
%     DeltaCL_alpha_clean = 0;
%     CLq_w_M0 = (CL_alpha_clean + DeltaCL_alpha_clean)*(0.5 + 2*chi);
%     int_CLq_w1 = (AR_w1 + 2*cos(Lambda_c4_w1))/(AR_w1*B_prandalt + 2*cos(Lambda_c4_w1));
%     CLq_w1 = ((S_w1_e*cmac_w1_e)/(S_ref*cmac_w1))*int_CLq_w1*CLq_w_M0;
    % W1 and HTP with no dynamic pressure correction
    if prop_wash_effect == 1
        eta_w1_afe_S_w1_afe_S_ref = afe.eta_w1_afe_S_w1_afe_S_ref;
        eta_w1_no_afe_S_w1_no_afe_S_ref = afe.eta_w1_no_afe_S_w1_no_afe_S_ref;
        CLq_w1_pw = eta_w1_afe_S_w1_afe_S_ref*CLq_w + eta_w1_no_afe_S_w1_no_afe_S_ref*CLq_w;
    else
        CLq_w1_pw = CLq_w;
    end   
else
    CLq_w1_pw = 0;
end

if HTP == 1
    % Reference wing
    S_HTP   = Geo_tier.S_HTP;
    S_HTP_e = Geo_tier.S_HTP_e;
    b_HTP   = Geo_tier.b_HTP;
    w_HTP = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_HTP);
    wingspan2bodydiam_HTP = b_HTP/w_HTP;
    conv_HTP = S_HTP/S_w1; %adimensionalizados con la S_HTP_total en flow5, no con la expuesta;
    
    x_xbar_HTP = Geo_tier.x_xbar_HTP;
    %         Vbar_h = ((x_xbar_HTP - x_XCG)/cmac_w1)*(S_HTP_e/S_ref);
    Vbar_h = ((x_xbar_HTP - x_XCG)/cmac_w1); %Sin *(S_HTP_e/S_ref) puesto que luego
    %multiplicas por eta_HTP_no_afe_S_HTP_no_afe_Sref y ya va
    %multiplicado ahí
    CL_alpha_HTP = Aero.CL_alpha_HTP_CR;
    
    if wingspan2bodydiam_HTP <= 2
        CLalpha_HTP_e     = CL_alpha_HTP*S_HTP_e/Geo_tier.S_HTP;
    else
        CLalpha_HTP_e     = CL_alpha_HTP;
    end
    
    x_xbar_HTP_e = Geo_tier.x_xbar_HTP_e;
    
    CLq_HTP = 2*CLalpha_HTP_e*Vbar_h; %eq 10.72 ROSKAM
    % W1 and HTP with no dynamic pressure correction
    if prop_wash_effect == 1
        eta_HTP_afe_S_HTP_afe_S_ref = afe.eta_HTP_afe_S_HTP_afe_S_ref;
        eta_HTP_no_afe_S_HTP_no_afe_S_ref = afe.eta_HTP_no_afe_S_HTP_no_afe_S_ref;
        CLq_HTP_pw = eta_HTP_afe_S_HTP_afe_S_ref*CLq_HTP + eta_HTP_no_afe_S_HTP_no_afe_S_ref*CLq_HTP;
    else
        CLq_HTP_pw = CLq_HTP*conv_HTP;
    end

else
    CLq_HTP_pw = 0;
end

if Vee == 1
    % Reference wing
    S_vee   = Geo_tier.S_vee;
    S_vee_e = Geo_tier.S_vee_e;
    b_vee   = Geo_tier.b_vee;
    
    w_vee = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee);
    wingspan2bodydiam_vee = b_vee/w_vee;
    conv_vee = S_vee/S_w1; %adimensionalizados con la S_vee_total en flow5, no con la expuesta;
    x_xbar_vee = Geo_tier.x_xbar_vee;
    x_xbar_vee_e = Geo_tier.x_xbar_vee_e;
    
    %         Vbar_vee = ((x_xbar_vee - x_XCG)/cmac_w1)*(S_vee_s/S_ref);
    Vbar_vee = ((x_xbar_vee - x_XCG)/cmac_w1); %Sin *(S_vee_e/S_ref) puesto que luego
    %multiplicas por eta_vee_no_afe_S_vee_no_afe_Sref y ya va
    %multiplicado ahí
    CL_alpha_wb_vee = Aero.CL_alpha_vee_CR;
    
    
    % Reference wing
    S_vee = Geo_tier.S_vee;
    S_vee_e = Geo_tier.S_vee_e;

    if wingspan2bodydiam_vee <= 2
        CLalpha_vee_e = CL_alpha_wb_vee*S_vee_e/Geo_tier.S_vee;
    else
        CLalpha_vee_e = CL_alpha_wb_vee;
    end    
    CLq_vee = 2*CLalpha_vee_e*Vbar_vee; %eq 10.72 ROSKAM
    % W1 and vee with no dynamic pressure correction
    if prop_wash_effect == 1
        eta_vee_afe_S_vee_afe_S_ref = afe.eta_vee_afe_S_vee_afe_S_ref;
        eta_vee_no_afe_S_vee_no_afe_S_ref = afe.eta_vee_no_afe_S_vee_no_afe_S_ref;
        CLq_vee_pw = eta_vee_afe_S_vee_afe_S_ref*CLq_vee + eta_vee_no_afe_S_vee_no_afe_S_ref*CLq_vee;
    else
        CLq_vee_pw = CLq_vee*conv_vee;
    end
else
    CLq_vee_pw = 0;
end

if Vee2 == 1
    % Reference wing
    S_vee2   = Geo_tier.S_vee2;
    S_vee2_e = Geo_tier.S_vee2_e;
    b_vee2   = Geo_tier.b_vee2;
    
    w_vee2 = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee2);
    wingspan2bodydiam_vee2 = b_vee2/w_vee2;
    conv_vee2 = S_vee2/S_w1; %adimensionalizados con la S_vee2_total en flow5, no con la expuesta;
    x_xbar_vee2 = Geo_tier.x_xbar_vee2;
    x_xbar_vee2_e = Geo_tier.x_xbar_vee2_e;
    
    %         Vbar_vee2 = ((x_xbar_vee2 - x_XCG)/cmac_w1)*(S_vee2_s/S_ref);
    Vbar_vee2 = ((x_xbar_vee2 - x_XCG)/cmac_w1); %Sin *(S_vee2_e/S_ref) puesto que luego
    %multiplicas por eta_vee2_no_afe_S_vee2_no_afe_Sref y ya va
    %multiplicado ahí
    CL_alpha_wb_vee2 = Aero.CL_alpha_vee2_CR;
    
    if wingspan2bodydiam_vee2 <= 2
        CLalpha_vee2_e = CL_alpha_vee2*S_vee2_e/Geo_tier.S_vee2;
    else
        CLalpha_vee2_e     = CL_alpha_wb_vee2;
    end    
    CLq_vee2 = 2*CLalpha_vee2_e*Vbar_vee2; %eq 10.72 ROSKAM
    % W1 and vee2 with no dynamic pressure correction
    if prop_wash_effect == 1
        eta_vee2_afe_S_vee2_afe_S_ref = afe.eta_vee2_afe_S_vee2_afe_S_ref;
        eta_vee2_no_afe_S_vee2_no_afe_S_ref = afe.eta_vee2_no_afe_S_vee2_no_afe_S_ref;
        CLq_vee2_pw = eta_vee2_afe_S_vee2_afe_S_ref*CLq_vee2 + eta_vee2_no_afe_S_vee2_no_afe_S_ref*CLq_vee2;
    else
        CLq_vee2_pw = CLq_vee2*conv_vee2;
    end
else
    CLq_vee2_pw = 0;
end

if Can == 1
    x_xbar_can = Geo_tier.x_xbar_can;
    % Reference wing
    S_can_e = Geo_tier.S_can_e;
    S_can   = Geo_tier.S_can;
    b_can   = Geo_tier.b_can;
    w_can = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_can);
    wingspan2bodydiam_can = b_can/w_can;
    %     conv_can = S_can_e/S_w1;
    conv_can = S_can/S_w1; %adimensionalizados con la S_can total en flow5, no con la expuesta;
    CLalpha_can = Aero.CL_alpha_can_CR;
    
    if wingspan2bodydiam_can <= 2
        CLalpha_can_e = CLalpha_can*S_can_e/S_can;
    else
        CLalpha_can_e  = CLalpha_can;
    end
    
%     Vbar_can = ((x_xbar_can + x_XCG)/cmac_w1)*(S_can_s/S_ref);
    Vbar_can = ((x_xbar_can - x_XCG)/cmac_w1);%Sin *(S_HTP_e/S_ref) puesto que luego
    %multiplicas por eta_can_no_afe_S_can_no_afe_Sref y ya va
    %multiplicado ahí
    CLq_can = -2*CLalpha_can_e*Vbar_can; %eq 10.72 ROSKAM
    % Can with no dynamic pressure correction
    if prop_wash_effect == 1
        eta_can_afe_S_can_afe_S_ref = afe.eta_can_afe_S_can_afe_S_ref;
        eta_can_no_afe_S_can_no_afe_S_ref = afe.eta_can_no_afe_S_can_no_afe_S_ref;
        CLq_can_pw = eta_can_afe_S_can_afe_S_ref*CLq_can + eta_can_no_afe_S_can_no_afe_S_ref*CLq_can;
    else
        CLq_can_pw = CLq_can*conv_can;
    end
        
else
    CLq_can_pw = 0;
end

% NAcelle contribution
if Nac== 1
    CLq_nac_pw = 0;
else
    CLq_nac_pw = 0;
end

% Total Derivative CLq
CLq_WB = CLq_w + CLq_fus;
CLq = CLq_WB + CLq_HTP_pw + CLq_vee_pw + CLq_vee2_pw + CLq_can_pw + CLq_nac_pw; %Pamadi

Stab_Der.CLq = CLq;
% Stores Derivatives per parts
Stab_Der_parts.CL_q = CLq;
Stab_Der_parts.CL_q_w1 = CLq_w1_pw;
Stab_Der_parts.CL_q_fus = CLq_fus;
Stab_Der_parts.CL_q_HTP = CLq_HTP_pw;
Stab_Der_parts.CL_q_vee = CLq_vee_pw;
Stab_Der_parts.CL_q_vee2 = CLq_vee2_pw;
Stab_Der_parts.CL_q_can = CLq_can_pw;
Stab_Der_parts.CL_q_nac = CLq_nac_pw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Cmq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuselage contribution


% Centroid fuselage
int_Sb1 = @(x) interp1(x_Area_body,Area_body,x).*(x);               %eq 4.514 PAMADI
% x_c_fus = (1/Vol_TOT)*quad(@(x) int_Sb1(x),0,length_fus);           %eq 4.514 PAMADI
x_c_fus = (1/Vol_TOT)*integral(@(x) int_Sb1(x),0,length_fus);           %eq 4.514 PAMADI
x_m = x_XCG;
xm_1 = x_m/length_fus;                      %eq 4.513 PAMADI
xc_1 = x_c_fus/length_fus;                  %eq 4.513 PAMADI
VB_1 = Vol_TOT/(Area_b_max*length_fus);     %eq 4.513 PAMADI

% Cálculo del centro aerodinamico del fuselaje

% Derivada del área:

% % Original calculation
% % dSdX = diff(Area_body)./diff(x_Area_body);
% dSdX = Body_Geo.dSdX;
% %  Calculates where fuselage ceases to be potential
% i=3;
% while dSdX(i) > 0
%     x_0_v = i;
%     i=i+1;
% end
% From Calc_Bpdy_Geometry_Dec218
% Volumen total geometría
% dSdX                = (Body_Geo.S_x(2) - Body_Geo.S_x(1))/(Body_Geo.x(2) - Body_Geo.x(1));
% dSdX                = [dSdX,(Body_Geo.S_x(3:end) - Body_Geo.S_x(1:(end-2)))./(Body_Geo.x(3:end) - Body_Geo.x(1:(end-2)))];
% dSdX                = [dSdX, (Body_Geo.S_x(end) - Body_Geo.S_x(end-1))/(Body_Geo.x(end) - Body_Geo.x(end-1))];
dSdX = Body_Geo.dSdX;
dSdX = Body_Geo.dSdX;

x_0_v = 2; %  does not include the first elements
while dSdX(x_0_v) > 0
    x_0_v=x_0_v+1;
end
x_0 = x_Area_body(x_0_v+1);

% int_Sb_CMq = @(x) interp1(x_Area_body(2:end),dSdX(2:end),x).*(x_m - x);    %eq 4.516 PAMADI
int_Sb_CMq = @(x) interp1(x_Area_body(1:end),dSdX(1:end),x).*(x_m - x);    %eq 4.516 PAMADI
CMalpha_B = (2*f_k2_k1/Vol_TOT)*quad(@(x) int_Sb_CMq(x),0,x_0);       %eq 4.516 PAMADI
CMalpha_B = (2*f_k2_k1/Vol_TOT)*integral(@(x) int_Sb_CMq(x),0,x_0);       %eq 4.516 PAMADI
CMalpha_p = CMalpha_B*VB_1;                                           %eq 4.515 PAMADI
CMq_B = 2*CMalpha_p*((1-xm_1)^2 - VB_1*(xc_1-xm_1))/(1-xm_1-VB_1);    %eq 4.512 PAMADI

% W1 and HTP with no dynamic pressure correction %eq 4.502 PAMADI
% Ojo revisar por la contribución de la presión dinámica

wingspan2bodydiam = b_w1/w_Area_b_max;
flagwingspan2bodydiam = OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam; 
Munk_fuselage_constribution = OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution;
%% Contribution of Fuselage
if   wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1 || Munk_fuselage_constribution == 1
    CMq_fus = CMq_B*(Area_b_max/S_w1)*((length_fus/cmac_w1)^2);                   %eq 4.502 PAMADI
else
    CMq_fus = 0;
end


if W1 == 1
    c1 = AR_w1_e^3*(tan(Lambda_c4_w1))^2;          %eq 4.505 PAMADI
    c2 = 3/B_prandalt;                           %eq 4.506 PAMADI
    c3 = AR_w1_e*B_prandalt + 6*cos(Lambda_c4_w1); %eq 4.507 PAMADI
    c4 = AR_w1_e + 6*cos(Lambda_c4_w1);            %eq 4.508 PAMADI
    c5 = AR_w1_e + 2*cos(Lambda_c4_w1);            %eq 5.509 PAMADI(length_fus/cmac_w1)^2
    % %         % Assumes that the CLalpha_w1 in 2D is the same as 3D
    % %         CLalpha_w1_2D = CLalpha_w1;
    
    %% Perfil revisado para el ala
    airfoil_data =  importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_w1);
    t_c = Geo_tier.t_c_ail;
    
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c); %Pamadi 3.17
    clalpha_HTP_M0 = 1.05*cla_clatheo*clatheo; %Pamadi 3.17
    CLalpha_w1_2D = clalpha_HTP_M0/sqrt(1-Mach^2); %Pamadi 3.17
    CMq_e_M02 = -0.7*CLalpha_w1_2D*cos(Lambda_c4_w1)*((AR_w1_e*(0.5*chi + 2*chi^2))/c5 +...
        (c1/(24*c4)) + 1/8);                    %eq 4.504 PAMADI
    CMq_e = CMq_e_M02*(c1/c3 + c2)/(c1/c4 + 3); %eq 4.503 PAMADI
    
    % W1 and HTP with no dynamic pressure correction %eq 4.502 PAMADI
    % Ojo revisar por la contribución de la presión dinámica
    KWB_w1 = 0.1714*(w_Area_b_max/b_w1)^2 + 0.8326*(w_Area_b_max/b_w1) + 0.9974; %eq 3.27 PAMADI
    KBW_w1 = 0.7810*(w_Area_b_max/b_w1)^2 + 1.1976*(w_Area_b_max/b_w1) + 0.0088; %eq 3.28 PAMADI
    CMq_w1 = (KWB_w1 + KBW_w1)*((cmac_w1_e/cmac_w1)^2)*CMq_e; %eq 4.502 PAMADI
    % Dynamic pressure correction
    if prop_wash_effect == 1
        eta_w1_afe_S_w1_afe_S_ref = afe.eta_w1_afe_S_w1_afe_S_ref;
        eta_w1_no_afe_S_w1_no_afe_S_ref = afe.eta_w1_no_afe_S_w1_no_afe_S_ref;
        CMq_w1_pw = eta_w1_afe_S_w1_afe_S_ref*CMq_w1 + eta_w1_no_afe_S_w1_no_afe_S_ref*CMq_w1;
    else
        CMq_w1_pw = CMq_w1;
    end       
else
    CMq_w1_pw = 0;
end

if HTP == 1
    %% Revisar OJO - contribuición de presión dinámica en ala w1
    CMq_HTP = - 2*CLalpha_HTP_e*Vbar_h*((x_xbar_HTP - x_XCG)/cmac_w1);              %eq 10.78 ROSKAM
    % Dynamic pressure correction
    if prop_wash_effect == 1
        eta_HTP_afe_S_HTP_afe_S_ref = afe.eta_HTP_afe_S_HTP_afe_S_ref;
        eta_HTP_no_afe_S_HTP_no_afe_S_ref = afe.eta_HTP_no_afe_S_HTP_no_afe_S_ref;
        CMq_HTP_pw = eta_HTP_afe_S_HTP_afe_S_ref*CMq_HTP + eta_HTP_no_afe_S_HTP_no_afe_S_ref*CMq_HTP;
    else
        CMq_HTP_pw = CMq_HTP*conv_HTP;
    end       
else
    CMq_HTP_pw = 0;
end

if Vee == 1
    %% Revisar OJO - contribuición de presiópn dinámica en ala w1
    CMq_vee = - 2*CLalpha_vee_e*Vbar_vee*((x_xbar_vee - x_XCG)/cmac_w1);              %eq 10.78 ROSKAM
    % Dynamic pressure correction
    if prop_wash_effect == 1
        eta_vee_afe_S_vee_afe_S_ref = afe.eta_vee_afe_S_vee_afe_S_ref;
        eta_vee_no_afe_S_vee_no_afe_S_ref = afe.eta_vee_no_afe_S_vee_no_afe_S_ref;
        CMq_vee_pw = eta_vee_afe_S_vee_afe_S_ref*CMq_vee + eta_vee_no_afe_S_vee_no_afe_S_ref*CMq_vee;
    else
        CMq_vee_pw = CMq_vee*conv_vee;
    end

else
    CMq_vee_pw = 0;
end


if Vee2 == 1
    %% Revisar OJO - contribuición de presiópn dinámica en ala w1
    CMq_vee2 = - 2*CLalpha_vee2_e*Vbar_vee2*((x_xbar_vee2 - x_XCG)/cmac_w1);              %eq 10.78 ROSKAM
    % Dynamic pressure correction
    if prop_wash_effect == 1
        eta_vee2_afe_S_vee2_afe_S_ref = afe.eta_vee2_afe_S_vee2_afe_S_ref;
        eta_vee2_no_afe_S_vee2_no_afe_S_ref = afe.eta_vee2_no_afe_S_vee2_no_afe_S_ref;
        CMq_vee2_pw = eta_vee2_afe_S_vee2_afe_S_ref*CMq_vee2 + eta_vee2_no_afe_S_vee2_no_afe_S_ref*CMq_vee2;
    else
        CMq_vee2_pw = CMq_vee2*conv_vee2;
    end

else
    CMq_vee2_pw = 0;
end

if Can == 1
    %% Revisar OJO - contribuición de presiópn dinámica en ala w1
    CMq_can = - 2*CLalpha_can_e*Vbar_can*((x_xbar_can - x_XCG)/cmac_w1);              %eq 10.78 ROSKAM
    % Dynamic pressure correction
    if prop_wash_effect == 1
        eta_can_afe_S_can_afe_S_ref = afe.eta_can_afe_S_can_afe_S_ref;
        eta_can_no_afe_S_can_no_afe_S_ref = afe.eta_can_no_afe_S_can_no_afe_S_ref;
        CMq_can_pw = eta_can_afe_S_can_afe_S_ref*CMq_can + eta_can_no_afe_S_can_no_afe_S_ref*CMq_can;
    else
        CMq_can_pw = CMq_can*conv_can;
    end
else
    CMq_can_pw = 0;
end

% NAcelle contribution
if Nac== 1
    CMq_nac_pw = 0;
else
    CMq_nac_pw = 0;
end

% Total Derivative CMq
CMq = CMq_w1_pw + CMq_fus + CMq_HTP_pw + CMq_vee_pw + CMq_vee2_pw + CMq_can_pw + CMq_nac_pw;

% Storing DATA
Stab_Der.CMq = CMq;
%     Stab_Der.CMq_w1_pw = CMq_w1_pw;
%     Stab_Der.CMq_HTP_pw = CMq_HTP_pw;
%     Stab_Der.CMq_can_pw = CMq_can_pw;
%     Stab_Der.CMq_fus = CMq_fus;
%

% Stores Derivatives per parts
Stab_Der_parts.CM_q = CMq;
Stab_Der_parts.CM_q_w1 = CMq_w1_pw;
Stab_Der_parts.CM_q_fus = CMq_fus;
Stab_Der_parts.CM_q_HTP = CMq_HTP_pw;
Stab_Der_parts.CM_q_vee = CMq_vee_pw;
Stab_Der_parts.CM_q_vee2 = CMq_vee2_pw;
Stab_Der_parts.CM_q_can = CMq_can_pw;
Stab_Der_parts.CM_q_nac = CMq_nac_pw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Cdq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The airplane drag-coefficient-due-to-pitch-rate derivative is negligible:

CDq_w1_pw = 0;
CDq_fus = 0;
CDq_HTP_pw = 0;
CDq_vee_pw = 0;
CDq_vee2_pw = 0;
CDq_can_pw = 0;
CDq_nac_pw = 0;
% Total Derivative CDq
CDq = CDq_w1_pw + CDq_fus + CDq_HTP_pw + CDq_vee_pw + CDq_vee2_pw + CDq_can_pw + CDq_nac_pw;

% Storing DATA
Stab_Der.CDq = CDq;

% Stores Derivatives per parts
Stab_Der_parts.CD_q = CDq;
Stab_Der_parts.CD_q_w1 = CDq_w1_pw;
Stab_Der_parts.CD_q_fus = CDq_fus;
Stab_Der_parts.CD_q_HTP = CDq_HTP_pw;
Stab_Der_parts.CD_q_vee = CDq_vee_pw;
Stab_Der_parts.CD_q_vee2 = CDq_vee2_pw;
Stab_Der_parts.CD_q_can = CDq_can_pw;
Stab_Der_parts.CD_q_nac = CDq_nac_pw;


CZq = -CLq;
CXq = 0;

Stab_Der.CXq = CXq;
Stab_Der.CZq = CZq;
    end