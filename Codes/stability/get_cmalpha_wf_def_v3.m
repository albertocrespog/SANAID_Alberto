%% Cálculo de CLalpha_wb de Roskam (Eq 3.35 Roskam Airplane Flight Dynamics pág. 116/608)

function [Trim_ITER,Stab_Der_parts,TRIM_RESULTS,Stab_Der] = get_cmalpha_wf_def_v3(AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,...
    Effects,afe,SM_des,OUTPUT_read_XLSX,Stab_Der,Trim_ITER)
% function [Cmalpha_wf, Cmalpha_pamadi] = get_cmalpha_wf_comparacion(AC_type,AC_CONFIGURATION,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,afe)

%SIN EFECTOS PROPULSIVOS!!
%% Input
W1      = AC_CONFIGURATION.W1;
HTP     = AC_CONFIGURATION.HTP;
VTP     = AC_CONFIGURATION.VTP;
Can     = AC_CONFIGURATION.Can;
Vee     = AC_CONFIGURATION.Vee;
Vee2     = AC_CONFIGURATION.Vee2;
Nac     = AC_CONFIGURATION.Nac;
AC_type = AC_CONFIGURATION.AC_type;


df                = Body_Geo.w_Area_b_max;
prop_wash_effect  = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;
x_cR_w1_LE        = Geo_tier.x_cR_w1_LE;
b_w1              = Geo_tier.b_w1;
S_w1              = Geo_tier.S_w1;
AR_w1_e           = Geo_tier.AR_w1_e;
AR_w1             = Geo_tier.AR_w1;
S_w1_e            = Geo_tier.S_w1_e;
S_ref             = Geo_tier.S_ref;
CLalpha_w         = Stab_Der_parts.CL_alpha_w1; %
x_w1_LE           = Geo_tier.x_w1_LE;
cmac_w1           = Geo_tier.cmac_w1;
length_x_position = Body_Geo.length_x_position;
width_x_position  = Body_Geo.width_x_position;
int_bf            = @(x) interp1(length_x_position,width_x_position,x);
cR_w1             = Geo_tier.cR_w1;
Lambda_c4_w1      = Geo_tier.Lambda_c4_w1;
Lambda_LE_w1      = Geo_tier.Lambda_LE_w1;
lambda_w1         = Geo_tier.lambda_w1;
lambda_w1_e       = Geo_tier.lambda_w1_e;
x_XCG             = conditions.x_XCG;
V                 = conditions.V;
rho               = Performance.rho;
a                 = Performance.a;
Mach              = V/a;

if prop_wash_effect == 1
    if W1 == 1
        eta_w1_afe_S_w1_afe_S_ref = afe.eta_w1_afe_S_w1_afe_S_ref;
        eta_w1_no_afe_S_w1_no_afe_S_ref = afe.eta_w1_no_afe_S_w1_no_afe_S_ref;
    end
    if HTP == 1
        eta_HTP_afe_S_HTP_afe_S_ref = afe.eta_HTP_afe_S_HTP_afe_S_ref;
        eta_HTP_no_afe_S_HTP_no_afe_S_ref = afe.eta_HTP_no_afe_S_HTP_no_afe_S_ref;
    end
    if Vee == 1
        eta_vee_afe_S_vee_afe_S_ref = afe.eta_vee_afe_S_vee_afe_S_ref;
        eta_vee_no_afe_S_vee_no_afe_S_ref = afe.eta_vee_no_afe_S_vee_no_afe_S_ref;
    end
    if Vee2 == 1
        eta_vee2_afe_S_vee2_afe_S_ref = afe.eta_vee2_afe_S_vee2_afe_S_ref;
        eta_vee2_no_afe_S_vee2_no_afe_S_ref = afe.eta_vee2_no_afe_S_vee2_no_afe_S_ref;
    end
    if Can == 1
        eta_can_afe_S_can_afe_S_ref = afe.eta_can_afe_S_can_afe_S_ref;
        eta_can_no_afe_S_can_no_afe_S_ref = afe.eta_can_no_afe_S_can_no_afe_S_ref;
    end    
end

if HTP == 1
    CL_alpha_HTP  = Stab_Der_parts.CL_alpha_HTP;
    downwash_HTP      = Effects.downwash_HTP;
    x_xbar_HTP     = Geo_tier.x_xbar_HTP;
    x_HTP_LE       = Geo_tier.x_HTP_LE;
    cmac_HTP       = Geo_tier.cmac_HTP;
    xac_HTP_bar    = (x_xbar_HTP - x_XCG)/cmac_w1;
    x_ac_HTP_bar_0 = x_xbar_HTP/cmac_w1;
end

if Vee == 1
    CL_alpha_vee  = Stab_Der_parts.CL_alpha_vee;
    downwash_vee      = Effects.downwash_vee;
    x_xbar_vee     = Geo_tier.x_xbar_vee;
    x_vee_LE       = Geo_tier.x_vee_LE;
    cmac_vee       = Geo_tier.cmac_vee;
    xac_vee_bar    = (x_xbar_vee - x_XCG)/cmac_w1;
    x_ac_vee_bar_0 = x_xbar_vee/cmac_w1;
end

if Vee2 == 1
    CL_alpha_vee2  = Stab_Der_parts.CL_alpha_vee2;
    downwash_vee2      = Effects.downwash_vee2;
    x_xbar_vee2     = Geo_tier.x_xbar_vee2;
    x_vee2_LE       = Geo_tier.x_vee2_LE;
    cmac_vee2       = Geo_tier.cmac_vee2;
    xac_vee2_bar    = (x_xbar_vee2 - x_XCG)/cmac_w1;
    x_ac_vee2_bar_0 = x_xbar_vee2/cmac_w1;
end

if Can == 1
    CL_alpha_can   = Stab_Der_parts.CL_alpha_can;
    upwash         = Effects.upwash;
    x_xbar_can     = Geo_tier.x_xbar_can;
    xac_can_bar    = (x_XCG - x_xbar_can)/cmac_w1;
    x_ac_can_bar_0 = x_xbar_can/cmac_w1;
end


q            = 0.5*rho*V^2;
cmac_w1_e    = Geo_tier.cmac_w1_e;
R2D          = conv_UNITS.R2D;
D2R          = conv_UNITS.D2R;
Area_body    = Body_Geo.Area_body;
x_Area_body  = Body_Geo.x_Area_body;
length_fus   = Body_Geo.l_fus;
w_Area_b_max = Body_Geo.w_Area_b_max;
x_Area_b_max = Body_Geo.x_Area_b_max;
Area_b_max   = Body_Geo.Area_b_max;
CL_alpha_ac  = Stab_Der_parts.CL_alpha_ac;

%% ROSKAM
% En este código se calcula CLalpha_wf*(xcg_bar - xac_wf_bar) (el resto de
% términos son idénticos en las formulaciones alternativas (Nelson y
% Pamadi)

% Cmalpha_roskam = CLalpha_wf*(xcg_bar - xac_wf_bar) -
% CLalpha_t*eta_t*St/S*(xact_bar-xcg_bar)*downwash +
% CLapha_c*eta_c*Sc/S*(xcg_bar-xac_c_bar)*upwash;

% Cmalpha_wf = CLalpha_wf*(xcg_bar - xac_wf_bar)
%donde:

k_wf = 1 + 0.025*(df/b_w1)-0.25*(df/b_w1)^2;
CLalpha_wf = k_wf*CLalpha_w;

xcg_bar = (x_XCG - x_w1_LE)/cmac_w1;

%Importante: Roskam define el xac_wf_bar = xac,w + Delta_xac_f
%donde Delta_xac_f es la variación del centro aerodinámico del conjunto
%ala-fuselaje debido al aumento en el ángulo de ataque que supone la
%circulación del flujo sobre el fuselaje.
% Se define en Roskam Part VI Sección 8.2.5.3. pág 357/582 del pdf. (325
% del libro)

% Nose section of Fuselage
% Digitalizacion Fig 3.9(b) Pamadi
% x_cre_1 = 0.8;
y_pdf_1 = 14.43; 
y_real_1 = 1;
x_d_epsu_d_alpha_1 = [0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4];
y_d_epsu_d_alpha_1 = (y_real_1/y_pdf_1)*[4.56,3.39,2.62,2.16,1.62,1.0,0.98,0.87,0.66,0.61];
% f_d_epsu_d_alpha_1  = interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,x_cre_1,'spline');

% Fuselage section right before wing
% Digitalizacion Fig 3.9(a) Pamadi
% x1bar_cre = 0.2;
y_pdf_2 = 57.88;
y_real_2 = 4;
x_d_epsu_d_alpha_2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
y_d_epsu_d_alpha_2 = (y_real_2/y_pdf_2)*[52.03,36.43,28.43,22.99,18.88,16.22,14.37,13.03,11.70,10.98];

f_d_epsu_d_alpha_f1  = @(x) interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,((x_w1_LE) - x)/cR_w1,'linear','extrap');
int_low_1 = 0;
int_up_1 = (x_w1_LE - cR_w1/3);

if int_up_1 > 0
    % if there are sections computes the integral
    int_Cmalpha_1 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f1(x)),int_low_1,int_up_1);
    int_Cmalpha_1 = integral(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f1(x)),int_low_1,int_up_1);
else
    % No sections
    int_Cmalpha_1 = 0;
end
dM_dalpha_f1 = (q/36.5)*(CLalpha_w*D2R/0.0785)*int_Cmalpha_1;
CM_alpha_f1 = (1/(36.5*S_w1*cmac_w1))*(CLalpha_w*D2R/0.0785)*int_Cmalpha_1;

f_d_epsu_d_alpha_f2  = @(x) interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,((x_w1_LE) - x)/cR_w1,'linear','extrap');

% Determines iing root chord - cR_w1
if int_up_1 > 0
    % if there are sections computes the integral with the limits
    int_low_2 = (x_w1_LE - cR_w1/3);
    int_up_2 = x_w1_LE;
    int_Cmalpha_2 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f2(x)),int_low_2,int_up_2);
    int_Cmalpha_2 = integral(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f2(x)),int_low_2,int_up_2);
else
    % if there are sections no sections farther away, computes if there are any portions right aft the wing less that 1 wing root chord from LE
    int_low_2 = 0;
    int_up_2 = x_w1_LE;
    if (int_up_2- int_low_2)>0
        int_Cmalpha_2 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f2(x)),int_low_2,int_up_2);
        int_Cmalpha_2 = integral(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f2(x)),int_low_2,int_up_2);
    else
        int_Cmalpha_2 = 0;
    end
end

dM_dalpha_f2 = (q/36.5)*(CLalpha_w*D2R/0.0785)*int_Cmalpha_2;
CM_alpha_f2 = (1/(36.5*S_w1*cmac_w1))*(CLalpha_w*D2R/0.0785)*int_Cmalpha_2;

int_low_3 = x_w1_LE + cR_w1; %¿No debería ser cR_w1 o x_w1_TE?
int_up_3 = length_fus;

% Determines if there are sections behind the wing.
% The limit is 1 wing root chord - cR_w1
if (int_up_3-int_low_3) > 0
    % if there are sections computes the integral
    if AC_type == 6 % No rear fuselage sections after wing
        int_Cmalpha_3 = 0;
    elseif HTP == 1 
        x_xbar_HTP = Geo_tier.x_xbar_HTP;
        x_HTP_LE = Geo_tier.x_HTP_LE;
        cmac_HTP = Geo_tier.cmac_HTP;
        xac_HTP_bar = (x_xbar_HTP - x_XCG)/cmac_w1;
        x_ac_HTP_bar_0 = x_xbar_HTP/cmac_w1;
        x_3 = (x_HTP_LE  + cmac_HTP/4) - (x_w1_LE + cR_w1);
        f_d_epsu_d_alpha_f3 = @(x) ((x-(x_w1_LE + cR_w1))/x_3).*downwash_HTP;
        int_Cmalpha_3 = quad(@(x) ((int_bf(x)).^2).*(f_d_epsu_d_alpha_f3(x)),int_low_3,int_up_3);
        int_Cmalpha_3 = integral(@(x) ((int_bf(x)).^2).*(f_d_epsu_d_alpha_f3(x)),int_low_3,int_up_3);
    elseif Vee == 1
        x_xbar_vee = Geo_tier.x_xbar_vee;
        x_vee_LE = Geo_tier.x_vee_LE;
        cmac_vee = Geo_tier.cmac_vee;
        xac_vee_bar = (x_xbar_vee - x_XCG)/cmac_w1;
        x_ac_vee_bar_0 = x_xbar_vee/cmac_w1;
        x_3 = (x_vee_LE  + cmac_vee/4) - (x_w1_LE + cR_w1);
        f_d_epsu_d_alpha_f3 = @(x) ((x-(x_w1_LE + cR_w1))/x_3).*downwash_vee;
        int_Cmalpha_3 = quad(@(x) ((int_bf(x)).^2).*(f_d_epsu_d_alpha_f3(x)),int_low_3,int_up_3);
        int_Cmalpha_3 = integral(@(x) ((int_bf(x)).^2).*(f_d_epsu_d_alpha_f3(x)),int_low_3,int_up_3);
     elseif Vee2 == 1
        x_xbar_vee2 = Geo_tier.x_xbar_vee2;
        x_vee2_LE = Geo_tier.x_vee2_LE;
        cmac_vee2 = Geo_tier.cmac_vee2;
        xac_vee2_bar = (x_xbar_vee2 - x_XCG)/cmac_w1;
        x_ac_vee2_bar_0 = x_xbar_vee2/cmac_w1;
        x_3 = (x_vee2_LE  + cmac_vee2/4) - (x_w1_LE + cR_w1);
        f_d_epsu_d_alpha_f3 = @(x) ((x-(x_w1_LE + cR_w1))/x_3).*downwash_vee2;
        int_Cmalpha_3 = quad(@(x) ((int_bf(x)).^2).*(f_d_epsu_d_alpha_f3(x)),int_low_3,int_up_3);
        int_Cmalpha_3 = integral(@(x) ((int_bf(x)).^2).*(f_d_epsu_d_alpha_f3(x)),int_low_3,int_up_3);
   else
        % No sections
        int_Cmalpha_3 = 0;
    end
else
    % Added for the case (int_up_3-int_low_3) < 0
    % No sections
    int_Cmalpha_3 = 0;
end

dM_dalpha_f3 = (q/36.5)*int_Cmalpha_3;
CM_alpha_f3 = (1/(36.5*S_w1*cmac_w1))*int_Cmalpha_3;

dM_dalpha = dM_dalpha_f1 + dM_dalpha_f2 + dM_dalpha_f3;

CM_alpha_fus_roskam = CM_alpha_f1 + CM_alpha_f2 + CM_alpha_f3;
% According to DATCOM all results are in per degree, hence correction
CM_alpha_fus_roskam = CM_alpha_fus_roskam*R2D;


% xac_w_bar = 0.25; %% Revisar, introducir la digitalización (+ adelante)
% Xac_cre_W_B = get_xac_cre_wing(lambda_w1, AR_w1,Lambda_LE_w1,Mach);
% xac_w_bar = Xac_cre_W_B*cR_w1/cmac_w1;
xac_w_bar = Geo_tier.xbar_w1/cmac_w1;
% Geo_tier.x_xbar_w1 = Xac_cre_W_B*cR_w1 + x_w1_LE;

% Geo_tier.xbar_w1 = Xac_cre_W_B*cR_w1;
Delta_xac_f_bar = -dM_dalpha/(q*S_w1*cmac_w1*CLalpha_w*D2R);
xac_wf_bar = xac_w_bar + Delta_xac_f_bar ;
Cmalpha_wf = CLalpha_wf*(xcg_bar - xac_wf_bar);

% Cálculo del centro aerodinámico:
x_ac_wf_bar_0 = xac_wf_bar + x_w1_LE/cmac_w1;
xac_w_bar_0 = xac_w_bar + x_w1_LE/cmac_w1;
%Eq 3.40
if AC_type == 1
    X_NP_bar = x_ac_wf_bar_0;
    X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 2
    X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_HTP/CLalpha_wf)*x_ac_HTP_bar_0*downwash_HTP)/(1 + (CL_alpha_HTP/CLalpha_wf)*downwash_HTP);
    X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 3
    X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_can/CLalpha_wf)*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_wf)*x_ac_HTP_bar_0*downwash_HTP)/(1 +(CL_alpha_can/CLalpha_wf)*upwash + (CL_alpha_HTP/CLalpha_wf)*downwash_HTP);
    X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 4
    X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_vee/CLalpha_wf)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_wf)*downwash_vee);
    X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 5
    X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_can/CLalpha_wf)*x_ac_can_bar_0*upwash + (CL_alpha_vee/CLalpha_wf)*x_ac_vee_bar_0*downwash_vee)/(1 +(CL_alpha_can/CLalpha_wf)*upwash + (CL_alpha_vee/CLalpha_wf)*downwash_vee);
    X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 6
    X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_can/CLalpha_wf)*x_ac_can_bar_0*upwash )/(1 +(CL_alpha_can/CLalpha_wf)*upwash);
    X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 7
    X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_vee/CLalpha_wf)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_wf)*downwash_vee);
    X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 8
    X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_vee/CLalpha_wf)*x_ac_vee_bar_0*downwash_vee + (CL_alpha_vee2/CLalpha_wf)*x_ac_vee2_bar_0*downwash_vee2)/(1 + (CL_alpha_vee/CLalpha_wf)*downwash_vee + (CL_alpha_vee2/CLalpha_wf)*downwash_vee2);
    X_NP = X_NP_bar*cmac_w1;
end

if AC_type == 1
    X_NP_no_fus_bar = xac_w_bar_0;
    X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 2
    X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_HTP/CLalpha_w)*x_ac_HTP_bar_0*downwash_HTP)/(1 + (CL_alpha_HTP/CLalpha_w)*downwash_HTP);
    X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 3
    X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*x_ac_HTP_bar_0*downwash_HTP)/(1 +(CL_alpha_can/CLalpha_w)*upwash + (CL_alpha_HTP/CLalpha_w)*downwash_HTP);
    X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 4
    X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee);
    X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 5
    X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 +(CL_alpha_can/CLalpha_w)*upwash + (CL_alpha_vee/CLalpha_w)*downwash_vee);
    X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 6
    X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash)/(1 +(CL_alpha_can/CLalpha_w)*upwash);
    X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 7
    X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee);
    X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 8
    X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee + (CL_alpha_vee2/CLalpha_w)*x_ac_vee2_bar_0*downwash_vee2)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee + (CL_alpha_vee2/CLalpha_w)*downwash_vee2);
    X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
end
 
Delta_X_NP_ac_fus = X_NP - X_NP_no_fus;
Percentage_inc_X_NP = ((Delta_X_NP_ac_fus)/X_NP_no_fus)*100;
SM_roskam = (X_NP - x_XCG)/cmac_w1;
SM_no_fus_roskam = (X_NP_no_fus - x_XCG)/cmac_w1;
CMalpha_ac_roskam = -CL_alpha_ac*SM_roskam;
CMalpha_ac_roskam_no_fus = -CL_alpha_ac*SM_no_fus_roskam;

%% Pamadi: Cmalpha = cmalpha_WB + cmalpha_t + cmalpha_c
%Versión WB para small wing-span-to-body-diameter ratios (pag 184 Pamadi,seccion 3.3.3)
% Pamadi (2004) notes that for configurations with relatively large values of the ratio of wing span to fuselage diameter ( b / d > 2) the mutual interference effects between the fuselage and the wing are small
% and can be neglected so that the effects of the body and the wing can be determined individually and summed. In general, experimental results suggest that the maximum lift coefficient for the wing alone and 
% for the wing-body combination typically differs by 5% or less. The lift curve slope is likewise affected little by the fuselage when b/d is large.
% Pag. 64, Sforza, Pasquale M.. Commercial Airplane Design Principles, Elsevier Science & Technology, 2014. ProQuest Ebook Central

wingspan2bodydiam = b_w1/w_Area_b_max;
% wingspan2bodydiam = 2.1;
%Cmalpha_WB = (xcg_bar - xac_WB)*CLalpha_WB; eq. (3.29)
if wingspan2bodydiam <= 2 %%Criterio Pamadi: Pasquale Sforza, in Commercial Airplane Design Principles, 2014, p.64
    %xcg_bar from LE, eq 3.30
    xcg_bar = (x_XCG - x_w1_LE)/cmac_w1;
   
    Fineness_Ratio = length_fus/w_Area_b_max;


    %digitaliazacion figura PAMADI CAP3 Fig 3.6
    x_f_k2_k1 = [4.,5.,6.,8.,10.,12.,14.,16.,18.,20.];
    y_f_k2_k1 = [.77,.825,.865,.91,.94,.955,.965,.97,.973,.975];
    f_k2_k1  = interp1(x_f_k2_k1,y_f_k2_k1,Fineness_Ratio,'spline');
    x1 = x_Area_b_max;
    % CLalpha_WB; eq. (3.24)
    %     C*********       KW(B) VS D/B  DATA FIGURE 4.3.1.2-10-A

    X10A=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    Y10A=[1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0];
    KWB = interp1(X10A, Y10A, w_Area_b_max/b_w1);

    % C ****         KB(W) VS D/B  DATA  FIGURE 4.3.1.2-10-B

    X10B=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    Y10B=[0.0,0.13,0.29,0.45,0.62,0.80,1.0,1.22,1.45,1.70,2.0];
    KBW = interp1(X10B, Y10B, w_Area_b_max/b_w1);
    CLalpha_N = 2*f_k2_k1*(Area_b_max/S_ref); %eq 3.26
    CLalpha_N = Body_Geo.CLa_fus;
    % CLalpha_N = 0.3227
    KN = (CLalpha_N/CLalpha_w)*(S_w1/S_w1_e); %eq 3.25
    CLalpha_WB =(KN + KWB + KBW)*CLalpha_w*(S_w1_e/S_w1);% eq. (3.24)
    % CLalpha_WB = Stab_Der_parts.a_WB_w1;
    CLalpha_W_B = CLalpha_w*(S_w1_e/S_w1)*KWB; %eq 3.33
    CLalpha_B_W = CLalpha_w*(S_w1_e/S_w1)*KBW;%eq 3.34

    % Aerodynamic center of fuselage
    % area derivative
    dSdX = diff(Area_body)./diff(x_Area_body);
    % From Calc_Bpdy_Geometry_Dec218
    dSdX = Body_Geo.dSdX;

    % Distance from the leading edge of the fuselage to the leading edge of the
    % exposed wing root
    l_N = x_cR_w1_LE;
    %  Calculates where fuselage ceases to be potential
    % i=1;
    % while dSdX(i) > 0
    %     x_0_v = i;
    %     i=i+1;
    % end
    % x_0 = x_Area_body(x_0_v+1);
    [max_Area_body,pos] = max(Area_body);
    x_0 = x_Area_body(pos);
    % Aerodynamic centers based in the leading edge
    int_Sb = @(x) interp1(x_Area_body(1:end),dSdX(1:end),x).*(l_N - x);     %eq 3.36 PAMADI
    Xac_cre_N = -(1/(cR_w1*Area_b_max))*quad(@(x) int_Sb(x),0,x_0); %eq 3.36 PAMADI
    Xac_cre_N = -(1/(cR_w1*Area_b_max))*integral(@(x) int_Sb(x),0,x_0); %eq 3.36 PAMADI

    % Digitalizacion figura 4.3.2.2-35 Datcom
    x_f_xi = [0.0,.025,.050,.075,.10,.15,.20,.30,.40,.50,.60,.70,.80];
    y_f_xi = [0.0,.056,.101,.130,.152,.190,.22,.266,.301,.33,.348,.365,.375];
    f_xi  = interp1(x_f_xi,y_f_xi,w_Area_b_max/b_w1,'spline');

    beta_f = sqrt(1 - Mach^2);
    Check_method = beta_f*AR_w1_e;

    if Check_method > 4
        Xac_cre_B_W = 1/4 + (b_w1 - w_Area_b_max)*f_xi*tan(Lambda_c4_w1)/(2*cR_w1); %eq 3.39 PAMADI
    else
        Warning = 'WARNING!!! Revise method for estimating Xac with WB - Code in PAUSE';

        %% Revisado Álvaro
        Xca_cre_BW4     = 1/4 + (b_w1 - w_Area_b_max)/(2*cR_w1)*f_xi*tan(Lambda_c4_w1);
        Xca_cre_BW0     = xac_cre_BW0_calc(AR_w1_e, lambda_w1, Lambda_LE_w1);
        Xca_cre_B_W      = interp1([0 4], [Xca_cre_BW0 Xca_cre_BW4], Check_method, 'linear');
    end

    %     Xac_cre_W_B = get_xac_cre_wing(lambda_w1, AR_w1,Lambda_LE_w1,Mach);
    % %     Xac_cre_W_B = (0.25*cmac_w1)/cR_w1; %eq 3.38 PAMADI % assumes the influence of the body on the location of the wing aerodynamic center is small; 25% de la cuerda hipótesis
    %     xac_w = Xac_cre_W_B*cR_w1/cmac_w1;
    xac_w = Geo_tier.xbar_w1/cmac_w1;
    xac_WB_cre = (Xac_cre_N*CLalpha_N + Xac_cre_W_B*CLalpha_W_B + Xac_cre_B_W*CLalpha_B_W)/CLalpha_WB;
    xac_WB_bar = xac_WB_cre*cR_w1/cmac_w1;
    Cmalpha_WB = (xcg_bar - xac_WB_bar)*CLalpha_WB;
    % Cálculo del centro aerodinámico:
    x_ac_WB_bar_0_pam = xac_WB_bar + x_w1_LE/cmac_w1;
    xac_w_bar_0_pam = xac_w + x_w1_LE/cmac_w1;
    %Eq 3.40
    if AC_type == 1
        X_NP_pam_bar = x_ac_WB_bar_0_pam
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 2
        X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_HTP/CLalpha_WB)*x_ac_HTP_bar_0*downwash_HTP)/(1 + (CL_alpha_HTP/CLalpha_WB)*downwash_HTP);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 3
        X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_can/CLalpha_WB)*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_WB)*x_ac_HTP_bar_0*downwash_HTP)/(1 +(CL_alpha_can/CLalpha_WB)*upwash + (CL_alpha_HTP/CLalpha_WB)*downwash_HTP);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 4
        X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_vee/CLalpha_WB)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_WB)*downwash_vee);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 5
        X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_can/CLalpha_WB)*x_ac_can_bar_0*upwash + (CL_alpha_vee/CLalpha_WB)*x_ac_vee_bar_0*downwash_vee)/(1 +(CL_alpha_can/CLalpha_WB)*upwash + (CL_alpha_vee/CLalpha_WB)*downwash_vee);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 6
        X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_can/CLalpha_WB)*x_ac_can_bar_0*upwash)/(1 +(CL_alpha_can/CLalpha_WB)*upwash);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 7
        X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_vee/CLalpha_WB)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_WB)*downwash_vee);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 8
        X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_vee/CLalpha_WB)*x_ac_vee_bar_0*downwash_vee + (CL_alpha_vee2/CLalpha_WB)*x_ac_vee2_bar_0*downwash_vee2)/(1 + (CL_alpha_vee/CLalpha_WB)*downwash_vee + (CL_alpha_vee2/CLalpha_WB)*downwash_vee2);
        X_NP_pam = X_NP_pam_bar*cmac_w1;

    end


    if AC_type == 1
        X_NP_nofus_pam_bar = xac_w_bar_0_pam;
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 2
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_HTP/CLalpha_w)*x_ac_HTP_bar_0*downwash_HTP)/(1 + (CL_alpha_HTP/CLalpha_w)*downwash_HTP);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 3
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*x_ac_HTP_bar_0*downwash_HTP)/(1 +(CL_alpha_can/CLalpha_w)*upwash + (CL_alpha_HTP/CLalpha_w)*downwash_HTP);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 4
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 5
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 +(CL_alpha_can/CLalpha_w)*upwash + (CL_alpha_vee/CLalpha_w)*downwash_vee);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 6
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash)/(1 +(CL_alpha_can/CLalpha_w)*upwash);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 7
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 8
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee + (CL_alpha_vee2/CLalpha_w)*x_ac_vee2_bar_0*downwash_vee2)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee + (CL_alpha_vee2/CLalpha_w)*downwash_vee2);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;

    end


    Delta_X_NP_ac_fus_pam = X_NP_pam - X_NP_nofus_pam;
    Percentage_inc_X_NP_pam = ((Delta_X_NP_ac_fus_pam)/X_NP_nofus_pam)*100;
    SM_pam = (X_NP_pam - x_XCG)/cmac_w1;
    SM_no_fus_pam = (X_NP_nofus_pam - x_XCG)/cmac_w1;
    CMalpha_ac_pam = -CL_alpha_ac*SM_pam;
    CMalpha_ac_pam_no_fus = -CL_alpha_ac*SM_no_fus_pam;

    CMalpha_fus = 0;

    %% Pamadi wing-span-to-body-diameter altos, se calcula el efecto del fuselaje y ala
    % por separado:
else
    %     Xac_cre_W_B = get_xac_cre_wing(lambda_w1, AR_w1,Lambda_LE_w1,Mach);
    % %     Xac_cre_W_B = (0.25*cmac_w1)/cR_w1; %eq 3.38 PAMADI % assumes the influence of the body on the location of the wing aerodynamic center is small; 25% de la cuerda hipótesis
    %     xac_w = Xac_cre_W_B*cR_w1/cmac_w1;
    %
    xac_w = Geo_tier.xbar_w1/cmac_w1;

    x_a_bar = xcg_bar - xac_w;
    CMalpha_wing = CLalpha_w*x_a_bar;

    wing_body_CLalpha = 0.0785;

    CMalpha_f1 = (CLalpha_w*D2R/(wing_body_CLalpha))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_1;

    CM_alpha_f1 = (1/(36.5*S_w1*cmac_w1))*(CLalpha_w*D2R/0.0785)*int_Cmalpha_1;

    CMalpha_f2 = (CLalpha_w*D2R/(wing_body_CLalpha))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_2;

    CMalpha_f3 = (pi/(2*S_w1*cmac_w1))*int_Cmalpha_3; %Both the curves in Figs. 3.9a and 3.9b are based on the wing-body lift-curve
    % slope of 0.0785/deg. For any other values of wing-body lift-curve slope, multiply
    % the values of deu/da obtained from Fig. 3.9a or 3.9b by the factor CLa_wb/0.0785



    CMalpha_fus = CMalpha_f1 + CMalpha_f2 + CMalpha_f3;

    CM_alpha_wing_fus = CMalpha_wing + CMalpha_fus;

    % Cálculo del centro aerodinámico:
    xac_w_bar_0_pam = xac_w + x_w1_LE/cmac_w1;
    %Eq 3.40
    if AC_type == 1
        X_NP_pam_bar = xac_w_bar_0_pam;
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 2
        X_NP_pam_bar = (xac_w_bar_0_pam - CMalpha_fus/CLalpha_w + (CL_alpha_HTP/CLalpha_w)*x_ac_HTP_bar_0*downwash_HTP)/(1 + (CL_alpha_HTP/CLalpha_w)*downwash_HTP);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 3
        X_NP_pam_bar = (xac_w_bar_0_pam - CMalpha_fus/CLalpha_w + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*x_ac_HTP_bar_0*downwash_HTP)/(1 +(CL_alpha_can/CLalpha_w)*upwash + (CL_alpha_HTP/CLalpha_w)*downwash_HTP);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 4
        X_NP_pam_bar = (xac_w_bar_0_pam - CMalpha_fus/CLalpha_w  + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 5
        X_NP_pam_bar = (xac_w_bar_0_pam - CMalpha_fus/CLalpha_w + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 +(CL_alpha_can/CLalpha_w)*upwash + (CL_alpha_vee/CLalpha_w)*downwash_vee);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 6
        X_NP_pam_bar = (xac_w_bar_0_pam - CMalpha_fus/CLalpha_w + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash)/(1 +(CL_alpha_can/CLalpha_w)*upwash);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 7
        X_NP_pam_bar = (xac_w_bar_0_pam - CMalpha_fus/CLalpha_w + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee);
        X_NP_pam = X_NP_pam_bar*cmac_w1;
    elseif AC_type == 8
        X_NP_pam_bar = (xac_w_bar_0_pam - CMalpha_fus/CLalpha_w + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee + (CL_alpha_vee2/CLalpha_w)*x_ac_vee2_bar_0*downwash_vee2)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee + (CL_alpha_vee2/CLalpha_w)*downwash_vee2);
        X_NP_pam = X_NP_pam_bar*cmac_w1;

    end

    if AC_type == 1
        X_NP_nofus_pam_bar = xac_w_bar_0_pam;
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 2
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_HTP/CLalpha_w)*x_ac_HTP_bar_0*downwash_HTP)/(1 + (CL_alpha_HTP/CLalpha_w)*downwash_HTP);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 3
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*x_ac_HTP_bar_0*downwash_HTP)/(1 +(CL_alpha_can/CLalpha_w)*upwash + (CL_alpha_HTP/CLalpha_w)*downwash_HTP);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 4
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 5
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 +(CL_alpha_can/CLalpha_w)*upwash + (CL_alpha_vee/CLalpha_w)*downwash_vee);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 6
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_can/CLalpha_w)*x_ac_can_bar_0*upwash)/(1 +(CL_alpha_can/CLalpha_w)*upwash);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 7
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    elseif AC_type == 8
        X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_vee/CLalpha_w)*x_ac_vee_bar_0*downwash_vee + (CL_alpha_vee2/CLalpha_w)*x_ac_vee2_bar_0*downwash_vee2)/(1 + (CL_alpha_vee/CLalpha_w)*downwash_vee + (CL_alpha_vee2/CLalpha_w)*downwash_vee2);
        X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
    end


    Delta_x_ac_fus = -CMalpha_fus/CLalpha_w;
    Delta_X_NP_ac_fus_pam = X_NP_pam - X_NP_nofus_pam;
    Percentage_inc_X_NP_pam = ((Delta_X_NP_ac_fus_pam)/X_NP_nofus_pam)*100;
    SM_pam = (X_NP_pam - x_XCG)/cmac_w1;
    SM_no_fus_pam = (X_NP_nofus_pam - x_XCG)/cmac_w1;
    CMalpha_ac_pam = -CL_alpha_ac*SM_pam;
    CMalpha_ac_pam_no_fus = -CL_alpha_ac*SM_no_fus_pam;

end

%% Contribution of Fuselage
wingspan2bodydiam = b_w1/w_Area_b_max;
flagwingspan2bodydiam = OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam; 
Munk_fuselage_constribution = OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution;
if   wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1 || Munk_fuselage_constribution == 1
    x_XCG_fus_des = X_NP_pam - SM_des*cmac_w1;
    CM_alpha_ac_des = -CL_alpha_ac*SM_des;
    x_XCG_des = x_XCG_fus_des;
    X_NP = X_NP_pam;

    % Calculates the CMalpha
    CMalpha_ac_pam = -CL_alpha_ac*SM_pam;
    SM_pam = SM_pam;
else
    x_XCG_No_fus_des = X_NP_nofus_pam - SM_des*cmac_w1;
    CM_alpha_ac_des = -CL_alpha_ac*SM_des;
    x_XCG_des = x_XCG_No_fus_des;
    X_NP = X_NP_nofus_pam;

    % Calculates the CMalpha
    CMalpha_ac_pam = - CL_alpha_ac*SM_no_fus_pam;
    SM_pam = SM_no_fus_pam;
end

Trim_ITER.X_NP_No_fus  = X_NP_nofus_pam;
Trim_ITER.SM_Excel_No_fus = SM_no_fus_pam;
Trim_ITER.CM_alpha_ac_No_fus = CMalpha_ac_pam_no_fus;
Trim_ITER.CM_alpha_ac_fus = CMalpha_ac_pam;
Trim_ITER.SM_Excel_fus = SM_pam;
Trim_ITER.X_NP_fus = X_NP_pam;
Trim_ITER.Delta_X_NP_ac_fus = Delta_X_NP_ac_fus_pam;
Trim_ITER.Percentage_inc_X_NP = Percentage_inc_X_NP_pam;



CM_alpha_ac= CMalpha_ac_pam;
TRIM_RESULTS.CM_alpha_ac = CMalpha_ac_pam; %CON EL XCG DE CONDITIONS
TRIM_RESULTS.X_NP = X_NP; %el calculado con pamadi, con o sin la contrib del fuselaje
TRIM_RESULTS.SM = SM_pam; %CON EL X_xCG DE CONDITIONS, EL INICIAL.
Stab_Der.CM_alpha_ac = CMalpha_ac_pam; %CON EL XCG DE CONDITIONS
Trim_ITER.CMalpha_f = CMalpha_fus; %pamadi
Stab_Der.CL_alpha_ac = CL_alpha_ac;

Stab_Der_parts.CL_alpha_ac = CL_alpha_ac;
Stab_Der_parts.CM_alpha_ac = CM_alpha_ac;
Stab_Der_parts.CM_alpha_fus = CMalpha_fus;

TRIM_RESULTS.CM_alpha_ac_des = CM_alpha_ac_des; %con el xcg DESEADO
TRIM_RESULTS.x_XCG_des = x_XCG_des; % el xcg deseado
TRIM_RESULTS.SM_des = SM_des; %con el xcg deseado
Stab_Der.CM_alpha_ac_des = CM_alpha_ac_des;%con el xcg deseado

% %----------------------------------------------------------------------
% CMalpha_ac_v =[CMalpha_ac_roskam CMalpha_ac_pam]
% CMalpha_ac_nofus_v=[CMalpha_ac_roskam_no_fus CMalpha_ac_pam_no_fus]
% SM_fus_v =[SM_roskam SM_pam ]
% SM_no_fus_v =[SM_no_fus_roskam SM_no_fus_pam]
% percentage_XNP = [Percentage_inc_X_NP Percentage_inc_X_NP_pam]
% X_NP_v = [X_NP X_NP_pam]
% X_NP_nofus_v = [X_NP_no_fus  X_NP_nofus_pam]
% Delta_x_ac_fus_v = [Delta_xac_f_bar Delta_x_ac_fus]
% pause
end
