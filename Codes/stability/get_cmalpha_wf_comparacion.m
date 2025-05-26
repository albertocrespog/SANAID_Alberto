%% Cálculo de CLalpha_wb de Roskam (Eq 3.35 Roskam Airplane Flight Dynamics pág. 116/608)

function [Cmalpha_wf, Cmalpha_WB,Cmalpha_wb_nelson] = get_cmalpha_wf_comparacion(AC_type,AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,Effects,afe)
% function [Cmalpha_wf, Cmalpha_WB,Cmalpha_wb_nelson] = get_cmalpha_wf_comparacion(AC_type,AC_CONFIGURATION,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,afe)

%SIN EFECTOS PROPULSIVOS!!

HTP = AC_CONFIGURATION.HTP;
Vee = AC_CONFIGURATION.Vee;
Can = AC_CONFIGURATION.Can;
df = Body_Geo.w_Area_b_max;
x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
b_w1 = Geo_tier.b_w1;
S_w1 = Geo_tier.S_w1; 
AR_w1_e = Geo_tier.AR_w1_e;
AR_w1 = Geo_tier.AR_w1;
S_w1_e = Geo_tier.S_w1_e;
S_ref = Geo_tier.S_ref;
CLalpha_w = Stab_Der_parts.CL_alpha_w1; % 
% CLalpha_w = 3.1154;
x_w1_LE = Geo_tier.x_w1_LE;
cmac_w1 = Geo_tier.cmac_w1;
length_x_position = Body_Geo.length_x_position;
width_x_position = Body_Geo.width_x_position;
int_bf = @(x) interp1(length_x_position,width_x_position,x);
cR_w1 = Geo_tier.cR_w1;
Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
lambda_w1 = Geo_tier.lambda_w1;
lambda_w1_e = Geo_tier.lambda_w1_e;
x_XCG = conditions.x_XCG;
% x_XCG = 0.9
V = conditions.V;
rho = Performance.rho;
a = Performance.a;
Mach = V/a;
eta_w1_afe_S_w1_afe_S_ref = afe.eta_w1_afe_S_w1_afe_S_ref;
eta_w1_no_afe_S_w1_no_afe_S_ref = afe.eta_w1_no_afe_S_w1_no_afe_S_ref;
if HTP == 1
CL_alpha_HTP = Stab_Der_parts.CL_alpha_htp;
downwash = Effects.downwash;
x_xbar_w2 = Geo_tier.x_xbar_w2;
x_w2_LE = Geo_tier.x_w2_LE;
cmac_w2 = Geo_tier.cmac_w2;
eta_w2_afe_S_w2_afe_S_ref = afe.eta_w2_afe_S_w2_afe_S_ref;
eta_w2_no_afe_S_w2_no_afe_S_ref = afe.eta_w2_no_afe_S_w2_no_afe_S_ref;
xac_w2_bar = (x_xbar_w2 - x_XCG)/cmac_w1;
x_ac_w2_bar_0 = x_xbar_w2/cmac_w1;
end
if Vee == 1
CL_alpha_Vee = Stab_Der_parts.CL_alpha_Vee;
downwash = Effects.downwash;
x_xbar_w2 = Geo_tier.x_xbar_w2;
 x_w2_LE = Geo_tier.x_w2_LE;
cmac_w2 = Geo_tier.cmac_w2;
eta_w2_afe_S_w2_afe_S_ref = afe.eta_w2_afe_S_w2_afe_S_ref;
eta_w2_no_afe_S_w2_no_afe_S_ref = afe.eta_w2_no_afe_S_w2_no_afe_S_ref;
xac_w2_bar = (x_xbar_w2 - x_XCG)/cmac_w1;
x_ac_w2_bar_0 = x_xbar_w2/cmac_w1;
end
if Can == 1
CL_alpha_can = Stab_Der_parts.CL_alpha_can;
upwash = Effects.upwash;
x_xbar_can = Geo_tier.x_xbar_can;
eta_can_afe_S_can_afe_S_ref = afe.eta_can_afe_S_can_afe_S_ref;
eta_can_no_afe_S_can_no_afe_S_ref = afe.eta_can_no_afe_S_can_no_afe_S_ref;
xac_can_bar = (x_XCG - x_xbar_can)/cmac_w1;
x_ac_can_bar_0 = x_xbar_can/cmac_w1;
end
q = 0.5*rho*V^2;
length_fus = Body_Geo.l_fus;
cmac_w1_e = Geo_tier.cmac_w1;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
% downwash = Effects.downwash;
% downwash = 0.7275;
Area_body = Body_Geo.Area_body;
x_Area_body = Body_Geo.x_Area_body;
length_fus = Body_Geo.l_fus;
w_Area_b_max = Body_Geo.w_Area_b_max;
x_Area_b_max = Body_Geo.x_Area_b_max;
Area_b_max = Body_Geo.Area_b_max;
CL_alpha_ac = Stab_Der_parts.CL_alpha_ac;

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
y_pdf_1 = 14.43; %qué son estos valores?
y_real_1 = 1;
x_d_epsu_d_alpha_1 = [0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4];
y_d_epsu_d_alpha_1 = (y_real_1/y_pdf_1)*[4.56,3.39,2.62,2.16,1.62,1.0,0.98,0.87,0.66,0.61];
% f_d_epsu_d_alpha_1  = interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,x_cre_1,'spline');

% Fuselage section right before wing
% Digitalizacion Fig 3.9(a) Pamadi
% x1bar_cre = 0.2;
y_pdf_2 = 57.88;%qué son estos valores?
y_real_2 = 4;
x_d_epsu_d_alpha_2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
y_d_epsu_d_alpha_2 = (y_real_2/y_pdf_2)*[52.03,36.43,28.43,22.99,18.88,16.22,14.37,13.03,11.70,10.98];

f_d_epsu_d_alpha_f1  = @(x) interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,((x_w1_LE) - x)/cR_w1,'linear','extrap'); 
int_low_1 = 0;
int_up_1 = (x_w1_LE - cR_w1/3);

if int_up_1 > 0
    % if there are sections computes the integral
    int_Cmalpha_1 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f1(x)),int_low_1,int_up_1);
else
    % No sections
    int_Cmalpha_1 = 0;
end
dM_dalpha_f1 = (q/36.5)*(CLalpha_w*D2R/0.0785)*int_Cmalpha_1
CM_alpha_f1 = (1/(36.5*S_w1*cmac_w1))*(CLalpha_w*D2R/0.0785)*int_Cmalpha_1

f_d_epsu_d_alpha_f2  = @(x) interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,((x_w1_LE) - x)/cR_w1,'linear','extrap');


% Determines iing root chord - cR_w1
if int_up_1 > 0
    % if there are sections computes the integral with the limits
    int_low_2 = (x_w1_LE - cR_w1/3);
    int_up_2 = x_w1_LE;
    int_Cmalpha_2 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f2(x)),int_low_2,int_up_2);
else
    % if there are sections no sections farther away, computes if there are any portions right aft the wing less that 1 wing root chord from LE
    int_low_2 = 0;
    int_up_2 = x_w1_LE;
    if (int_up_2- int_low_2)>0
    int_Cmalpha_2 = quad(@(x) ((int_bf(x)).^2).*(1+f_d_epsu_d_alpha_f2(x)),int_low_2,int_up_2);
    else
        int_Cmalpha_2 = 0;
    end
end

dM_dalpha_f2 = (q/36.5)*(CLalpha_w*D2R/0.0785)*int_Cmalpha_2
CM_alpha_f2 = (1/(36.5*S_w1*cmac_w1))*(CLalpha_w*D2R/0.0785)*int_Cmalpha_2


int_low_3 = x_w1_LE + cR_w1; %¿No debería ser cR_w1 o x_w1_TE?
int_up_3 = length_fus;
% Determines if there are sections behind the wing.
% The limit is 1 wing root chord - cR_w1
if (int_up_3-int_low_3) > 0
    % if there are sections computes the integral
    x_3 = (x_w2_LE  + cmac_w2/4) - (x_w1_LE + cR_w1);
    f_d_epsu_d_alpha_f3 = @(x) ((x-(x_w1_LE + cR_w1))/x_3).*downwash;
    int_Cmalpha_3 = quad(@(x) ((int_bf(x)).^2).*(f_d_epsu_d_alpha_f3(x)),int_low_3,int_up_3);
else
    % No sections
    int_Cmalpha_3 = 0;
end
dM_dalpha_f3 = (q/36.5)*int_Cmalpha_3
CM_alpha_f3 = (1/(36.5*S_w1*cmac_w1))*int_Cmalpha_3

dM_dalpha = dM_dalpha_f1 + dM_dalpha_f2 + dM_dalpha_f3;
CM_alpha_fus_roskam = CM_alpha_f1 + CM_alpha_f2 + CM_alpha_f3
% xac_w_bar = 0.25; %% Revisar, introducir la digitalización (+ adelante)
Xac_cre_W_B = get_xac_cre_wing(lambda_w1, AR_w1,Lambda_LE_w1,Mach)
xac_w_bar = Xac_cre_W_B*cR_w1/cmac_w1;
Delta_xac_f_bar = -dM_dalpha/(q*S_w1*cmac_w1*CLalpha_w*D2R)
xac_wf_bar = xac_w_bar + Delta_xac_f_bar ;
Cmalpha_wf = CLalpha_wf*(xcg_bar - xac_wf_bar);

% Cálculo del centro aerodinámico:
x_ac_wf_bar_0 = xac_wf_bar + x_w1_LE/cmac_w1;
xac_w_bar_0 = xac_w_bar + x_w1_LE/cmac_w1;
%Eq 3.40 
if AC_type == 1
X_NP_bar = x_ac_wf_bar_0
X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 2
X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_HTP/CLalpha_wf)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_wf)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 3
X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_can/CLalpha_wf)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_wf)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_wf)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_wf)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 4
X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_Vee/CLalpha_wf)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_wf)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP = X_NP_bar*cmac_w1;
elseif AC_type == 5
X_NP_bar = (x_ac_wf_bar_0 + (CL_alpha_can/CLalpha_wf)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_wf)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_wf)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_wf)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP = X_NP_bar*cmac_w1;
end

if AC_type == 1
X_NP_no_fus_bar = xac_w_bar_0;
X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 2
X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 3
X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 4
X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
elseif AC_type == 5
X_NP_no_fus_bar = (xac_w_bar_0 + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_no_fus = X_NP_no_fus_bar*cmac_w1;
end

Delta_X_NP_ac_fus = X_NP - X_NP_no_fus
Percentage_inc_X_NP = ((Delta_X_NP_ac_fus)/X_NP_no_fus)*100
SM_roskam = (X_NP - x_XCG)/cmac_w1
SM_no_fus_roskam = (X_NP_no_fus - x_XCG)/cmac_w1
CMalpha_ac_roskam = -CL_alpha_ac*SM_roskam 
CMalpha_ac_roskam_no_fus = -CL_alpha_ac*SM_no_fus_roskam

%% Pamadi: Cmalpha = cmalpha_WB + cmalpha_t + cmalpha_c 
%Versión WB para small wing-span-to-body-diameter ratios (pag 184 Pamadi,seccion 3.3.3)

%Cmalpha_WB = (xcg_bar - xac_WB)*CLalpha_WB; eq. (3.29)

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
CLalpha_N = 2*f_k2_k1*(Area_b_max/S_ref) %eq 3.26
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
int_Sb = @(x) interp1(x_Area_body(1:end-1),dSdX,x).*(l_N - x);     %eq 3.36 PAMADI
Xac_cre_N = -(1/(cR_w1*Area_b_max))*quad(@(x) int_Sb(x),0,x_0); %eq 3.36 PAMADI

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
    Xca_cre_B_W      = interp1([0 4], [Xca_cre_BW0 Xca_cre_BW4], Check_method, 'linear')
end

Xac_cre_W_B = get_xac_cre_wing(lambda_w1, AR_w1,Lambda_LE_w1,Mach)
% Xac_cre_W_B = (0.25*cmac_w1)/cR_w1; %eq 3.38 PAMADI % assumes the influence of the body on the location of the wing aerodynamic center is small; 25% de la cuerda hipótesis
xac_w = Xac_cre_W_B*cR_w1/cmac_w1
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
X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_HTP/CLalpha_WB)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_WB)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_pam = X_NP_pam_bar*cmac_w1;
elseif AC_type == 3
X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_can/CLalpha_WB)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_WB)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_WB)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_WB)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_pam = X_NP_pam_bar*cmac_w1;
elseif AC_type == 4
X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_Vee/CLalpha_WB)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_WB)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_pam = X_NP_pam_bar*cmac_w1;
elseif AC_type == 5
X_NP_pam_bar = (x_ac_WB_bar_0_pam + (CL_alpha_can/CLalpha_WB)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_WB)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_WB)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_WB)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_pam = X_NP_pam_bar*cmac_w1;
end

if AC_type == 1
X_NP_nofus_pam_bar = xac_w_bar_0_pam;
X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
elseif AC_type == 2
X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
elseif AC_type == 3
X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
elseif AC_type == 4
X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
elseif AC_type == 5
X_NP_nofus_pam_bar = (xac_w_bar_0_pam + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_nofus_pam = X_NP_nofus_pam_bar*cmac_w1;
end


Delta_X_NP_ac_fus_pam = X_NP_pam - X_NP_nofus_pam
Percentage_inc_X_NP_pam = ((Delta_X_NP_ac_fus_pam)/X_NP_nofus_pam)*100
SM_pam = (X_NP_pam - x_XCG)/cmac_w1
SM_no_fus_pam = (X_NP_nofus_pam - x_XCG)/cmac_w1
CMalpha_ac_pam = -CL_alpha_ac*SM_pam
CMalpha_ac_pam_no_fus = -CL_alpha_ac*SM_no_fus_pam


%% Pamadi wing-span-to-body-diameter altos, se calcula el efecto del fuselaje y ala
% por separado:
Xac_cre_W_B = get_xac_cre_wing(lambda_w1, AR_w1,Lambda_LE_w1,Mach)
% Xac_cre_W_B = (0.25*cmac_w1)/cR_w1; %eq 3.38 PAMADI % assumes the influence of the body on the location of the wing aerodynamic center is small; 25% de la cuerda hipótesis
xac_w = Xac_cre_W_B*cR_w1/cmac_w1
x_a_bar = xcg_bar - xac_w;
CMalpha_wing = CLalpha_w*x_a_bar;
CMalpha_fus = Stab_Der_parts.CM_alpha_fus
CM_alpha_wing_fus = CMalpha_wing + CMalpha_fus

% Cálculo del centro aerodinámico:
xac_w_bar_0_pam_2 = xac_w + x_w1_LE/cmac_w1;
%Eq 3.40 
if AC_type == 1
X_NP_pam_bar_2 = xac_w_bar_0_pam_2
X_NP_pam_2 = X_NP_pam_bar_2*cmac_w1;
elseif AC_type == 2
X_NP_pam_bar_2 = (xac_w_bar_0_pam_2 - CMalpha_fus/CLalpha_w + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_pam_2 = X_NP_pam_bar_2*cmac_w1;
elseif AC_type == 3
X_NP_pam_bar_2 = (xac_w_bar_0_pam_2 - CMalpha_fus/CLalpha_w + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_pam_2 = X_NP_pam_bar_2*cmac_w1;
elseif AC_type == 4
X_NP_pam_bar_2 = (xac_w_bar_0_pam_2 - CMalpha_fus/CLalpha_w  + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_pam_2 = X_NP_pam_bar_2*cmac_w1;
elseif AC_type == 5
X_NP_pam_bar_2 = (xac_w_bar_0_pam_2 - CMalpha_fus/CLalpha_w + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_pam_2 = X_NP_pam_bar_2*cmac_w1;
end

if AC_type == 1
X_NP_nofus_pam_bar_2 = xac_w_bar_0_pam_2;
X_NP_nofus_pam_2 = X_NP_nofus_pam_bar_2*cmac_w1;
elseif AC_type == 2
X_NP_nofus_pam_bar_2 = (xac_w_bar_0_pam_2 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_nofus_pam_2 = X_NP_nofus_pam_bar_2*cmac_w1;
elseif AC_type == 3
X_NP_nofus_pam_bar_2 = (xac_w_bar_0_pam_2 + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_nofus_pam_2 = X_NP_nofus_pam_bar_2*cmac_w1;
elseif AC_type == 4
X_NP_nofus_pam_bar_2 = (xac_w_bar_0_pam_2 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_nofus_pam_2 = X_NP_nofus_pam_bar_2*cmac_w1;
elseif AC_type == 5
X_NP_nofus_pam_bar_2 = (xac_w_bar_0_pam_2 + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
X_NP_nofus_pam_2 = X_NP_nofus_pam_bar_2*cmac_w1;
end



Delta_X_NP_ac_fus_pam_2 = X_NP_pam_2 - X_NP_nofus_pam_2
Percentage_inc_X_NP_pam_2 = ((Delta_X_NP_ac_fus_pam_2)/X_NP_nofus_pam_2)*100
SM_pam_2 = (X_NP_pam_2 - x_XCG)/cmac_w1
SM_no_fus_pam_2 = (X_NP_nofus_pam_2 - x_XCG)/cmac_w1
CMalpha_ac_pam2 = -CL_alpha_ac*SM_pam_2
CMalpha_ac_pam2_no_fus = -CL_alpha_ac*SM_no_fus_pam_2


%% Pamadi con la fórmula de Roskam/Nelson CMalpha_fus = 1/(36.5*S*c)*int_0^length_fus (bf^2*dE_dalpha*dx)
% wing_body_CLalpha = 0.0785; 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cm_alpha%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Nose section of Fuselage
% % Digitalizacion Fig 3.9(b) Pamadi
% % x_cre_1 = 0.8;
% y_pdf_1 = 14.43; %qué son estos valores?
% y_real_1 = 1;
% x_d_epsu_d_alpha_1 = [0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4];
% y_d_epsu_d_alpha_1 = (y_real_1/y_pdf_1)*[4.56,3.39,2.62,2.16,1.62,1.0,0.98,0.87,0.66,0.61]*(CLalpha_w*D2R/(wing_body_CLalpha));
% % f_d_epsu_d_alpha_1  = interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,x_cre_1,'spline');
% 
% % Fuselage section right before wing
% % Digitalizacion Fig 3.9(a) Pamadi
% % x1bar_cre = 0.2;
% y_pdf_2 = 57.88;%qué son estos valores?
% y_real_2 = 4;
% x_d_epsu_d_alpha_2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
% y_d_epsu_d_alpha_2 = (y_real_2/y_pdf_2)*[52.03,36.43,28.43,22.99,18.88,16.22,14.37,13.03,11.70,10.98]*(CLalpha_w*D2R/(wing_body_CLalpha));
% % f_d_epsu_d_alpha_2  = interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,x1bar_cre,'spline');
% 
% % selection of segments
% % x1bar_sec = x_w1_LE - cmac_w_e;
% % x1bar_sec = x_w1_LE - cR_w1;  %% x1bar_sec = cR_w1, no?
% 
% % wing_body_CLalpha = 0.0785; %% IMP: PER DEGREE!
% 
% % int_bf = @(x) interp1(length_x_positi
% % Section from tip-nose to front (from graphic)
% % f_d_epsu_d_alpha_f1  = @(x) interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,x,'linear','extrap').*((x_w1_LE) - x);
% % Método III eq 3.7 PAMADI
% f_d_epsu_d_alpha_f1  = @(x) interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,((x_w1_LE) - x)/cR_w1,'linear','extrap'); 
% int_low_1 = 0;
% int_up_1 = (x_w1_LE - cR_w1/3);
% % Determines if there are sections farther away from the LE of the wing.
% % The limit is 1 wing root chord - cR_w1
% if int_up_1 > 0
%     % if there are sections computes the integral
%     int_Cmalpha_1_sinuno = quad(@(x) ((int_bf(x)).^2).*((f_d_epsu_d_alpha_f1(x))),int_low_1,int_up_1);
% else
%     % No sections
%     int_Cmalpha_1_sinuno = 0;
% end
% %   CMalpha_f1 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_1;
% CMalpha_f1_sinuno = (1/(36.5*S_w1*cmac_w1))*int_Cmalpha_1_sinuno
% 
% 
% 
% 
% % Section from right before wing (from graphic)
% % f_d_epsu_d_alpha_f2  = @(x) interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,x,'linear','extrap').*((x_w1_LE) - x);
% f_d_epsu_d_alpha_f2  = @(x) interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,((x_w1_LE) - x)/cR_w1,'linear','extrap');
% 
% 
% % Determines iing root chord - cR_w1
% if int_up_1 > 0
%     % if there are sections computes the integral with the limits
%     int_low_2 = (x_w1_LE - cR_w1/3);
%     int_up_2 = x_w1_LE;
%     int_Cmalpha_2_sinuno = quad(@(x) ((int_bf(x)).^2).*((f_d_epsu_d_alpha_f2(x))),int_low_2,int_up_2);
% else
%     % if there are sections no sections farther away, computes if there are any portions right aft the wing less that 1 wing root chord from LE
%     int_low_2 = 0;
%     int_up_2 = x_w1_LE;
%     if (int_up_2- int_low_2)>0
%     int_Cmalpha_2_sinuno = quad(@(x) ((int_bf(x)).^2).*((f_d_epsu_d_alpha_f2(x))),int_low_2,int_up_2);
%     else
%         int_Cmalpha_2_sinuno = 0;
%     end
% end
% 
% %   CMalpha_f2 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_2;
% CMalpha_f2_sinuno = (1/(36.5*S_w1*cmac_w1))*int_Cmalpha_2_sinuno
% 
% 
% 
% % Section from behind the wing (from eq 3.10 PAMADI)
% int_low_3 = x_w1_LE + cR_w1; %¿No debería ser cR_w1 o x_w1_TE?
% int_up_3 = length_fus;
% % Determines if there are sections behind the wing.
% % The limit is 1 wing root chord - cR_w1
% if (int_up_3-int_low_3) > 0
%     % if there are sections computes the integral
%     x_3 = (x_w2_LE  + cmac_w2/4) - (x_w1_LE + cR_w1);
%     f_d_epsu_d_alpha_f3 = @(x) ((x-(x_w1_LE + cR_w1))/x_3).*downwash;
%     int_Cmalpha_3_sinuno = quad(@(x) ((int_bf(x)).^2).*(f_d_epsu_d_alpha_f3(x)),int_low_3,int_up_3);
% else
%     % No sections
%     int_Cmalpha_3_sinuno = 0;
% end
% 
% %   CMalpha_f3 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_3;
% % CMalpha_f3 = (a_WB_w1/(wing_body_CLalpha))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_3;
% CMalpha_f3_sinuno = (1/(36.5*S_w1*cmac_w1))*int_Cmalpha_3_sinuno %Both the curves in Figs. 3.9a and 3.9b are based on the wing-body lift-curve
%                                                   % slope of 0.0785/deg. For any other values of wing-body lift-curve slope, multiply
%                                                   % the values of deu/da obtained from Fig. 3.9a or 3.9b by the factor CLa_wb/0.0785
% 
% 
% 
% C_Malpha_fus_sinuno = (CMalpha_f1_sinuno + CMalpha_f2_sinuno + CMalpha_f3_sinuno)*R2D
% % Stab_Der_parts.CM_alpha_fus = C_Malpha_fus;%f there are sections right aft the LE of the wing.
% % The limit is 1 w
% 
% % Cálculo del centro aerodinámico:
% xac_w_bar_0_pam_3 = xac_w + x_w1_LE/cmac_w1;
% %Eq 3.40 
% if AC_type == 1
% X_NP_pam_bar_3 = xac_w_bar_0_pam_3
% X_NP_pam_3 = X_NP_pam_bar_3*cmac_w1;
% elseif AC_type == 2
% X_NP_pam_bar_3 = (xac_w_bar_0_pam_3 - C_Malpha_fus_sinuno/CLalpha_w + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% X_NP_pam_3 = X_NP_pam_bar_3*cmac_w1;
% elseif AC_type == 3
% X_NP_pam_bar_3 = (xac_w_bar_0_pam_3 - C_Malpha_fus_sinuno/CLalpha_w + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% X_NP_pam_3 = X_NP_pam_bar_3*cmac_w1;
% elseif AC_type == 4
% X_NP_pam_bar_3 = (xac_w_bar_0_pam_3 - C_Malpha_fus_sinuno/CLalpha_w  + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% X_NP_pam_3 = X_NP_pam_bar_3*cmac_w1;
% elseif AC_type == 5
% X_NP_pam_bar_3 = (xac_w_bar_0_pam_3 - C_Malpha_fus_sinuno/CLalpha_w + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% X_NP_pam_3 = X_NP_pam_bar_3*cmac_w1;
% end
% 
% if AC_type == 1
% X_NP_nofus_pam_bar_3 = xac_w_bar_0_pam_3;
% X_NP_nofus_pam_3 = X_NP_nofus_pam_bar_3*cmac_w1;
% elseif AC_type == 2
% X_NP_nofus_pam_bar_3 = (xac_w_bar_0_pam_3 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% X_NP_nofus_pam_3 = X_NP_nofus_pam_bar_3*cmac_w1;
% elseif AC_type == 3
% X_NP_nofus_pam_bar_3 = (xac_w_bar_0_pam_3 + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% X_NP_nofus_pam_3 = X_NP_nofus_pam_bar_3*cmac_w1;
% elseif AC_type == 4
% X_NP_nofus_pam_bar_3 = (xac_w_bar_0_pam_3 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% X_NP_nofus_pam_3 = X_NP_nofus_pam_bar_3*cmac_w1;
% elseif AC_type == 5
% X_NP_nofus_pam_bar_3 = (xac_w_bar_0_pam_3 + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% X_NP_nofus_pam_3 = X_NP_nofus_pam_bar_3*cmac_w1;
% end
% 
% 
% 
% Delta_X_NP_ac_fus_pam_3 = X_NP_pam_3 - X_NP_nofus_pam_3
% Percentage_inc_X_NP_pam_3 = ((Delta_X_NP_ac_fus_pam_3)/X_NP_nofus_pam_3)*100
% SM_pam_3 = (X_NP_pam_3 - x_XCG)/cmac_w1
% SM_no_fus_pam_3 = (X_NP_nofus_pam_3 - x_XCG)/cmac_w1
% CMalpha_ac_pam3 = -CL_alpha_ac*SM_pam_3
% CMalpha_ac_pam3_no_fus = -CL_alpha_ac*SM_no_fus_pam_3
% % % % % %% Nelson: 
% % % % % 
% % % % % 
% % % % % CM_dalpha_f1 = ((CLalpha_w*D2R/0.0785)/(36.5*S_w1*cmac_w1))*int_Cmalpha_1;
% % % % % CM_dalpha_f2 = ((CLalpha_w*D2R/0.0785)/(36.5*S_w1*cmac_w1))*int_Cmalpha_2;
% % % % % CM_dalpha_f3 = ((1)/(36.5*S_w1*cmac_w1))*int_Cmalpha_3;
% % % % % Cmalpha_fus_nelson = CM_dalpha_f1 + CM_dalpha_f2 + CM_dalpha_f3
% % % % % % Cmalpha_fus_nelson = Cmalpha_fus_nelson
% % % % % % CLalpha_w = Stab_Der_parts.CL_alpha_w1; % 
% % % % % % CLalpha_w = 3.1154;
% % % % % xcg_bar = (x_XCG - x_w1_LE)/cmac_w1;
% % % % % % xac_bar = 0.25;
% % % % % Xac_cre_W_B = get_xac_cre_wing(lambda_w1, AR_w1,Lambda_LE_w1,Mach)
% % % % % xac_bar = Xac_cre_W_B*cR_w1/cmac_w1;
% % % % % Cmalpha_wing_nelson = CLalpha_w*(xcg_bar - xac_bar);
% % % % % Cmalpha_wb_nelson = CLalpha_w*(xcg_bar - xac_bar) + Cmalpha_fus_nelson ;
% % % % % 
% % % % % % Cálculo del centro aerodinámico:
% % % % % xac_w_bar_0_nelson = xac_bar + x_w1_LE/cmac_w1;
% % % % % %Eq 3.40 
% % % % % if AC_type == 1
% % % % % X_NP_nelson_bar = xac_w_bar_0_nelson
% % % % % X_NP_nelson = X_NP_nelson_bar*cmac_w1;
% % % % % elseif AC_type == 2
% % % % % X_NP_nelson_bar = (xac_w_bar_0_nelson - Cmalpha_fus_nelson/CLalpha_w + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% % % % % X_NP_nelson = X_NP_nelson_bar*cmac_w1;
% % % % % elseif AC_type == 3
% % % % % X_NP_nelson_bar = (xac_w_bar_0_nelson - Cmalpha_fus_nelson/CLalpha_w + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% % % % % X_NP_nelson = X_NP_nelson_bar*cmac_w1;
% % % % % elseif AC_type == 4
% % % % % X_NP_nelson_bar = (xac_w_bar_0_nelson - Cmalpha_fus_nelson/CLalpha_w  + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% % % % % X_NP_nelson = X_NP_nelson_bar*cmac_w1;
% % % % % elseif AC_type == 5
% % % % % X_NP_nelson_bar = (xac_w_bar_0_nelson - Cmalpha_fus_nelson/CLalpha_w + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% % % % % X_NP_nelson = X_NP_nelson_bar*cmac_w1;
% % % % % end
% % % % % 
% % % % % if AC_type == 1
% % % % % X_NP_nofus_nelson_bar = xac_w_bar_0_nelson;
% % % % % X_NP_no_fus_nelson = X_NP_nofus_nelson_bar*cmac_w1;
% % % % % elseif AC_type == 2
% % % % % X_NP_nofus_nelson_bar = (xac_w_bar_0_nelson + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% % % % % X_NP_no_fus_nelson = X_NP_nofus_nelson_bar*cmac_w1;
% % % % % elseif AC_type == 3
% % % % % X_NP_nofus_nelson_bar = (xac_w_bar_0_nelson + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_HTP/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% % % % % X_NP_no_fus_nelson = X_NP_nofus_nelson_bar*cmac_w1;
% % % % % elseif AC_type == 4
% % % % % X_NP_nofus_nelson_bar = (xac_w_bar_0_nelson + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% % % % % X_NP_no_fus_nelson = X_NP_nofus_nelson_bar*cmac_w1;
% % % % % elseif AC_type == 5
% % % % % X_NP_nofus_nelson_bar = (xac_w_bar_0_nelson + (CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*x_ac_can_bar_0*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*x_ac_w2_bar_0*downwash)/(1 +(CL_alpha_can/CLalpha_w)*eta_can_no_afe_S_can_no_afe_S_ref*upwash + (CL_alpha_Vee/CLalpha_w)*eta_w2_no_afe_S_w2_no_afe_S_ref*downwash)
% % % % % X_NP_no_fus_nelson = X_NP_nofus_nelson_bar*cmac_w1;
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % Delta_X_NP_ac_fus_nelson = X_NP_nelson - X_NP_no_fus_nelson
% % % % % Percentage_inc_X_NP_nelson = ((Delta_X_NP_ac_fus_nelson)/X_NP_no_fus_nelson)*100
% % % % % SM_nelson = (X_NP_nelson - x_XCG)/cmac_w1
% % % % % SM_no_fus_nelson = (X_NP_no_fus_nelson - x_XCG)/cmac_w1
% % % % % CMalpha_ac_nelson = -CL_alpha_ac*SM_nelson
% % % % % CMalpha_ac_nelson_no_fus = -CL_alpha_ac*SM_no_fus_nelson


% CMalpha_ac_v =[CMalpha_ac_roskam CMalpha_ac_pam CMalpha_ac_pam2 CMalpha_ac_pam3 CMalpha_ac_nelson]
% CMalpha_ac_nofus_v=[CMalpha_ac_roskam_no_fus CMalpha_ac_pam_no_fus CMalpha_ac_pam2_no_fus CMalpha_ac_pam3_no_fus CMalpha_ac_nelson_no_fus]
% SM_fus_v =[SM_roskam SM_pam SM_pam_2 SM_pam_3 SM_nelson]
% SM_no_fus_v =[SM_no_fus_roskam SM_no_fus_pam SM_no_fus_pam_2 SM_no_fus_pam_3 SM_no_fus_nelson]
% percentage_XNP = [Percentage_inc_X_NP Percentage_inc_X_NP_pam Percentage_inc_X_NP_pam_2 Percentage_inc_X_NP_pam_3 Percentage_inc_X_NP_nelson ]
% X_NP_v = [X_NP X_NP_pam X_NP_pam_2 X_NP_pam_3 X_NP_nelson]
% X_NP_nofus_v = [X_NP_no_fus  X_NP_nofus_pam X_NP_nofus_pam_2 X_NP_nofus_pam_3 X_NP_no_fus_nelson]


CMalpha_ac_v =[CMalpha_ac_roskam CMalpha_ac_pam CMalpha_ac_pam2 CMalpha_ac_nelson]
CMalpha_ac_nofus_v=[CMalpha_ac_roskam_no_fus CMalpha_ac_pam_no_fus CMalpha_ac_pam2_no_fus CMalpha_ac_nelson_no_fus]
SM_fus_v =[SM_roskam SM_pam SM_pam_2 SM_nelson]
SM_no_fus_v =[SM_no_fus_roskam SM_no_fus_pam SM_no_fus_pam_2 SM_no_fus_nelson]
percentage_XNP = [Percentage_inc_X_NP Percentage_inc_X_NP_pam Percentage_inc_X_NP_pam_2 Percentage_inc_X_NP_nelson ]
X_NP_v = [X_NP X_NP_pam X_NP_pam_2 X_NP_nelson]
X_NP_nofus_v = [X_NP_no_fus  X_NP_nofus_pam X_NP_nofus_pam_2 X_NP_no_fus_nelson]
end
 
