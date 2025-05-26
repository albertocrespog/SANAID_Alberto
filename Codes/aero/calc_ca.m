%%%%%%%%%%%%%%%%%CALCULO CENTRO AERODINAMICO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_NP,XAC_WB_LE_adim]=calc_ca(AC_type,conditions,Performance,Geo_tier,Body_Geo,Stab_Der_parts, Effects,Design_criteria,AC_CONFIGURATION)

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Nac = AC_CONFIGURATION.Nac;


V = conditions.V;
a = Performance.a;
Mach = V/a;

AR_w1 = Geo_tier.AR_w1;
AR_w1_e = Geo_tier.AR_w1_e;
b_w1 = Geo_tier.b_w1;
Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
Area_body = Body_Geo.Area_body;
x_Area_body = Body_Geo.x_Area_body;
w_Area_b_max = Body_Geo.w_Area_b_max;
Area_b_max = Body_Geo.Area_b_max; 
xbar_w1 = Geo_tier.xbar_w1;
xbar_w1_e = Geo_tier.xbar_w1; % effective
x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
x_w1_LE = Geo_tier.x_w1_LE;
cmac_w1 = Geo_tier.cmac_w1;
cmac_w1_e = Geo_tier.cmac_w1;
cR_w1 = Geo_tier.cR_w1;

Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
lambda_w1 = Geo_tier.lambda_w1;
lambda_w1_e = Geo_tier.lambda_w1_e;

CLa_fus = Body_Geo.CLa_fus;
aN_w1 = CLa_fus;
a_wb_w1 = Stab_Der_parts.a_WB_w1;
a_bw_w1 = Stab_Der_parts.a_bw_w1;
a_WB_w1 = Stab_Der_parts.a_WB_w1;

CL_alpha_w1 = Stab_Der_parts.CL_alpha_w1;
CL0_w1_e_corrected = Stab_Der_parts.CL0_w1_e_corrected;
i_w1 = Design_criteria.i_w1;

if HTP == 1
CL_alpha_HTP = Stab_Der_parts.CL_alpha_htp;
downwash = Effects.downwash;
x_xbar_w2 = Geo_tier.x_xbar_w2;
end
if Vee == 1
CL_alpha_Vee = Stab_Der_parts.CL_alpha_Vee;
downwash = Effects.downwash;
x_xbar_w2 = Geo_tier.x_xbar_w2;
end
if Can == 1
CL_alpha_can = Stab_Der_parts.CL_alpha_can;
upwash = Effects.upwash;
x_xbar_can = Geo_tier.x_xbar_can;
end

beta_f = sqrt(1 - Mach^2);
Check_method = beta_f*AR_w1_e;

% Digitalizacion figura 4.3.2.2-35 Datcom
x_f_xi = [0.0,.025,.050,.075,.10,.15,.20,.30,.40,.50,.60,.70,.80];
y_f_xi = [0.0,.056,.101,.130,.152,.190,.22,.266,.301,.33,.348,.365,.375];
f_xi  = interp1(x_f_xi,y_f_xi,w_Area_b_max/b_w1,'spline');
% Revised version Álvaro
chi1    = chi_calc(b_w1, w_Area_b_max);% Performance, Stability, Dynamics and Control; Pamadi, Bandu. Pag. 193,  Figure 3.20
% Aerodynamic center of fuselage
% area derivative
dSdX = diff(Area_body)./diff(x_Area_body);
% Distance from the leading edge of the fuselage to the leading edge of the
% exposed wing root
l_N = x_cR_w1_LE;
%  Calculates where fuselage ceases to be potential
i=1;
while dSdX(i) > 0
    x_0_v = i;
    i=i+1;
end
x_0 = x_Area_body(x_0_v+1);

% Aerodynamic centers based in the leading edge
int_Sb = @(x) interp1(x_Area_body(1:end-1),dSdX,x).*(l_N - x);     %eq 3.36 PAMADI
Xac_cre_N = -(1/(cR_w1*Area_b_max))*quad(@(x) int_Sb(x),0,x_0) %eq 3.36 PAMADI

if Check_method > 4
    Xac_cre_BW = 1/4 + (b_w1 - w_Area_b_max)/(2*cR_w1)*f_xi*tan(Lambda_c4_w1) %eq 3.39 PAMADI
else
    Warning = 'WARNING!!! Revise method for estimating Xac with WB - Code in PAUSE';
    %     disp(Warning)
    %     pause
    % From Fig 3.20 obtain Xac_cre_BW for beta*Ae = 4
%     Xac_cre_BW_4 = 1/4 + (b_w1 - w_Area_b_max)/(2*cmac_w1_e)*f_xi*tan(Lambda_c4_w1); %eq 3.39 PAMADI
%     beta_Ae_4 = 4;
    % From Fig 3.21 obtain Xac_cre_BW for beta*Ae = 0
    % x-axis of Fig 3.21 (TFM fig 3.6)
%     x_axis_F_3_21 = 1/4*(AR_w1_e*(1 + lambda_w1_e)*tan(Lambda_LE_w1));
%     if x_axis_F_3_21 > 1
%         Xac_cre_BW_0 = 0.5;
%     else
%         Xac_cre_BW_0 = 0.5*x_axis_F_3_21;
%     end
%     Xac_cre_BW_0 = 0.0; % Fig 3.21 with above results
%     beta_Ae_0 = 0;
    % Linear interpolation
%     slope_Xac_cre_BW_0 = (Xac_cre_BW_4 - Xac_cre_BW_0)/(beta_Ae_4 - beta_Ae_0); % Fig 3.21
%     Xac_cre_BW = Xac_cre_BW_0 + slope_Xac_cre_BW_0*beta_f*AR_w1_e;
    %% Revisado Álvaro
    Xca_cre_BW4     = 1/4 + (b_w1 - w_Area_b_max)/(2*cR_w1)*chi1*tan(Lambda_c4_w1);
    Xca_cre_BW0     = xac_cre_BW0_calc(AR_w1_e, lambda_w1, Lambda_LE_w1);
    Xca_cre_BW      = interp1([0 4], [Xca_cre_BW0 Xca_cre_BW4], Check_method, 'linear')
end
Xac_cre_WB = get_xac_cre_wing(lambda_w1, AR_w1,Lambda_LE_w1,Mach)
Xac_cre_WB = xbar_w1/cR_w1  %eq 3.38 PAMADI % assumes the influence of the body on the location of the wing aerodynamic center is small; 25% de la cuerda hipótesis

 
XAC_WB_cre = ((Xac_cre_N*aN_w1 + Xac_cre_WB*a_wb_w1 + Xac_cre_BW*a_bw_w1)/a_WB_w1) %eq 3.32 PAMADI
% Xace at the appex (w1 croot LE)
XAC_WB_LE = XAC_WB_cre*(cR_w1)  %eq 3.36 PAMADI
XAC_WB_LE_adim = XAC_WB_LE/cmac_w1_e
% Convertin Xac appex to reference origin
XAC_WB = XAC_WB_LE + x_w1_LE                                                    %eq 3.36 PAMADI
x_xbar_wb_w1 = XAC_WB;


% corregir con cre y para dimensionalizar mult por cre 

if AC_type == 1
    CL_alpha_ac = CL_alpha_w1;
    X_NP = (CL_alpha_w1*XAC_WB)/CL_alpha_ac;
     
elseif AC_type == 2
    CL_alpha_ac = CL_alpha_w1 + CL_alpha_HTP*(downwash);
    X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_HTP*(downwash)*x_xbar_w2)/CL_alpha_ac;
  
elseif AC_type == 3
    CL_alpha_ac = CL_alpha_w1 + CL_alpha_can*(upwash) + CL_alpha_HTP*(downwash);
    X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can + CL_alpha_HTP*(downwash)*x_xbar_w2)/CL_alpha_ac;
    
elseif AC_type == 4
    CL_alpha_ac = CL_alpha_w1 + CL_alpha_Vee*(downwash);
    X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_Vee*(downwash)*x_xbar_w2)/CL_alpha_ac;
    
elseif AC_type == 5
    CL_alpha_ac = CL_alpha_w1 + CL_alpha_can*(upwash) + CL_alpha_Vee*(downwash);
    X_NP = (CL_alpha_w1*XAC_WB + CL_alpha_can*(upwash)*x_xbar_can + CL_alpha_Vee*(downwash)*x_xbar_w2)/CL_alpha_ac;
  
end
