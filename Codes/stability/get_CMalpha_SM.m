function [Stab_Der_parts,Trim_ITER,Stab_Der,TRIM_RESULTS] = get_CMalpha_SM(AC_CONFIGURATION,Geo_tier,Body_Geo,Stab_Der,Stab_Der_parts,conv_UNITS,Effects,OUTPUT_read_XLSX,SM_des,X_NP,conditions)
% Nota: Me pareció que salían valores altos comparando con la literatura,
% no obstante en el archivo ejemplo_calculo_cmalphafus_3_1_pamadi.m he
% comprobado que la estructura y la programación da valores del mismo orden
% que el ejemplo, o sea que en el caso de que haya algo mal, se comprobará
% más adelante identificando que todos los términos que contribuyen al
% cálculo de CMalpha_fus sean correctos:Sale tan alto por tener una S_w1
% tan pequeña.

% Nota: Esto tiene que ir despues de las propulsivas:

HTP = AC_CONFIGURATION.HTP;
Vee = AC_CONFIGURATION.Vee;
x_XCG = conditions.x_XCG;
x_w1_LE = Geo_tier.x_w1_LE;
cR_w1 = Geo_tier.cR_w1;
length_x_position = Body_Geo.length_x_position;
width_x_position = Body_Geo.width_x_position;
int_bf = @(x) interp1(length_x_position,width_x_position,x);
downwash = Effects.downwash;
% CL_alpha_wb_w1 = Stab_Der_parts.CL_alpha_wb_w1;
% a_WB_w1 = CL_alpha_wb_w1;
R2D = conv_UNITS.R2D;
S_w1 = Geo_tier.S_w1; 
cmac_w1 = Geo_tier.cmac_w1;
cmac_w1_e = Geo_tier.cmac_w1;
if HTP == 1 | Vee == 1
x_w2_LE = Geo_tier.x_w2_LE;
cmac_w2 = Geo_tier.cmac_w2;
end
length_fus = Body_Geo.l_fus;
CLalpha_w = Stab_Der_parts.CL_alpha_w1;
CMTalpha = Stab_Der.CMtalpha;
CL_alpha_ac = Stab_Der_parts.CL_alpha_ac;
wing_body_CLalpha = 0.0785; 

D2R=conv_UNITS.D2R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cm_alpha%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% f_d_epsu_d_alpha_2  = interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,x1bar_cre,'spline');

% selection of segments
% x1bar_sec = x_w1_LE - cmac_w_e;
% x1bar_sec = x_w1_LE - cR_w1;  %% x1bar_sec = cR_w1, no?

% wing_body_CLalpha = 0.0785; %% IMP: PER DEGREE!

% int_bf = @(x) interp1(length_x_positi
% Section from tip-nose to front (from graphic)
% f_d_epsu_d_alpha_f1  = @(x) interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,x,'linear','extrap').*((x_w1_LE) - x);
% Método III eq 3.7 PAMADI
f_d_epsu_d_alpha_f1  = @(x) interp1(x_d_epsu_d_alpha_1,y_d_epsu_d_alpha_1,((x_w1_LE) - x)/cR_w1,'linear','extrap'); 
int_low_1 = 0;
int_up_1 = (x_w1_LE - cR_w1/3);
% Determines if there are sections farther away from the LE of the wing.
% The limit is 1 wing root chord - cR_w1
if int_up_1 > 0
    % if there are sections computes the integral
    int_Cmalpha_1 = quad(@(x) ((int_bf(x)).^2).*(1+(f_d_epsu_d_alpha_f1(x))),int_low_1,int_up_1)
else
    % No sections
    int_Cmalpha_1 = 0;
end
%   CMalpha_f1 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_1;
CMalpha_f1 = (CLalpha_w*D2R/(wing_body_CLalpha))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_1




% Section from right before wing (from graphic)
% f_d_epsu_d_alpha_f2  = @(x) interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,x,'linear','extrap').*((x_w1_LE) - x);
f_d_epsu_d_alpha_f2  = @(x) interp1(x_d_epsu_d_alpha_2,y_d_epsu_d_alpha_2,((x_w1_LE) - x)/cR_w1,'linear','extrap');


% Determines iing root chord - cR_w1
if int_up_1 > 0
    % if there are sections computes the integral with the limits
    int_low_2 = (x_w1_LE - cR_w1/3);
    int_up_2 = x_w1_LE;
    int_Cmalpha_2 = quad(@(x) ((int_bf(x)).^2).*(1+(f_d_epsu_d_alpha_f2(x))),int_low_2,int_up_2)
else
    % if there are sections no sections farther away, computes if there are any portions right aft the wing less that 1 wing root chord from LE
    int_low_2 = 0;
    int_up_2 = x_w1_LE;
    if (int_up_2- int_low_2)>0
    int_Cmalpha_2 = quad(@(x) ((int_bf(x)).^2).*(1+(f_d_epsu_d_alpha_f2(x))),int_low_2,int_up_2)
    else
        int_Cmalpha_2 = 0;
    end
end

%   CMalpha_f2 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_2;
CMalpha_f2 = (CLalpha_w*D2R/(wing_body_CLalpha))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_2



% Section from behind the wing (from eq 3.10 PAMADI)
int_low_3 = x_w1_LE + cR_w1; %¿No debería ser cR_w1 o x_w1_TE?
int_up_3 = length_fus;
% Determines if there are sections behind the wing.
% The limit is 1 wing root chord - cR_w1
if (int_up_3-int_low_3) > 0
    % if there are sections computes the integral
    x_3 = (x_w2_LE  + cmac_w2/4) - (x_w1_LE + cR_w1);
    f_d_epsu_d_alpha_f3 = @(x) ((x-(x_w1_LE + cR_w1))/x_3).*downwash;
    int_Cmalpha_3 = quad(@(x) ((int_bf(x)).^2).*(f_d_epsu_d_alpha_f3(x)),int_low_3,int_up_3)
else
    % No sections
    int_Cmalpha_3 = 0;
end

%   CMalpha_f3 = (a_WB/(wing_body_CLalpha*R2D))*(pi*f_k2_k1/(2*S_w*cmac_w))*int_Cmalpha_3;
% CMalpha_f3 = (a_WB_w1/(wing_body_CLalpha))*(pi/(2*S_w1*cmac_w1))*int_Cmalpha_3;
CMalpha_f3 = (pi/(2*S_w1*cmac_w1))*int_Cmalpha_3 %Both the curves in Figs. 3.9a and 3.9b are based on the wing-body lift-curve
                                                  % slope of 0.0785/deg. For any other values of wing-body lift-curve slope, multiply
                                                  % the values of deu/da obtained from Fig. 3.9a or 3.9b by the factor CLa_wb/0.0785



C_Malpha_fus = CMalpha_f1 + CMalpha_f2 + CMalpha_f3
Stab_Der_parts.CM_alpha_fus = C_Malpha_fus;%f there are sections right aft the LE of the wing.
% The limit is 1 w






%% NOTE IMPORTANT
%% Estimation of SM with only CLalpha contribution
SM_Excel_No_fus = (X_NP - x_XCG)/cmac_w1
x_XCG_Excel = x_XCG;
Trim_ITER.x_XCG_Excel = x_XCG_Excel;
Trim_ITER.X_NP = X_NP;
Trim_ITER.SM_Excel_No_fus = SM_Excel_No_fus;



%% NOTE Corrections of x_XCG taking into account fuselage CMalpha contribution
% CM_alpha CONTRIBUTION ONLY TAKING INTO ACCOUNT CLalpha

CM_alpha_ac_No_fus = -CL_alpha_ac*SM_Excel_No_fus
CM_alpha_ac_fus = CM_alpha_ac_No_fus + CMTalpha + C_Malpha_fus % Total CMalpha with fuselage and propulsive contribution - Uses Munk's contribution
SM_Excel_fus = -CM_alpha_ac_fus/CL_alpha_ac % Real SM with fuselage contribution

% Shift in wing+fuselage aerodynamic center from wing aerodynamic center
% caused by the so called Munk effect
X_NP_fus = (x_XCG - cmac_w1*CM_alpha_ac_fus/CL_alpha_ac);
Delta_X_NP_ac_fus = X_NP_fus - X_NP;
Percentage_inc_X_NP = ((Delta_X_NP_ac_fus)/X_NP)*100;

Trim_ITER.CM_alpha_ac_No_fus = CM_alpha_ac_No_fus;
Trim_ITER.CM_alpha_ac_fus = CM_alpha_ac_fus;
Trim_ITER.SM_Excel_fus = SM_Excel_fus
Trim_ITER.X_NP_fus = X_NP_fus;
Trim_ITER.Delta_X_NP_ac_fus = Delta_X_NP_ac_fus;
Trim_ITER.Percentage_inc_X_NP = Percentage_inc_X_NP;

% Adjusting x_XCG to ensure longitudinal stability due to the fuselage
% contribution

% Gathers all the flags
Munk_fuselage_constribution = OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution;
if Munk_fuselage_constribution == 1;
    x_XCG_fus = X_NP_fus - SM_des*cmac_w1;
    CM_alpha_ac = CM_alpha_ac_fus;
    x_XCG = x_XCG_fus;
    X_NP = X_NP + Delta_X_NP_ac_fus;
    SM = (X_NP - x_XCG)/cmac_w1;
else
    x_XCG_No_fus = X_NP - SM_des*cmac_w1;
    CM_alpha_ac = CM_alpha_ac_No_fus;
    x_XCG = x_XCG_No_fus;
    X_NP = X_NP;
    SM = (X_NP - x_XCG)/cmac_w1;
end

TRIM_RESULTS.CM_alpha_ac = CM_alpha_ac;
TRIM_RESULTS.x_XCG = x_XCG;
TRIM_RESULTS.X_NP = X_NP;
TRIM_RESULTS.SM = SM;

% % x_XCG with no fuselage contribution
% % Revised values
% study_var_xcg = conditions.study_var_xcg;
% if study_var_xcg == 1
%     CM_alpha_ac = -CL_alpha_ac*((X_NP-x_XCG_No_fus)/cmac_w1) + CMTalpha + C_Malpha_fus;
%     x_XCG = x_XCG_No_fus;
% else
%     CM_alpha_ac = -CL_alpha_ac*((X_NP-x_XCG_fus)/cmac_w1) + CMTalpha + C_Malpha_fus;
%     x_XCG = x_XCG_fus;
% end

Stab_Der.CL_alpha_ac = CL_alpha_ac;
Stab_Der.CM_alpha_ac = CM_alpha_ac;

% SM_actual = -CM_alpha_ac/CL_alpha_ac
% SM_real = (X_NP - x_XCG_fus)/cmac_w1
% pause

% Store DATA
% Trim_ITER.CM_alpha_ac_CLalpha = CM_alpha_ac_CLalpha;
% Trim_ITER.CM_alpha_ac_fus = CM_alpha_ac_fus;
Trim_ITER.CMTalpha = CMTalpha;
Trim_ITER.CMalpha_f = C_Malpha_fus;
% Trim_ITER.SM_Excel_fus = SM_Excel_fus;
% Trim_ITER.CM_alpha_ac = CM_alpha_ac;
% Trim_ITER.SM_actual = SM_actual;
% Trim_ITER.SM_real = SM_real;


end