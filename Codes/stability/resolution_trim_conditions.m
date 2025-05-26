function [TRIM_RESULTS,Trim_ITER,Stab_Der] = resolution_trim_conditions(AC_CONFIGURATION,Stab_Der,Stab_Der_parts,Design_criteria,Effects,...
    Geo_tier,conv_UNITS,conditions,Performance,Aero_TH,Aero,TRIM_RESULTS,Propulsion,Trim_ITER)

AC_type = AC_CONFIGURATION.AC_type;
HTP     = AC_CONFIGURATION.HTP;
Vee     = AC_CONFIGURATION.Vee;
Can     = AC_CONFIGURATION.Can;

CL_delta_e = Stab_Der.CL_delta_e;
CD_delta_e = Stab_Der.CD_delta_e;
CM_delta_e = Stab_Der.CM_delta_e; 

CL_delta_can = Stab_Der.CL_delta_can;
CD_delta_can = Stab_Der.CD_delta_can;
CM_delta_can = Stab_Der.CM_delta_can; 

CL_delta_rv = Stab_Der.CL_delta_rv;
CD_delta_rv = Stab_Der.CD_delta_rv;
CM_delta_rv = Stab_Der.CM_delta_rv; 

if Can == 1
    CL_alpha_can_e = Stab_Der_parts.CLalpha_can_e_pw;
    i_can = Design_criteria.i_w3;
    eps_can = Effects.eps_can;
    upwash = Effects.upwash;
    CL0_can_e_corrected = Stab_Der_parts.CL0_can_e_corrected;
end
if Vee == 1
    delta_rudvtr_min = Geo_tier.delta_rudvtr_min;
    CL_alpha_Vee_e = Stab_Der_parts.CLalpha_Vee_e_pw;
    CL0_Vee_e_corrected = Stab_Der_parts.CL0_Vee_e_corrected;
    i_w2 = Design_criteria.i_w2;
    eps_w2 = Effects.eps_w2;
    downwash = Effects.downwash;
end

if HTP == 1
    CL0_HTP_e_corrected = Stab_Der_parts.CL0_HTP_e_corrected;
    CL_alpha_HTP_e = Stab_Der_parts.CLalpha_HTP_e_pw;
    i_w2 = Design_criteria.i_w2;
    eps_w2 = Effects.eps_w2;
    downwash = Effects.downwash;
    delta_ele_min = Geo_tier.delta_ele_min;
end

if AC_type == 1 % AC_type = 1 - flying wing
    delta_ele_min = Geo_tier.delta_elevon_min;
end

g = conv_UNITS.g;
m_TOW = conditions.m_TOW;
w_T0 = m_TOW*g;
V = conditions.V;
rho = Performance.rho;
q_inf = 0.5*rho*V^2;
S_ref = Geo_tier.S_ref;
CL0_ac = Stab_Der_parts.CL0_ac;
CM0_ac = Stab_Der.CM0_ac;
CL_alpha_ac = Stab_Der_parts.CL_alpha_ac;
CM_alpha_ac_des = Stab_Der.CM_alpha_ac_des;
CM_alpha_ac = Stab_Der.CM_alpha_ac;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
% C_D0 = Aero_TH.CD0;
% C_D1 = Aero_TH.CD1;
% C_D2 = Aero_TH.CD2;
C_D0 = Aero.Polar.C_D0;
C_D1 = Aero.Polar.C_D1;
C_D2 = Aero.Polar.C_D2;
CL0_w1_e_corrected = Stab_Der_parts.CL0_w1_e_corrected;
CL_alpha_w1_e = Stab_Der_parts.CLalpha_w1_e_pw;
i_w1 = Design_criteria.i_w1;
alpha_max_w1_ope = Performance.alpha_max_w1_ope;
cmac_w1 = Geo_tier.cmac_w1;
X_NP = TRIM_RESULTS.X_NP;

%% Aircraft type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
switch AC_type
    case 1 % AC_type = 1 - flying wing
        CL_delta = CL_delta_e;
        CD_delta = CD_delta_e;
        CM_delta = CM_delta_e;
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        CL_delta = CL_delta_e;
        CD_delta = CD_delta_e;
        CM_delta = CM_delta_e;
    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        CL_delta = CL_delta_can + CL_delta_e;
        CD_delta = CD_delta_can + CD_delta_e;
        CM_delta = CM_delta_can + CM_delta_e;
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        CL_delta = CL_delta_rv;
        CD_delta = CD_delta_rv;
        CM_delta = CM_delta_rv;
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        CL_delta = CL_delta_can + CL_delta_rv;
        CD_delta = CD_delta_can + CD_delta_rv;
        CM_delta = CM_delta_can + CM_delta_rv;
    case 6 % AC_type = 6 - 3 surface: cannard + wing + VTP
        CL_delta = CL_delta_can;
        CD_delta = CD_delta_can;
        CM_delta = CM_delta_can;
end

% A_mat = [CL_alpha_ac, CL_delta; CM_alpha_ac, CM_delta];
% B_mat = [w_T0/(q_inf*S_ref) - CL0_ac; - CM0_ac]
% trim = B_mat\A_mat
% trim_deg = trim*R2D
% pause
%% Resolution of Trim Conditions
% num_delta_e = w_T0/(q_i32nf*S_ref) - CL0_ac + CL_alpha_ac*CM0_ac/CM_alpha_ac;
% den_delta_e = -CL_alpha_ac*CM_delta/CM_alpha_ac + CL_delta;
% trim_delta_e = num_delta_e/den_delta_e;
% trim_delta_e_deg = trim_delta_e*R2D;
% trim_alpha = - (CM_delta_e*trim_delta_e + CM0_ac)/CM_alpha_ac;
% trim_alpha_deg = trim_alpha*R2D;

%% Calcula el trimado con el XCG inicial
n = conditions.n;
% n = 3.8;
% n = 1;
phi_n = n-1;
q = (g/V)*phi_n;

q_ref = (2*V/0.469);

% Pitch Derivatives for Maneu ver static stability
CLq = Stab_Der.CLq;
CMq = Stab_Der.CMq;

CL_needed = n*w_T0/(q_inf*S_ref) - CLq*q/q_ref;
CM0_needed = CM0_ac + CMq*q/q_ref;

trim_alpha = ((CL_needed - CL0_ac)*CM_delta + CL_delta*CM0_needed)/(CL_alpha_ac*CM_delta - CL_delta*CM_alpha_ac);
trim_delta_e = -((CL_needed - CL0_ac)*CM_alpha_ac + CL_alpha_ac*CM0_needed)/(CL_alpha_ac*CM_delta - CL_delta*CM_alpha_ac);
trim_alpha_deg = trim_alpha*R2D;
trim_delta_e_deg = trim_delta_e*R2D;

x0 = [0, 0];
trim_fun = @(x) trim_eqs(x, CL_needed, CL0_ac, CL_alpha_ac, CL_delta, CM0_needed, CM_alpha_ac, CM_delta);
options = optimoptions('fsolve','Display','off');
 x = fsolve(trim_fun, x0, options);

TRIM_RESULTS.V = conditions.V;
TRIM_RESULTS.m_TOW = conditions.m_TOW;
TRIM_RESULTS.x_XCG = conditions.x_XCG;
TRIM_RESULTS.z_XCG = conditions.z_XCG;
TRIM_RESULTS.CL0_ac = CL0_ac;
TRIM_RESULTS.CL_alpha_ac = CL_alpha_ac;
TRIM_RESULTS.CL_delta = CL_delta;
TRIM_RESULTS.CM0_ac = CM0_ac;
TRIM_RESULTS.CM_delta = CM_delta;
TRIM_RESULTS.trim_alpha = trim_alpha;
TRIM_RESULTS.trim_delta_e = trim_delta_e;
TRIM_RESULTS.trim_alpha_deg = trim_alpha_deg;
TRIM_RESULTS.trim_delta_e_deg = trim_delta_e_deg;
TRIM_RESULTS.delta_T = Propulsion.delta_T;


% Euler angle
theta1 = trim_alpha;

% Calculo de valores totales
% Total Lift coefficient
CL_needed = w_T0/(q_inf*S_ref);
CL = CL0_ac + CL_alpha_ac*trim_alpha + CL_delta*trim_delta_e + CLq*q;
CD_trim_delta = abs(CD_delta*trim_delta_e);
CD_Total = C_D0 + C_D1*CL + C_D2*CL^2 + CD_trim_delta;  %polar aeronave


%% Aircraft Type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
% AC_type = 6 - 2 surface: cannard + wing + VTP
switch AC_type
    case 1 % AC_type = 1 - flying wing
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        CL_Total = CL_w1;
        Trim_ITER.CL_w1 = CL_w1;
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
        
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        CL_HTP = CL0_HTP_e_corrected + CL_alpha_HTP_e*(i_w2 - eps_w2) + CL_alpha_HTP_e*(downwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_HTP;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_HTP = CL_HTP;
        
        % angle of attack of HTP
        alpha_HTP = ((i_w2 - eps_w2) + (downwash)*trim_alpha);
        alpha_HTP_deg = alpha_HTP*R2D;
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;

    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        CL_HTP = CL0_HTP_e_corrected + CL_alpha_HTP_e*(i_w2 - eps_w2) + CL_alpha_HTP_e*(downwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_HTP + CL_can;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_HTP = CL_HTP;
        Trim_ITER.CL_can = CL_can;
        
        % angle of attack of HTP
        alpha_HTP = ((i_w2 - eps_w2) + (downwash)*trim_alpha);
        alpha_HTP_deg = alpha_HTP*R2D;
        % angle of attack of Can
        alpha_can = ((i_can + eps_can) + (upwash)*trim_alpha);
        alpha_can_deg = alpha_HTP*R2D;
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
        
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        CL_Vee = CL0_Vee_e_corrected + CL_alpha_Vee_e*(i_w2 - eps_w2) + CL_alpha_Vee_e*(downwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_Vee;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_Vee = CL_Vee;
        % angle of attack of HTP
        alpha_Vee = ((i_w2 - eps_w2) + (downwash)*trim_alpha);
        alpha_Vee_deg = alpha_Vee*R2D;
        
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_rudvtr_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
        
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        %         CL_HTP = CL0_Vee_e_corrected + CL_alpha_wb_Vee*(i_w2 - eps_w2) + CLalpha_Vee_e_pw*(downwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Vee = CL0_Vee_e_corrected + CL_alpha_Vee_e*(i_w2 - eps_w2) + CL_alpha_Vee_e*(downwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_Vee + CL_can;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_Vee = CL_Vee;
        Trim_ITER.CL_can = CL_can;
        
        % angle of attack of HTP
        alpha_Vee = ((i_w2 - eps_w2) + (downwash)*trim_alpha);
        alpha_Vee_deg = alpha_Vee*R2D;
        % angle of attack of Can
        alpha_can = ((i_can + eps_can) + (upwash)*trim_alpha);
        alpha_can_deg = alpha_Vee*R2D;
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_rudvtr_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_can;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_can = CL_can;
        
        % angle of attack of Can
        alpha_can = ((i_can + eps_can) + (upwash)*trim_alpha);
        alpha_can_deg = alpha_Vee*R2D;
        %% Most Forward XCG
        x_XCG_fwd =  -(delta_rudvtr_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;

end

% Store DATA
Trim_ITER.CL = CL;
Trim_ITER.CL_needed = CL_needed;
Trim_ITER.CD_trim_delta = CD_trim_delta;
Trim_ITER.CL_Total = CL_Total;
Trim_ITER.CD_Total = CD_Total;


% Total CM - Checking that is zero

CM = CM0_ac + CM_alpha_ac*trim_alpha + CM_delta*trim_delta_e;
TRIM_RESULTS.CM = CM;
Stab_Der.CM = CM;

x_offset_engine = Geo_tier.x_eng_xbar - conditions.x_XCG;

Trim_ITER.x_offset_engine = x_offset_engine;




%% Most Rearward XCG
SM_min = 0.10;
x_XCG_rwd  = X_NP - SM_min*cmac_w1; 
% Storing DATA
TRIM_RESULTS.x_XCG_fwd = x_XCG_fwd;
TRIM_RESULTS.x_XCG_rwd = x_XCG_rwd;

Stab_Der.CL_Total = CL_Total;


C_D0 = Aero.Polar.C_D0;
C_D1 = Aero.Polar.C_D1;
C_D2 = Aero.Polar.C_D2;

CL_alpha_ac = TRIM_RESULTS.CL_alpha_ac;
trim_alpha = TRIM_RESULTS.trim_alpha;
CM_alpha_ac = TRIM_RESULTS.CM_alpha_ac_des;
CD = Stab_Der.CD;
CD_alpha = (C_D1 + 2*C_D2*CL_needed)*CL_alpha_ac;

CX_alfa = - CD_alpha + CL_alpha_ac*trim_alpha + CL_needed; % Roskam 3.128
CX_alfa = - CD_alpha  + CL; % Steady State Flight condition Roskam 3.128
CX_alpha = - CD_alpha + CL; % Pamadi 4.447

% CD_trim_delta = abs(CD_delta*trim_delta_e);
%     CD = C_D0 + C_D1*CL + C_D2*CL^2 + CD_trim_delta;  %polar aeronave
CZ_alpha = - CL_alpha_ac -CD_alpha*trim_alpha - CD; % Roskam 3.131
CZ_alpha = - CL_alpha_ac - CD; % Steady State Flight condition Roskam 3.131
CM_alpha = CM_alpha_ac;

Stab_Der.CX_alpha = CX_alpha;
Stab_Der.CZ_alpha = CZ_alpha;
Stab_Der.CM_alpha = CM_alpha;

end
