function [TRIM_RESULTS,Trim_ITER,Stab_Der] = resolution_trim_conditions_v3(OUTPUT_read_XLSX,AC_CONFIGURATION,Stab_Der,Stab_Der_parts,Design_criteria,Effects,...
    Geo_tier,conv_UNITS,conditions,Performance,Aero_TH,Aero,TRIM_RESULTS,Propulsion,Trim_ITER,Prop_data,Weight_tier)

AC_type = AC_CONFIGURATION.AC_type;
%% Input
W1  = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

% Available control surfaces
d_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_ail; %Definition of available control surface - aileron
d_ele = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_ele; %Definition of available control surface - elevator
d_elevon = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_elevon; %Definition of available control surface - elevon
d_flap = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_flap; %Definition of available control surface - flap
d_rudder = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudder; %Definition of available control surface - rudder
d_rudvtr = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudvtr; %Definition of available control surface - ruddervator
d_rudvtr2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudvtr2; %Definition of available control surface - ruddervator
d_can = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_can; %Definition of available control surface - canard

CL_delta_e = Stab_Der.CL_delta_e;
CD_delta_e = Stab_Der.CD_delta_e;
CM_delta_e = Stab_Der.CM_delta_e; 

CL_delta_elevon = Stab_Der.CL_delta_elevon;
CD_delta_elevon = Stab_Der.CD_delta_elevon;
CM_delta_elevon = Stab_Der.CM_delta_elevon; 

CL_delta_can = Stab_Der.CL_delta_can;
CD_delta_can = Stab_Der.CD_delta_can;
CM_delta_can = Stab_Der.CM_delta_can; 

CL_delta_rv = Stab_Der.CL_delta_rv;
CD_delta_rv = Stab_Der.CD_delta_rv;
CM_delta_rv = Stab_Der.CM_delta_rv; 

CL_delta_rv2 = Stab_Der.CL_delta_rv2;
CD_delta_rv2 = Stab_Der.CD_delta_rv2;
CM_delta_rv2 = Stab_Der.CM_delta_rv2; 

if Can == 1
    CL_alpha_can_e = Stab_Der_parts.CLalpha_can_e_pw;
    i_can = Design_criteria.i_can;
    eps_can = Effects.eps_can;
    upwash = Effects.upwash;
    CL0_can_e_corrected = Stab_Der_parts.CL0_can_e_corrected;
end

if Vee == 1
    delta_rudvtr_min = Geo_tier.delta_rudvtr_min;
    CL_alpha_vee_e = Stab_Der_parts.CLalpha_vee_e_pw;
    CL0_vee_e_corrected = Stab_Der_parts.CL0_vee_e_corrected;
    i_vee = Design_criteria.i_vee;
    eps_vee = Effects.eps_vee;
    downwash_vee = Effects.downwash_vee;
end

if Vee2 == 1
    delta_rudvtr2_min = Geo_tier.delta_rudvtr_min;
    CL_alpha_vee2_e = Stab_Der_parts.CLalpha_vee2_e_pw;
    CL0_vee2_e_corrected = Stab_Der_parts.CL0_vee2_e_corrected;
    i_vee2 = Design_criteria.i_vee2;
    eps_vee2 = Effects.eps_vee2;
    downwash_vee2 = Effects.downwash_vee2;
end

if HTP == 1
    CL0_HTP_e_corrected = Stab_Der_parts.CL0_HTP_e_corrected;
    CL_alpha_HTP_e = Stab_Der_parts.CLalpha_HTP_e_pw;
    i_HTP = Design_criteria.i_HTP;
    eps_HTP = Effects.eps_HTP;
    downwash_HTP = Effects.downwash_HTP;
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

        delta_ele_min = Geo_tier.delta_elevon_min;

    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        CL_delta = CL_delta_e;
        CD_delta = CD_delta_e;
        CM_delta = CM_delta_e;

        delta_ele_min = Geo_tier.delta_ele_min;

    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        % logic that determine s whcih controkl is used

        if d_ele == 1
            CL_delta = CL_delta_e;
            CD_delta = CD_delta_e;
            CM_delta = CM_delta_e;

            delta_ele_min = Geo_tier.delta_ele_min;

        end

        if d_can == 1
            CL_delta = CL_delta_can;
            CD_delta = CD_delta_can;
            CM_delta = CM_delta_can;

            delta_ele_min = Geo_tier.delta_can_min;
        end

        if  d_ele == 1 && d_can == 1
            CL_delta = CL_delta_e;
            CD_delta = CD_delta_e;
            CM_delta = CM_delta_e;

            delta_ele_min = Geo_tier.delta_ele_min;

            Warning = 'WARNING!!! 2 Control Signals are considered for Longitudinal Control: Canard and HTP deflections. Is this version only considered HTP-elvator. If a different one please modiffy the control selection - Code in PAUSE';
            disp(Warning)
            pause
        end

        % CL_delta = CL_delta_can + CL_delta_e;
        % CD_delta = CD_delta_can + CD_delta_e;
        % CM_delta = CM_delta_can + CM_delta_e;
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        CL_delta = CL_delta_rv;
        CD_delta = CD_delta_rv;
        CM_delta = CM_delta_rv;

        delta_ele_min = Geo_tier.delta_rudvtr_min;

    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        % CL_delta = CL_delta_can + CL_delta_rv;
        % CD_delta = CD_delta_can + CD_delta_rv;
        % CM_delta = CM_delta_can + CM_delta_rv;

        if d_rudvtr == 1
            CL_delta = CL_delta_rv;
            CD_delta = CD_delta_rv;
            CM_delta = CM_delta_rv;

            delta_ele_min = Geo_tier.delta_rudvtr_min;
        end

        if d_can == 1
            CL_delta = CL_delta_can;
            CD_delta = CD_delta_can;
            CM_delta = CM_delta_can;

            delta_ele_min = Geo_tier.delta_can_min;
        end

        if  d_rudvtr == 1 && d_can == 1
            CL_delta = CL_delta_rv;
            CD_delta = CD_delta_rv;
            CM_delta = CM_delta_rv;

            delta_ele_min = Geo_tier.delta_rudvtr_min;

            Warning = 'WARNING!!! 2 Control Signals are considered for Longitudinal Control: Canard and Vtail deflections. Is this version only considered Vtail ruddervator. If a different one please modiffy the control selection - Code in PAUSE';
            disp(Warning)
            pause
        end
   
    case 6 % AC_type = 6 - 3 surface: cannard + wing + VTP

        if d_elevon == 1
            CL_delta = CL_delta_elevon;
            CD_delta = CD_delta_elevon;
            CM_delta = CM_delta_elevon;
            delta_ele_min = Geo_tier.delta_elevon_min;
        end

        if d_can == 1
            CL_delta = CL_delta_can;
            CD_delta = CD_delta_can;
            CM_delta = CM_delta_can;
            delta_ele_min = Geo_tier.delta_can_min;
        end

        if  d_elevon == 1 && d_can == 1
            CL_delta = CL_delta_elevon;
            CD_delta = CD_delta_elevon;
            CM_delta = CM_delta_elevon;
            delta_ele_min = Geo_tier.delta_elevon_min;

            Warning = 'WARNING!!! 2 Control Signals are considered for Longitudinal Control: Canard and Wing deflections. Is this version only considered Wing-elevon. If a different one please modiffy the control selection - Code in PAUSE';
            disp(Warning)
            pause
        end

        

     case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        CL_delta = CL_delta_rv;
        CD_delta = CD_delta_rv;
        CM_delta = CM_delta_rv;

        delta_ele_min = Geo_tier.delta_rudvtr_min;

    case 8 % AC_type = 8 - 2 surface: wing + V-tail + V-tail
   
        if d_rudvtr == 1
            CL_delta = CL_delta_rv;
            CD_delta = CD_delta_rv;
            CM_delta = CM_delta_rv;

            delta_ele_min = Geo_tier.delta_rudvtr_min;
        end

        if d_rudvtr2 == 1
            CL_delta = CL_delta_rv2;
            CD_delta = CD_delta_rv2;
            CM_delta = CM_delta_rv2;

            delta_ele_min = Geo_tier.delta_rudvtr_min;
        end

        if  d_rudvtr == 1 && d_rudvtr2 == 1
            CL_delta = CL_delta_rv + CL_delta_rv2;
            CD_delta = CD_delta_rv + CD_delta_rv2;
            CM_delta = CM_delta_rv + CM_delta_rv2;

            delta_ele_min = Geo_tier.delta_rudvtr_min;
        end

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
TRIM_RESULTS.delta_T = Propulsion.delta_T;

%% Solves using FSOLVE
% x0 = [0, 0];
% trim_fun = @(x) trim_eqs(x, CL_needed, CL0_ac, CL_alpha_ac, CL_delta, CM0_needed, CM_alpha_ac, CM_delta);
% options = optimoptions('fsolve','Display','off');
% x1 = fsolve(trim_fun, x0, options);
% trim_alpha = x1(1);
% trim_delta_e = x1(2);

%% Solves 3 equations
% Selects the source of the polar model
C_D0 = Aero.Polar.C_D0;
C_D1 = Aero.Polar.C_D1;
C_D2 = Aero.Polar.C_D2;

% CD_alpha = Stab_Der.CD_alpha;

% z_d_T = Geo_tier.z_d_T;
% x_d_T = Geo_tier.x_d_T; % Positive for engine behind Xcg : x_d_T = x_eng_xbar - x_XCG; 
% y_d_T = Geo_tier.y_d_T;

% Flight Path Angle
% gamma = 0;
% Bank Angle
phi = 0;
beta_ang = 0;

% conditions.gamma = gamma; 
conditions.phi = phi;
conditions.beta_ang = beta_ang;

% Resolver el sistema de ecuaciones con fsolve
x0B = [0, 0, 0, 0];
options = optimoptions('fsolve', 'Display', 'none', 'TolFun', 1e-12, 'TolX', 1e-12);
% options = optimoptions('fmincon', 'Display', 'none', 'TolFun', 1e-8, 'TolX', 1e-8);
if conditions.V_n == 1;
    [x, fval, exitflag] = fsolve(@(x) trim_eqs_SS_PULL_UP(x,CD_delta,CL_delta,CM_delta,AC_CONFIGURATION,Aero_TH,Aero,Prop_data,Geo_tier,...
        conditions, Performance, conv_UNITS, OUTPUT_read_XLSX,Stab_Der_parts,Stab_Der,Weight_tier), x0B, options);
    test=1;
else
    [x, fval, exitflag] = fsolve(@(x) trim_eqs_SS_STRAIGHT_LINE_FLIGHT(x,CD_delta,CL_delta,CM_delta,AC_CONFIGURATION,Aero_TH,Aero,Prop_data,Geo_tier,...
        conditions, Performance, conv_UNITS, OUTPUT_read_XLSX,Stab_Der_parts,Stab_Der,Weight_tier), x0B, options);
end

% Resultados
trim_alpha = x(1); % Ángulo de ataque
trim_delta_e = x(2); % Deflexión del control
trim_alpha_deg = trim_alpha*R2D;
trim_delta_e_deg = trim_delta_e*R2D;

T = x(3);     % Empuje
delta_T = x(4); % Palanca de ajuste del empuje
TRIM_RESULTS.delta_T = delta_T;

% disp('Resultados:');
% disp(['Alpha (rad): ', num2str(alpha)]);
% disp(['Delta (rad): ', num2str(delta)]);
% disp(['Thrust (T): ', num2str(T)]);
% disp(['Throttle Position (deltaT): ', num2str(delta_T)]);

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


% Euler angle
theta1 = trim_alpha;

CD_q = Stab_Der.CDq;
CL_q = Stab_Der.CLq;
CM_q = Stab_Der.CMq;

% Calculo de valores totales
% Total Lift coefficient
CL_needed = w_T0/(q_inf*S_ref);
CL = CL0_ac + CL_alpha_ac*trim_alpha + CL_delta*trim_delta_e + CL_q*(g*cmac_w1*(n-1)*cos(phi))/(2*V^2);
CD_trim_delta = abs(CD_delta*trim_delta_e);
CD_Total = C_D0 + C_D1*CL + C_D2*CL^2 + CD_trim_delta + CD_q*(g*cmac_w1*(n-1)*cos(phi))/(2*V^2);  %polar aeronave


%% Aircraft Type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
% AC_type = 6 - 2 surface: cannard + wing + VTP

TRIM_RESULTS.alpha_HTP_deg = 0;
TRIM_RESULTS.alpha_can_deg = 0;
TRIM_RESULTS.alpha_vee_deg = 0;
TRIM_RESULTS.alpha_vee2_deg = 0;

switch AC_type
    case 1 % AC_type = 1 - flying wing
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1;
        Trim_ITER.CL_w1 = CL_w1;
        % Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
        
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        CL_HTP = CL0_HTP_e_corrected + CL_alpha_HTP_e*(i_HTP - eps_HTP) + CL_alpha_HTP_e*(downwash_HTP)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_HTP;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_HTP = CL_HTP;
        
        % angle of attack of HTP
        alpha_HTP = ((i_HTP - eps_HTP) + (downwash_HTP)*trim_alpha);
        alpha_HTP_deg = alpha_HTP*R2D;
        % Storing
        TRIM_RESULTS.alpha_HTP_deg = alpha_HTP_deg;

        % Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;

    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;

        if d_ele == 1
            CL_HTP = CL0_HTP_e_corrected + CL_alpha_HTP_e*(i_HTP - eps_HTP) + CL_alpha_HTP_e*(downwash_HTP)*trim_alpha + CL_delta*trim_delta_e;
            CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha;
        end

        if d_can == 1
            CL_HTP = CL0_HTP_e_corrected + CL_alpha_HTP_e*(i_HTP - eps_HTP) + CL_alpha_HTP_e*(downwash_HTP)*trim_alpha;
            CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha + CL_delta*trim_delta_e;
        end

        if  d_ele == 1 && d_can == 1
            CL_HTP = CL0_HTP_e_corrected + CL_alpha_HTP_e*(i_HTP - eps_HTP) + CL_alpha_HTP_e*(downwash_HTP)*trim_alpha + CL_delta*trim_delta_e;
            CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha;

            Warning = 'WARNING!!! 2 Control Signals are considered for Longitudinal Control: Canard and HTP deflections. Is this version only considered HTP-elvator. If a different one please modiffy the control selection - Code in PAUSE';
            disp(Warning)
            pause
        end

        % CL_HTP = CL0_HTP_e_corrected + CL_alpha_HTP_e*(i_HTP - eps_HTP) + CL_alpha_HTP_e*(downwash_HTP)*trim_alpha + CL_delta*trim_delta_e;
        % CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha + CL_delta*trim_delta_e;
        
        CL_Total = CL_w1 + CL_HTP + CL_can;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_HTP = CL_HTP;
        Trim_ITER.CL_can = CL_can;
        
        % angle of attack of HTP
        alpha_HTP = ((i_HTP - eps_HTP) + (downwash_HTP)*trim_alpha);
        alpha_HTP_deg = alpha_HTP*R2D;
        % angle of attack of Can
        alpha_can = ((i_can + eps_can) + (upwash)*trim_alpha);
        alpha_can_deg = alpha_can*R2D;
        
        % Storing
        TRIM_RESULTS.alpha_HTP_deg = alpha_HTP_deg;
        TRIM_RESULTS.alpha_can_deg = alpha_can_deg;
        
        % Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
        
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        CL_vee = CL0_vee_e_corrected + CL_alpha_vee_e*(i_vee - eps_vee) + CL_alpha_vee_e*(downwash_vee)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_vee;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_vee = CL_vee;
        % angle of attack of HTP
        alpha_vee = ((i_vee - eps_vee) + (downwash_vee)*trim_alpha);
        alpha_vee_deg = alpha_vee*R2D;
        
        %Storing
        TRIM_RESULTS.alpha_vee_deg = alpha_vee_deg;

        % Most Forward XCG
        x_XCG_fwd =  -(delta_rudvtr_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;
        
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        %         CL_HTP = CL0_vee_e_corrected + CL_alpha_wb_vee*(i_vee - eps_vee) + CLalpha_vee_e_pw*(downwash_vee)*trim_alpha + CL_delta*trim_delta_e;
        
        if d_rudvtr == 1
            CL_vee = CL0_vee_e_corrected + CL_alpha_vee_e*(i_vee - eps_vee) + CL_alpha_vee_e*(downwash_vee)*trim_alpha + CL_delta*trim_delta_e;
            CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha;
        end

        if d_can == 1
            CL_vee = CL0_vee_e_corrected + CL_alpha_vee_e*(i_vee - eps_vee) + CL_alpha_vee_e*(downwash_vee)*trim_alpha;
            CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha + CL_delta*trim_delta_e;
        end

        if  d_rudvtr == 1 && d_can == 1
            CL_vee = CL0_vee_e_corrected + CL_alpha_vee_e*(i_vee - eps_vee) + CL_alpha_vee_e*(downwash_vee)*trim_alpha + CL_delta*trim_delta_e;
            CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha;

            Warning = 'WARNING!!! 2 Control Signals are considered for Longitudinal Control: Canard and Vtail deflections. Is this version only considered Vtail ruddervator. If a different one please modiffy the control selection - Code in PAUSE';
            disp(Warning)
            pause
        end
        
        CL_Total = CL_w1 + CL_vee + CL_can;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_vee = CL_vee;
        Trim_ITER.CL_can = CL_can;
        
        % angle of attack of HTP
        alpha_vee = ((i_vee - eps_vee) + (downwash_vee)*trim_alpha);
        alpha_vee_deg = alpha_vee*R2D;
        % angle of attack of Can
        alpha_can = ((i_can + eps_can) + (upwash)*trim_alpha);
        alpha_can_deg = alpha_can*R2D;

        %Storing
        TRIM_RESULTS.alpha_can_deg = alpha_can_deg;
        TRIM_RESULTS.alpha_vee_deg = alpha_vee_deg;
        
        % Most Forward XCG
        x_XCG_fwd =  -(delta_rudvtr_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;

    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        % CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        % CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha + CL_delta*trim_delta_e;
        % 
        if d_elevon == 1
            CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha + CL_delta*trim_delta_e;
            CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha;
        end

        if d_can == 1
            CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
            CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha + CL_delta*trim_delta_e;
        end

        if  d_elevon == 1 && d_can == 1
            CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha + CL_delta*trim_delta_e;
            CL_can = CL0_can_e_corrected + CL_alpha_can_e*(i_can + eps_can) + CL_alpha_can_e*(upwash)*trim_alpha;

            Warning = 'WARNING!!! 2 Control Signals are considered for Longitudinal Control: Canard and Wing deflections. Is this version only considered Wing-elevon. If a different one please modiffy the control selection - Code in PAUSE';
            disp(Warning)
            pause
        end

        CL_Total = CL_w1 + CL_can; 
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_can = CL_can;
        
        % angle of attack of Can
        alpha_can = ((i_can + eps_can) + (upwash)*trim_alpha);
        alpha_can_deg = alpha_can*R2D;
        
        %Storing
        TRIM_RESULTS.alpha_can_deg = alpha_can_deg;
      
        % Most Forward XCG
        x_XCG_fwd =  -(delta_ele_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;

    case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;
        CL_vee = CL0_vee_e_corrected + CL_alpha_vee_e*(i_vee - eps_vee) + CL_alpha_vee_e*(downwash_vee)*trim_alpha + CL_delta*trim_delta_e;
        CL_Total = CL_w1 + CL_vee;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_vee = CL_vee;
        % angle of attack of HTP
        alpha_vee = ((i_vee - eps_vee) + (downwash_vee)*trim_alpha);
        alpha_vee_deg = alpha_vee*R2D;
        
        %Storing
        TRIM_RESULTS.alpha_vee_deg = alpha_vee_deg;
        
        % Most Forward XCG
        x_XCG_fwd =  -(delta_rudvtr_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;

    case 8 % AC_type = 8 - 2 surface: wing + V-tail + V-tail
        CL_w1 = CL0_w1_e_corrected + CL_alpha_w1_e*i_w1 + CL_alpha_w1_e*trim_alpha;

        CL_vee = CL0_vee_e_corrected + CL_alpha_vee_e*(i_vee - eps_vee) + CL_alpha_vee_e*(downwash_vee)*trim_alpha + CL_delta*trim_delta_e;
        CL_vee2 = CL0_vee2_e_corrected + CL_alpha_vee2_e*(i_vee2 - eps_vee2) + CL_alpha_vee2_e*(downwash_vee2)*trim_alpha + CL_delta*trim_delta_e;
        
        CL_Total = CL_w1 + CL_vee;
        Trim_ITER.CL_w1 = CL_w1;
        Trim_ITER.CL_vee = CL_vee;
        Trim_ITER.CL_vee2 = CL_vee2;
        % angle of attack of HTP
        alpha_vee = ((i_vee - eps_vee) + (downwash_vee)*trim_alpha);
        alpha_vee2 = ((i_vee2 - eps_vee2) + (downwash_vee2)*trim_alpha);
        alpha_vee_deg = alpha_vee*R2D;
        alpha_vee2_deg = alpha_vee2*R2D;
        
        %Storing
        TRIM_RESULTS.alpha_vee_deg = alpha_vee_deg;
        TRIM_RESULTS.alpha_vee2_deg = alpha_vee2_deg;

        % Most Forward XCG
        x_XCG_fwd =  -(delta_rudvtr_min + (CM0_ac/CM_delta))/(alpha_max_w1_ope*D2R*CL_alpha_ac/(CM_delta*cmac_w1)) + X_NP;

end

% Store DATA
Trim_ITER.CL = CL;
Trim_ITER.CL_needed = CL_needed;
Trim_ITER.CD_trim_delta = CD_trim_delta;
Trim_ITER.CL_Total = CL_Total;
Trim_ITER.CD_Total = CD_Total;

% Total CM - Checking that is zero

CM = CM0_ac + CM_alpha_ac*trim_alpha + CM_delta*trim_delta_e + CM_q*(g*cmac_w1*(n-1)*cos(phi))/(2*V^2);
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
CZ_alpha = - CL_alpha_ac - CD_alpha*trim_alpha - CD; % Roskam 3.131
CZ_alpha = - CL_alpha_ac - CD; % Steady State Flight condition Roskam 3.131
CM_alpha = CM_alpha_ac;

Stab_Der.CX_alpha = CX_alpha;
Stab_Der.CZ_alpha = CZ_alpha;
Stab_Der.CM_alpha = CM_alpha;

end
