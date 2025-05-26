function  [Stab_Der,Trim_ITER] = get_propulsive_derivatives(AC_CONFIGURATION, Propulsion, Aero_TH,Aero,Prop_data,Geo_tier,Stab_Der_parts,afe,conv_UNITS,conditions,Performance,Trim_ITER,Stab_Der)

propul = AC_CONFIGURATION.propulsion;
type_engine = propul(1);
bypass_ratio = propul(5); % By-pass Ratio
type_prop = AC_CONFIGURATION.type_prop; 

C_D0 = Aero.Polar.C_D0;
C_D1 = Aero.Polar.C_D1;
C_D2 = Aero.Polar.C_D2;

S_w1 = Geo_tier.S_w1;
S_ref = Geo_tier.S_ref;
cmac_w1 = Geo_tier.cmac_w1;
n_eng = Prop_data.n_eng;
T_eng = Propulsion.Ti_eng;
J = Propulsion.J;
CL_alpha_ac = Stab_Der_parts.CL_alpha_ac;
depsu_dalpha = afe.depsu_dalpha;
D_prop = Prop_data.D_prop;
R_prop = D_prop/2;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;
m22ft2 = conv_UNITS.m22ft2;
m2ft = conv_UNITS.m2ft;
W2hp = conv_UNITS.W2hp;
g = conv_UNITS.g;
rho_SI2rho_IMP = conv_UNITS.rho_SI2rho_IMP;

m_TOW = conditions.m_TOW;

z_d_T = Geo_tier.z_d_T;
x_d_T = Geo_tier.x_d_T; % Positive for engine behind Xcg : x_d_T = x_eng_xbar - x_XCG; 
y_d_T = Geo_tier.y_d_T;
phi_T = 0;
d_T_d_V = Propulsion.d_T_d_V;
d_etap_d_J = Propulsion.d_etap_d_J;
d_CP_d_V = Propulsion.d_CP_d_V;
etha_emp = Propulsion.etha_emp;
d_T = x_d_T*sin(phi_T) + z_d_T*cos(phi_T);
V = conditions.V;
rho = Performance.rho;
q_inf = 0.5*rho*V^2; 
w_T0 = m_TOW*g;
CL = w_T0/(q_inf*S_ref); % assume equilibry steady state flight
%% Propulsive Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%DERIVADA PROPULSIVA LONGITUDINAL%%%%%%%%%%%%%%%%
CD_Total_prel = C_D0 + C_D1*CL + C_D2*CL^2;  %polar aeronave
Drag = q_inf*S_ref*CD_Total_prel;

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CTx1 %%%%%%%%%%%%%%%%%%%%%%%%%
phiT = 0; %	is the thrust line inclination angle.
T_set = Propulsion.Ti;
CTx1 = CD_Total_prel;
% CTx1 = T_set*cos(phiT+trim_alpha)/q_inf*S_ref;
Stab_Der.T = T_set;

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CMT1 %%%%%%%%%%%%%%%%%%%%%%%%%
%% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
switch type_engine
    case 1 % TIPO DE MOTOR --> 1_TURBOFAN 
        % Increment of Pitching Moment Coefficient due to the Lift Component of the Propeller Normal Force Aproximada a 0
        Delta_CM_Nprop = 0;
        Delta_CM_Tprop = -(d_T/cmac_w1)*CTx1;
        CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;
    case 2 % TIPO DE MOTOR --> 2_TURBOPROP
        % Increment of Pitching Moment Coefficient due to the Lift Component of the Propeller Normal Force Aproximada a 0
        Delta_CM_Nprop = 0;
        Delta_CM_Tprop = -(d_T/cmac_w1)*CTx1;
        CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;
    case 3 % TIPO DE MOTOR --> 3_PISTON 
        Delta_CM_Nprop = 0;
        CMT1 = - (d_T/cmac_w1)*CTx1;
        CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;
    case 4 % TIPO DE MOTOR --> 4_ELECTRIC_PROP
        % Increment of Pitching Moment Coefficient due to the Lift Component of the Propeller Normal Force Aproximada a 0
        Delta_CM_Nprop = 0;
        Delta_CM_Tprop = -(d_T/cmac_w1)*CTx1;
        CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;
    case 5 % TIPO DE MOTOR --> 5_PISTO_CUSTOM
        % Increment of Pitching Moment Coefficient due to the Lift Component of the Propeller Normal Force Aproximada a 0
        Delta_CM_Nprop = 0;
        Delta_CM_Tprop = -(d_T/cmac_w1)*CTx1;
        CMT1 = Delta_CM_Nprop + Delta_CM_Tprop;
end

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CMTalpha %%%%%%%%%%%%%%%%%%%%%%%%%
%% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
switch type_engine
    case 1 % TIPO DE MOTOR --> 1_TURBOFAN
        n_j = n_eng; % number of jet engines
        switch bypass_ratio
            case 1 % BPR from 0 to 1
                k_gas = 0.0003;
            case 2 % BPR from 1 to 2
                k_gas = 0.0007;
            case 3 % BPR from 2 to 4
                k_gas = 0.0009;
            case 4 % BPR from 4 to 6
                k_gas = 0.0011;
        end
        slugs_s2kg_s = 14.593903; % Convertsfrom slugs/s to kg/s
        m_dot_gass = T_eng*k_gas*slugs_s2kg_s;
        m_dot_cool = 0.006*m_dot_gass;
        m_dot = m_dot_gass + m_dot_cool;
        % is the inlet cross sectional area for Max Thrust
        % effect of thrustline offset on longitudinal stability
        dCM_dCL_TL = 0 ; % Negligible for Jet engines
        % effect of propeller or inlet normal force on longitudinal stability
        dCM_dCL_N = 0.035*m_dot*(-x_d_T)*depsu_dalpha/(S_w1*cmac_w1*rho*V*CL_alpha_ac);
        Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;        
        CMTalpha = Delta_CM_CL_T*CL_alpha_ac;
%         CMTalpha = 0;
        % Data from NACA REport 640 - properties of airfoil
        w_R_30 = 0.0525*2;
        w_R_60 = 0.073*2;
        w_R_90 = 0.045*2;        

    case 2 % TIPO DE MOTOR --> 2_TURBOPROP
        Nprop = 1; % number of props
        % Data from NACA REport 640 - properties of airfoil
        w_R_30 = 0.0525*2;
        w_R_60 = 0.073*2;
        w_R_90 = 0.045*2;        
        % Relation between propeller geometry - pitch and beta
        pitch =12*D2R; % Pitch of propeller 22x12W
        beta_3_4R = atan(pitch/(2*pi*(2/3)*R_prop));
        % DATCOM - PROPELLER INFLOW FACTOR - FIGURE 4.6.1-25B
        % Nominal blade angle at 0.75R Radius in deg
        Beta_blade =[15.,20.,25.,30.,35.,40.,60.];
        % Data for CNa
        beta_blade = beta_3_4R*R2D;
        % N_blades = [2.,3.,4.,6.];
        n_blades=2;
        if n_blades==2
            CNa_n = [.08,.10,.115,.126,.140,.150,.192];
        elseif n_blades==3
            CNa_n = [.11,.139,.160,.182,.20,.216,.275];
        elseif n_blades==5
            CNa_n = [.136,.172,.20,.226,.25,.272,.35];
        elseif n_blades==6
            CNa_n = [.196,.237,.275,.315,.35,.382,.5];
        end
        C_Na_807  = interp1(Beta_blade,CNa_n,beta_blade,'spline');
        % Airplane Design; Roskam, Jan; Part VI, pag 342 (2032 PDF)
        % depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop, cr_we, XLE_w, l_fus, downw, CLa_WB*S_w/Sref); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
        % n_blades = 2;
        % beta_blade = 20; % beta_blade in degrees!!
        % CNa_807         = CNa_807_calc(n_blades, beta_blade)
        KN              = 262*(2*w_R_30/D_prop) + 262*(2*w_R_60/D_prop) + 135*(2*w_R_90/D_prop);
        dCN_dalpha  = C_Na_807*(1 + 0.8*(KN/80.7 - 1));
        
        % Distance from the end of wing to the end of fuselage
        % x_1 = x_m_propeller - (x_w_LE + cmac_w1);
        % l_h = l_ht - (3*cmac_w1/4);
        % depsu_dalpha = (1-deps_dalpha)*(x1/l_h)-1;
        
        % % Location of the prop disk (source of thrust)
        % x_prop_cF = Geo_tier.x_prop_cF;
        % AR_w1 = Geo_tier.AR_w1;
        % depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);
        
        % rho_met2imp = conv_UNITS.rho_met2imp;
        % W2hp = conv_UNITS.W2hp;
        
        T_set_eng = Propulsion.Ti_eng;
        CT_p_eng = T_set_eng/(q_inf*S_ref);
        % P_SHP = P_SET*W2hp;
        Pi_eng = Propulsion.Pi_eng; % power per engine
        P_SHP = Pi_eng*W2hp;
        % rho_met2imp = 0.0017645;
        W_current_lb = m_TOW*2.20462;
        % The intermediate calculation parameter is given by:
        dTc_dCL = (3/2)*550*P_SHP*sqrt(rho*rho_SI2rho_IMP)*sqrt(CL)/...
            (sqrt((2*W_current_lb/(S_w1*m22ft2))^3)*(D_prop*m2ft)^2);
        % effect of thrustline offset on longitudinal stability
        dCM_dCL_TL = (Nprop*(2*((D_prop)^2)*d_T)/(S_ref*cmac_w1))*dTc_dCL;
        % DATCOM -
        % ----PROPELLER INFLOW FACTOR
        %      ----FIGURE 4.6.1-25B
        f_inflow_input_vec = [0.,1.,2.,3.,4.,6.,8.,10.,14.,19.,22.];
        f_inflow_input = S_ref*(CT_p_eng)/(8*(R_prop^2));
        f_inflow_vec = [1.0,1.55,1.94,2.20,2.40,2.75,3.05,3.30,3.75,4.25,4.54];
        finflow  = interp1(f_inflow_input_vec,f_inflow_vec,f_inflow_input,'spline');
        % effect of propeller or inlet normal force on longitudinal stability is given by:
        l_prop = (-1)*x_d_T*cos(phi_T) - (-1)*y_d_T*sin(phi_T);
        dCM_dCL_N = (pi/4)*finflow*Nprop*l_prop*((D_prop)^2)...
            *dCN_dalpha*(1 - depsu_dalpha)/(S_w1*cmac_w1*CL_alpha_ac);
        Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;
        CMTalpha = Delta_CM_CL_T*CL_alpha_ac;
    case 3 % TIPO DE MOTOR --> 3_PISTON
        Nprop = 1; % number of props
        % Data from NACA REport 640 - properties of airfoil
        w_R_30 = 0.0525*2;
        w_R_60 = 0.073*2;
        w_R_90 = 0.045*2;       
        % Relation between propeller geometry - pitch and beta
        pitch =12*D2R; % Pitch of propeller 22x12W
        beta_3_4R = atan(pitch/(2*pi*(2/3)*R_prop));
        % DATCOM - PROPELLER INFLOW FACTOR - FIGURE 4.6.1-25B
        % Nominal blade angle at 0.75R Radius in deg
        Beta_blade =[15.,20.,25.,30.,35.,40.,60.];
        % Data for CNa
        beta_blade = beta_3_4R*R2D;
        % N_blades = [2.,3.,4.,6.];
        n_blades=2;
        if n_blades==2
            CNa_n = [.08,.10,.115,.126,.140,.150,.192];
        elseif n_blades==3
            CNa_n = [.11,.139,.160,.182,.20,.216,.275];
        elseif n_blades==5
            CNa_n = [.136,.172,.20,.226,.25,.272,.35];
        elseif n_blades==6
            CNa_n = [.196,.237,.275,.315,.35,.382,.5];
        end
        C_Na_807  = interp1(Beta_blade,CNa_n,beta_blade,'spline');
        % Airplane Design; Roskam, Jan; Part VI, pag 342 (2032 PDF)
        % depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop, cr_we, XLE_w, l_fus, downw, CLa_WB*S_w/Sref); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
        % n_blades = 2;
        % beta_blade = 20; % beta_blade in degrees!!
        % CNa_807         = CNa_807_calc(n_blades, beta_blade)
        KN              = 262*(2*w_R_30/D_prop) + 262*(2*w_R_60/D_prop) + 135*(2*w_R_90/D_prop);
        dCN_dalpha  = C_Na_807*(1 + 0.8*(KN/80.7 - 1));
        
        % Distance from the end of wing to the end of fuselage
        % x_1 = x_m_propeller - (x_w_LE + cmac_w1);
        % l_h = l_ht - (3*cmac_w1/4);
        % depsu_dalpha = (1-deps_dalpha)*(x1/l_h)-1;
        
        % % Location of the prop disk (source of thrust)
        % x_prop_cF = Geo_tier.x_prop_cF;
        % AR_w1 = Geo_tier.AR_w1;
        % depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);
        
        % rho_met2imp = conv_UNITS.rho_met2imp;
        % W2hp = conv_UNITS.W2hp;
        
        T_set_eng = Propulsion.Ti_eng;
        CT_p_eng = T_set_eng/(q_inf*S_ref);
        % P_SHP = P_SET*W2hp;
        Pi_eng = Propulsion.Pi_eng; % power per engine
        P_SHP = Pi_eng*W2hp;
        % rho_met2imp = 0.0017645;
        W_current_lb = m_TOW*2.20462;
        % The intermediate calculation parameter is given by:
        dTc_dCL = (3/2)*550*P_SHP*sqrt(rho*rho_SI2rho_IMP)*sqrt(CL)/...
            (sqrt((2*W_current_lb/(S_w1*m22ft2))^3)*(D_prop*m2ft)^2);
        % effect of thrustline offset on longitudinal stability
        dCM_dCL_TL = (Nprop*(2*((D_prop)^2)*d_T)/(S_ref*cmac_w1))*dTc_dCL;
        % DATCOM -
        % ----PROPELLER INFLOW FACTOR
        %      ----FIGURE 4.6.1-25B
        f_inflow_input_vec = [0.,1.,2.,3.,4.,6.,8.,10.,14.,19.,22.];
        f_inflow_input = S_ref*(CT_p_eng)/(8*(R_prop^2));
        f_inflow_vec = [1.0,1.55,1.94,2.20,2.40,2.75,3.05,3.30,3.75,4.25,4.54];
        finflow  = interp1(f_inflow_input_vec,f_inflow_vec,f_inflow_input,'spline');
        % effect of propeller or inlet normal force on longitudinal stability is given by:
        l_prop = (-1)*x_d_T*cos(phi_T) - (-1)*y_d_T*sin(phi_T);
        dCM_dCL_N = (pi/4)*finflow*Nprop*l_prop*((D_prop)^2)...
            *dCN_dalpha*(1 - depsu_dalpha)/(S_w1*cmac_w1*CL_alpha_ac);
        Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;
        CMTalpha = Delta_CM_CL_T*CL_alpha_ac;
    case 4 % TIPO DE MOTOR --> 4_ELECTRIC_PROP
        Nprop = 1; % number of props
        % Data from NACA REport 640 - properties of airfoil
        w_R_30 = 0.0525*2;
        w_R_60 = 0.073*2;
        w_R_90 = 0.045*2;        
        % Relation between propeller geometry - pitch and beta
        pitch =12*D2R; % Pitch of propeller 22x12W
        beta_3_4R = atan(pitch/(2*pi*(2/3)*R_prop));
        % DATCOM - PROPELLER INFLOW FACTOR - FIGURE 4.6.1-25B
        % Nominal blade angle at 0.75R Radius in deg
        Beta_blade =[15.,20.,25.,30.,35.,40.,60.];
        % Data for CNa
        beta_blade = beta_3_4R*R2D;
        % N_blades = [2.,3.,4.,6.];
        n_blades=2;
        if n_blades==2
            CNa_n = [.08,.10,.115,.126,.140,.150,.192];
        elseif n_blades==3
            CNa_n = [.11,.139,.160,.182,.20,.216,.275];
        elseif n_blades==5
            CNa_n = [.136,.172,.20,.226,.25,.272,.35];
        elseif n_blades==6
            CNa_n = [.196,.237,.275,.315,.35,.382,.5];
        end
        C_Na_807  = interp1(Beta_blade,CNa_n,beta_blade,'spline');
        
        % Airplane Design; Roskam, Jan; Part VI, pag 342 (2032 PDF)
        % depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop, cr_we, XLE_w, l_fus, downw, CLa_WB*S_w/Sref); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
        % n_blades = 2;
        % beta_blade = 20; % beta_blade in degrees!!
        % CNa_807         = CNa_807_calc(n_blades, beta_blade)
        KN              = 262*(2*w_R_30/D_prop) + 262*(2*w_R_60/D_prop) + 135*(2*w_R_90/D_prop);
        dCN_dalpha  = C_Na_807*(1 + 0.8*(KN/80.7 - 1));
        
        % Distance from the end of wing to the end of fuselage
        % x_1 = x_m_propeller - (x_w_LE + cmac_w1);
        % l_h = l_ht - (3*cmac_w1/4);
        % depsu_dalpha = (1-deps_dalpha)*(x1/l_h)-1;
        
        % % Location of the prop disk (source of thrust)
        % x_prop_cF = Geo_tier.x_prop_cF;
        % AR_w1 = Geo_tier.AR_w1;
        % depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);
        
        % rho_met2imp = conv_UNITS.rho_met2imp;
        % W2hp = conv_UNITS.W2hp;
        
        T_set_eng = Propulsion.Ti_eng;
        CT_p_eng = T_set_eng/(q_inf*S_ref);
        % P_SHP = P_SET*W2hp;
        Pi_eng = Propulsion.Pi_eng; % power per engine
        P_SHP = Pi_eng*W2hp;
        % rho_met2imp = 0.0017645;
        W_current_lb = m_TOW*2.20462;
        % The intermediate calculation parameter is given by:
        dTc_dCL = (3/2)*550*P_SHP*sqrt(rho*rho_SI2rho_IMP)*sqrt(CL)/...
            (sqrt((2*W_current_lb/(S_w1*m22ft2))^3)*(D_prop*m2ft)^2);
        % effect of thrustline offset on longitudinal stability
        dCM_dCL_TL = (Nprop*(2*((D_prop)^2)*d_T)/(S_ref*cmac_w1))*dTc_dCL;
        % DATCOM -
        % ----PROPELLER INFLOW FACTOR
        %      ----FIGURE 4.6.1-25B
        f_inflow_input_vec = [0.,1.,2.,3.,4.,6.,8.,10.,14.,19.,22.];
        f_inflow_input = S_ref*(CT_p_eng)/(8*(R_prop^2));
        f_inflow_vec = [1.0,1.55,1.94,2.20,2.40,2.75,3.05,3.30,3.75,4.25,4.54];
        finflow  = interp1(f_inflow_input_vec,f_inflow_vec,f_inflow_input,'spline');
        % effect of propeller or inlet normal force on longitudinal stability is given by:
        l_prop = (-1)*x_d_T*cos(phi_T) - (-1)*y_d_T*sin(phi_T);
        dCM_dCL_N = (pi/4)*finflow*Nprop*l_prop*((D_prop)^2)...
            *dCN_dalpha*(1 - depsu_dalpha)/(S_w1*cmac_w1*CL_alpha_ac);
        Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;
        CMTalpha = Delta_CM_CL_T*CL_alpha_ac;
        
    case 5 % TIPO DE MOTOR --> 3_PISTON_CUSTOM
        Nprop = 1; % number of props
        % Data from NACA REport 640 - properties of airfoil
        w_R_30 = 0.0525*2;
        w_R_60 = 0.073*2;
        w_R_90 = 0.045*2;       
        % Relation between propeller geometry - pitch and beta
        pitch =12*D2R; % Pitch of propeller 22x12W
        beta_3_4R = atan(pitch/(2*pi*(2/3)*R_prop));
        % DATCOM - PROPELLER INFLOW FACTOR - FIGURE 4.6.1-25B
        % Nominal blade angle at 0.75R Radius in deg
        Beta_blade =[15.,20.,25.,30.,35.,40.,60.];
        % Data for CNa
        beta_blade = beta_3_4R*R2D;
        % N_blades = [2.,3.,4.,6.];
        n_blades=2;
        if n_blades==2
            CNa_n = [.08,.10,.115,.126,.140,.150,.192];
        elseif n_blades==3
            CNa_n = [.11,.139,.160,.182,.20,.216,.275];
        elseif n_blades==5
            CNa_n = [.136,.172,.20,.226,.25,.272,.35];
        elseif n_blades==6
            CNa_n = [.196,.237,.275,.315,.35,.382,.5];
        end
        C_Na_807  = interp1(Beta_blade,CNa_n,beta_blade,'spline');
        % Airplane Design; Roskam, Jan; Part VI, pag 342 (2032 PDF)
        % depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop, cr_we, XLE_w, l_fus, downw, CLa_WB*S_w/Sref); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
        % n_blades = 2;
        % beta_blade = 20; % beta_blade in degrees!!
        % CNa_807         = CNa_807_calc(n_blades, beta_blade)
        KN              = 262*(2*w_R_30/D_prop) + 262*(2*w_R_60/D_prop) + 135*(2*w_R_90/D_prop);
        dCN_dalpha  = C_Na_807*(1 + 0.8*(KN/80.7 - 1));
        
        % Distance from the end of wing to the end of fuselage
        % x_1 = x_m_propeller - (x_w_LE + cmac_w1);
        % l_h = l_ht - (3*cmac_w1/4);
        % depsu_dalpha = (1-deps_dalpha)*(x1/l_h)-1;
        
        % % Location of the prop disk (source of thrust)
        % x_prop_cF = Geo_tier.x_prop_cF;
        % AR_w1 = Geo_tier.AR_w1;
        % depsu_dalpha = -upwash_calc(AR_w1, 0.25+x_xbar_w1, x_prop_cF, cmac_w1);
        
        % rho_met2imp = conv_UNITS.rho_met2imp;
        % W2hp = conv_UNITS.W2hp;
        
        T_set_eng = Propulsion.Ti_eng;
        CT_p_eng = T_set_eng/(q_inf*S_ref);
        % P_SHP = P_SET*W2hp;
        Pi_eng = Propulsion.Pi_eng; % power per engine
        P_SHP = Pi_eng*W2hp;
        % rho_met2imp = 0.0017645;
        W_current_lb = m_TOW*2.20462;
        % The intermediate calculation parameter is given by:
        dTc_dCL = (3/2)*550*P_SHP*sqrt(rho*rho_SI2rho_IMP)*sqrt(CL)/...
            (sqrt((2*W_current_lb/(S_w1*m22ft2))^3)*(D_prop*m2ft)^2);
        % effect of thrustline offset on longitudinal stability
        dCM_dCL_TL = (Nprop*(2*((D_prop)^2)*d_T)/(S_ref*cmac_w1))*dTc_dCL;
        % DATCOM -
        % ----PROPELLER INFLOW FACTOR
        %      ----FIGURE 4.6.1-25B
        f_inflow_input_vec = [0.,1.,2.,3.,4.,6.,8.,10.,14.,19.,22.];
        f_inflow_input = S_ref*(CT_p_eng)/(8*(R_prop^2));
        f_inflow_vec = [1.0,1.55,1.94,2.20,2.40,2.75,3.05,3.30,3.75,4.25,4.54];
        finflow  = interp1(f_inflow_input_vec,f_inflow_vec,f_inflow_input,'spline');
        % effect of propeller or inlet normal force on longitudinal stability is given by:
        l_prop = (-1)*x_d_T*cos(phi_T) - (-1)*y_d_T*sin(phi_T);
        dCM_dCL_N = (pi/4)*finflow*Nprop*l_prop*((D_prop)^2)...
            *dCN_dalpha*(1 - depsu_dalpha)/(S_w1*cmac_w1*CL_alpha_ac);
        Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;
        CMTalpha = Delta_CM_CL_T*CL_alpha_ac;
end
%%
%% REVIEW
% Cancel out as not sure if its working
CMTalpha = 0;
% Stab_Der.CMTalpha = CMTalpha;

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas%%%%%%%%%%%%%%%%%%%%%%%%%
% Coeffcicients of second order Power model
% P = A_power*V^2 + B_power*V + C_power

%conversiones de unidades imperiales
% qmet2qimp = 10.76391*(3.28084^2)*(1.2250/0.0023772);
% W2pftsec = 0.7375621;
% m22ft2 = 10.76391;
% rho_SI2rho_IMP = (0.0023772/1.2250);
% qmet2qimp = (3.28084^2)*(1.2250/0.0023772);
% qmet2qimp = m22ft2*rho_SI2rho_IMP;
% N2lbf = 0.2248089;
% m2ft = 3.28083;

% For future versions
alpha = 0;
beta = 0;
% RPS = Propulsion.RPS;
% adimensional_P_Heli =rho*(RPS^3)*(D_prop^5);
% adimensional_P_Aero =q_inf*S_ref;
% adimensional_conv = adimensional_P_Aero/adimensional_P_Heli;
% dcP_du = d_CP_d_V*adimensional_conv;
% ctxu = (1/(q_inf*qmet2qimp*S_w1*m22ft2))*(2*A_power*V + B_power)*W2pftsec;

%% Type engine
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CTxu %%%%%%%%%%%%%%%%%%%%%%%%%
switch type_engine
    case 1 % TIPO DE MOTOR --> 1_TURBOFAN 
                CTxu = (V/(q_inf*S_ref))*d_T_d_V - 2*CTx1;
    case 2 % TIPO DE MOTOR --> 2_TURBOPROP 
        switch type_prop
            case 1
                %         CTxu = (1/(q_inf*S_ref))*d_P_d_V - 2*CTx1  %por helice a paso fijo
                CTxu = d_CP_d_V - 2*CTx1;  %por helice a paso fijo
                CTxu = - 3*CTx1 + (CTx1*J*d_etap_d_J/etha_emp);   %por helice a paso fijo
            case 2
                CTxu = - 3*CTx1;  %por helice a paso variable
        end
    case 3 % TIPO DE MOTOR --> 3_PISTON
        switch type_prop
            case 1
                %         CTxu = (1/(q_inf*S_ref))*d_P_d_V - 2*CTx1  %por helice a paso fijo
                CTxu = d_CP_d_V - 2*CTx1;  %por helice a paso fijo
                CTxu = - 3*CTx1 + (CTx1*J*d_etap_d_J/etha_emp);   %por helice a paso fijo
            case 2
                CTxu = - 3*CTx1;  %por helice a paso variable
        end
    case 4% TIPO DE MOTOR --> 4_ELECTRIC_PROP
        switch type_prop
            case 1
%                 CTxu = (1/(q_inf*S_ref))*d_P_d_V - 2*CTx1;  %por helice a paso fijo
                CTxu = d_CP_d_V - 2*CTx1;  %por helice a paso fijo
                CTxu = - 3*CTx1 + (CTx1*J*d_etap_d_J/etha_emp);   %por helice a paso fijo
            case 2
                CTxu = - 3*CTx1;  %por helice a paso variable
        end
    case 5 % TIPO DE MOTOR --> 5_PISTON_CUSTOM
        switch type_prop
            case 1
                %         CTxu = (1/(q_inf*S_ref))*d_P_d_V - 2*CTx1  %por helice a paso fijo
                CTxu = d_CP_d_V - 2*CTx1;  %por helice a paso fijo
                CTxu = - 3*CTx1 + (CTx1*J*d_etap_d_J/etha_emp);   %por helice a paso fijo
            case 2
                CTxu = - 3*CTx1;  %por helice a paso variable
        end
end

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CTxalpha %%%%%%%%%%%%%%%%%%%%%%%%%
CTxalpha = 0;

%%%%%%%%%%%%%%%%%%%%%%%derivadas propulsivas - CMtu %%%%%%%%%%%%%%%%%%%%%%%%%
CMtu = - (d_T/cmac_w1)*CTxu;
Trim_ITER.w_R_30 = w_R_30;
Trim_ITER.w_R_60 = w_R_60;
Trim_ITER.w_R_90 = w_R_90;

Stab_Der.CTx1 = CTx1;
Stab_Der.CTxu = CTxu;
Stab_Der.CTxalpha = CTxalpha;
Stab_Der.CMt1 = CMT1;
Trim_ITER.CMT1 = CMT1;
Trim_ITER.CMTalpha = CMTalpha;
Stab_Der.CMtu = CMtu;
Stab_Der.CMtalpha = CMTalpha;
Stab_Der.CL = CL;
% Stab_Der.CD = CD;
end