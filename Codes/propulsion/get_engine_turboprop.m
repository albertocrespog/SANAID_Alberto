function Propulsion = get_engine_turboprop(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX) % 1: TIPO DE MOTOR -->  2_TURBOPROP 

h_SL = 0;
[Temp_SL,rho_SL,p_SL,a_SL]=atmos_inter_mio(h_SL);
[Temp,rho,p,a]=atmos_inter_mio(h);

% Mach angle at Sea Level and at altitude
M_SL = V/a_SL;
M = V/a;

S_ref = Geo_tier.S_ref;
cmac = Geo_tier.cmac_w1;

qbar = 0.5*rho*V^2;

D_prop = Prop_data.D_prop;
A_prop = Prop_data.A_prop;
R_heli = D_prop/2;
D = Prop_data.D_prop;

RPM_min = OUTPUT_read_XLSX.Propulsive_flags.RPM_min; %
RPM_max = OUTPUT_read_XLSX.Propulsive_flags.RPM_max; %

propul = AC_CONFIGURATION.propulsion;

% propul(1) = 3; % Type of engine
% propul(2) = n_eng; % Number of engines
% propul(3) = 5; % Thrust (lbf) or Power (shp) per engine
% propul(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
% propul(5) = 0; % By-pass
% propul(6) = 0.82; % Prop efficiency
% propul(7) = 1; % Normativa
% propul(8) = 0; % Capacidad calorifica del combustible
% propul(9) = 28*2.54/100; % Diámetro de la hélice

% PROPULSION
tipo_motor = propul(1);
n_eng = propul(2);
T_SL = propul(3);
c_SL = propul(4);
civil_mil = propul(5);
eta_p = propul(6);
derivacion = propul(7);

RPM_min = OUTPUT_read_XLSX.Propulsive_flags.RPM_min; %
RPM_max = OUTPUT_read_XLSX.Propulsive_flags.RPM_max; %
beta_pitch = OUTPUT_read_XLSX.Propulsive_flags.beta_pitch; %
beta_variable_pitch = OUTPUT_read_XLSX.Propulsive_flags.beta_variable_pitch; %
b_prop = OUTPUT_read_XLSX.Propulsive_flags.b_prop;

% Calculate Propuilsion model based on Propeller
T_eng = Fdes/n_eng;

if beta_variable_pitch==1;
    
    % Cargar los modelos polinómicos
    load('modelos_ct.mat', 'modelos_ct');
    load('modelos_cp.mat', 'modelos_cp');
    load('modelos_ct_red.mat', 'modelos_ct_red');
    load('modelos_cp_red.mat', 'modelos_cp_red');

    % Ángulos de hélice disponibles
    angulos_helice = [15, 20, 25, 30, 35, 40, 45];

    % Encontrar el índice correspondiente al ángulo de hélice
    indice_curva = find(angulos_helice == beta_pitch);
    if isempty(indice_curva)
        error('Ángulo de hélice no válido. Debe ser uno de los siguientes: [15, 20, 25, 30, 35, 40, 45]');
    end

    % Obtener el modelo polinómico correspondiente
    p_ct = modelos_ct{indice_curva};
    p_cp = modelos_cp{indice_curva};
    p_ct_red = modelos_ct_red{indice_curva};
    p_cp_red = modelos_ct_red{indice_curva};

    %n = get_propeller_models(V,Prop_data,T_eng,rho,p_ct,p_ct_red);
    n = Selection_J_CT_F_des_v2(V,Prop_data,T_eng,rho,p_ct,p_ct_red,OUTPUT_read_XLSX);

    RPM = n*60;
    Delta_RPM = RPM_max - RPM_min;
    delta_T = n*60/Delta_RPM;

    Omega = RPM*(2*pi/60);
    % n = Selection_J_CT_F_des(V,Prop_data,T_eng,rho)
    check_RPMMAX = n*60; % changes to RPM

    eta_gear = 1;
    eta_m = Prop_data.eta_m;
    eta_esc = 1;
    eta_dist = Prop_data.eta_dist;

    
    N_order_CT = length(p_ct)-1;
    N_order_CP = length(p_cp)-1;

    CT_Polyfit = p_ct;
    CP_Polyfit = p_cp;

    % Calculates the Advanced PArameter Ratio
    J = V/(n*D);
    % Determines CT as a funcion of approximation polynomial
    Ct_total =  CT_Polyfit(N_order_CT+1);
    for j=1:N_order_CT
        Ct_intermediate = CT_Polyfit(j)*J.^(N_order_CT+1-j);
        Ct_total = Ct_total + Ct_intermediate;
    end
    % CT = CT0 + CT1*J + CT2*J^2;
    CT = Ct_total;

    % Determines CP as a funcion of approximation polynomial
    Cp_total =  CP_Polyfit(N_order_CP+1);
    for j=1:N_order_CP
        Cp_intermediate = CP_Polyfit(j)*J.^(N_order_CP+1-j);
        Cp_total = Cp_total + Cp_intermediate;
    end
    % CP = CP0 + CP1*J + CP2*J^2 + CP3*J^3;
    CP = Cp_total;

    % % Determines CP as a funcion of approximation polynomial
    % etamp_total =  etamp_Polyfit(N_order_etamp+1);
    % for j=1:N_order_etamp
    %     etamp_intermediate = etamp_Polyfit(j)*J.^(N_order_etamp+1-j);
    %     etamp_total = etamp_total + etamp_intermediate;
    % end
    % % ethamp = etha_mp0 + etha_mp1*J + etha_mp2*J^2 + etha_mp3*J^3
    % ethamp = etamp_total;


    % Determines CQ as a funcion of approximation polynomial
    % Cq_total =  CQ_Polyfit(N_order_CQ+1);
    % for j=1:N_order_CQ
    %     Cq_intermediate = CQ_Polyfit(j)*J.^(N_order_CQ+1-j);
    %     Cq_total = Cq_total + Cq_intermediate;
    % end
    % % CQ = CQ0 + CQ1*J + CQ2*J^2 + CQ3*J^3;
    % CQ = Cq_total
    
    % Determines etha as a funcion of approximation polynomial
    % etamp_total =  etamp_Polyfit(N_order_etamp+1);
    % for j=1:N_order_etamp
    %     etamp_intermediate = etamp_Polyfit(j)*J.^(N_order_etamp+1-j);
    %     etamp_total = etamp_total + etamp_intermediate;
    % end
    % % ethamp = etha_mp0 + etha_mp1*J + etha_mp2*J^2 + etha_mp3*J^3
    % ethamp = etamp_total;
    ethamp = J*CT/CP;
    CQ = Cp_total/(J*ethamp);

    Nprop_correction = b_prop/2;
    Ti_eng = Nprop_correction*CT*rho*(n^2)*(D^4);
    Ti = Ti_eng*n_eng;
    Pi_eng = Nprop_correction*CP*rho*(n^3)*(D^5);
    Pi = Pi_eng*n_eng;
    Qi_eng = Nprop_correction*CQ*rho*(n^2)*(D^5);
    etha_emp = ethamp*eta_m*eta_gear*eta_esc*eta_dist;
    Pe = Pi/etha_emp;
    Pe_eng = Pe/n_eng;

    % Determines d CT/d J
    CT_dJ_total =  CT_Polyfit(N_order_CT);
    for j=1:(N_order_CT-1)
        CT_dJ_intermediate = ((N_order_CT+1-j))*CT_Polyfit(j)*J.^(N_order_CT-j);
        CT_dJ_total = CT_dJ_total + CT_dJ_intermediate;
    end
    d_CT_d_J = CT_dJ_total;
    d_T_d_J = n_eng*d_CT_d_J*rho*(n^2)*(D^4);
    d_Teng_d_J = d_CT_d_J*rho*(n^2)*(D^4);

    % Determines d CT/d V
    CT_dV_total =  CT_Polyfit(N_order_CT)/(n*D);
    for j=1:(N_order_CT-1)
        CT_dV_intermediate = ((N_order_CT+1-j))*CT_Polyfit(j)*J.^(N_order_CT-j);
        CT_dV_total = CT_dV_total + CT_dV_intermediate;
    end
    d_CT_d_V = CT_dV_total;
    d_T_d_V = n_eng*d_CT_d_V*rho*(n^2)*(D^4);
    d_Teng_d_V = d_CT_d_V*rho*(n^2)*(D^4);

    % Determines d CP/d J
    CP_dJ_total =  CP_Polyfit(N_order_CP);
    for j=1:(N_order_CP-1)
        CP_dJ_intermediate = ((N_order_CP+1-j))*CP_Polyfit(j)*J.^(N_order_CP-j);
        CP_dJ_total = CP_dJ_total + CP_dJ_intermediate;
    end
    d_CP_d_J = CP_dJ_total;
    d_P_d_J = n_eng*d_CP_d_J*rho*(n^3)*(D^5);
    d_Peng_d_J = d_CP_d_J*rho*(n^3)*(D^5);

    % Determines d CP/d V
    CP_dV_total =  CP_Polyfit(N_order_CP)/(n*D);
    for j=1:(N_order_CP-1)
        CP_dV_intermediate = ((N_order_CP+1-j))*CP_Polyfit(j)*J.^(N_order_CP-j);
        CP_dV_total = CP_dV_total + CP_dV_intermediate;
    end
    d_CP_d_V = CP_dV_total;
    d_P_d_V = n_eng*d_CP_d_V*rho*(n^3)*(D^5);
    d_Peng_d_V = d_CP_d_V*rho*(n^3)*(D^5);

    % Determines deta/dJ
    % etamp_dJ_total =  etamp_Polyfit(N_order_etamp);
    % for j=1:(N_order_etamp-1)
    %     etamp_dJ_intermediate = ((N_order_etamp+1-j))*etamp_Polyfit(j)*J.^(N_order_etamp-j);
    %     etamp_dJ_total = etamp_dJ_total + etamp_dJ_intermediate;
    % end
    d_etap_d_J = CT/CP + J*d_CT_d_J/CP - J*CT*d_CP_d_J/CP^2;

    % Determines deta/dV
    % etamp_dV_total =  etamp_Polyfit(N_order_etamp)/(n*D);
    % for j=1:(N_order_etamp-1)
    %     etamp_dV_intermediate = ((N_order_etamp+1-j))*etamp_Polyfit(j)*J.^(N_order_etamp-j);
    %     etamp_dV_total = etamp_dV_total + etamp_dV_intermediate;
    % end
    % d_etap_d_V = etamp_dV_total;
    d_J_d_V = (1/(n*D));
    d_etap_d_V = d_J_d_V*CT/CP + J*d_CT_d_J*d_J_d_V/CP - J*CT*d_CP_d_J*d_J_d_V/CP^2;

    % d_etha_d_V = 3*V^2*etha_mp3/(n^3*D^3)+2*V*etha_mp2/(n^2*D^2)+etha_mp1/(n*D)
    % d_CP_d_V = CP3*V^3/(n^3*D^3) + CP2*V^2/(n^2*D^2) + CP1*V/(n*D) + CP0;

    % Estimation of Prop data
    A_prop = pi*(D/2)^2;
    v_i = -(1/2)*V + sqrt(1/4*V^2 + Ti_eng/(2*rho*A_prop)); % Induced velocity at prop disk
    R_inf = (D/2)*sqrt((V + v_i)/(V + 2*v_i));
    v_inf = 2*v_i;

    %% Variables for Prop engine Model
    Propulsion.delta_T = delta_T;
    % Thrust
    Propulsion.Ti = Ti;
    Propulsion.Ti_eng = Ti_eng;
    % Thrust Derivatives
    Propulsion.d_T_d_V = d_T_d_V;
    Propulsion.d_Teng_d_V = d_Teng_d_V;
    Propulsion.d_CT_d_V = d_CT_d_V;
    Propulsion.d_T_d_M = d_T_d_V*a;
    Propulsion.d_CT_d_M = d_CT_d_V;
    % Power
    Propulsion.Pi_eng = Pi_eng;
    Propulsion.Pi = Pi;
    Propulsion.Pe = Pe;
    Propulsion.Pe_eng = Pe_eng;
    % power Derivatives
    Propulsion.d_P_d_V = d_P_d_V;
    Propulsion.d_Peng_d_V = d_Peng_d_V;
    Propulsion.d_CP_d_V = d_CP_d_V;
    Propulsion.d_P_d_M = d_P_d_V*a;
    Propulsion.d_CP_d_M = d_CP_d_V*a;
    % Induced Velocity
    Propulsion.v_i = v_i; % Induced velocity at prop disk
    Propulsion.R_inf = R_inf; % Radius of prop wash at infinity
    Propulsion.v_inf = v_inf; % induced velocity at propwash at infinity
    Propulsion.delta_T = delta_T;

    Propulsion.Ti = Ti;
    Propulsion.RPM_max = RPM_max;
    Propulsion.RPM = RPM;
    Propulsion.RPS = n;
    Propulsion.Omega = Omega;
    Propulsion.Qi_eng = Qi_eng;

    Propulsion.J = J;
    Propulsion.CT = CT;
    Propulsion.CP = CP;
    Propulsion.CQ = CQ;

    Propulsion.ethamp = ethamp;
    Propulsion.etha_emp = etha_emp;

    Propulsion.d_T_d_J = d_T_d_J;
    Propulsion.d_CT_d_J = d_CT_d_J;
    Propulsion.d_P_d_J = d_P_d_J;
    Propulsion.d_Peng_d_J = d_Peng_d_J;
    Propulsion.d_CP_d_J = d_CP_d_J;

    Propulsion.d_etap_d_V = d_etap_d_V;
    Propulsion.d_etap_d_J = d_etap_d_J;

    % Dummy for Electric Engine
    Propulsion.C = 0;
    Propulsion.correccion = 0;

else

    % Coefficients correctores consumo específico
    a_1=3.559957437510763;
    a_2=-10.739698199171459;
    a_3= 11.989635150373475;
    a_4=-5.869876557884609;
    a_5=2.059994459180667;

    delta_T = (Fdes/n_eng)/(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
    correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
    C = c_SL.*correccion.*(1+1.44*M).*sqrt(Temp/Temp_SL).* V./eta_p;

    RPM_nominal = (RPM_max-RPM_min)*delta_T;

    Ti_eng = delta_T.*(T_SL.*eta_p./V.*(1+0.2.*M.^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
    Ti = Ti_eng*n_eng;
    Pi_eng = Ti_eng*V;
    Pi = Pi_eng*n_eng;
    Pe = Pi/eta_p;
    Pe_eng = Pe/n_eng;
    CT = Ti_eng / (qbar*S_ref);
    CP = Pi_eng / (qbar*S_ref*cmac);
    CQ = CP/(2*pi);
    Q_eng = CQ;

    % Propulsive derivatives as a funcion of Speed
    dT_dV = 1.400000000*delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^2.500000000*rho*Temp/(rho_SL*Temp_SL*a^2)-...
        delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^3.500000000*rho*Temp/(V^2*rho_SL*Temp_SL);
    dCT_dV = 2.800000000*delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^2.500000000*Temp/(V^2*rho_SL*Temp_SL*S_ref*a^2)-...
        6.000000000*delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^3.500000000*Temp/(V^4*rho_SL*Temp_SL*S_ref);
    dP_dV = 1.400000000*delta_T*T_SL*(1+.2*V^2/a^2)^2.500000000*rho*Temp*V/(rho_SL*Temp_SL*a^2);
    dCP_dV = 2.800000000*delta_T*T_SL*(1+.2*V^2/a^2)^2.500000000*Temp/(V*rho_SL*Temp_SL*S_ref*cmac*a^2)-...
        4.000000000*delta_T*T_SL*(1+.2*V^2/a^2)^3.500000000*Temp/(V^3*rho_SL*Temp_SL*S_ref*cmac);

    % Propulsive derivatives as a funcion of Mach
    dT_dM = 1.400000000*delta_T*T_SL*eta_p*(1+.2*M^2)^2.500000000*rho*Temp/(a*rho_SL*Temp_SL)-...
        delta_T*T_SL*eta_p*(1+.2*M^2)^3.500000000*rho*Temp/(M^2*a*rho_SL*Temp_SL);
    dCT_dM = 2.800000000*delta_T*T_SL*eta_p*(1+.2*M^2)^2.500000000*Temp/(M^2*a^3*rho_SL*Temp_SL*S_ref)-...
        6.000000000*delta_T*T_SL*eta_p*(1+.2*M^2)^3.500000000*Temp/(M^4*a^3*rho_SL*Temp_SL*S_ref);
    dP_dM = 1.400000000*delta_T*T_SL*(1+.2*M^2)^2.500000000*rho*Temp*M/(rho_SL*Temp_SL);
    dCP_dM = 2.800000000*delta_T*T_SL*(1+.2*M^2)^2.500000000*Temp/(M*a^2*rho_SL*Temp_SL*S_ref*cmac)-...
        4.000000000*delta_T*T_SL*(1+.2*M^2)^3.500000000*Temp/(M^3*a^2*rho_SL*Temp_SL*S_ref*cmac);

    % Estimation of Prop data
    A_prop = pi*(D/2)^2;
    v_i = -(1/2)*V + sqrt(1/4*V^2 + Ti_eng/(2*rho*A_prop)); % Induced velocity at prop disk
    R_inf = (D/2)*sqrt((V + v_i)/(V + 2*v_i));
    v_inf = 2*v_i;

    %% Variables for Prop engine Model
    Propulsion.delta_T = delta_T;
    Propulsion.C = C;
    Propulsion.correccion = correccion;
    % Thrust
    Propulsion.Ti = Ti;
    Propulsion.Ti_eng = Ti_eng;
    % Thrust Derivatives
    Propulsion.d_T_d_V = dT_dV;
    Propulsion.d_CT_d_V = dCT_dV;
    Propulsion.d_T_d_M = dT_dM;
    Propulsion.d_CT_d_M = dCT_dM;
    % Power
    Propulsion.Pi_eng = Pi_eng;
    Propulsion.Pi = Pi;
    Propulsion.Pe = Pe;
    Propulsion.Pe_eng = Pe_eng;
    % power Derivatives
    Propulsion.d_P_d_V = dP_dV;
    Propulsion.d_CP_d_V = dCP_dV;
    Propulsion.d_P_d_M = dP_dM;
    Propulsion.d_CP_d_M = dCP_dM;
    % Induced Velocity
    Propulsion.v_i = v_i; % Induced velocity at prop disk
    Propulsion.R_inf = R_inf; % Radius of prop wash at infinity
    Propulsion.v_inf = v_inf; % induced velocity at propwash at infinity
    Propulsion.delta_T = delta_T;

    n = RPM_nominal/60;
    Omega = 2*pi*n;

    % Dummy variables for propulsion
    Propulsion.RPM_max = RPM_max;
    Propulsion.RPM = RPM_nominal;
    Propulsion.RPS = n;
    Propulsion.Omega = Omega;
    Propulsion.Qi_eng = Q_eng;

    J = V/(n*D);

    Propulsion.J = J;
    Propulsion.CP = CP;
    Propulsion.CQ = CQ;

    Propulsion.ethamp = 0.82;
    Propulsion.etha_emp = 0;

    Propulsion.d_T_d_J = 0;
    Propulsion.d_CT_d_J = 0;
    Propulsion.d_P_d_J = 0;
    Propulsion.d_Peng_d_J = 0;
    Propulsion.d_CP_d_J = 0;

    Propulsion.d_etap_d_V = 0;
    Propulsion.d_etap_d_J = 0;
end


