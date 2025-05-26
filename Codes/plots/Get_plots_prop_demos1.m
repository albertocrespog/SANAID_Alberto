

% [Propulsion] = get_EngineProperties_v3(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,n_eng)
% [Propulsion] =get_EngineProperties_AP(h,V,alpha,beta,Geo_tier,Fdes,AC_CONFIGURATION)


g = conv_UNITS.g;
W2hp = conv_UNITS.W2hp;

S_w = Geo_calc_ITER.S_w;
C_D0_CR = Aero_Polar_PI.C_D0_CR_1;
C_D1_CR = Aero_Polar_PI.C_D1_CR_1;
C_D2_CR = Aero_Polar_PI.C_D2_CR_1;

% V_min = Performance_calc_ITER_h.V_min;
% V_opt_E = Performance_calc_ITER_h.V_opt_E;
V_stall_CR = Performance_calc_ITER_h.V_stall_CR;
% V_stall_TO = Performance_calc_ITER_h.V_stall_TO;
% V_stall_LD = Performance_calc_ITER_h.V_stall_LD;

h = 3000;

% Calculo CL vuelo
[Temp,rho,p,a]=atmos_inter_mio(h);
q = 0.5*rho*V_cruise^2;

m_TO_1 = m_no_fuel + m_fuel;

Delta_m = 0.050;
M_vec = (m_TO_1:-Delta_m:(m_TO_1-m_fuel));

for i=1:length(M_vec)
    W_TO_1(i) = g*M_vec(i);
    CL_CR(i) = W_TO_1(i)/(S_w*q);
    CD_CR(i) = C_D0_CR + C_D1_CR*CL_CR(i) + C_D2_CR*CL_CR(i)^2;
    D(i) =  S_ref*q*CD_CR(i);
    L(i) =  S_ref*q*CL_CR(i);
    T(i) = D(i);
    P(i) = T(i)*V_cruise;
    T_calc = 0;
    delta_p = 0.1;
    while T_calc < T(i)
        [T_calc,c,Vnmax] = propulsive_model_v3(V_cruise,h,delta_p);
        delta_p = delta_p + 0.01;
    end
    delta_p_calc(i) = delta_p;
end


switch tipo_motor
    case 1
        [Propulsion] = get_EngineProperties_AP(h,V,alpha,beta,Geo_tier,Fdes,AC_CONFIGURATION);
        % Storing Data
        Propulsion.Ti = Ti;
        Propulsion.Ti_eng = Ti_eng;
        Propulsion.dT_dV = dT_dV;
        Propulsion.dCT_dV = dCT_dV;
        Propulsion.T = T;
        Propulsion.dT_dV = dT_dV;
        Propulsion.dCT_dV = dCT_dV;
    case 2
        [Propulsion] = get_EngineProperties_AP(h,V,alpha,beta,Geo_tier,Fdes,AC_CONFIGURATION);
        % Storing Data
        Propulsion.Ti = Ti;
        Propulsion.Ti_eng = Ti_eng;
        Propulsion.dT_dV = dT_dV;
        Propulsion.dCT_dV = dCT_dV;
        Propulsion.Pi_eng = Pi_eng;
        Propulsion.Pi = Pi;
        Propulsion.Pe = Pe;
        Propulsion.Pe_eng = Pe_eng;
        Propulsion.dP_dV = dP_dV;
        Propulsion.dCP_dV = dCP_dV;
    case 3
        [Propulsion] = get_EngineProperties_AP(h,V,alpha,beta,Geo_tier,Fdes,AC_CONFIGURATION);
        % Storing Data
        Propulsion.Ti = Ti;
        Propulsion.Ti_eng = Ti_eng;
        Propulsion.dT_dV = dT_dV;
        Propulsion.dCT_dV = dCT_dV;
        Propulsion.Pi_eng = Pi_eng;
        Propulsion.Pi = Pi;
        Propulsion.Pe = Pe;
        Propulsion.Pe_eng = Pe_eng;
        Propulsion.dP_dV = dP_dV;
        Propulsion.dCP_dV = dCP_dV;
    case 4
        [Propulsion] = get_EngineProperties_v3(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,n_eng);
        % Storing Data
        Propulsion.Ti = Ti;
        Propulsion.Ti_eng = Ti_eng;
        Propulsion.dT_dV = dT_dV;
        Propulsion.dCT_dV = dCT_dV;
        Propulsion.Pi_eng = Pi_eng;
        Propulsion.Pi = Pi;
        Propulsion.Pe = Pe;
        Propulsion.Pe_eng = Pe_eng;
        Propulsion.dP_dV = dP_dV;
        Propulsion.dCP_dV = dCP_dV;
end

