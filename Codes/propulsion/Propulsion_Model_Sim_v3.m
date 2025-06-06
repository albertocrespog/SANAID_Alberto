function P_poly = Propulsion_Model_Sim_v3(m_no_fuel,m_fuel,V_cruise,Performance_calc_ITER_h,...
    a,Aero_Polar_PI,conv_UNITS,Geo_calc_ITER,CG,h)

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
    D(i) =  S_w*q*CD_CR(i);
    L(i) =  S_w*q*CL_CR(i);
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

delta_p_medio = (max(delta_p_calc) + min(delta_p_calc))/2;
V_max = V_cruise/delta_p_medio;
V_max = 50;
V_vec = linspace(V_stall_CR,V_max, 200)

delta_p_medio = 0.10;

for i=1:length(V_vec)
    [T_calc(i),c(i),Vnmax] = propulsive_model_v3(V_vec(i),h,delta_p_medio);
end

%Determinación de la penalización por altitud: Modelo de Raymer.
delta_rho = (1-6.5e-3*h/288.15)^4.25193; %Troposfera ISA.
rho       = 1.225*delta_rho;
delta_h   = delta_rho-(1-delta_rho)/7.55;

for i=1:length(T_calc)
    P_calc(i) = T_calc(i)*V_vec(i);
end

P_fit_coef = polyfit(V_vec,P_calc,2);
P_fit_c = polyfit(V_vec,c,3);

for i=1:length(V_vec)
    P_poly(i) = P_fit_coef(3) + P_fit_coef(2)*V_vec(i) + P_fit_coef(1)*V_vec(i)^2;
    c_poly(i) = P_fit_c(4) + P_fit_c(3)*V_vec(i) + P_fit_c(2)*V_vec(i)^2 + P_fit_c(1)*V_vec(i)^3;
end

SAVE_FIGS=0;
LS = 2;
TS = 10;

Figures = 1;
if Figures == 1
    figure(1)
    plot(M_vec,T,'LineWidth', LS)
    title('T vs. mass')
    xlabel('M (kg)')
    ylabel('Thrust (N)')
    grid on
    
    figure(2)
    plot(M_vec,CL_CR,'LineWidth', LS)
    title('C_L vs. mass')
    xlabel('M (kg)')
    ylabel('C_L')
    grid on
    
    figure(3)
    plot(M_vec,CD_CR,'LineWidth', LS)
    title('C_D vs. mass')
    xlabel('M (kg)')
    ylabel('C_D')
    grid on
    
    figure(4)
    plot(M_vec,L,'LineWidth', LS)
    title('L vs. mass')
    xlabel('M (kg)')
    ylabel('Lift (N)')
    grid on
    
    figure(5)
    plot(M_vec,delta_p_calc,'LineWidth', LS)
    title('\delta_p vs. mass')
    xlabel('M (kg)')
    ylabel('\delta_p')
    grid on
    
    figure(6)
    plot(M_vec,P*W2hp,'LineWidth', LS)
    title('P vs. mass')
    xlabel('M (kg)')
    ylabel('Power (hp)')
    grid on
    
    figure(7)
    plot(V_vec,T_calc,'LineWidth', LS)
    title('T vs. V')
    xlabel('V (m/s)')
    ylabel('T (N)')
    grid on
    
    figure(8)
    plot(V_vec,P_calc,'LineWidth', LS)
    hold on
    plot(V_vec,P_poly,'r','LineWidth', LS)
    title('P vs. V')
    xlabel('V (m/s)')
    ylabel('P (W)')
    grid on

    figure(9)
    plot(V_vec,c,'LineWidth', LS)
    hold on
    plot(V_vec,c_poly,'r','LineWidth', LS)
    title('c vs. V')
    xlabel('V (m/s)')
    ylabel('c (W)')
    grid on

end

if PI == 0
    save Data_P_poly.mat P_poly
elseif PI==1
    save Data_P_poly_PI.mat P_poly
end
