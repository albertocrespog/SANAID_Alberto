function PROP_var = get_Propulsion_v1(m_TO_1,V,conv_UNITS,S_w,...
    h,C_D0,C_D1,C_D2,rho,prop_NEW)
% (m_no_fuel,m_fuel,V_cruise,Performance_calc_ITER_h,...
%     a,Aero_Polar_PI,conv_UNITS,Geo_calc_ITER,CG,h

g = conv_UNITS.g;
W2hp = conv_UNITS.W2hp;

 
% % Calculo CL vuelo
% [Temp,rho,p,a]=atmos_inter_mio(h);
q = 0.5*rho*V^2;

Delta_m = 0.050;
Delta_p = 0.010;

W_TO_1 = g*m_TO_1;
CL = W_TO_1/(S_w*q);
CD = C_D0 + C_D1*CL + C_D2*CL^2;
D = S_w*q*CD;
L = S_w*q*CL;
T = D;
P = T*V;
T_calc = 0;
delta_p = 0.1;
while T_calc < T

    if prop_NEW == 1
        [T_calc,c,Vnmax,delta_act] = propulsive_model_v5(V,h,delta_p); % Actualizacion Oct 2015
        else
            [T_calc,c,Vnmax] = propulsive_model_v3(V,h,delta_p);
    end
        
    delta_p = delta_p + Delta_p;
end
delta_p_calc = delta_p;

% delta_p_medio = (max(delta_p_calc) + min(delta_p_calc))/2;
% V_max = V/delta_p_calc;
V_max = 50;
% V_max = Vnmax;

V_vec = linspace(0.80*V,V_max,200);


for i=1:length(V_vec)
        if prop_NEW == 1
            [T_calc(i),c(i),Vnmax,delta_act] = propulsive_model_v5(V_vec(i),h,delta_p_calc); % Revisión Oct 2015
        else
            [T_calc(i),c(i),Vnmax] = propulsive_model_v3(V_vec(i),h,delta_p_calc);
        end
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

PROP_var.P_fit_coef = P_fit_coef;
PROP_var.P_poly = P_poly;
PROP_var.c_poly = c_poly;
PROP_var.delta_p_calc = delta_p_calc;

SAVE_FIGS=0;
LS = 2;
TS = 10;

Figures = 1;
if Figures == 1
%     figure(1)
%     plot(M_vec,T,'LineWidth', LS)
%     title('T vs. mass')
%     xlabel('M (kg)')
%     ylabel('Thrust (N)')
%     grid on
%     
%     figure(2)
%     plot(M_vec,CL,'LineWidth', LS)
%     title('C_L vs. mass')
%     xlabel('M (kg)')
%     ylabel('C_L')
%     grid on
%     
%     figure(3)
%     plot(M_vec,CD,'LineWidth', LS)
%     title('C_D vs. mass')
%     xlabel('M (kg)')
%     ylabel('C_D')
%     grid on
%     
%     figure(4)
%     plot(M_vec,L,'LineWidth', LS)
%     title('L vs. mass')
%     xlabel('M (kg)')
%     ylabel('Lift (N)')
%     grid on
%     
%     figure(5)
%     plot(M_vec,delta_p_calc,'LineWidth', LS)
%     title('\delta_p vs. mass')
%     xlabel('M (kg)')
%     ylabel('\delta_p')
%     grid on
%     
%     figure(6)
%     plot(M_vec,P*W2hp,'LineWidth', LS)
%     title('P vs. mass')
%     xlabel('M (kg)')
%     ylabel('Power (hp)')
%     grid on
    
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

% if PI == 0
%     save Data_P_poly.mat P_poly
% elseif PI==1
save Data_P_poly_v5.mat P_poly
% end