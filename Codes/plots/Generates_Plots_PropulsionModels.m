function [Fig] = Generates_Plots_PropulsionModels(Prop_data,Data_Trim,V,alpha,D,Plot_Options,Fig)

% n = Prop_data.n;
% Density
rho = Data_Trim.rho;
% Parameter advanced ratio

eta_gear = Prop_data.eta_gear;
eta_m = Prop_data.eta_m;
eta_esc = Prop_data.eta_esc;
eta_dist = Prop_data.eta_dist;

N_order_CT = Prop_data.N_order_CT;
N_order_CP = Prop_data.N_order_CP;
N_order_CQ = Prop_data.N_order_CQ;
N_order_etamp = Prop_data.N_order_etamp;

CT_Polyfit = Prop_data.CT_Polyfit;
CP_Polyfit = Prop_data.CP_Polyfit;
CQ_Polyfit = Prop_data.CQ_Polyfit;
etamp_Polyfit = Prop_data.etamp_Polyfit;

% Calculates the Advanced PArameter Ratio
J_vec = linspace(0,1,100);
% Determines CT as a funcion of approximation polynomial
Ct_total =  CT_Polyfit(N_order_CT+1);
for j=1:N_order_CT
    Ct_intermediate = CT_Polyfit(j)*J.^(N_order_CT+1-j);
    Ct_total = Ct_total + Ct_intermediate;
end
% CT = CT0 + CT1*J + CT2*J^2;
CT_vec = Ct_total;
% Limits J to the Jmax for CT = 0
J_max = interp1(CT_vec,J_vec,0,'spline');
J = linspace(0,J_max,100);

% Determines CT as a funcion of approximation polynomial
Ct_total =  CT_Polyfit(N_order_CT+1);
for j=1:N_order_CT
    Ct_intermediate = CT_Polyfit(j)*J.^(N_order_CT+1-j);
    Ct_total = Ct_total + Ct_intermediate;
end
% CT = CT0 + CT1*J + CT2*J^2;
CT_vec = Ct_total;

% Determines CP as a funcion of approximation polynomial
Cp_total =  CP_Polyfit(N_order_CP+1);
for j=1:N_order_CP
    Cp_intermediate = CP_Polyfit(j)*J.^(N_order_CP+1-j);
    Cp_total = Cp_total + Cp_intermediate;
end
% CP = CP0 + CP1*J + CP2*J^2 + CP3*J^3;
CP_vec = Cp_total;

% Determines CP as a funcion of approximation polynomial
Cq_total =  CQ_Polyfit(N_order_CQ+1);
for j=1:N_order_CQ
    Cq_intermediate = CQ_Polyfit(j)*J.^(N_order_CQ+1-j);
    Cq_total = Cq_total + Cq_intermediate;
end
% CQ = CQ0 + CQ1*J + CQ2*J^2 + CQ3*J^3;
CQ_vec = Cq_total;

% Determines CP as a funcion of approximation polynomial
etamp_total =  etamp_Polyfit(N_order_etamp+1);
for j=1:N_order_etamp
    etamp_intermediate = etamp_Polyfit(j)*J.^(N_order_etamp+1-j);
    etamp_total = etamp_total + etamp_intermediate;
end
% ethamp = etha_mp0 + etha_mp1*J + etha_mp2*J^2 + etha_mp3*J^3
% Propulsive efficiency
etha_ep_vec = etamp_total;
% Motor efficiency
etha_emp_vec = eta_m*eta_gear*eta_esc*eta_dist;

PRINT_PLOTS_CT_Model = Plot_Options.PRINT_PLOTS_CT_Model;
% % Polyfit Coefficients for APC 22x10
% % CT_Polyfit{k} = CT{k}(3) + CT{k}(2).*XData_CT{k} + CT{k}(1).*XData_CT{k}.^2;
% CT0 = Prop_data.CT0;
% CT1 = Prop_data.CT1;
% CT2 = Prop_data.CT2;
% 
% CP0 = Prop_data.CT0;
% CP1 = Prop_data.CP1;
% CP2 = Prop_data.CP2;
% CP3 = Prop_data.CP3;
% 
% CQ0 = Prop_data.CQ0;
% CQ1 = Prop_data.CQ1;
% CQ2 = Prop_data.CQ2;
% CQ3 = Prop_data.CQ3;
% 
% etha_mp3 = Prop_data.etha_mp3;
% etha_mp2 = Prop_data.etha_mp2;
% etha_mp1 = Prop_data.etha_mp1;
% etha_mp0 = Prop_data.etha_mp0;
% 
% CT_vec = CT0 + CT1*J_vec + CT2*(J_vec.^2);
% % Limits J to the Jmax
% J_max = interp1(CT_vec,J_vec,0,'spline');
% J_vec = linspace(0,J_max,100);
% CT_vec = CT0 + CT1*J_vec + CT2*(J_vec.^2);
% CP_vec = CP0 + CP1*J_vec + CP2*(J_vec.^2) + CP3*(J_vec.^3);
% CQ_vec = CQ0 + CQ1*J_vec + CQ2*(J_vec.^2) + CQ3*(J_vec.^3);
% etha_vec = etha_mp0 + etha_mp1*J_vec + etha_mp2*(J_vec.^2) + etha_mp3*(J_vec.^3);

%  PRINT_PLOTS_PROPULSION
% if PRINT_PLOTS_CT_Model ==1
    Fig = Fig + 1;
    figure(Fig)
    plot(J_vec,CT_vec)
    title('C_T vs. J')
    xlabel('J')
    ylabel('C_T')
    grid on
    
    Fig = Fig + 1;
    figure(Fig)
    plot(J_vec,CP_vec)
    title('C_P vs. J')
    xlabel('J')
    ylabel('C_P')
    grid on
    
    Fig = Fig + 1;
    figure(Fig)
    plot(J_vec,CQ_vec)
    title('C_Q vs. J')
    xlabel('J')
    ylabel('C_P')
    grid on
    
    Fig = Fig + 1;
    figure(Fig)
    plot(J_vec,etha_emp_vec)
    title('\eta_{MP} vs. J')
    xlabel('J')
    ylabel('\eta_{MP}')
    grid on
    
    Fig = Fig + 1;
    figure(Fig)
    plot(J_vec,etha_ep_vec)
    title('\eta_{P} vs. J')
    xlabel('J')
    ylabel('\eta_{P}')
    grid on
% end