function n = Selection_J_CT_F_des_v2(V,Prop_data,T_eng,rho,p_ct,p_ct_red,OUTPUT_read_XLSX)

% Estimacion iniciali apara inicialización
D = Prop_data.D_prop;
n_eng = Prop_data.n_eng;
CT0 = p_ct_red(3);
CT1 = p_ct_red(2);
CT2 = p_ct_red(1);

b_prop = OUTPUT_read_XLSX.Propulsive_flags.b_prop;
Nprop_correction = b_prop/2;

% D = Prop_data.D_prop;
% n_eng = Prop_data.n_eng;
% CT0 = Prop_data.CT_Polyfit(3);
% CT1 = Prop_data.CT_Polyfit(2);
% CT2 = Prop_data.CT_Polyfit(1);

% N_order_CT = Prop_data.N_order_CT;
% Ajuste polinómico 
N_order_CT = length(p_ct)-1;
CT_Polyfit = p_ct;

% CT_Polyfit = Prop_data.CT_Polyfit;
% CP_Polyfit = Prop_data.CP_Polyfit;
% CQ_Polyfit = Prop_data.CQ_Polyfit;
% etamp_Polyfit = Prop_data.etamp_Polyfit;

% Calculates the Advanced Parameter Ratio
J = @(n) V/(n*D);
% Determines CT as a funcion of approximation polynomial
Ct_total =  @(n) CT_Polyfit(N_order_CT+1);
for j=1:N_order_CT
    Ct_intermediate = @(n) CT_Polyfit(j)*J(n).^(N_order_CT+1-j);
    Ct_total = @(n) (Ct_total(n) + Ct_intermediate(n));
end
CT = @(n) Ct_total(n);

% initialization
n_sel = -(1/2)*(CT1*V*n_eng*rho*D-sqrt(-4*CT0*CT2*D^2*V^2*n_eng^2*rho^2 + CT1^2*D^2*V^2*n_eng^2*rho^2 + ...
    4*CT0*T_eng*n_eng*rho))/(CT0*n_eng*rho*D^2);

% J = @(n) V/(n*D);
% CT = @(n) CT0 + CT1*J(n) + CT2*J(n)^2;
Ti_eng = @(n) (Nprop_correction*n_eng*CT(n)*rho*(n^2)*(D^4));  
n = fzero(@(n) (T_eng - Ti_eng(n)),n_sel);

% J = @(n) V/(n*D);
% CT = @(n) CT0 + CT1*J(n) + CT2*J(n)^2;
% Ti_eng = @(n) n_eng*CT(n)*rho*(n^2)*(D^4);  
% n = fzero(@(n) T_eng - Ti_eng(n),n_sel);