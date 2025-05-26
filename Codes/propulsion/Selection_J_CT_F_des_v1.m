function n = Selection_J_CT_F_des_v1(V,Prop_data,T_eng,rho,OUTPUT_read_XLSX)

b_prop = OUTPUT_read_XLSX.Propulsive_flags.b_prop;
% Corrects for the number of prop blades
if b_prop< 2
    Warning = 'WARNING!!! The number of prop blades cannot be smaller than 2. Please Check the Input data b_prop and correct.';
                disp(Warning)
                pause
end

Nprop_correction = b_prop/2;

D = Prop_data.D_prop;
n_eng = Prop_data.n_eng;
CT0 = Prop_data.CT_Polyfit(3);
CT1 = Prop_data.CT_Polyfit(2);
CT2 = Prop_data.CT_Polyfit(1);

N_order_CT = Prop_data.N_order_CT;
N_order_CP = Prop_data.N_order_CP;
N_order_CQ = Prop_data.N_order_CQ;
N_order_etamp = Prop_data.N_order_etamp;

CT_Polyfit = Prop_data.CT_Polyfit;
CP_Polyfit = Prop_data.CP_Polyfit;
CQ_Polyfit = Prop_data.CQ_Polyfit;
etamp_Polyfit = Prop_data.etamp_Polyfit;

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

% To avoid using initial value imaginary when the equation out of range
n_sel = real(n_sel);

% J = @(n) V/(n*D);
% CT = @(n) CT0 + CT1*J(n) + CT2*J(n)^2;
Ti_eng = @(n) (Nprop_correction*n_eng*CT(n)*rho*(n^2)*(D^4));  
n = fzero(@(n) (T_eng - Ti_eng(n)),n_sel);
dummy = 0;
% J = @(n) V/(n*D);
% CT = @(n) CT0 + CT1*J(n) + CT2*J(n)^2;
% Ti_eng = @(n) n_eng*CT(n)*rho*(n^2)*(D^4);  
% n = fzero(@(n) T_eng - Ti_eng(n),n_sel);
