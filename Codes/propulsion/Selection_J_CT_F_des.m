function n = Selection_J_CT_F_des(V,Prop_data,T_eng,rho)

D = Prop_data.D_prop;
n_eng = Prop_data.n_eng;
CT0 = Prop_data.CT_Polyfit(3);
CT1 = Prop_data.CT_Polyfit(2);
CT2 = Prop_data.CT_Polyfit(1);

n_sel = -(1/2)*(CT1*V*n_eng*rho*D-sqrt(-4*CT0*CT2*D^2*V^2*n_eng^2*rho^2 + CT1^2*D^2*V^2*n_eng^2*rho^2 + ...
    4*CT0*T_eng*n_eng*rho))/(CT0*n_eng*rho*D^2);

J = @(n) V/(n*D);
CT = @(n) CT0 + CT1*J(n) + CT2*J(n)^2;
Ti_eng = @(n) n_eng*CT(n)*rho*(n^2)*(D^4);  
n = fzero(@(n) T_eng - Ti_eng(n),n_sel);


