function  D_p = determination_D_prop_horizontal(Aero_TH,pp,RPMMAX_APC,Vh,Prop_data,Weight_tier,conv_UNITS,Performance_preliminar,Geo_tier,AC_CONFIGURATION,Drag)

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
g = conv_UNITS.g;

m_TOW = Weight_tier.m_TOW;
W = m_TOW*g;

S_ref = Geo_tier.S_ref;

rho = Performance_preliminar.rho;

n_eng = AC_CONFIGURATION.n_eng;

% q_inf = 0.5*rho*Vh^2;
% CL = (m_TOW*g)/(q_inf*S_ref);
% CD0 = Aero_TH.CD0;
% CD1 = Aero_TH.CD1;
% CD2 = Aero_TH.CD2;
% CD = CD0 + CD1*CL + CD2*CL^2;
% Drag = CD*q_inf*S_ref;

% % CT_Polyfit{k} = CT{k}(3) + CT{k}(2).*XData_CT{k} + CT{k}(1).*XData_CT{k}.^2;
% CT0 = Prop_data.CT0;
% CT1 = Prop_data.CT1;
% CT2 = Prop_data.CT2;
% % CP_Polyfit{k} = CP{k}(4) + CP{k}(3).*XData_CP{k} + CP{k}(2).*XData_CP{k}.^2 + CP{k}(1).*XData_CP{k}.^3;         
% CP0 = Prop_data.CP0;
% CP1 = Prop_data.CP1;
% CP2 = Prop_data.CP2;
% CP3 = Prop_data.CP3;
% % CQ_Polyfit{k} = CQ{k}(4) + CQ{k}(3).*XData_CQ{k} + CQ{k}(2).*XData_CQ{k}.^2 + CQ{k}(1).*XData_CQ{k}.^3;
% CQ0 = Prop_data.CQ0;
% CQ1 = Prop_data.CQ1;
% CQ2 = Prop_data.CQ2;
% CQ3 = Prop_data.CQ3;
% % ethamp_Polyfit{k} = ethamp{k}(4) + ethamp{k}(3).*XData_ethamp{k} + ethamp{k}(2).*XData_ethamp{k}.^2 + ethamp{k}(1).*XData_ethamp{k}.^3;
% etha_mp0 = Prop_data.etha_mp0;
% etha_mp1 = Prop_data.etha_mp1;
% etha_mp2 = Prop_data.etha_mp2;
% etha_mp3 = Prop_data.etha_mp3;

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

n_MAX = @(D_p) (pp*RPMMAX_APC/(100*D_p/2.54))/60; 
J = @(D_p) Vh/(n_MAX(D_p)*D_p);

% Determines CT as a funcion of approximation polynomial
Ct_total =  CT_Polyfit(N_order_CT+1);
for j=1:N_order_CT
    Ct_intermediate = @(D_p) CT_Polyfit(j)*J(D_p).^(N_order_CT+1-j);
    Ct_total = @(D_p) Ct_total + Ct_intermediate(D_p);
end
CT = @(D_p) Ct_total(D_p);

% n_MAX = @(D_p) (pp*RPMMAX_APC/(100*D_p/2.54))/60; 
% J = @(D_p) Vh/(n_MAX(D_p)*D_p);
CT0 = 0.0891;   
CT1 = -0.0267;  
CT2 = -0.1626;
CT = @(D_p) CT0 + CT1*J(D_p) + CT2*J(D_p)^2;
% CP = CP0 + CP1*J(D_prop) + CP2*J(D_prop)^2 + CP3*J(D_prop)^3;
% CQ = CQ0 + CQ1*J(D_prop) + CQ2*J(D_prop)^2 + CQ3*J(D_prop)^3;
Ti_eng = @(D_p) n_eng*CT(D_p)*rho*(n_MAX(D_p)^2)*(D_p^4);  
% ethamp = etha_mp0 + etha_mp1*J + etha_mp2*J^2 + etha_mp3*J^3;
D_p = fzero(@(D_p) Drag - Ti_eng(D_p),22*2.54/100);
if isnan(D_p)
    Warning = 'WARNING!!! fzero not producing solution for the Horizontal Speed interval - Code in PAUSE';
    disp(Warning)
end

