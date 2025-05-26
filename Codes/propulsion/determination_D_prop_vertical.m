function  D_p = determination_D_prop_vertical(Aero_TH,pp,RPMMAX_APC,Vv,Prop_data,Weight_tier,conv_UNITS,Performance_preliminar,Geo_tier,AC_CONFIGURATION)

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
g = conv_UNITS.g;

m_TOW = Weight_tier.m_TOW;
W = m_TOW*g;

S_ref = Geo_tier.S_ref;

rho = Performance_preliminar.rho;
n_eng = AC_CONFIGURATION.n_eng;

C_DV = 0.5;

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
J = @(D_p) Vv/(n_MAX(D_p)*D_p);

% Determines CT as a funcion of approximation polynomial
Ct_total =  CT_Polyfit(N_order_CT+1);
for j=1:N_order_CT
    Ct_intermediate = @(D_p) CT_Polyfit(j)*J(D_p).^(N_order_CT+1-j);
    Ct_total = @(D_p) Ct_total + Ct_intermediate(D_p);
end
% CT = CT0 + CT1*J + CT2*J^2;
CT = @(D_p) Ct_total(D_p);

CT0 = 0.0891;   
CT1 = -0.0267;  
CT2 = -0.1626;
CT = @(D_p) CT0 + CT1*J(D_p) + CT2*J(D_p)^2;

% CT = @(D_p) CT0 + CT1*J(D_p) + CT2*J(D_p)^2;
% CP = CP0 + CP1*J(D_prop) + CP2*J(D_prop)^2 + CP3*J(D_prop)^3;
% CQ = CQ0 + CQ1*J(D_prop) + CQ2*J(D_prop)^2 + CQ3*J(D_prop)^3;
Ti_eng = @(D_p) n_eng*CT(D_p)*rho*(n_MAX(D_p)^2)*(D_p^4);  
% ethamp = etha_mp0 + etha_mp1*J + etha_mp2*J^2 + etha_mp3*J^3;
T_V = 0.5*rho*Vv^2*S_ref*C_DV + W;
D_p = fzero(@(D_p) T_V - Ti_eng(D_p),22*2.54/100);
if isnan(D_p)
    Warning = 'WARNING!!! fzero not producing solution for the Vertical Speed interval - Code in PAUSE';
    disp(Warning)
    pause
end
