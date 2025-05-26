function Stab_Der = getCnTbeta(modelo,Stab_Der)
n_eng        = modelo.propulsion.engines;
n_blades     = modelo.propulsion.blades;
beta_blade   = modelo.propulsion.beta;
X_Prop       = modelo.propulsion.X;
w_R_30       = modelo.propulsion.wR30;
w_R_60       = modelo.propulsion.wR60;
w_R_90       = modelo.propulsion.wR90;
D_prop       = modelo.propulsion.D;

Sref         = modelo.general.Sref;
CLa_WB       = modelo.general.CLa_WB;
downw        = modelo.general.downwash;
Xcg          = modelo.general.Xcg;

cr_we        = modelo.ala.MAC_e;
XLE_w        = modelo.ala.XLE;
S_w          = modelo.ala.S;
b_w          = modelo.ala.b;

l_fus        = modelo.fuselaje.l;

m22ft2       = modelo.conversion.m22ft2;
m2ft         = modelo.conversion.m2ft;

%% CnT_beta Airplane Design (2088 PDF)
for jj = 1:n_eng
    depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop, cr_we, XLE_w, l_fus, downw, CLa_WB*S_w/Sref); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
    CNa_807         = CNa_807_calc(n_blades, beta_blade);
    KN              = 262*(2*w_R_30/D_prop) + 262*(2*w_R_60/D_prop) + 135*(2*w_R_90/D_prop);
    dCN_dalpha(jj)  = CNa_807*(1 + 0.8*(KN/80.7 - 1)); 
end

CnT_beta    = (pi/4)/(Sref*m22ft2)/(b_w*m2ft)*sum(((D_prop*m2ft).^2).*(X_Prop - Xcg).*dCN_dalpha);

Stab_Der.CNTb = CnT_beta;
end