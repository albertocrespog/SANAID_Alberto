function [Trim_ITER_LAT_Viraje,Fig] = Calculo_Trim_ITER_LAT_Viraje_v3(conv_UNITS,conditions_TRIM_turning,...
    Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,Storing_STABILITY_DATA_1,conditions_TRIM_lat)

Geo_tier = Storing_GEO_DATA.Geo_tier;
Performance = Storing_AERO_DATA.Performance;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;
Stab_Der = Storing_STABILITY_DATA_1.Stab_Der;

Ixx = Weight_tier.Ixx;
Iyy = Weight_tier.Iyy;
Izz = Weight_tier.Izz;
Ixz = Weight_tier.Ixz;

% Study conditions
phi = conditions_TRIM_turning.phi;
phi_vec = conditions_TRIM_turning.phi_vec;
n_viraje = conditions_TRIM_turning.n_viraje;

% Conversion units and constans
g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
ft2m = conv_UNITS.ft2m;
m2ft = conv_UNITS.m2ft;
W2hp = conv_UNITS.W2hp;
mps2ftps = conv_UNITS.mps2ftps;
kg2lb = conv_UNITS.kg2lb;
lb2kg = conv_UNITS.lb2kg;
ftpm2mps = conv_UNITS.ftpm2mps;
in2m = conv_UNITS.in2m;
W2pftsec = conv_UNITS.W2pftsec;
m22ft2 = conv_UNITS.m22ft2;
rho_SI2rho_IMP = conv_UNITS.rho_SI2rho_IMP;
qmet2qimp = conv_UNITS.qmet2qimp;
N2lbf = conv_UNITS.N2lbf;

% Performance
V = Performance.V;
rho = conditions_TRIM_turning.rho;

% Geometric data
S_ref = Geo_tier.S_ref;
S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
cmac_w1 = Geo_tier.cmac_w1;
m_TOW = Weight_tier.m_TOW;

% Constants
g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
q_inf = Performance.q_inf;
         
S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;

T = Stab_Der.T;
CD = Stab_Der.CD;
D = q_inf*S_w1*CD;
gamma = asin((T-D)/(m_TOW*g));
gamma_deg = gamma*R2D;
Trim_ITER_LAT_Viraje.gamma_deg = gamma_deg;

Cyb = Stab_Der.Cyb;
Clb = Stab_Der.Clb;
Cnb = Stab_Der.Cnb;

Cyr = Stab_Der.Cyr;
Clr = Stab_Der.Clr;
Cnr = Stab_Der.Cnr;

Cydeltaa = Stab_Der.Cydeltaa;
Cldeltaa = Stab_Der.Cldeltaa;
Cndeltaa = Stab_Der.Cndeltaa;

Cydeltar = Stab_Der.Cydeltar;
Cldeltar = Stab_Der.Cldeltar;
Cndeltar = Stab_Der.Cndeltar;

%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL TRIMADO Virage
%%%%%%%%%%%%%%%%%%%%%Stab_DerERAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_viraje = conditions_TRIM_turning.n_viraje;
phi_viraje = acos(1/n_viraje);
phi_grados_viraje = phi_viraje*R2D;
n_viraje = 1/cos(phi_viraje);
psi_dot_viraje = g*tan(phi_viraje)/V;
R_t_viraje = V^2/(g*tan(phi_viraje));

Stab_Der.phi_viraje = phi_viraje;
Stab_Der.phi_grados_viraje = phi_grados_viraje;
Stab_Der.n_viraje = n_viraje;
Stab_Der.psi_dot_viraje = psi_dot_viraje;
Stab_Der.R_t_viraje = R_t_viraje;

b1 = -Cyr*b_w1*g*sin(phi_viraje)/(2*(V^2));
b2 = (Izz-Iyy)*(g^2)*((sin(phi_viraje))^3)/(q_inf*S_w1*b_w1*(V^2)*cos(phi_viraje)) - ...
    Clr*b_w1*g*sin(phi_viraje)/(2*(V^2));
b3 = Ixz*(g^2)*((sin(phi_viraje))^3)/(q_inf*S_w1*b_w1*(V^2)*cos(phi_viraje)) - ...
    Cnr*b_w1*g*sin(phi_viraje)/(2*(V^2));

X1_viraje   =   [b1;b2;b3];
XX1_viraje  = [Cyb Cydeltaa Cydeltar; Clb Cldeltaa Cldeltar; Cnb Cndeltaa Cndeltar];

trim=inv(XX1_viraje)*X1_viraje;         %el primer valor es alpha y el segundo delta

beta_viraje = trim(1,1);
deltaa_viraje = trim(2,1);
deltar_viraje = trim(3,1);

beta_deg_viraje = beta_viraje*R2D;
deltaa_deg_viraje = deltaa_viraje*R2D;
deltar_deg_viraje = deltar_viraje*R2D;

Trim_ITER_LAT_Viraje.beta_viraje = beta_viraje;
Trim_ITER_LAT_Viraje.deltaa_viraje = deltaa_viraje;
Trim_ITER_LAT_Viraje.deltar_viraje = deltar_viraje;

Trim_ITER_LAT_Viraje.beta_deg_viraje = beta_deg_viraje;
Trim_ITER_LAT_Viraje.deltaa_deg_viraje = deltaa_deg_viraje;
Trim_ITER_LAT_Viraje.deltar_deg_viraje = deltar_deg_viraje;

% Variable Beta study
for i=1:length(phi_vec)

    n_viraje_var(i) = 1/cos(phi_vec(i));
    psi_dot_viraje_var(i) = g*tan(phi_vec(i))/V;
    R_t_viraje_var(i) = V^2/(g*tan(phi_vec(i)));

    q_inf = 0.5*rho*V^2;
    b1 = -Cyr*b_w1*g*sin(phi_vec(i))/(2*(V^2));
    b2 = (Izz-Iyy)*(g^2)*((sin(phi_vec(i)))^3)/(q_inf*S_w1*b_w1*(V^2)*cos(phi_vec(i))) -...
        Clr*b_w1*g*sin(phi_vec(i))/(2*(V^2));
    b3 = Ixz*(g^2)*((sin(phi_vec(i)))^3)/(q_inf*S_w1*b_w1*(V^2)*cos(phi_vec(i))) -...
        Cnr*b_w1*g*sin(phi_vec(i))/(2*(V^2));
   
    X1_viraje   = [b1;b2;b3];
    XX1_viraje  = [Cyb Cydeltaa Cydeltar; Clb Cldeltaa Cldeltar; Cnb Cndeltaa Cndeltar];    
    
    trim=inv(XX1_viraje)*X1_viraje;         %el primer valor es alpha y el segundo delta
    trim=XX1_viraje\X1_viraje;         %el primer valor es alpha y el segundo delta
    
    beta_viraje_var(i) = trim(1,1);
    deltaa_viraje_var(i) = trim(2,1);
    deltar_viraje_var(i) = trim(3,1);
    
    beta_deg_viraje_var(i) = beta_viraje_var(i)*R2D;
    deltaa_deg_viraje_var(i) = deltaa_viraje_var(i)*R2D;
    deltar_deg_viraje_var(i) = deltar_viraje_var(i)*R2D;
end

Trim_ITER_LAT_Viraje.n_viraje_var = n_viraje_var;
Trim_ITER_LAT_Viraje.psi_dot_viraje_var = psi_dot_viraje_var;
Trim_ITER_LAT_Viraje.R_t_viraje_var = R_t_viraje_var;

Trim_ITER_LAT_Viraje.beta_viraje_var = beta_viraje_var;
Trim_ITER_LAT_Viraje.deltaa_viraje_var = deltaa_viraje_var;
Trim_ITER_LAT_Viraje.deltar_viraje_var = deltar_viraje_var;

Trim_ITER_LAT_Viraje.beta_deg_viraje_var = beta_deg_viraje_var;
Trim_ITER_LAT_Viraje.deltaa_deg_viraje_var = deltaa_deg_viraje_var;
Trim_ITER_LAT_Viraje.deltar_deg_viraje_var = deltar_deg_viraje_var;