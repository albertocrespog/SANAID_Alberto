function [Trim_ITER_LAT] = Calculo_Trim_ITER_LAT_accelerations(conv_UNITS,conditions_TRIM_lat,...
    Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,Storing_STABILITY_DATA_1,conditions,case_AC,OUTPUT_read_XLSX)

Geo_tier = Storing_GEO_DATA.Geo_tier;
Performance = Storing_AERO_DATA.Performance;
Aero = Storing_AERO_DATA.Aero;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;
Stab_Der = Storing_STABILITY_DATA_1.Stab_Der;

% Study conditions
beta = conditions_TRIM_lat.beta;
beta_vec = conditions_TRIM_lat.beta_vec;

% Inertias
Ixx = Weight_tier.Ixx;
Iyy = Weight_tier.Iyy;
Izz = Weight_tier.Izz;
Ixz = Weight_tier.Ixz;

% Geometric data
S_ref = Geo_tier.S_ref;
S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
cmac_w1 = Geo_tier.cmac_w1;
% m_TOW = Weight_tier.m_TOW;

% Constants
% Constants
g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

Cydeltaa = Stab_Der.Cydeltaa;
Cldeltaa = Stab_Der.Cldeltaa;
Cndeltaa = Stab_Der.Cndeltaa;

Cydeltar = Stab_Der.Cydeltar;
Cldeltar = Stab_Der.Cldeltar;
Cndeltar = Stab_Der.Cndeltar;

I3 = Izz/(Ixx*Izz - Ixz^2);
I7 = (Ixz)/(Ixx*Izz-Ixz*Ixz);
I8 = (Ixx)/(Ixx*Izz-Ixz*Ixz); 

V_min = 22;
V_max = 52;
V_vec = linspace(V_min,V_max,50);

deltaa_min = 0*D2R;
deltaa_max = 25*D2R;
deltaa_vec = linspace(deltaa_min,deltaa_max,50);

deltar_min = 0*D2R;
deltar_max = 25*D2R;
deltar_vec = linspace(deltar_min,deltar_max,50);

rho = conditions_TRIM_lat.rho;

% Variable Beta study
for i=1:length(V_vec)
    V = V_vec(i);
    for j=1:length(deltaa_vec)
        q_inf(i) = 0.5*rho*V^2;
        pdot_vec(i,j) = I3*(Cldeltaa*deltaa_vec(j))*(q_inf(i)*S_w1*b_w1);
    end
    for k=1:length(deltar_vec)
        q_inf(k) = 0.5*rho*V^2;
        deltaa_const = 0;
        rdot_vec(i,k) = I7*(Cldeltaa*deltaa_const + Cldeltar*deltar_vec(k))*(q_inf(i)*S_w1*b_w1) +...
            I8*(Cndeltaa*deltaa_const +  + Cndeltar*deltar_vec(k))*(q_inf(i)*S_w1*b_w1);    
    end
end
Trim_ITER_LAT.V_vec = V_vec;
Trim_ITER_LAT.deltaa_vec = deltaa_vec;
Trim_ITER_LAT.deltar_vec = deltar_vec;
Trim_ITER_LAT.pdot_vec = pdot_vec;
Trim_ITER_LAT.rdot_vec = rdot_vec;
