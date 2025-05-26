function [Aero,DATA_PL] = Generate_Aero_Data(DATA_Ae,Design_criteria,Performance,Geo_tier,Weight_tier,conv_UNITS)

% Reference wing
S_w1 = Geo_tier.S_w1;
S_w2 = Geo_tier.S_w2;
S_ref = Geo_tier.S_w1;

% Reference wing
S_w1_s = Geo_tier.S_w1_s;
S_w2_s = Geo_tier.S_w2_s;
S_w1_e = Geo_tier.S_w1_e;
S_w2_e = Geo_tier.S_w2_e;

conv_w1 = S_w1_s/S_w1_e;
conv_w2 = S_w2_s/S_w2_e;

% incidence of the FW, RW and the overall configuration
i_w1 = Design_criteria.i_w1; % incidence of Front Wing 
i_w2 = Design_criteria.i_w2; % incidence of Rear Wing

% Flight Safe Margin to calculate aerodynamic properties
Flight_SF = Design_criteria.Flight_SF;% Stall Safe Margin Conditions

g = conv_UNITS.g;
h = Performance.h;
Mach = Performance.Mach;
Temp = Performance.Temp;
rho = Performance.rho;
a = Performance.a;

V = Performance.V;
V_max = Performance.V_max;

% Weights to estimate aerodynamic properties
m_TOW = Weight_tier.m_TOW;

% Select txt files associated with each aerodynamic study
% Make sure to check with case asssigned in read_aero_files_Aug2018.m
index_w1 = Design_criteria.index_w1;
index_w2 = Design_criteria.index_w2;

% Allows the user to select the min AoA used to determine Curve Lift Slope
alpha_selected_w1 = Design_criteria.alpha_selected_w1;
alpha_selected_w2 = Design_criteria.alpha_selected_w2;

Design_criteria.alpha_selected_w1 = alpha_selected_w1;
Design_criteria.alpha_selected_w2 = alpha_selected_w2;

% CLmax 
CL_max_w1_CR = max(DATA_Ae(index_w1).CL*conv_w1);
CL_max_w2_CR = max(DATA_Ae(index_w2).CL*conv_w2);

Aero.C_L_max_w1_CR = CL_max_w1_CR;
Aero.C_L_max_w2_CR = CL_max_w2_CR;

% Determination of Stall Speed (min allowable speed)
V_stall_CR = sqrt(2*m_TOW*g/(rho*CL_max_w1_CR*S_ref));
V_min_CR = V_stall_CR;

% Min and Max CL (Flight Conditions associated to Min and max speeds)
CL_min_V_max_CR = m_TOW*g/(0.5*rho*(V_max^2)*S_ref);
CL_max_V_min_CR = m_TOW*g/(0.5*rho*(V_min_CR^2)*S_ref);

% Defines limits on CL for all aerodynamic surfaces
CL_w1_limit_min = m_TOW*g/(0.5*rho*(V_max^2)*S_w1);
CL_w1_limit_max = m_TOW*g/(0.5*rho*(V_min_CR^2)*S_w1);
CL_w2_limit_min = m_TOW*g/(0.5*rho*(V_max^2)*S_w2);
CL_w2_limit_max = m_TOW*g/(0.5*rho*(V_min_CR^2)*S_w2);

Performance.V_min_CR = V_min_CR;
Performance.V_stall_CR = V_stall_CR;

Aero.CL_min_V_max_CR = CL_min_V_max_CR;
Aero.CL_max_V_min_CR = CL_max_V_min_CR;

Aero.CL_w1_limit_min = CL_w1_limit_min;
Aero.CL_w1_limit_max = CL_w1_limit_max;
Aero.CL_w2_limit_min = CL_w2_limit_min;
Aero.CL_w2_limit_max = CL_w2_limit_max;

% Alpha at max CL
alpha_max_w1_CR = interp1(DATA_Ae(index_w1).CL*conv_w1,DATA_Ae(index_w1).alpha,CL_max_w1_CR,'spline');
alpha_max_w2_CR = interp1(DATA_Ae(index_w2).CL*conv_w2,DATA_Ae(index_w2).alpha,CL_max_w2_CR,'spline');

Aero.alpha_max_w1_CR = alpha_max_w1_CR;
Aero.alpha_max_w2_CR = alpha_max_w2_CR;

% CL at 0 AoA
alpha_CL_0_w1_CR = interp1(DATA_Ae(index_w1).CL*conv_w1,DATA_Ae(index_w1).alpha,0,'spline');
alpha_CL_0_w2_CR = interp1(DATA_Ae(index_w2).CL*conv_w2,DATA_Ae(index_w2).alpha,0,'spline');

Aero.alpha_CL_0_w1_CR = alpha_CL_0_w1_CR;
Aero.alpha_CL_0_w2_CR = alpha_CL_0_w2_CR;

CL_0_w1_CR = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CL*conv_w1,0,'spline');
CL_0_w2_CR = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).CL*conv_w2,0,'spline');

Aero.CL_0_w1_CR = CL_0_w1_CR;
Aero.CL_0_w2_CR = CL_0_w2_CR;

% Lift coefficient during zero AoA (fixed i_w1, and i_w2)
% VLM results
CL_w1_CR_iw = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CL*conv_w1,i_w1,'spline');
CL_w2_CR_iw = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).CL*conv_w2,i_w2,'spline');

Aero.CL_w1_CR_iw = CL_w1_CR_iw;
Aero.CL_w2_CR_iw = CL_w2_CR_iw;

% Moment at alpha=0
CM_0_w1_CR = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).Cm*conv_w1,0,'spline');
CM_0_w2_CR = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).Cm*conv_w2,0,'spline');

Aero.CM_0_w1_CR = CM_0_w1_CR;
Aero.CM_0_w2_CR = CM_0_w2_CR;

% Alpha at max CL with Fligth Safe Margin - operative, where operative
% equals to V_ope = Flight_SF*V_stall
alpha_w1_CR_ope  = interp1(DATA_Ae(index_w1).CL*conv_w1,DATA_Ae(index_w1).alpha,CL_max_w1_CR/Flight_SF^2,'spline');
alpha_w2_CR_ope  = interp1(DATA_Ae(index_w2).CL*conv_w2,DATA_Ae(index_w2).alpha,CL_max_w2_CR/Flight_SF^2,'spline');

CL_w1_CR_ope = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CL*conv_w1,alpha_w1_CR_ope,'spline');
CL_w2_CR_ope = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).CL*conv_w2,alpha_w2_CR_ope,'spline');

Aero.alpha_w1_CR_ope = alpha_w1_CR_ope;
Aero.alpha_w2_CR_ope = alpha_w2_CR_ope;

Aero.CL_w1_CR_ope = CL_w1_CR_ope;
Aero.CL_w2_CR_ope = CL_w2_CR_ope;

% Claculates Curve Lift Slope According to the min AoA selected (alpha
% selected) and the max CL (CLmax/1.2^2)
CL_alpha_w1_CR = (CL_w1_CR_ope - DATA_Ae(index_w1).CL(DATA_Ae(index_w1).alpha==alpha_selected_w1)*conv_w1)/...
    (alpha_w1_CR_ope - DATA_Ae(index_w1).alpha(DATA_Ae(index_w1).alpha==alpha_selected_w1)*conv_w1)*180/pi;
CL_alpha_w2_CR = (CL_w2_CR_ope - DATA_Ae(index_w2).CL(DATA_Ae(index_w2).alpha==alpha_selected_w2)*conv_w2)/...
    (alpha_w2_CR_ope - DATA_Ae(index_w2).alpha(DATA_Ae(index_w2).alpha==alpha_selected_w2)*conv_w2)*180/pi;

Aero.CL_alpha_w1_CR = CL_alpha_w1_CR;
Aero.CL_alpha_w2_CR = CL_alpha_w2_CR;

% Polar DATA
CD_w1_CR_ope = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CD*conv_w1,alpha_w1_CR_ope,'spline');
CD_w2_CR_ope = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).CD*conv_w2,alpha_w2_CR_ope,'spline');

% Efficiency at operative condition
E_w1_CR_ope = CL_w1_CR_ope/CD_w1_CR_ope;
E_w2_CR_ope = CL_w2_CR_ope/CD_w2_CR_ope;

Aero.CD_w1_CR_ope = CD_w1_CR_ope;
Aero.CD_w2_CR_ope = CD_w2_CR_ope;

Aero.E_w1_CR_ope = E_w1_CR_ope;
Aero.E_w1_CR_ope = E_w2_CR_ope;

% Efficiency vectors
% Max Range Jet
E_w1_max_R_jet_vec = (DATA_Ae(index_w1).CL*conv_w1)./DATA_Ae(index_w1).CD*conv_w1.^(3/2);
E_w2_max_R_jet_vec = (DATA_Ae(index_w2).CL*conv_w2)./DATA_Ae(index_w2).CD*conv_w2.^(3/2);

E_w1_max_R_jet = max(E_w1_max_R_jet_vec);
E_w2_max_R_jet = max(E_w2_max_R_jet_vec);

Aero.E_w1_max_R_jet_vec = E_w1_max_R_jet_vec;
Aero.E_w2_max_R_jet_vec = E_w2_max_R_jet_vec;

Aero.E_w1_max_R_jet = E_w1_max_R_jet;
Aero.E_w2_max_R_jet = E_w2_max_R_jet;

% Max Range Prop
E_w1_max_R_prop_vec = (DATA_Ae(index_w1).CL*conv_w1)./DATA_Ae(index_w1).CD*conv_w1;
E_w2_max_R_prop_vec = (DATA_Ae(index_w2).CL*conv_w2)./DATA_Ae(index_w2).CD*conv_w2;

E_w1_max_R_prop = max(E_w1_max_R_prop_vec);
E_w2_max_R_prop = max(E_w2_max_R_prop_vec);

Aero.E_w1_max_R_prop_vec = E_w1_max_R_prop_vec;
Aero.E_w2_max_R_prop_vec = E_w2_max_R_prop_vec;

Aero.E_w1_max_R_prop = E_w1_max_R_prop;
Aero.E_w2_max_R_prop = E_w2_max_R_prop;

% Max Endurance Jet
E_w1_max_E_jet_vec = (DATA_Ae(index_w1).CL*conv_w1)./DATA_Ae(index_w1).CD*conv_w1;
E_w2_max_E_jet_vec = (DATA_Ae(index_w2).CL*conv_w2)./DATA_Ae(index_w2).CD*conv_w2;

E_w1_max_E_jet = max(E_w1_max_E_jet_vec);
E_w2_max_E_jet = max(E_w2_max_E_jet_vec);

Aero.E_w1_max_E_jet_vec = E_w1_max_E_jet_vec;
Aero.E_w2_max_E_jet_vec = E_w2_max_E_jet_vec;

Aero.E_w1_max_E_jet = E_w1_max_E_jet;
Aero.E_w2_max_E_jet = E_w2_max_E_jet;

% Max Endurance Prop
E_w1_max_E_prop_vec = sqrt(3)*E_w1_max_R_prop;
E_w2_max_E_prop_vec = sqrt(3)*E_w2_max_R_prop;

E_w1_max_E_prop = max(E_w1_max_E_prop_vec);
E_w2_max_E_prop = max(E_w2_max_E_prop_vec);

Aero.E_w1_max_E_prop_vec = E_w1_max_E_prop_vec;
Aero.E_w2_max_E_prop_vec = E_w2_max_E_prop_vec;

Aero.E_w1_max_E_prop = E_w1_max_E_prop;
Aero.E_w2_max_E_prop = E_w2_max_E_prop;

CL_ref_limit_min = CL_min_V_max_CR;
CL_ref_limit_max = CL_max_V_min_CR;

% Polyfit Polar
PL_w1 = polyfit(DATA_Ae(index_w1).CL(DATA_Ae(index_w1).CL>= CL_w1_limit_min & DATA_Ae(index_w1).CL<=CL_w1_limit_max)*conv_w1,...
    DATA_Ae(index_w1).CD(DATA_Ae(index_w1).CL>=CL_w1_limit_min & DATA_Ae(index_w1).CL<= CL_w1_limit_max)*conv_w1,2);
PL_w2 = polyfit(DATA_Ae(index_w2).CL(DATA_Ae(index_w2).CL>= CL_w2_limit_min & DATA_Ae(index_w2).CL<=CL_w2_limit_max)*conv_w2,...
    DATA_Ae(index_w2).CD(DATA_Ae(index_w2).CL>=CL_w2_limit_min & DATA_Ae(index_w2).CL<= CL_w2_limit_max)*conv_w2,2);

Aero.PL_w1 = PL_w1;
Aero.PL_w2 = PL_w2;

CD0_w1 = PL_w1(3);
CD1_w1 = PL_w1(2);
CD2_w1 = PL_w1(1);

Aero.CD0_w1 = CD0_w1;
Aero.CD1_w1 = CD1_w1;
Aero.CD2_w1 = CD2_w1;

CD0_w2 = PL_w2(3);
CD1_w2 = PL_w2(2);
CD2_w2 = PL_w2(1);

Aero.CD0_w2 = CD0_w2;
Aero.CD1_w2 = CD1_w2;
Aero.CD2_w2 = CD2_w2;

% Generates Vector with polyfit approximation
Aero.CD_poly_w1 = PL_w1(3) + PL_w1(2).*DATA_Ae(index_w1).CL*conv_w1 + PL_w1(1).*DATA_Ae(index_w1).CL*conv_w1.^2;
Aero.CD_poly_w2 = PL_w2(3) + PL_w2(2).*DATA_Ae(index_w2).CL*conv_w2 + PL_w2(1).*DATA_Ae(index_w2).CL*conv_w2.^2;

DATA_PL(1).CD_poly = Aero.CD_poly_w1;
DATA_PL(2).CD_poly = Aero.CD_poly_w2;