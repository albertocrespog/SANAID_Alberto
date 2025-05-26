function [Aero,DATA_PL,Performance] = Generate_Aero_Data_2021_v1(DATA_Ae,Design_criteria,Performance,Geo_tier,...
    Weight_tier,conv_UNITS,AC_CONFIGURATION,Body_Geo,OUTPUT_read_XLSX)

%% identifies the aerodynamic surfaces being used
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;

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

Performance.h = h;
Performance.Mach = Mach;
Performance.Temp = Temp;
Performance.rho = rho;
Performance.a = a;
Performance.V = V;
Performance.V_max = V_max;

q_inf = 0.5*rho*V^2;
Performance.q_inf = q_inf;

% Weights to estimate aerodynamic properties
m_TOW = Weight_tier.m_TOW;

S_ref = Geo_tier.S_ref;

if W1 == 1
    % Reference wing
    S_w1 = Geo_tier.S_w1;
    S_ref = Geo_tier.S_w1;
    
    % Reference wing
    S_w1_s = Geo_tier.S_w1_s;
    S_w1_e = Geo_tier.S_w1_e;
    
  
    % Assume no surface conversion: aerodynamic data is as original from
    % Aerodynamic Software
    conv_w1 = 1;
    
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_w1 = Design_criteria.index_w1;
    
    % incidence of the FW, RW and the overall configuration
    i_w1 = Design_criteria.i_w1; % incidence of Front Wing
    
    % Allows the user to select the min AoA used to determine Curve Lift Slope
    alpha_selected_w1 = Design_criteria.alpha_selected_w1;
    
    % CLmax
    [CL_max_w1_CR,max_index] = max(DATA_Ae(index_w1).CL*conv_w1);
    Aero.CL_max_w1_CR = CL_max_w1_CR;
    alpha_max_w1_CR = DATA_Ae(index_w1).alpha(max_index);
    Aero.alpha_max_w1_CR = alpha_max_w1_CR;
    
    % generates the CL and alpha such that avoids the problem of "Sample
    % points must be unique and sorted in ascending order." that it is
    % encountered
    CL_new = DATA_Ae(index_w1).CL(1:max_index)*conv_w1;
    alpha_new = DATA_Ae(index_w1).alpha(1:max_index);
    
    %% Determination of Stall conditions
    CL_max_w1_ope = CL_max_w1_CR/(Flight_SF^2);
    V_stall_w1 = sqrt((2*m_TOW*g)/(rho*S_ref*CL_max_w1_CR));
    V_min_w1 = Flight_SF*V_stall_w1;
    %     alpha_max_w1_ope = interp1(DATA_Ae(index_w1).CL,DATA_Ae(index_w1).alpha,CL_max_w1_ope,'spline');
    alpha_max_w1_ope = interp1(CL_new,alpha_new,CL_max_w1_ope,'spline');
    
    Performance.V_stall_w1 = V_stall_w1;
    Performance.V_min_w1 = V_min_w1;
    Performance.CL_max_w1_ope = CL_max_w1_ope;
    Performance.alpha_max_w1_ope = alpha_max_w1_ope;
    Performance.alpha_max = alpha_max_w1_CR;
    Performance.CL_max_w1 = CL_max_w1_CR;
    
    q_inf = 0.5*rho*V^2;
    Performance.q_inf = q_inf;
    
    % Min and Max CL (Flight Conditions associated to Min and max speeds)
    CL_min_V_max_CR = m_TOW*g/(0.5*rho*(V_max^2)*S_ref);
    CL_max_V_min_CR = m_TOW*g/(0.5*rho*(V_min_w1^2)*S_ref);
    
    % Defines limits on CL for all aerodynamic surfaces
    CL_w1_limit_min = m_TOW*g/(0.5*rho*(V_max^2)*S_w1);
    CL_w1_limit_max = m_TOW*g/(0.5*rho*(V_min_w1^2)*S_w1);
    
    Aero.CL_min_V_max_CR = CL_min_V_max_CR;
    Aero.CL_max_V_min_CR = CL_max_V_min_CR;
    
    Aero.CL_w1_limit_min = CL_w1_limit_min;
    Aero.CL_w1_limit_max = CL_w1_limit_max;
    
    % CL at 0 AoA
    %     alpha_CL_0_w1_CR = interp1(DATA_Ae(index_w1).CL*conv_w1,DATA_Ae(index_w1).alpha,0,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    alpha_CL_0_w1_CR = interp1(CL_new,alpha_new,0,'spline');
    Aero.alpha_CL_0_w1_CR = alpha_CL_0_w1_CR;
    % CL_0_w1_CR = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CL*conv_w1,0,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    CL_0_w1_CR = interp1(alpha_new,CL_new,0,'spline');
    Aero.CL_0_w1_CR = CL_0_w1_CR;
    
    % Lift coefficient during zero AoA (fixed i_w1, and i_w2)
    % VLM results
    % CL_w1_CR_iw = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CL*conv_w1,i_w1,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    CL_w1_CR_iw = interp1(alpha_new,CL_new,i_w1,'spline');
    
    Aero.CL_w1_CR_iw = CL_w1_CR_iw;
    
    % Moment at alpha=0
    CM_0_w1_CR = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).Cm*conv_w1,0,'spline');
    Aero.CM_0_w1_CR = CM_0_w1_CR;
    
    % Alpha at max CL with Fligth Safe Margin - operative, where operative
    % equals to V_ope = Flight_SF*V_stall
    %     alpha_w1_CR_ope  = interp1(DATA_Ae(index_w1).CL*conv_w1,DATA_Ae(index_w1).alpha,CL_max_w1_CR/Flight_SF^2,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    alpha_w1_CR_ope = interp1(CL_new,alpha_new,CL_max_w1_CR/Flight_SF^2,'spline');
    
    %     CL_w1_CR_ope = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CL*conv_w1,alpha_w1_CR_ope,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    CL_w1_CR_ope = interp1(alpha_new,CL_new,alpha_w1_CR_ope,'spline');
    Aero.alpha_w1_CR_ope = alpha_w1_CR_ope;
    Aero.CL_w1_CR_ope = CL_w1_CR_ope;
    
    % Claculates Curve Lift Slope According to the min AoA selected (alpha
    % selected) and the max CL (CLmax/1.2^2)
    CL_alpha_w1_CR = (CL_w1_CR_ope - DATA_Ae(index_w1).CL(DATA_Ae(index_w1).alpha==alpha_selected_w1)*conv_w1)/...
        (alpha_w1_CR_ope - DATA_Ae(index_w1).alpha(DATA_Ae(index_w1).alpha==alpha_selected_w1)*conv_w1)*180/pi;
    Aero.CL_alpha_w1_CR = CL_alpha_w1_CR;
    
    % Polar DATA
    CD_w1_CR_ope = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CD*conv_w1,alpha_w1_CR_ope,'spline');
    
    % Efficiency at operative condition
    E_w1_CR_ope = CL_w1_CR_ope/CD_w1_CR_ope;
    Aero.CD_w1_CR_ope = CD_w1_CR_ope;
    
    Aero.E_w1_CR_ope = E_w1_CR_ope;
    
    % Efficiency vectors
    % Max Range Jet
    E_w1_max_R_jet_vec = (DATA_Ae(index_w1).CL*conv_w1)./DATA_Ae(index_w1).CD*conv_w1.^(3/2);
    E_w1_max_R_jet = max(E_w1_max_R_jet_vec);
    Aero.E_w1_max_R_jet_vec = E_w1_max_R_jet_vec;
    Aero.E_w1_max_R_jet = E_w1_max_R_jet;
    % Max Range Prop
    E_w1_max_R_prop_vec = (DATA_Ae(index_w1).CL*conv_w1)./DATA_Ae(index_w1).CD*conv_w1;
    E_w1_max_R_prop = max(E_w1_max_R_prop_vec);
    Aero.E_w1_max_R_prop_vec = E_w1_max_R_prop_vec;
    % Max Endurance Jet
    E_w1_max_E_jet_vec = (DATA_Ae(index_w1).CL*conv_w1)./DATA_Ae(index_w1).CD*conv_w1;
    E_w1_max_E_jet = max(E_w1_max_E_jet_vec);
    Aero.E_w1_max_E_jet_vec = E_w1_max_E_jet_vec;
    Aero.E_w1_max_E_jet = E_w1_max_E_jet;
    % Max Endurance Prop
    E_w1_max_E_prop_vec = sqrt(3)*E_w1_max_R_prop;
    E_w1_max_E_prop = max(E_w1_max_E_prop_vec);
    Aero.E_w1_max_E_prop_vec = E_w1_max_E_prop_vec;
    Aero.E_w1_max_E_prop = E_w1_max_E_prop;
    CL_ref_limit_min = CL_min_V_max_CR;
    CL_ref_limit_max = CL_max_V_min_CR;
    % Polyfit Polar
    PL_w1 = polyfit(DATA_Ae(index_w1).CL(DATA_Ae(index_w1).CL>= CL_w1_limit_min & DATA_Ae(index_w1).CL<=CL_w1_limit_max)*conv_w1,...
        DATA_Ae(index_w1).CD(DATA_Ae(index_w1).CL>=CL_w1_limit_min & DATA_Ae(index_w1).CL<= CL_w1_limit_max)*conv_w1,2);
    Aero.PL_w1 = PL_w1;
    CD0_w1 = PL_w1(3);
    CD1_w1 = PL_w1(2);
    CD2_w1 = PL_w1(1);
    Aero.CD0_w1 = CD0_w1;
    Aero.CD1_w1 = CD1_w1;
    Aero.CD2_w1 = CD2_w1;
    
    % Generates Vector with polyfit approximation
    Aero.CD_poly_w1 = CD0_w1 + CD1_w1.*DATA_Ae(index_w1).CL*conv_w1 + CD2_w1.*(conv_w1*DATA_Ae(index_w1).CL.^2);
    DATA_PL.CD_poly = Aero.CD_poly_w1;
else
    Aero.CD0_w1 = 0;
    Aero.CD1_w1 = 0;
    Aero.CD2_w1 = 0;
    S_w1 = 0;
end

if Can == 1
    % Reference wing
    S_can = Geo_tier.S_can;
    
    % Reference wing
    S_can_s = Geo_tier.S_can_s;
    S_can_e = Geo_tier.S_can_e;
    
    conv_can = S_can_s/S_w1_e;
    % Assume no surface conversion: aerodynamic data is as original from
    % Aerodynamic Software
    conv_can = 1;
    
    
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_w3 = Design_criteria.index_w3;
    
    % incidence of the FW, RW and the overall configuration
    i_w3 = Design_criteria.i_w3; % incidence of Front Wing
    
    % Allows the user to select the min AoA used to determine Curve Lift Slope
    alpha_selected_w3 = Design_criteria.alpha_selected_w3;
    
    % CLmax
    [CL_max_can_CR,max_index] = max(DATA_Ae(index_w3).CL*conv_can);
    Aero.CL_max_can_CR = CL_max_can_CR;
    alpha_max_can_CR = DATA_Ae(index_w3).alpha(max_index);
    Aero.alpha_max_can_CR = alpha_max_can_CR;
    
    % generates the CL and alpha such that avoids the problem of "Sample
    % points must be unique and sorted in ascending order." that it is
    % encountered
    CL_new = DATA_Ae(index_w3).CL(1:max_index)*conv_can;
    alpha_new = DATA_Ae(index_w3).alpha(1:max_index);
    
    %% Determination of Stall conditions
    CL_max_can_ope = CL_max_can_CR/(Flight_SF^2);
    V_stall_can = sqrt((2*m_TOW*g)/(rho*S_ref*CL_max_can_CR));
    V_min_can = Flight_SF*V_stall_can;
    %     alpha_max_can_ope = interp1(DATA_Ae(index_w3).CL,DATA_Ae(index_w3).alpha,CL_max_can_ope,'spline');
    alpha_max_can_ope = interp1(CL_new,alpha_new,CL_max_can_ope,'spline');
    
    Performance.V_stall_can = V_stall_can;
    Performance.V_min_can = V_min_can;
    Performance.CL_max_can_ope = CL_max_can_ope;
    Performance.alpha_max_can_ope = alpha_max_can_ope;
    Performance.alpha_max = alpha_max_can_CR;
    Performance.CL_max_can = CL_max_can_CR;
    
    q_inf = 0.5*rho*V^2;
    Performance.q_inf = q_inf;
    
    % Min and Max CL (Flight Conditions associated to Min and max speeds)
    CL_min_V_max_CR = m_TOW*g/(0.5*rho*(V_max^2)*S_ref);
    CL_max_V_min_CR = m_TOW*g/(0.5*rho*(V_min_can^2)*S_ref);
    
    % Defines limits on CL for all aerodynamic surfaces
    CL_can_limit_min = m_TOW*g/(0.5*rho*(V_max^2)*S_can);
    CL_can_limit_max = m_TOW*g/(0.5*rho*(V_min_can^2)*S_can);
    
    Aero.CL_min_V_max_CR = CL_min_V_max_CR;
    Aero.CL_max_V_min_CR = CL_max_V_min_CR;
    
    Aero.CL_can_limit_min = CL_can_limit_min;
    Aero.CL_can_limit_max = CL_can_limit_max;
    
    % CL at 0 AoA
    %     alpha_CL_0_can_CR = interp1(DATA_Ae(index_w3).CL*conv_can,DATA_Ae(index_w3).alpha,0,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    alpha_CL_0_can_CR = interp1(CL_new,alpha_new,0,'spline');
    Aero.alpha_CL_0_can_CR = alpha_CL_0_can_CR;
    %     CL_0_can_CR = interp1(DATA_Ae(index_w3).alpha,DATA_Ae(index_w3).CL*conv_can,0,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    CL_0_can_CR = interp1(alpha_new,CL_new,0,'spline');
    Aero.CL_0_can_CR = CL_0_can_CR;
    
    % Lift coefficient during zero AoA (fixed i_w3, and i_w2)
    % VLM results
    %     CL_can_CR_iw = interp1(DATA_Ae(index_w3).alpha,DATA_Ae(index_w3).CL*conv_can,i_w3,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    CL_can_CR_iw = interp1(alpha_new,CL_new,i_w3,'spline');
    
    Aero.CL_can_CR_iw = CL_can_CR_iw;
    
    % Moment at alpha=0
    CM_0_can_CR = interp1(DATA_Ae(index_w3).alpha,DATA_Ae(index_w3).Cm*conv_can,0,'spline');
    Aero.CM_0_can_CR = CM_0_can_CR;
    
    % Alpha at max CL with Fligth Safe Margin - operative, where operative
    % equals to V_ope = Flight_SF*V_stall
    %     alpha_can_CR_ope  = interp1(DATA_Ae(index_w3).CL*conv_can,DATA_Ae(index_w3).alpha,CL_max_can_CR/Flight_SF^2,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    alpha_can_CR_ope = interp1(CL_new,alpha_new,CL_max_can_CR/Flight_SF^2,'spline');
    
    %     CL_can_CR_ope = interp1(DATA_Ae(index_w3).alpha,DATA_Ae(index_w3).CL*conv_can,alpha_can_CR_ope,'spline');
    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    CL_can_CR_ope = interp1(alpha_new,CL_new,alpha_can_CR_ope,'spline');
    Aero.alpha_can_CR_ope = alpha_can_CR_ope;
    Aero.CL_can_CR_ope = CL_can_CR_ope;
    
    % Claculates Curve Lift Slope According to the min AoA selected (alpha
    % selected) and the max CL (CLmax/1.2^2)
    CL_alpha_can_CR = (CL_can_CR_ope - DATA_Ae(index_w3).CL(DATA_Ae(index_w3).alpha==alpha_selected_w3)*conv_can)/...
        (alpha_can_CR_ope - DATA_Ae(index_w3).alpha(DATA_Ae(index_w3).alpha==alpha_selected_w3)*conv_can)*180/pi;
    Aero.CL_alpha_can_CR = CL_alpha_can_CR;
    
    % Polar DATA
    CD_can_CR_ope = interp1(DATA_Ae(index_w3).alpha,DATA_Ae(index_w3).CD*conv_can,alpha_can_CR_ope,'spline');
    
    % Efficiency at operative condition
    E_can_CR_ope = CL_can_CR_ope/CD_can_CR_ope;
    Aero.CD_can_CR_ope = CD_can_CR_ope;
    
    Aero.E_can_CR_ope = E_can_CR_ope;
    
    % Efficiency vectors
    % Max Range Jet
    E_can_max_R_jet_vec = (DATA_Ae(index_w3).CL*conv_can)./DATA_Ae(index_w3).CD*conv_can.^(3/2);
    E_can_max_R_jet = max(E_can_max_R_jet_vec);
    Aero.E_can_max_R_jet_vec = E_can_max_R_jet_vec;
    Aero.E_can_max_R_jet = E_can_max_R_jet;
    % Max Range Prop
    E_can_max_R_prop_vec = (DATA_Ae(index_w3).CL*conv_can)./DATA_Ae(index_w3).CD*conv_can;
    E_can_max_R_prop = max(E_can_max_R_prop_vec);
    Aero.E_can_max_R_prop_vec = E_can_max_R_prop_vec;
    % Max Endurance Jet
    E_can_max_E_jet_vec = (DATA_Ae(index_w3).CL*conv_can)./DATA_Ae(index_w3).CD*conv_can;
    E_can_max_E_jet = max(E_can_max_E_jet_vec);
    Aero.E_can_max_E_jet_vec = E_can_max_E_jet_vec;
    Aero.E_can_max_E_jet = E_can_max_E_jet;
    % Max Endurance Prop
    E_can_max_E_prop_vec = sqrt(3)*E_can_max_R_prop;
    E_can_max_E_prop = max(E_can_max_E_prop_vec);
    Aero.E_can_max_E_prop_vec = E_can_max_E_prop_vec;
    Aero.E_can_max_E_prop = E_can_max_E_prop;
    CL_ref_limit_min = CL_min_V_max_CR;
    CL_ref_limit_max = CL_max_V_min_CR;
    % Polyfit Polar
    PL_can = polyfit(DATA_Ae(index_w3).CL(DATA_Ae(index_w3).CL>= CL_can_limit_min & DATA_Ae(index_w3).CL<=CL_can_limit_max)*conv_can,...
        DATA_Ae(index_w3).CD(DATA_Ae(index_w3).CL>=CL_can_limit_min & DATA_Ae(index_w3).CL<= CL_can_limit_max)*conv_can,2);
    Aero.PL_can = PL_can;
    CD0_can = PL_can(3);
    CD1_can = PL_can(2);
    CD2_can = PL_can(1);
    Aero.CD0_can = CD0_can;
    Aero.CD1_can = CD1_can;
    Aero.CD2_can = CD2_can;
    
    % Generates Vector with polyfit approximation
    Aero.CD_poly_can = CD0_can + CD1_can.*DATA_Ae(index_w3).CL*conv_can + CD2_can.*(conv_can*DATA_Ae(index_w3).CL.^2);
    DATA_PL.CD_poly = Aero.CD_poly_can;
else
    Aero.CD0_can = 0;
    Aero.CD1_can = 0;
    Aero.CD2_can = 0;
    S_can = 0;
end

if HTP == 1
    % Reference wing
    S_w2 = Geo_tier.S_w2;
    % Reference wing
    S_w2_s = Geo_tier.S_w2_s;
    S_w2_e = Geo_tier.S_w2_e;
    S_w1_e = Geo_tier.S_w1_e;
    
    conv_w2 = S_w2_s/S_w1_e;
    % Assume no surface conversion: aerodynamic data is as original from
    % Aerodynamic Software
    conv_w2 = 1;
    
    % incidence of the FW, RW and the overall configuration
    i_w2 = Design_criteria.i_w2; % incidence of Rear Wing
    
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_w2 = Design_criteria.index_w2;
    % incidence of the FW, RW and the overall configuration
    i_w2 = Design_criteria.i_w2; % incidence of Rear Wing
    
    % Allows the user to select the min AoA used to determine Curve Lift Slope
    alpha_selected_w2 = Design_criteria.alpha_selected_w2;
    % CLmax
    CL_max_w2_CR = max(DATA_Ae(index_w2).CL*conv_w2);
    Aero.CL_max_w2_CR = CL_max_w2_CR;
    % Alpha at max CL
    alpha_max_w2_CR = interp1(DATA_Ae(index_w2).CL*conv_w2,DATA_Ae(index_w2).alpha,CL_max_w2_CR,'spline');
    Aero.alpha_max_w2_CR = alpha_max_w2_CR;
    
    %% Determination of Stall conditions
    CL_max_w2_ope = CL_max_w2_CR/(Flight_SF^2);
    alpha_max_w2_ope = interp1(DATA_Ae(index_w2).CL,DATA_Ae(index_w2).alpha,CL_max_w2_ope,'spline');
    % Defines limits on CL for all aerodynamic surfaces
    % CL at 0 AoA
    alpha_CL_0_w2_CR = interp1(DATA_Ae(index_w2).CL*conv_w2,DATA_Ae(index_w2).alpha,0,'spline');
    Aero.alpha_CL_0_w2_CR = alpha_CL_0_w2_CR;
    CL_0_w2_CR = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).CL*conv_w2,0,'spline');
    Aero.CL_0_w2_CR = CL_0_w2_CR;
    
    % Lift coefficient during zero AoA (fixed i_w1, and i_w2)
    % VLM results
    CL_w2_CR_iw = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).CL*conv_w2,i_w2,'spline');
    Aero.CL_w2_CR_iw = CL_w2_CR_iw;
    % Moment at alpha=0
    CM_0_w2_CR = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).Cm*conv_w2,0,'spline');
    Aero.CM_0_w2_CR = CM_0_w2_CR;
    
    % Alpha at max CL with Fligth Safe Margin - operative, where operative
    % equals to V_ope = Flight_SF*V_stall
    alpha_w2_CR_ope  = interp1(DATA_Ae(index_w2).CL*conv_w2,DATA_Ae(index_w2).alpha,CL_max_w2_CR/Flight_SF^2,'spline');
    CL_w2_CR_ope = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).CL*conv_w2,alpha_w2_CR_ope,'spline');
    Aero.alpha_w2_CR_ope = alpha_w2_CR_ope;
    Aero.CL_w2_CR_ope = CL_w2_CR_ope;
    
    % Claculates Curve Lift Slope According to the min AoA selected (alpha
    % selected) and the max CL (CLmax/1.2^2)
    CL_alpha_w2_CR = (CL_w2_CR_ope - DATA_Ae(index_w2).CL(DATA_Ae(index_w2).alpha==alpha_selected_w2)*conv_w2)/...
        (alpha_w2_CR_ope - DATA_Ae(index_w2).alpha(DATA_Ae(index_w2).alpha==alpha_selected_w2)*conv_w2)*180/pi;
    Aero.CL_alpha_w2_CR = CL_alpha_w2_CR;
    
    % Polar DATA
    CD_w2_CR_ope = interp1(DATA_Ae(index_w2).alpha,DATA_Ae(index_w2).CD*conv_w2,alpha_w2_CR_ope,'spline');
    % Efficiency at operative condition
    Aero.CD_w2_CR_ope = CD_w2_CR_ope;
    
    % Efficiency at operative condition
    E_w2_CR_ope = CL_w2_CR_ope/CD_w2_CR_ope;
    Aero.CD_w2_CR_ope = CD_w2_CR_ope;
    
    Aero.E_w2_CR_ope = E_w2_CR_ope;
    
    % Efficiency vectors
    % Max Range Jet
    E_w2_max_R_jet_vec = (DATA_Ae(index_w2).CL*conv_w2)./DATA_Ae(index_w2).CD*conv_w2.^(3/2);
    E_w2_max_R_jet = max(E_w2_max_R_jet_vec);
    Aero.E_w2_max_R_jet_vec = E_w2_max_R_jet_vec;
    Aero.E_w2_max_R_jet = E_w2_max_R_jet;
    % Max Range Prop
    E_w2_max_R_prop_vec = (DATA_Ae(index_w2).CL*conv_w2)./DATA_Ae(index_w2).CD*conv_w2;
    E_w2_max_R_prop = max(E_w2_max_R_prop_vec);
    Aero.E_w2_max_R_prop_vec = E_w2_max_R_prop_vec;
    % Max Endurance Jet
    E_w2_max_E_jet_vec = (DATA_Ae(index_w2).CL*conv_w2)./DATA_Ae(index_w2).CD*conv_w2;
    E_w2_max_E_jet = max(E_w2_max_E_jet_vec);
    Aero.E_w2_max_E_jet_vec = E_w2_max_E_jet_vec;
    Aero.E_w2_max_E_jet = E_w2_max_E_jet;
    % Max Endurance Prop
    E_w2_max_E_prop_vec = sqrt(3)*E_w2_max_R_prop;
    E_w2_max_E_prop = max(E_w2_max_E_prop_vec);
    Aero.E_w2_max_E_prop_vec = E_w2_max_E_prop_vec;
    Aero.E_w2_max_E_prop = E_w2_max_E_prop;
    CL_ref_limit_min = CL_min_V_max_CR;
    CL_ref_limit_max = CL_max_V_min_CR;
    % Polyfit Polar

    CL_w2_limit_max = max(DATA_Ae(index_w2).CL);
    CL_w2_limit_min = min(DATA_Ae(index_w2).CL);

    % PL_w2 = polyfit(DATA_Ae(index_w2).CL(DATA_Ae(index_w2).CL>= CL_w1_limit_min & DATA_Ae(index_w1).CL<=CL_w1_limit_max)*conv_w2,...
    %     DATA_Ae(index_w2).CD(DATA_Ae(index_w2).CL>=CL_w2_limit_min & DATA_Ae(index_w2).CL<= CL_w2_limit_max)*conv_w2,2);
     
    % PL_w2 = polyfit(DATA_Ae(index_w2).CL(DATA_Ae(index_w1).CL>= CL_w1_limit_min & DATA_Ae(index_w1).CL<=CL_w1_limit_max)*conv_w2,...
    %     DATA_Ae(index_w2).CD(DATA_Ae(index_w1).CL>=CL_w1_limit_min & DATA_Ae(index_w1).CL<= CL_w1_limit_max)*conv_w2,2);
    
    % PL_w2 = polyfit(DATA_Ae(index_w2).CL(DATA_Ae(index_w2).CL>= CL_w2_limit_min & DATA_Ae(index_w2).CL<=CL_w2_limit_max)*conv_w2,...
    %     DATA_Ae(index_w2).CD(DATA_Ae(index_w2).CL>=CL_w2_limit_min & DATA_Ae(index_w2).CL<= CL_w2_limit_max)*conv_w2,2);

     PL_w2 = polyfit(DATA_Ae(index_w2).CL(DATA_Ae(index_w2).CL>= CL_w2_limit_min & DATA_Ae(index_w2).CL<=CL_w2_limit_max)*conv_w2,...
        DATA_Ae(index_w2).CD(DATA_Ae(index_w2).CL>=CL_w2_limit_min & DATA_Ae(index_w2).CL<= CL_w2_limit_max)*conv_w2,2);
    
    Aero.PL_w2 = PL_w2;
    CD0_w2 = PL_w2(3);
    CD1_w2 = PL_w2(2);
    CD2_w2 = PL_w2(1);
    Aero.CD0_w2 = CD0_w2;
    Aero.CD1_w2 = CD1_w2;
    Aero.CD2_w2 = CD2_w2;
    % Generates Vector with polyfit approximation
    Aero.CD_poly_w2 = CD0_w2 + CD1_w2.*DATA_Ae(index_w2).CL*conv_w2 + CD2_w2.*(conv_w2*DATA_Ae(index_w2).CL.^2);
    
    DATA_PL.CD_poly = Aero.CD_poly_w1;
else
    Aero.CD0_w2 = 0;
    Aero.CD1_w2 = 0;
    Aero.CD2_w2 = 0;
    S_w2 = 0;
end

if Vee == 1
    % Reference wing
    S_vee = Geo_tier.S_vee;
    % Reference wing
    S_vee_s = Geo_tier.S_vee_s;
    S_vee_e = Geo_tier.S_vee_e;
    S_w1_e = Geo_tier.S_w1_e;
    
    % conv_w2 = S_w2_s/S_w1_e;
    % Assume no surface conversion: aerodynamic data is as original from
    % Aerodynamic Software no needed 
    conv_vee = 1;
    
    % incidence of the FW, RW and the overall configuration
    i_vee = Design_criteria.i_vee; % incidence of Rear Wing
    index_vee_cy = Design_criteria.index_vee_cy; %
    
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_vee = Design_criteria.index_vee;
    index_vee_cy = Design_criteria.index_vee_cy; %
%     index_vee_cy = 7; %

    % incidence of the FW, RW and the overall configuration
    i_vee = Design_criteria.i_vee; % incidence of Rear Wing
    
    % Allows the user to select the min AoA used to determine Curve Lift Slope
    alpha_selected_vee = Design_criteria.alpha_selected_vee;
    
    % CLmax
    CL_max_vee_CR = max(DATA_Ae(index_vee).CL*conv_vee);
    Aero.CL_max_vee_CR = CL_max_vee_CR;
    % Alpha at max CL
    alpha_max_vee_CR = interp1(DATA_Ae(index_vee).CL*conv_vee,DATA_Ae(index_vee).alpha,CL_max_vee_CR,'spline');
    Aero.alpha_max_vee_CR = alpha_max_vee_CR;
    
    %% Determination of Stall conditions
    CL_max_vee_ope = CL_max_vee_CR/(Flight_SF^2);
    V_stall_vee = sqrt((2*m_TOW*g)/(rho*S_ref*CL_max_vee_ope));
    V_min_vee = Flight_SF*V_stall_vee;
    
    alpha_max_vee_ope = interp1(DATA_Ae(index_vee).CL,DATA_Ae(index_vee).alpha,CL_max_vee_ope,'spline');
    % Defines limits on CL for all aerodynamic surfaces
    CL_vee_limit_min = m_TOW*g/(0.5*rho*(V_max^2)*S_vee);
    CL_vee_limit_max = m_TOW*g/(0.5*rho*(V_min_vee^2)*S_vee);
    Aero.CL_vee_limit_min = CL_vee_limit_min;
    Aero.CL_vee_limit_max = CL_vee_limit_max;
    % CL at 0 AoA
    alpha_CL_0_vee_CR = interp1(DATA_Ae(index_vee).CL*conv_vee,DATA_Ae(index_vee).alpha,0,'spline');
    Aero.alpha_CL_0_vee_CR = alpha_CL_0_vee_CR;
    CL_0_vee_CR = interp1(DATA_Ae(index_vee).alpha,DATA_Ae(index_vee).CL*conv_vee,0,'spline');
    Aero.CL_0_vee_CR = CL_0_vee_CR;
    
    % Lift coefficient during zero AoA (fixed i_w1, and i_vee)
    % VLM results
    CL_vee_CR_iw = interp1(DATA_Ae(index_vee).alpha,DATA_Ae(index_vee).CL*conv_vee,i_vee,'spline');
    Aero.CL_vee_CR_iw = CL_vee_CR_iw;
    % Moment at alpha=0
    CM_0_vee_CR = interp1(DATA_Ae(index_vee).alpha,DATA_Ae(index_vee).Cm*conv_vee,0,'spline');
    Aero.CM_0_vee_CR = CM_0_vee_CR;
    
    % Alpha at max CL with Fligth Safe Margin - operative, where operative
    % equals to V_ope = Flight_SF*V_stall
    alpha_vee_CR_ope  = interp1(DATA_Ae(index_vee).CL*conv_vee,DATA_Ae(index_vee).alpha,CL_max_vee_CR/Flight_SF^2,'spline');
    CL_vee_CR_ope = interp1(DATA_Ae(index_vee).alpha,DATA_Ae(index_vee).CL*conv_vee,alpha_vee_CR_ope,'spline');
    Aero.alpha_vee_CR_ope = alpha_vee_CR_ope;
    Aero.CL_vee_CR_ope = CL_vee_CR_ope;
    
    % Claculates Curve Lift Slope According to the min AoA selected (alpha
    % selected) and the max CL (CLmax/1.2^2)
    CL_alpha_vee_CR = (CL_vee_CR_ope - DATA_Ae(index_vee).CL(DATA_Ae(index_vee).alpha==alpha_selected_vee)*conv_vee)/...
        (alpha_vee_CR_ope - DATA_Ae(index_vee).alpha(DATA_Ae(index_vee).alpha==alpha_selected_vee)*conv_vee)*180/pi;
    Aero.CL_alpha_vee_CR = CL_alpha_vee_CR;
    
    K_nacareport = get_k_823(Geo_tier.AR_vee_e,Geo_tier.lambda_vee);
    [KVB, KBV] = get_Vtail_body_interference(Body_Geo, Geo_tier.x_vee_LE, Geo_tier.b_vee_e);
    
    datos_sin_diedro = 0; % proyectado, si proviene de datos de wind tunnel o xflr5 con diedro datos_sin_diedro = 0;
    if datos_sin_diedro == 1
        CL_alpha_wb_Vee_0 = Aero.CL_alpha_vee_CR;
        %ExpresiÃ³n obtenida de Lift Curve Slope-V-Tail (by DARcorporation)
        %en Help
        CLalpha_vee = CL_alpha_wb_Vee_0*cos((Geo_tier.dihedral_vee))^2*(Geo_tier.S_vee_e/Geo_tier.S_vee)*(KVB + KBV); %multiplicado por cos(diedro)^2 y con la interferencia del fuselaje.
        CYbeta_vee = -K_nacareport*CL_alpha_wb_Vee_0*sin((Geo_tier.dihedral_vee))^2*(Geo_tier.S_vee_e/Geo_tier.S_vee)*(KVB + KBV);%multiplicado por sin(diedro)^2 y con la interferencia del fuselaje.
        
    else
        if index_vee_cy ~= 0
            CLalpha_vee = Aero.CL_alpha_vee_CR;
            %              DATA_Ae(index_w1).CL(DATA_Ae(index_w1).alpha==alpha_selected_w1
            CY_max_vee_CR = max(DATA_Ae(index_vee_cy).CY);
            
            beta_vee_CR_ope  = interp1(DATA_Ae(index_vee_cy).CY,DATA_Ae(index_vee_cy).beta,CY_max_vee_CR/Flight_SF^2,'spline');
            CY_vee_CR_ope = interp1(DATA_Ae(index_vee_cy).beta,DATA_Ae(index_vee_cy).CY,beta_vee_CR_ope,'spline');
%             CL_vee_CR_ope = interp1(DATA_Ae(index_vee).alpha,DATA_Ae(index_vee).CL*conv_vee,alpha_vee_CR_ope,'spline');

            CYbeta_vee = (CY_vee_CR_ope - DATA_Ae(index_vee_cy).CY(DATA_Ae(index_vee_cy).beta==0))/...
                (beta_vee_CR_ope - DATA_Ae(index_vee_cy).beta(DATA_Ae(index_vee_cy).beta==0))*180/pi;
            
            [CY_max_Vee_CR,max_index] = max(DATA_Ae(index_vee_cy).CY);
            [CY_min_Vee_CR,min_index] = min(DATA_Ae(index_vee_cy).CY);
            
            beta_max_Vee_CR = DATA_Ae(index_vee_cy).beta(max_index);
            beta_min_Vee_CR = DATA_Ae(index_vee_cy).beta(min_index);
            % CLmax
            CYbeta_vee2 = -(CY_max_Vee_CR -CY_min_Vee_CR)*conv_vee/...
                abs(beta_max_Vee_CR - beta_min_Vee_CR)*180/pi;
            
            %         CYbeta_vee = -K_nacareport*CLalpha_vee*tan((Geo_tier.dihedral_vee))^2; %o windtunnel directamente
            %         CYbeta_vee = windtunnel;
        else
            CLalpha_vee = Aero.CL_alpha_vee_CR;
            CYbeta_vee = -K_nacareport*CLalpha_vee*tan(Geo_tier.dihedral_vee)^2;
            
        end       
        
        
    end
    Aero.CLalpha_vee = CLalpha_vee;
    Aero.CYbeta_vee = CYbeta_vee;
    
    % Efficiency vectors
    % Max Range Jet
    E_vee_max_R_jet_vec = (DATA_Ae(index_vee).CL*conv_vee)./DATA_Ae(index_vee).CD*conv_vee.^(3/2);
    E_vee_max_R_jet = max(E_vee_max_R_jet_vec);
    Aero.E_vee_max_R_jet_vec = E_vee_max_R_jet_vec;
    Aero.E_vee_max_R_jet = E_vee_max_R_jet;
    % Max Range Prop
    E_vee_max_R_prop_vec = (DATA_Ae(index_vee).CL*conv_vee)./DATA_Ae(index_vee).CD*conv_vee;
    E_vee_max_R_prop = max(E_vee_max_R_prop_vec);
    Aero.E_vee_max_R_prop_vec = E_vee_max_R_prop_vec;
    % Max Endurance Jet
    E_vee_max_E_jet_vec = (DATA_Ae(index_vee).CL*conv_vee)./DATA_Ae(index_vee).CD*conv_vee;
    E_vee_max_E_jet = max(E_vee_max_E_jet_vec);
    Aero.E_vee_max_E_jet_vec = E_vee_max_E_jet_vec;
    Aero.E_vee_max_E_jet = E_vee_max_E_jet;
    % Max Endurance Prop
    E_vee_max_E_prop_vec = sqrt(3)*E_vee_max_R_prop;
    E_vee_max_E_prop = max(E_vee_max_E_prop_vec);
    Aero.E_vee_max_E_prop_vec = E_vee_max_E_prop_vec;
    Aero.E_vee_max_E_prop = E_vee_max_E_prop;
    CL_ref_limit_min = CL_min_V_max_CR;
    CL_ref_limit_max = CL_max_V_min_CR;
    % Polyfit Polar
%     index_w1
%     index_vee
%     DATA_Ae
%     DATA_Ae(index_w1)
%     DATA_Ae(index_vee)
%     CL_w1_limit_min
%     pause
    PL_vee = polyfit(DATA_Ae(index_vee).CL(DATA_Ae(index_vee).CL>= CL_w1_limit_min & DATA_Ae(index_vee).CL<=CL_w1_limit_max)*conv_vee,...
        DATA_Ae(index_vee).CD(DATA_Ae(index_vee).CL>=CL_w1_limit_min & DATA_Ae(index_vee).CL<= CL_w1_limit_max)*conv_vee,2);
%     PL_vee = PL_w1;
    Aero.PL_vee = PL_vee;
    CD0_vee = PL_vee(3);
    CD1_vee = PL_vee(2);
    CD2_vee = PL_vee(1);
    Aero.CD0_vee = CD0_vee;
    Aero.CD1_vee = CD1_vee;
    Aero.CD2_vee = CD2_vee;
    % Generates Vector with polyfit approximation
    Aero.CD_poly_vee = CD0_vee + CD1_vee.*DATA_Ae(index_vee).CL*conv_vee + CD2_vee.*(conv_vee*DATA_Ae(index_vee).CL.^2);
% else
%     Aero.CD0_vee = 0;
%     Aero.CD1_vee = 0;
%     Aero.CD2_vee = 0;
%     S_vee = 0;
end

if Vee2 == 1
    % Reference wing
    S_vee2 = Geo_tier.S_vee2;
    % Reference wing
    S_vee2_s = Geo_tier.S_vee2_s;
    S_vee2_e = Geo_tier.S_vee2_e;
    S_w1_e = Geo_tier.S_w1_e;
    
    % conv_w2 = S_w2_s/S_w1_e;
    % Assume no surface conversion: aerodynamic data is as original from
    % Aerodynamic Software no needed 
    conv_vee2 = 1;
    
    % incidence of the FW, RW and the overall configuration
    i_vee2 = Design_criteria.i_vee2; % incidence of Rear Wing
    index_vee2_cy = Design_criteria.index_vee2_cy; %
    
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_vee2 = Design_criteria.index_vee2;
    index_vee2_cy = Design_criteria.index_vee2_cy; %
%     index_vee2_cy = 7; %

    % incidence of the FW, RW and the overall configuration
    i_vee2 = Design_criteria.i_vee2; % incidence of Rear Wing
    
    % Allows the user to select the min AoA used to determine Curve Lift Slope
    alpha_selected_vee2 = Design_criteria.alpha_selected_vee2;
    
    % CLmax
    CL_max_vee2_CR = max(DATA_Ae(index_vee2).CL*conv_vee2);
    Aero.CL_max_vee2_CR = CL_max_vee2_CR;
    % Alpha at max CL
    alpha_max_vee2_CR = interp1(DATA_Ae(index_vee2).CL*conv_vee2,DATA_Ae(index_vee2).alpha,CL_max_vee2_CR,'spline');
    Aero.alpha_max_vee2_CR = alpha_max_vee2_CR;
    
    %% Determination of Stall conditions
    CL_max_vee2_ope = CL_max_vee2_CR/(Flight_SF^2);
    V_stall_vee2 = sqrt((2*m_TOW*g)/(rho*S_ref*CL_max_vee2_ope));
    V_min_vee2 = Flight_SF*V_stall_vee2;
    
    alpha_max_vee2_ope = interp1(DATA_Ae(index_vee2).CL,DATA_Ae(index_vee2).alpha,CL_max_vee2_ope,'spline');
    % Defines limits on CL for all aerodynamic surfaces
    CL_vee2_limit_min = m_TOW*g/(0.5*rho*(V_max^2)*S_vee2);
    CL_vee2_limit_max = m_TOW*g/(0.5*rho*(V_min_vee2^2)*S_vee2);
    Aero.CL_vee2_limit_min = CL_vee2_limit_min;
    Aero.CL_vee2_limit_max = CL_vee2_limit_max;
    % CL at 0 AoA
    alpha_CL_0_vee2_CR = interp1(DATA_Ae(index_vee2).CL*conv_vee2,DATA_Ae(index_vee2).alpha,0,'spline');
    Aero.alpha_CL_0_vee2_CR = alpha_CL_0_vee2_CR;
    CL_0_vee2_CR = interp1(DATA_Ae(index_vee2).alpha,DATA_Ae(index_vee2).CL*conv_vee2,0,'spline');
    Aero.CL_0_vee2_CR = CL_0_vee2_CR;
    
    % Lift coefficient during zero AoA (fixed i_w1, and i_vee2)
    % VLM results
    CL_vee2_CR_iw = interp1(DATA_Ae(index_vee2).alpha,DATA_Ae(index_vee2).CL*conv_vee2,i_vee2,'spline');
    Aero.CL_vee2_CR_iw = CL_vee2_CR_iw;
    % Moment at alpha=0
    CM_0_vee2_CR = interp1(DATA_Ae(index_vee2).alpha,DATA_Ae(index_vee2).Cm*conv_vee2,0,'spline');
    Aero.CM_0_vee2_CR = CM_0_vee2_CR;
    
    % Alpha at max CL with Fligth Safe Margin - operative, where operative
    % equals to V_ope = Flight_SF*V_stall
    alpha_vee2_CR_ope  = interp1(DATA_Ae(index_vee2).CL*conv_vee2,DATA_Ae(index_vee2).alpha,CL_max_vee2_CR/Flight_SF^2,'spline');
    CL_vee2_CR_ope = interp1(DATA_Ae(index_vee2).alpha,DATA_Ae(index_vee2).CL*conv_vee2,alpha_vee2_CR_ope,'spline');
    Aero.alpha_vee2_CR_ope = alpha_vee2_CR_ope;
    Aero.CL_vee2_CR_ope = CL_vee2_CR_ope;
    
    % Claculates Curve Lift Slope According to the min AoA selected (alpha
    % selected) and the max CL (CLmax/1.2^2)
    CL_alpha_vee2_CR = (CL_vee2_CR_ope - DATA_Ae(index_vee2).CL(DATA_Ae(index_vee2).alpha==alpha_selected_vee2)*conv_vee2)/...
        (alpha_vee2_CR_ope - DATA_Ae(index_vee2).alpha(DATA_Ae(index_vee2).alpha==alpha_selected_vee2)*conv_vee2)*180/pi;
    Aero.CL_alpha_vee2_CR = CL_alpha_vee2_CR;
    
    K_nacareport = get_k_823(Geo_tier.AR_vee2_e,Geo_tier.lambda_vee2);
    [KVB, KBV] = get_Vtail_body_interference(Body_Geo, Geo_tier.x_vee2_LE, Geo_tier.b_vee2_e);
    
    datos_sin_diedro = 0; % proyectado, si proviene de datos de wind tunnel o xflr5 con diedro datos_sin_diedro = 0;
    if datos_sin_diedro == 1
        CL_alpha_wb_vee2_0 = Aero.CL_alpha_vee2_CR;
        %ExpresiÃ³n obtenida de Lift Curve Slope-V-Tail (by DARcorporation)
        %en Help
        CLalpha_vee2 = CL_alpha_wb_vee2_0*cos((Geo_tier.dihedral_vee2))^2*(Geo_tier.S_vee2_e/Geo_tier.S_vee2)*(KVB + KBV); %multiplicado por cos(diedro)^2 y con la interferencia del fuselaje.
        CYbeta_vee2 = -K_nacareport*CL_alpha_wb_vee2_0*sin((Geo_tier.dihedral_vee2))^2*(Geo_tier.S_vee2_e/Geo_tier.S_vee2)*(KVB + KBV);%multiplicado por sin(diedro)^2 y con la interferencia del fuselaje.
        
    else
        if index_vee2_cy ~= 0
            CLalpha_vee2 = Aero.CL_alpha_vee2_CR;
            %              DATA_Ae(index_w1).CL(DATA_Ae(index_w1).alpha==alpha_selected_w1
            CY_max_vee2_CR = max(DATA_Ae(index_vee2_cy).CY);
            
            beta_vee2_CR_ope  = interp1(DATA_Ae(index_vee2_cy).CY,DATA_Ae(index_vee2_cy).beta,CY_max_vee2_CR/Flight_SF^2,'spline');
            CY_vee2_CR_ope = interp1(DATA_Ae(index_vee2_cy).beta,DATA_Ae(index_vee2_cy).CY,beta_vee2_CR_ope,'spline');
%             CL_vee2_CR_ope = interp1(DATA_Ae(index_vee2).alpha,DATA_Ae(index_vee2).CL*conv_vee2,alpha_vee2_CR_ope,'spline');

            CYbeta_vee2 = (CY_vee2_CR_ope - DATA_Ae(index_vee2_cy).CY(DATA_Ae(index_vee2_cy).beta==0))/...
                (beta_vee2_CR_ope - DATA_Ae(index_vee2_cy).beta(DATA_Ae(index_vee2_cy).beta==0))*180/pi;
            
            [CY_max_vee2_CR,max_index] = max(DATA_Ae(index_vee2_cy).CY);
            [CY_min_vee2_CR,min_index] = min(DATA_Ae(index_vee2_cy).CY);
            
            beta_max_vee2_CR = DATA_Ae(index_vee2_cy).beta(max_index);
            beta_min_vee2_CR = DATA_Ae(index_vee2_cy).beta(min_index);
            % CLmax
            CYbeta_vee22 = -(CY_max_vee2_CR -CY_min_vee2_CR)*conv_vee2/...
                abs(beta_max_vee2_CR - beta_min_vee2_CR)*180/pi;
            
            %         CYbeta_vee2 = -K_nacareport*CLalpha_vee2*tan((Geo_tier.dihedral_vee2))^2; %o windtunnel directamente
            %         CYbeta_vee2 = windtunnel;
        else
            CLalpha_vee2 = Aero.CL_alpha_vee2_CR;
            CYbeta_vee2 = -K_nacareport*CLalpha_vee2*tan(Geo_tier.dihedral_vee2)^2;
            
        end       
        
        
    end
    Aero.CLalpha_vee2 = CLalpha_vee2;
    Aero.CYbeta_vee2 = CYbeta_vee2;
    
    % Efficiency vectors
    % Max Range Jet
    E_vee2_max_R_jet_vec = (DATA_Ae(index_vee2).CL*conv_vee2)./DATA_Ae(index_vee2).CD*conv_vee2.^(3/2);
    E_vee2_max_R_jet = max(E_vee2_max_R_jet_vec);
    Aero.E_vee2_max_R_jet_vec = E_vee2_max_R_jet_vec;
    Aero.E_vee2_max_R_jet = E_vee2_max_R_jet;
    % Max Range Prop
    E_vee2_max_R_prop_vec = (DATA_Ae(index_vee2).CL*conv_vee2)./DATA_Ae(index_vee2).CD*conv_vee2;
    E_vee2_max_R_prop = max(E_vee2_max_R_prop_vec);
    Aero.E_vee2_max_R_prop_vec = E_vee2_max_R_prop_vec;
    % Max Endurance Jet
    E_vee2_max_E_jet_vec = (DATA_Ae(index_vee2).CL*conv_vee2)./DATA_Ae(index_vee2).CD*conv_vee2;
    E_vee2_max_E_jet = max(E_vee2_max_E_jet_vec);
    Aero.E_vee2_max_E_jet_vec = E_vee2_max_E_jet_vec;
    Aero.E_vee2_max_E_jet = E_vee2_max_E_jet;
    % Max Endurance Prop
    E_vee2_max_E_prop_vec = sqrt(3)*E_vee2_max_R_prop;
    E_vee2_max_E_prop = max(E_vee2_max_E_prop_vec);
    Aero.E_vee2_max_E_prop_vec = E_vee2_max_E_prop_vec;
    Aero.E_vee2_max_E_prop = E_vee2_max_E_prop;
    CL_ref_limit_min = CL_min_V_max_CR;
    CL_ref_limit_max = CL_max_V_min_CR;
    % Polyfit Polar
%     index_w1
%     index_vee2
%     DATA_Ae
%     DATA_Ae(index_w1)
%     DATA_Ae(index_vee2)
%     CL_w1_limit_min
%     pause
    PL_vee2 = polyfit(DATA_Ae(index_vee2).CL(DATA_Ae(index_vee2).CL>= CL_w1_limit_min & DATA_Ae(index_vee2).CL<=CL_w1_limit_max)*conv_vee2,...
        DATA_Ae(index_vee2).CD(DATA_Ae(index_vee2).CL>=CL_w1_limit_min & DATA_Ae(index_vee2).CL<= CL_w1_limit_max)*conv_vee2,2);
%     PL_vee2 = PL_w1;
    Aero.PL_vee2 = PL_vee2;
    CD0_vee2 = PL_vee2(3);
    CD1_vee2 = PL_vee2(2);
    CD2_vee2 = PL_vee2(1);
    Aero.CD0_vee2 = CD0_vee2;
    Aero.CD1_vee2 = CD1_vee2;
    Aero.CD2_vee2 = CD2_vee2;
    % Generates Vector with polyfit approximation
    Aero.CD_poly_vee2 = CD0_vee2 + CD1_vee2.*DATA_Ae(index_vee2).CL*conv_vee2 + CD2_vee2.*(conv_vee2*DATA_Ae(index_vee2).CL.^2);
% else
%     Aero.CD0_vee2 = 0;
%     Aero.CD1_vee2 = 0;
%     Aero.CD2_vee2 = 0;
%     S_vee2 = 0;
end

if VTP == 1
    % Reference wing
    S_VTP = Geo_tier.S_VTP;
    % Reference wing
    S_VTP_s = Geo_tier.S_VTP_s;
    S_VTP_e = Geo_tier.S_VTP_e;
    conv_VTP = S_VTP_s/S_VTP_e;
    
    % AC_type = 1 - flying wing
    % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
    % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
    % AC_type = 4 - 2 surface: wing + V-tail
    % AC_type = 5 - 3 surface: cannard + wing + V-tail
    % AC_type = 6 - 2 surface: cannard + wing + VTP
    % AC_type = 7 - 3 surface: wing + V-tail + VTP
    % AC_type = 8 - 2 surface: wing + V-tail + V-tail


    AC_type = OUTPUT_read_XLSX.AC_Data_flags.AC_type;
    switch AC_type
        case 1
            % Aerodynamic results conducted with the wing
            S_w1_e = Geo_tier.S_w1_e;
            % conv_VTP = S_w1_e/S_w1_e;
            conv_VTP = S_w1_e/S_VTP_e;
            if AC_CONFIGURATION.twin_VTP == 1
                conv_VTP = S_w1_e/S_VTP_e;
            end

        case 2
            % If there is HTP the Aerodynamic studies in FLOW have been conducted
            % HTP + VTP hence the area used for adimentionalizing is the one from the HTP
            if HTP == 1
                S_w1_e = Geo_tier.S_w1_e;
                S_HTP_e = Geo_tier.S_w2_e;
                % conv_VTP = S_HTP_e/S_w1_e;
                conv_VTP = S_HTP_e/S_VTP_e;
            end
        case 3
            % If there is HTP the Aerodynamic studies in FLOW have been conducted
            % HTP + VTP hence the area used for adimentionalizing is the one from the HTP
            if HTP == 1
                S_w1_e = Geo_tier.S_w1_e;
                S_HTP_e = Geo_tier.S_w2_e;
                % conv_VTP = S_HTP_e/S_w1_e;
                conv_VTP = S_HTP_e/S_VTP_e;
            end
        case 4
            % No VTP
        case 5
            % No VTP
        case 6
            % Aerodynamic results conducted with the wing
            S_w1_e = Geo_tier.S_w1_e;
            conv_VTP = S_w1_e/S_w1_e;
            conv_VTP = S_w1_e/S_VTP_e;
        case 7
            % Aerodynamic results conducted with the wing
            S_w1_e = Geo_tier.S_w1_e;
            conv_VTP = S_w1_e/S_w1_e;
            conv_VTP = S_w1_e/S_VTP_e;
        case 8
            % No VTP    
    end
    
    % Assume no surface conversion: aerodynamic data is as original from
    % Aerodynamic Software
    % conv_VTP = 1;
    
    % Select txt files associated with each aerodynamic study
    % Make sure to check with case asssigned in read_aero_files_Aug2018.m
    index_VTP = Design_criteria.index_VTP;
    
    % Allows the user to select the min AoA used to determine Curve Lift Slope
    % alpha_selected_VTP = Design_criteria.alpha_selected_VTP;
    % incidence of the FW, RW and the overall configuration
    i_VTP = Design_criteria.i_w1; % incidence of Front Wing
    
    % CYmax
    %     CY_max_VTP_CR = max(abs(DATA_Ae(index_VTP).CY*conv_VTP));
    %     CY_min_VTP_CR = min(abs(DATA_Ae(index_VTP).CY*conv_VTP));
    %     beta_max_VTP_CR = interp1(DATA_Ae(index_VTP).CY*conv_VTP,DATA_Ae(index_VTP).beta,CY_max_VTP_CR,'spline');
    %     beta_min_VTP_CR = interp1(DATA_Ae(index_VTP).CY*conv_VTP,DATA_Ae(index_VTP).beta,CY_min_VTP_CR,'spline');
    %     beta_max_VTP_CR = max(DATA_Ae(index_VTP).beta);
    %     beta_min_VTP_CR = max(DATA_Ae(index_VTP).beta);
    
    [CY_max_VTP_CR,max_index] = max(DATA_Ae(index_VTP).CY);
    [CY_min_VTP_CR,min_index] = min(DATA_Ae(index_VTP).CY);
    
    beta_max_VTP_CR = DATA_Ae(index_VTP).beta(max_index);
    beta_min_VTP_CR = DATA_Ae(index_VTP).beta(min_index);
    % CLmax
  
    CL_alpha_VTP_CR = conv_VTP*(abs(CY_max_VTP_CR - CY_min_VTP_CR)/(abs(beta_max_VTP_CR - beta_min_VTP_CR)*pi/180));

    % CL_alpha_VTP_CR = conv_VTP*((abs(CY_max_VTP_CR -CY_min_VTP_CR)/...
    %     abs(beta_max_VTP_CR - beta_min_VTP_CR)*pi/180))
       
    Aero.CL_alpha_VTP_CR = CL_alpha_VTP_CR;

    % generates the CL and alpha such that avoids the problem of "Sample
    % points must be unique and sorted in ascending order." that it is
    % encountered
    % DATA_Ae(index_VTP).CY
    % DATA_Ae(index_VTP).beta
    % CY_new = DATA_Ae(index_VTP).CY(1:max_index)
    % beta_new = DATA_Ae(index_VTP).beta(1:max_index)

    % Correction to avoid "Sample points must be unique and sorted in ascending order."
    CY_0_VTP_CR = conv_VTP*interp1(DATA_Ae(index_VTP).beta,DATA_Ae(index_VTP).CY,0,'spline');
    % CY_0_VTP_CR = conv_VTP*interp1(beta_new,CY_new,0,'spline');


    Aero.CY_0_VTP_CR = CY_0_VTP_CR;

    %     switch AC_type
%         case 1
%             % Aerodynamic results conducted with the wing
%             S_w1_e = Geo_tier.S_w1_e;
%             conv_VTP = S_w1_e/S_w1_e;
%             CL_alpha_VTP_CR = CL_alpha_VTP_CR*S_w1_e/S_VTP_e;
%             
%         case 2
%             % If there is HTP the Aerodynamic studies in FLOW have been conducted
%             % HTP + VTP hence the area used for adimentionalizing is the one from the HTP
%             if HTP == 1
%                 S_w1_e = Geo_tier.S_w1_e;
%                 S_HTP_e = Geo_tier.S_w2_e;
%                 conv_VTP = S_HTP_e/S_w1_e;
%                 CL_alpha_VTP_CR = CL_alpha_VTP_CR*S_w2_e/S_VTP_e;
%                 
%             end
%         case 3
%             % If there is HTP the Aerodynamic studies in FLOW have been conducted
%             % HTP + VTP hence the area used for adimentionalizing is the one from the HTP
%             if HTP == 1
%                 S_w1_e = Geo_tier.S_w1_e;
%                 S_HTP_e = Geo_tier.S_w2_e;
%                 conv_VTP = S_HTP_e/S_w1_e;
%                 CL_alpha_VTP_CR = CL_alpha_VTP_CR*S_w2_e/S_VTP_e;
%             end
%         case 4
%             % No VTP
%         case 5
%             % No VTP
%         case 6
%             % Aerodynamic results conducted with the wing
%             S_w1_e = Geo_tier.S_w1_e;
%             conv_VTP = S_w1_e/S_w1_e;
%             CL_alpha_VTP_CR = CL_alpha_VTP_CR*S_w1_e/S_VTP_e;
%     end
%     Aero.CL_alpha_VTP_CR = CL_alpha_VTP_CR;
    
end

% Fuse polar model
Aero.CD0_ac =  Aero.CD0_w1*S_w1/S_ref + Aero.CD0_w2*S_w2/S_ref + Aero.CD0_can*S_can/S_ref;
Aero.CD1_ac =  Aero.CD1_w1*S_w1/S_ref + Aero.CD1_w2*S_w2/S_ref + Aero.CD1_can*S_can/S_ref;
Aero.CD2_ac =  Aero.CD2_w1*S_w1/S_ref + Aero.CD2_w2*S_w2/S_ref + Aero.CD2_can*S_can/S_ref;

% Calculates the total CLmax associated to the differente control surfaces
switch OUTPUT_read_XLSX.AC_Data_flags.AC_type
    case 1 % AC_type = 1 - flying wing + VTP
        CL_max_ac_ope = CL_max_w1_ope*S_w1/S_ref;
        CL_max_ac = CL_max_w1_CR*S_w1/S_ref;        
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        CL_max_ac_ope = CL_max_w1_ope*S_w1/S_ref;
        CL_max_ac = CL_max_w1_CR*S_w1/S_ref;        
    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        CL_max_ac_ope = CL_max_w1_ope*S_w1/S_ref + CL_max_can_ope*S_can/S_ref;
        CL_max_ac = CL_max_w1_CR*S_w1/S_ref + CL_max_can_CR*S_can/S_ref;
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        CL_max_ac_ope = CL_max_w1_ope*S_w1/S_ref;
        CL_max_ac = CL_max_w1_CR*S_w1/S_ref;        
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        CL_max_ac_ope = CL_max_w1_ope*S_w1/S_ref + CL_max_can_ope*S_can/S_ref;
        CL_max_ac = CL_max_w1_CR*S_w1/S_ref + CL_max_can_CR*S_can/S_ref;
    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        CL_max_ac_ope = CL_max_w1_ope*S_w1/S_ref + CL_max_can_ope*S_can/S_ref;
        CL_max_ac = CL_max_w1_CR*S_w1/S_ref + CL_max_can_CR*S_can/S_ref;
    case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        CL_max_ac_ope = CL_max_w1_ope*S_w1/S_ref;
        CL_max_ac = CL_max_w1_CR*S_w1/S_ref ;
    case 8 % AC_type = 8 - 2 surface: wing + V-tail + V-tail
        CL_max_ac_ope = CL_max_w1_ope*S_w1/S_ref;
        CL_max_ac = CL_max_w1_CR*S_w1/S_ref ;
end

Aero.CL_max_ac_ope = CL_max_ac_ope;
Aero.CL_max_ac = CL_max_ac;

V_stall_ac = sqrt((2*m_TOW*g)/(rho*S_ref*CL_max_ac));
V_min_ac = Flight_SF*V_stall_ac;

Performance.V_stall_ac = V_stall_ac;
Performance.V_min_ac = V_min_ac;

%% High Lift Devices
Flap_type = OUTPUT_read_XLSX.Aerodynamic_Data_flags.Flap_type; %% High Lift Devices  - Leading Edge
LED_type = OUTPUT_read_XLSX.Aerodynamic_Data_flags.LED_type; %%High Lift Devices  - Trailing  Edge (Flaps)	

% Estimation of DeltaCLmax associated to Raymer
% Flaps Figure 12.6 Rymer and Table 12.2 Raymer
switch Flap_type
    case 0 % No Flap
        Delta_Clmax_TE = 0.0;  
    case 1 % Plain and split
        Delta_Clmax_TE = 0.9;  
    case 2 % Slotted
        Delta_Clmax_TE = 1.3;  
    case 3 % Fowler
        cprime = 53.77;
        c = 48.21;
        cprime_c = cprime/c; % approximated
        Delta_Clmax_TE = 1.3*cprime_c;  
    case 4 % Doble slotted
        cprime = 50.62;
        c = 48.21;
        cprime_c = cprime/c; % approximated
        Delta_Clmax_TE = 1.6*cprime_c;  
    case 5 % Triple slotted
        cprime = 52.30;
        c = 48.21;
        cprime_c = cprime/c; % approximated
        Delta_Clmax_TE = 1.9*cprime_c;  
end
% HL
K_y1_flap_w1 = Geo_tier.K_y1_flap_w1;
K_y2_flap_w1 = Geo_tier.K_y2_flap_w1;
K_y1y2_flap_w1 = K_y2_flap_w1 - K_y1_flap_w1;
Lambda_c34_w1 = Geo_tier.Lambda_c34_w1;
Delta_CLmax_TE = Delta_Clmax_TE*K_y1y2_flap_w1*cos(Lambda_c34_w1);

% Estimation of DeltaCLmax associated to Raymer
% Leading Edge Devices Figure 12.6 Rymer and Table 12.2 Raymer
switch LED_type
    case 0 % No leading Edge Devices
        Delta_Clmax_LE = 0.0;  
    case 1 % Fixed Slot
        Delta_Clmax_LE = 0.2;  
    case 2 % Leading Edge Flap
        Delta_Clmax_LE = 0.3;  
    case 3 % Krugger Flap
        Delta_Clmax_LE = 0.3;  
    case 4 % Slat
        cprime = 52.90;
        c = 49.44;
        cprime_c = cprime/c; % approximated
        Delta_Clmax_LE = 0.4*cprime_c;
end
Delta_CLmax_LE = Delta_Clmax_LE*K_y1y2_flap_w1*cos(Lambda_c34_w1);
Aero.Delta_CLmax_TE = Delta_CLmax_TE;
Aero.Delta_CLmax_LE = Delta_CLmax_LE;
Delta_CLmax = Delta_CLmax_LE + Delta_CLmax_TE;
Aero.Delta_CLmax = Delta_CLmax;

%% Calculates the total CLmax associated to the differente control surfaces
% for the use of High Lift Devices
switch OUTPUT_read_XLSX.AC_Data_flags.AC_type
    case 1 % AC_type = 1 - flying wing + VTP
        CL_max_ac_ope_TO = (CL_max_w1_ope + Delta_CLmax)*S_w1/S_ref;
        CL_max_ac_TO = (CL_max_w1_CR + Delta_CLmax)*S_w1/S_ref;        
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        CL_max_ac_ope_TO = (CL_max_w1_ope + Delta_CLmax)*S_w1/S_ref;
        CL_max_ac_TO = (CL_max_w1_CR + Delta_CLmax)*S_w1/S_ref;        
    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        CL_max_ac_ope_TO = (CL_max_w1_ope + Delta_CLmax)*S_w1/S_ref + CL_max_can_ope*S_can/S_ref;
        CL_max_ac_TO = (CL_max_w1_CR + Delta_CLmax)*S_w1/S_ref + CL_max_can_CR*S_can/S_ref;
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        CL_max_ac_ope_TO = (CL_max_w1_ope + Delta_CLmax)*S_w1/S_ref;
        CL_max_ac_TO = (CL_max_w1_CR + Delta_CLmax)*S_w1/S_ref;        
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        CL_max_ac_ope_TO = (CL_max_w1_ope + Delta_CLmax)*S_w1/S_ref + CL_max_can_ope*S_can/S_ref;
        CL_max_ac_TO = (CL_max_w1_CR + Delta_CLmax)*S_w1/S_ref + CL_max_can_CR*S_can/S_ref;
    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        CL_max_ac_ope_TO = (CL_max_w1_ope + Delta_CLmax)*S_w1/S_ref + CL_max_can_ope*S_can/S_ref;
        CL_max_ac_TO = (CL_max_w1_CR + Delta_CLmax)*S_w1/S_ref + CL_max_can_CR*S_can/S_ref;
    case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        CL_max_ac_ope_TO = (CL_max_w1_ope + Delta_CLmax)*S_w1/S_ref;
        CL_max_ac_TO = (CL_max_w1_CR + Delta_CLmax)*S_w1/S_ref;        
    case 8 % AC_type = 8 - 2 surface: wing + V-tail + V-tail
        CL_max_ac_ope_TO = (CL_max_w1_ope + Delta_CLmax)*S_w1/S_ref;
        CL_max_ac_TO = (CL_max_w1_CR + Delta_CLmax)*S_w1/S_ref;        
end

Aero.CL_max_ac_ope_TO = CL_max_ac_ope_TO;
Aero.CL_max_ac_TO = CL_max_ac_TO;

V_stall_ac_TO = sqrt((2*m_TOW*g)/(rho*S_ref*CL_max_ac_TO));
V_min_ac_TO = Flight_SF*V_stall_ac_TO;

Performance.V_stall_ac_TO = V_stall_ac_TO;
Performance.V_min_ac_TO = V_min_ac_TO;

