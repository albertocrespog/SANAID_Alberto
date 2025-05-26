function [Stab_Der_parts Stab_Der] = get_hinge_derivatives(AC_CONFIGURATION,modelo,TRIM_RESULTS,Stab_Der,Stab_Der_parts,Trim_ITER,OUTPUT_read_XLSX,Geo_tier,conditions,Performance)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

% Available control surfaces
d_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_ail; %Definition of available control surface - aileron
d_ele = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_ele; %Definition of available control surface - elevator
d_elevon = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_elevon; %Definition of available control surface - elevon
d_flap = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_flap; %Definition of available control surface - flap
d_rudder = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudder; %Definition of available control surface - rudder
d_rudvtr = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudvtr; %Definition of available control surface - ruddervator
d_rudvtr2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudvtr2; %Definition of available control surface - ruddervator
d_can = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_can; %Definition of available control surface - canard

d_TT_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_TT_ail;
d_TT_ele = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_TT_ele;
d_TT_rudder = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_TT_rudder;
d_TT_sp = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_TT_sp;

% Geo_tier = Storing_GEO_DATA.Geo_tier;
% Body_Geo = Storing_GEO_DATA.Body_Geo;
% Performance = Storing_AERO_DATA.Performance;
% Aero = Storing_AERO_DATA.Aero;
% Weight_tier = Storing_WEIGHT_DATA.Weight_tier;

prop_wash_effect = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;
Mach = Performance.Mach;

% Geometry aileron (Wing)
if W1 == 1
    K_y1_ail_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_ail_w1;
    K_y2_ail_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_ail_w1;
    Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
    Lambda_c34_w1 = Geo_tier.Lambda_c34_w1;
    Lambda_TE_w1 = Geo_tier.Lambda_TE_w1;
    Lambda_HL_aileron = Lambda_c34_w1; % approximate to 3/4 of chord
    AR_w1 = Geo_tier.AR_w1;
end

% Geometry rudder (VTP)
if VTP == 1
    Lambda_c4_VTP = Geo_tier.Lambda_c4_VTP;
    Lambda_c34_VTP = Geo_tier.Lambda_c34_VTP;
    Lambda_TE_VTP = Geo_tier.Lambda_TE_VTP;
    Lambda_HL_rudder = Lambda_c34_VTP; % approximate to 3/4 of chord
    AR_VTP = Geo_tier.AR_VTP;
    K_y1_rudder_VTP = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_rudder_VTP;
    K_y2_rudder_VTP = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_rudder_VTP;
end

if d_ail == 1
    % Geometry ailerontab (Wing)
    K_y1_TT_ail_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_TT_ail_w1;
    K_y2_TT_ail_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_TT_ail_w1;
    Lambda_HL_TT_aileron = Lambda_c34_w1; % approximate to 3/4 of chord

    %% CH delta_aileron
    % Chdelta_a = (Chdela_theory_uncorrected + DeltaChdelta_horn + DeltaChdelta_PartSpan)*cos(Sweep_c_4_w1)*cos(sweep_HL)
    sweep_factors = cos(Lambda_c4_w1)*cos(Lambda_HL_aileron);
    % aileron has a ct/cf = 0.29

    Chdelta_prima_Chdelta_theory = 0.9; % Fig 10.69a
    Chdelta_theory = -0.825; % Fig 10.69b
    Chdelta_prima = Chdelta_prima_Chdelta_theory*Chdelta_theory;
    Chdelta_2prima = Chdelta_prima;
    % Chdelta_bal_Chdelta_2prima = 0.4; % Fig 10.71 for round nose NACA 0015
    Chdelta_bal_Chdelta_2prima = 0.7; % Fig 10.71 for round nose NACA 0015
    Chdelta_bal = Chdelta_2prima*Chdelta_bal_Chdelta_2prima; % EQ 10.137
    Chdelta_M = Chdelta_bal/sqrt(1-Mach^2); % Eq 10.133
    fact1 =Chdelta_M;

    % Second term
    alpha_delta = 0.5; %Figure 8.17 for cf/c 0.30 and deflection of 20deg

    Chalpha_theory = -0.5; % figure 10.63
    Chalpha_prima_Chalpha_theory = 0.775; % figure 10.63
    Chalpha_prima = Chalpha_prima_Chalpha_theory*Chalpha_theory;
    Chalpha_2prima = Chalpha_prima;
    Chalpha_bal_Chalpha_2prima = 0.8; % figure 10.65 - for round nose NACA 0015
    Chalpha_bal = Chalpha_2prima*Chalpha_bal_Chalpha_2prima;
    Chalpha_M = Chalpha_bal/sqrt(1-Mach^2);
    fact2 = alpha_delta*Chalpha_M*(2*Lambda_c4_w1)/(AR_w1 + 2*Lambda_c4_w1);

    % Third term
    DeltaChdelta_adim = 0.05; % Fig 10.78a for AR = 9.7 Cessna 208
    Cldelta = 0.45;% Fig 8.1
    B2 = 1; % Fig 10.77b
    Kdelta_i = 1.75; % Fig 10.78
    Kdelta_o = 4.3;
    Kdelta = Kdelta_i*(1-K_y1_ail_w1) - Kdelta_o*(1-K_y2_ail_w1)/(K_y2_ail_w1-K_y1_ail_w1);
    DeltaChdelta = DeltaChdelta_adim*Cldelta*B2*Kdelta*sweep_factors;
    fact3 = DeltaChdelta;
    Ch_delta_a = sweep_factors*(fact1 + fact2 + fact3);
else
    Ch_delta_a = 0;
end

if d_rudder == 1
    %% CH delta_rudder
    % Chdelta_a = (Chdela_theory_uncorrected + DeltaChdelta_horn + DeltaChdelta_PartSpan)*cos(Sweep_c_4_w1)*cos(sweep_HL)
    sweep_factors = cos(Lambda_c4_VTP)*cos(Lambda_HL_rudder);
    % rudder has a ct/cf = 0.29

    Chdelta_prima_Chdelta_theory = 0.9; % Fig 10.69a
    Chdelta_theory = -0.825; % Fig 10.69b
    Chdelta_prima = Chdelta_prima_Chdelta_theory*Chdelta_theory;
    Chdelta_2prima = Chdelta_prima;
    Chdelta_bal_Chdelta_2prima = 0.4; % Fig 10.71 for round nose NACA 0015
    Chdelta_bal = Chdelta_2prima*Chdelta_bal_Chdelta_2prima; % EQ 10.137
    Chdelta_M = Chdelta_bal/sqrt(1-Mach^2); % Eq 10.133
    fact1 =Chdelta_M;

    % Second term
    alpha_delta = 0.5; %Figure 8.17 for cf/c 0.30 and deflection of 20deg

    Chalpha_theory = -0.5; % figure 10.63
    Chalpha_prima_Chalpha_theory = 0.775; % figure 10.63
    Chalpha_prima = Chalpha_prima_Chalpha_theory*Chalpha_theory;
    Chalpha_2prima = Chalpha_prima;
    Chalpha_bal_Chalpha_2prima = 0.8; % figure 10.65 - for round nose NACA 0015
    Chalpha_bal = Chalpha_2prima*Chalpha_bal_Chalpha_2prima;
    Chalpha_M = Chalpha_bal/sqrt(1-Mach^2);
    fact2 = alpha_delta*Chalpha_M*(2*Lambda_c4_VTP)/(AR_w1 + 2*Lambda_c4_VTP);

    % Third term
    DeltaChdelta_adim = 0.05; % Fig 10.78a for AR = 9.7 Cessna 208
    Cldelta = 0.45;% Fig 8.14
    B2 = 1; % Fig 10.77b
    Kdelta_i = 1.75; % Fig 10.78
    Kdelta_o = 4.3;
    Kdelta = Kdelta_i*(1-K_y1_rudder_VTP) - Kdelta_o*(1-K_y2_rudder_VTP)/(K_y2_rudder_VTP-K_y1_rudder_VTP);
    DeltaChdelta = DeltaChdelta_adim*Cldelta*B2*Kdelta*sweep_factors;
    fact3 = DeltaChdelta;
    Ch_delta_r = sweep_factors*(fact1 + fact2 + fact3);
else
    Ch_delta_r = 0;
end


%% CH delta_aileron_TT
if d_TT_ail == 1
    % Chdelta_a = (Chdela_theory_uncorrected + DeltaChdelta_horn + DeltaChdelta_PartSpan)*cos(Sweep_c_4_w1)*cos(sweep_HL)
    sweep_factors = cos(Lambda_HL_TT_aileron)*cos(Lambda_HL_TT_aileron);

    % Trim tab has a ct/cf = 0.29 and a ct/c 0.08
    Chdelta_prima_Chdelta_theory = 0.9; % Fig 10.69a
    Chdelta_theory = -0.825; % Fig 10.69b
    Chdelta_theory = -0.64; % Fig 10.69b
    Chdelta_prima = Chdelta_prima_Chdelta_theory*Chdelta_theory;
    Chdelta_2prima = Chdelta_prima;
    % Chdelta_bal_Chdelta_2prima = 0.4; % Fig 10.71 for round nose NACA 0015
    Chdelta_bal_Chdelta_2prima = 0.7; % Fig 10.71 for round nose NACA 0015
    Chdelta_bal = Chdelta_2prima*Chdelta_bal_Chdelta_2prima; % EQ 10.137
    Chdelta_M = Chdelta_bal/sqrt(1-Mach^2); % Eq 10.133
    fact1 =Chdelta_M;

    % Second term
    alpha_delta = 0.5; %Figure 8.17 for cf/c 0.30 and deflection of 20deg
    alpha_delta = 0.4; %Figure 8.17 for cf/c 0.30 and deflection of 20deg

    Chalpha_theory = -0.5; % figure 10.63
    Chalpha_theory = -0.23; % figure 10.63
    Chalpha_prima_Chalpha_theory = 0.775; % figure 10.63
    Chalpha_prima_Chalpha_theory = 0.7; % figure 10.63
    Chalpha_prima = Chalpha_prima_Chalpha_theory*Chalpha_theory;
    Chalpha_2prima = Chalpha_prima;
    Chalpha_bal_Chalpha_2prima = 0.8; % figure 10.65 - for round nose NACA 0015
    Chalpha_bal = Chalpha_2prima*Chalpha_bal_Chalpha_2prima;
    Chalpha_M = Chalpha_bal/sqrt(1-Mach^2);
    fact2 = alpha_delta*Chalpha_M*(2*Lambda_c4_w1)/(AR_w1 + 2*Lambda_c4_w1);

    % Third term
    DeltaChdelta_adim = 0.05; % Fig 10.78a for AR = 9.7 Cessna 208
    Cldelta = 0.45;% Fig 8.14
    B2 = 1; % Fig 10.77b
    B2 = 0.6; % Fig 10.77b
    Kdelta_i = 3.5; % Fig 10.78
    Kdelta_o = 4.3;
    Kdelta = Kdelta_i*(1-K_y1_TT_ail_w1) - Kdelta_o*(1-K_y2_TT_ail_w1)/(K_y2_TT_ail_w1-K_y1_TT_ail_w1);
    DeltaChdelta = DeltaChdelta_adim*Cldelta*B2*Kdelta*sweep_factors;
    fact3 = DeltaChdelta;
    Ch_delta_at = sweep_factors*(fact1 + fact2 + fact3);

    % method2 - 2D
    Chdelta_at_cldelta = -0.014*180/pi; % Fig 10.72 converted from deg^-1 to rad^-1
    Chcl_delta_at = -0.055; % Fig 10.73
    CLalpha_2D_wing = 0.104*180/pi; % From Table 8.1a NACA 23018
    alpha_delta_at_cl = -0.3 ; % Fig 10.74 for ct/c = 0.08
    Ch_delta_at2 = Chdelta_at_cldelta - (Chcl_delta_at*CLalpha_2D_wing*alpha_delta_at_cl);
else
    Ch_delta_at = 0;
end


if d_TT_rudder == 1

    %% CH delta_rudder_TT
    %% CH delta_rudder
    % Chdelta_a = (Chdela_theory_uncorrected + DeltaChdelta_horn + DeltaChdelta_PartSpan)*cos(Sweep_c_4_w1)*cos(sweep_HL)
    Lambda_HL_TT_rudder = Lambda_c4_VTP;
    Lambda_HL_TT_rudder = Lambda_HL_TT_rudder;
    sweep_factors = cos(Lambda_HL_TT_rudder)*cos(Lambda_HL_TT_rudder);

    % rudder has a ct/cf = 0.29

    Chdelta_prima_Chdelta_theory = 0.9; % Fig 10.69a
    Chdelta_theory = -0.825; % Fig 10.69b
    Chdelta_prima = Chdelta_prima_Chdelta_theory*Chdelta_theory;
    Chdelta_2prima = Chdelta_prima;
    Chdelta_bal_Chdelta_2prima = 0.4; % Fig 10.71 for round nose NACA 0015
    Chdelta_bal = Chdelta_2prima*Chdelta_bal_Chdelta_2prima; % EQ 10.137
    Chdelta_M = Chdelta_bal/sqrt(1-Mach^2); % Eq 10.133
    fact1 =Chdelta_M;

    % Second term
    alpha_delta = 0.5; %Figure 8.17 for cf/c 0.30 and deflection of 20deg

    Chalpha_theory = -0.5; % figure 10.63
    Chalpha_prima_Chalpha_theory = 0.775; % figure 10.63
    Chalpha_prima = Chalpha_prima_Chalpha_theory*Chalpha_theory;
    Chalpha_2prima = Chalpha_prima;
    Chalpha_bal_Chalpha_2prima = 0.8; % figure 10.65 - for round nose NACA 0015
    Chalpha_bal = Chalpha_2prima*Chalpha_bal_Chalpha_2prima;
    Chalpha_M = Chalpha_bal/sqrt(1-Mach^2);
    fact2 = alpha_delta*Chalpha_M*(2*Lambda_c4_VTP)/(AR_w1 + 2*Lambda_c4_VTP);

    % Third term
    DeltaChdelta_adim = 0.05; % Fig 10.78a for AR = 9.7 Cessna 208
    Cldelta = 0.45;% Fig 8.14
    B2 = 1; % Fig 10.77b
    Kdelta_i = 1.75; % Fig 10.78
    Kdelta_o = 4.3;
    Kdelta = Kdelta_i*(1-K_y1_rudder_VTP) - Kdelta_o*(1-K_y2_rudder_VTP)/(K_y2_rudder_VTP-K_y1_rudder_VTP);
    DeltaChdelta = DeltaChdelta_adim*Cldelta*B2*Kdelta*sweep_factors;
    fact3 = DeltaChdelta;
    Ch_delta_rt = sweep_factors*(fact1 + fact2 + fact3);
else
    Ch_delta_rt = 0;
end

Stab_Der.Ch_delta_a = Ch_delta_a;
Stab_Der.Ch_delta_r = Ch_delta_r;
Stab_Der.Ch_delta_at = Ch_delta_at;
Stab_Der.Ch_delta_rt = Ch_delta_rt;

Stab_Der_parts.Ch_delta_a = Ch_delta_a;
Stab_Stab_Der_partsDer.Ch_delta_r = Ch_delta_r;
Stab_Der_parts.Ch_delta_at = Ch_delta_at;
Stab_Der_parts.Ch_delta_rt = Ch_delta_rt;
