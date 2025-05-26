%% Function that writes into the excel the selected structure of data
function [Geo,Aero,Stab,Perfo,Misc,name_geo,name_aero,name_stab,name_perfo,name_misc] = Write_DATA_complete_AC_sheet1(Storing_DATA,...
    M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filename,Sheet_AC1,prefix)

DUMMY = 1;

Storing_GEO_DATA_1 = Storing_DATA.Storing_GEO_DATA_1;
Storing_WEIGHT_DATA_1 = Storing_DATA.Storing_WEIGHT_DATA_1;
Storing_AERO_DATA_1 = Storing_DATA.Storing_AERO_DATA_1;
Storing_PROPULSION_DATA_1 = Storing_DATA.Storing_PROPULSION_DATA_1;
% Performance
Storing_PERFORMANCE_DATA_1 = Storing_DATA.Storing_PERFORMANCE_DATA_1;
Storing_PERFORMANCE_DATA_21 = Storing_DATA.Storing_PERFORMANCE_DATA_21;
Storing_PERFORMANCE_DATA_22 = Storing_DATA.Storing_PERFORMANCE_DATA_22;
Storing_PERFORMANCE_DATA_23 = Storing_DATA.Storing_PERFORMANCE_DATA_23;
Storing_PERFORMANCE_DATA_24 = Storing_DATA.Storing_PERFORMANCE_DATA_24;
Storing_PERFORMANCE_DATA_25 = Storing_DATA.Storing_PERFORMANCE_DATA_25;
Storing_PERFORMANCE_DATA_26 = Storing_DATA.Storing_PERFORMANCE_DATA_26;
Storing_PERFORMANCE_DATA_27 = Storing_DATA.Storing_PERFORMANCE_DATA_27;
% Stability
Storing_STABILITY_DATA_1 = Storing_DATA.Storing_STABILITY_DATA_1;
Storing_STABILITY_DATA_2 = Storing_DATA.Storing_STABILITY_DATA_2;
Storing_STABILITY_DATA_2B = Storing_DATA.Storing_STABILITY_DATA_2B;
Storing_STABILITY_DATA_3 = Storing_DATA.Storing_STABILITY_DATA_3;
Storing_STABILITY_DATA_4A = Storing_DATA.Storing_STABILITY_DATA_4A;
Storing_STABILITY_DATA_4B = Storing_DATA.Storing_STABILITY_DATA_4B;
Storing_STABILITY_DATA_4C = Storing_DATA.Storing_STABILITY_DATA_4C;
Storing_STABILITY_DATA_4D = Storing_DATA.Storing_STABILITY_DATA_4D;
Storing_STABILITY_DATA_5 = Storing_DATA.Storing_STABILITY_DATA_5;

AC_type = OUTPUT_read_XLSX.AC_Data_flags.AC_type;

% Selects the data according to the type of aircraft
switch AC_type
    case 1 % AC_type = 1 - flying wing
        AR_w2 = 0;
        b_w2 = 0;
        cR_w2 = 0;
        cT_w2 = 0;
        S_w2 = 0;
        cmac_w2 = 0;
        xbar_w2 = 0;
        ybar_w2 = 0;
        x_xbar_w2 = 0;
        l_arm_w1w2 = 0;
        AR_w3 = 0;
        b_w3 = 0;
        cR_w3 = 0;
        cT_w3 = 0;
        S_w3 = 0;
        cmac_w3 = 0;
        xbar_w3 = 0;
        ybar_w3 = 0;
        x_xbar_w3 = 0;
        l_arm_w1w3 = 0;
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        AR_w2 = Storing_GEO_DATA_1.Geo_tier.AR_HTP;
        b_w2 = Storing_GEO_DATA_1.Geo_tier.b_HTP;
        cR_w2 = Storing_GEO_DATA_1.Geo_tier.cR_HTP;
        cT_w2 = Storing_GEO_DATA_1.Geo_tier.cT_HTP;
        S_w2 = Storing_GEO_DATA_1.Geo_tier.S_HTP;
        cmac_w2 = Storing_GEO_DATA_1.Geo_tier.cmac_HTP;
        xbar_w2 = Storing_GEO_DATA_1.Geo_tier.xbar_HTP;
        ybar_w2 = Storing_GEO_DATA_1.Geo_tier.ybar_HTP;
        x_xbar_w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_HTP;
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_HTP - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;
        AR_w3 = 0;
        b_w3 = 0;
        cR_w3 = 0;
        cT_w3 = 0;
        S_w3 = 0;
        cmac_w3 = 0;
        xbar_w3 = 0;
        ybar_w3 = 0;
        x_xbar_w3 = 0;
        l_arm_w1w3 = 0;
    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        AR_w2 = Storing_GEO_DATA_1.Geo_tier.AR_can;
        b_w2 = Storing_GEO_DATA_1.Geo_tier.b_can;
        cR_w2 = Storing_GEO_DATA_1.Geo_tier.cR_can;
        cT_w2 = Storing_GEO_DATA_1.Geo_tier.cT_can;
        S_w2 = Storing_GEO_DATA_1.Geo_tier.S_can;
        cmac_w2 = Storing_GEO_DATA_1.Geo_tier.cmac_can;
        xbar_w2 = Storing_GEO_DATA_1.Geo_tier.xbar_can;
        ybar_w2 = Storing_GEO_DATA_1.Geo_tier.ybar_can;
        x_xbar_w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_can;
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w1 - Storing_GEO_DATA_1.Geo_tier.x_xbar_HTP;
        AR_w3 = Storing_GEO_DATA_1.Geo_tier.AR_HTP;
        b_w3 = Storing_GEO_DATA_1.Geo_tier.b_HTP;
        cR_w3 = Storing_GEO_DATA_1.Geo_tier.cR_HTP;
        cT_w3 = Storing_GEO_DATA_1.Geo_tier.cT_HTP;
        S_w3 = Storing_GEO_DATA_1.Geo_tier.S_HTP;
        cmac_w3 = Storing_GEO_DATA_1.Geo_tier.cmac_HTP;
        xbar_w3 = Storing_GEO_DATA_1.Geo_tier.xbar_HTP;
        ybar_w3 = Storing_GEO_DATA_1.Geo_tier.ybar_HTP;
        x_xbar_w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_HTP;
        l_arm_w1w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_HTP - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        AR_w2 = Storing_GEO_DATA_1.Geo_tier.AR_vee;
        b_w2 = Storing_GEO_DATA_1.Geo_tier.b_vee;
        cR_w2 = Storing_GEO_DATA_1.Geo_tier.cR_vee;
        cT_w2 = Storing_GEO_DATA_1.Geo_tier.cT_vee;
        S_w2 = Storing_GEO_DATA_1.Geo_tier.S_vee;
        cmac_w2 = Storing_GEO_DATA_1.Geo_tier.cmac_vee;
        xbar_w2 = Storing_GEO_DATA_1.Geo_tier.xbar_vee;
        ybar_w2 = Storing_GEO_DATA_1.Geo_tier.ybar_vee;
        x_xbar_w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee;
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;
        AR_w3 = 0;
        b_w3 = 0;
        cR_w3 = 0;
        cT_w3 = 0;
        S_w3 = 0;
        cmac_w3 = 0;
        xbar_w3 = 0;
        ybar_w3 = 0;
        x_xbar_w3 = 0;
        l_arm_w1w3 = 0;
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        AR_w2 = Storing_GEO_DATA_1.Geo_tier.AR_can;
        b_w2 = Storing_GEO_DATA_1.Geo_tier.b_can;
        cR_w2 = Storing_GEO_DATA_1.Geo_tier.cR_can;
        cT_w2 = Storing_GEO_DATA_1.Geo_tier.cT_can;
        S_w2 = Storing_GEO_DATA_1.Geo_tier.S_can;
        cmac_w2 = Storing_GEO_DATA_1.Geo_tier.cmac_can;
        xbar_w2 = Storing_GEO_DATA_1.Geo_tier.xbar_can;
        ybar_w2 = Storing_GEO_DATA_1.Geo_tier.ybar_can;
        x_xbar_w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_can;

        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w1 - Storing_GEO_DATA_1.Geo_tier.x_xbar_vee;
        AR_w3 = Storing_GEO_DATA_1.Geo_tier.AR_vee;
        b_w3 = Storing_GEO_DATA_1.Geo_tier.b_vee;
        cR_w3 = Storing_GEO_DATA_1.Geo_tier.cR_vee;
        cT_w3 = Storing_GEO_DATA_1.Geo_tier.cT_vee;
        S_w3 = Storing_GEO_DATA_1.Geo_tier.S_vee;
        cmac_w3 = Storing_GEO_DATA_1.Geo_tier.cmac_vee;
        xbar_w3 = Storing_GEO_DATA_1.Geo_tier.xbar_vee;
        ybar_w3 = Storing_GEO_DATA_1.Geo_tier.ybar_vee;
        x_xbar_w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee;
        l_arm_w1w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;

    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        AR_w2 = Storing_GEO_DATA_1.Geo_tier.AR_can;
        b_w2 = Storing_GEO_DATA_1.Geo_tier.b_can;
        cR_w2 = Storing_GEO_DATA_1.Geo_tier.cR_can;
        cT_w2 = Storing_GEO_DATA_1.Geo_tier.cT_can;
        S_w2 = Storing_GEO_DATA_1.Geo_tier.S_can;
        cmac_w2 = Storing_GEO_DATA_1.Geo_tier.cmac_can;
        xbar_w2 = Storing_GEO_DATA_1.Geo_tier.xbar_can;
        ybar_w2 = Storing_GEO_DATA_1.Geo_tier.ybar_can;
        x_xbar_w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_can;
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w1 - Storing_GEO_DATA_1.Geo_tier.x_xbar_can;
        AR_w3 = 0;
        b_w3 = 0;
        cR_w3 = 0;
        cT_w3 = 0;
        S_w3 = 0;
        cmac_w3 = 0;
        xbar_w3 = 0;
        ybar_w3 = 0;
        x_xbar_w3 = 0;
        l_arm_w1w3 = 0;
    case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        AR_w2 = Storing_GEO_DATA_1.Geo_tier.AR_vee;
        b_w2 = Storing_GEO_DATA_1.Geo_tier.b_vee;
        cR_w2 = Storing_GEO_DATA_1.Geo_tier.cR_vee;
        cT_w2 = Storing_GEO_DATA_1.Geo_tier.cT_vee;
        S_w2 = Storing_GEO_DATA_1.Geo_tier.S_vee;
        cmac_w2 = Storing_GEO_DATA_1.Geo_tier.cmac_vee;
        xbar_w2 = Storing_GEO_DATA_1.Geo_tier.xbar_vee;
        ybar_w2 = Storing_GEO_DATA_1.Geo_tier.ybar_vee;
        x_xbar_w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee;
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;
        AR_w3 = 0;
        b_w3 = 0;
        cR_w3 = 0;
        cT_w3 = 0;
        S_w3 = 0;
        cmac_w3 = 0;
        xbar_w3 = 0;
        ybar_w3 = 0;
        x_xbar_w3 = 0;
        l_arm_w1w3 = 0;

    case 8 % AC_type = 8 - 2 surface: wing + V-tail + V-tail
        AR_w2 = Storing_GEO_DATA_1.Geo_tier.AR_vee;
        b_w2 = Storing_GEO_DATA_1.Geo_tier.b_vee;
        cR_w2 = Storing_GEO_DATA_1.Geo_tier.cR_vee;
        cT_w2 = Storing_GEO_DATA_1.Geo_tier.cT_vee;
        S_w2 = Storing_GEO_DATA_1.Geo_tier.S_vee;
        cmac_w2 = Storing_GEO_DATA_1.Geo_tier.cmac_vee;
        xbar_w2 = Storing_GEO_DATA_1.Geo_tier.xbar_vee;
        ybar_w2 = Storing_GEO_DATA_1.Geo_tier.ybar_vee;
        x_xbar_w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee;
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;

        AR_w3 = Storing_GEO_DATA_1.Geo_tier.AR_vee2;
        b_w3 = Storing_GEO_DATA_1.Geo_tier.b_vee2;
        cR_w3 = Storing_GEO_DATA_1.Geo_tier.cR_vee2;
        cT_w3 = Storing_GEO_DATA_1.Geo_tier.cT_vee2;
        S_w3 = Storing_GEO_DATA_1.Geo_tier.S_vee2;
        cmac_w3 = Storing_GEO_DATA_1.Geo_tier.cmac_vee2;
        xbar_w3 = Storing_GEO_DATA_1.Geo_tier.xbar_vee2;
        ybar_w3 = Storing_GEO_DATA_1.Geo_tier.ybar_vee2;
        x_xbar_w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee2;
        l_arm_w1w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_vee2 - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;
end

Geo = [Storing_GEO_DATA_1.Geo_tier.AR_w1;...
    Storing_GEO_DATA_1.Geo_tier.b_w1;...
    Storing_GEO_DATA_1.Geo_tier.cR_w1;...
    Storing_GEO_DATA_1.Geo_tier.cT_w1;
    Storing_GEO_DATA_1.Geo_tier.S_w1;...
    Storing_GEO_DATA_1.Geo_tier.cmac_w1;...
    Storing_GEO_DATA_1.Geo_tier.xbar_w1;...
    Storing_GEO_DATA_1.Geo_tier.ybar_w1;...
    Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;...
    AR_w2;...
    b_w2;...
    cR_w2;...
    cT_w2;...
    S_w2;...
    cmac_w2;...
    xbar_w2;...
    ybar_w2;...
    x_xbar_w2;...
    l_arm_w1w2;...
    AR_w3;...
    b_w3;...
    cR_w3;...
    cT_w3;...
    S_w3;...
    cmac_w3;...
    xbar_w3;...
    ybar_w3;...
    x_xbar_w3;...
    l_arm_w1w3] ;

switch AC_type
    case 1 % AC_type = 1 - flying wing
        CL_0_w2_CR = 0;
        CL_alpha_w2_CR = 0;
        CM_0_w2_CR = 0;
        CL_max_w2_CR = 0;
        alpha_max_w2_CR = 0;
        E_w2_max_E_prop = 0;

        CL_0_w3_CR = 0;
        CL_alpha_w3_CR = 0;
        CM_0_w3_CR = 0;
        CL_max_w3_CR = 0;
        alpha_max_w3_CR = 0;
        E_w3_max_E_prop = 0;
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        CL_0_w2_CR = Storing_AERO_DATA_1.Aero.CL_0_HTP_CR;
        CL_alpha_w2_CR = Storing_AERO_DATA_1.Aero.CL_alpha_HTP_CR;
        CM_0_w2_CR = Storing_AERO_DATA_1.Aero.CM_0_HTP_CR;
        CL_max_w2_CR = Storing_AERO_DATA_1.Aero.CL_max_HTP_CR;
        alpha_max_w2_CR = Storing_AERO_DATA_1.Aero.alpha_max_HTP_CR;
        E_w2_max_E_prop = Storing_AERO_DATA_1.Aero.E_HTP_max_E_prop;

        CL_0_w3_CR = 0;
        CL_alpha_w3_CR = 0;
        CM_0_w3_CR = 0;
        CL_max_w3_CR = 0;
        alpha_max_w3_CR = 0;
        E_w3_max_E_prop = 0;
    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        CL_0_w2_CR = Storing_AERO_DATA_1.Aero.CL_0_can_CR;
        CL_alpha_w2_CR = Storing_AERO_DATA_1.Aero.CL_alpha_can_CR;
        CM_0_w2_CR = Storing_AERO_DATA_1.Aero.CM_0_can_CR;
        CL_max_w2_CR = Storing_AERO_DATA_1.Aero.CL_max_can_CR;
        alpha_max_w2_CR = Storing_AERO_DATA_1.Aero.alpha_max_can_CR;
        E_w2_max_E_prop = Storing_AERO_DATA_1.Aero.E_can_max_E_prop;
        CL_0_w3_CR = Storing_AERO_DATA_1.Aero.CL_0_HTP_CR;
        CL_alpha_w3_CR = Storing_AERO_DATA_1.Aero.CL_alpha_HTP_CR;
        CM_0_w3_CR = Storing_AERO_DATA_1.Aero.CM_0_HTP_CR;
        CL_max_w3_CR = Storing_AERO_DATA_1.Aero.CL_max_HTP_CR;
        alpha_max_w3_CR = Storing_AERO_DATA_1.Aero.alpha_max_HTP_CR;
        E_w3_max_E_prop = Storing_AERO_DATA_1.Aero.E_HTP_max_E_prop;
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        CL_0_w2_CR = Storing_AERO_DATA_1.Aero.CL_0_vee_CR;
        CL_alpha_w2_CR = Storing_AERO_DATA_1.Aero.CL_alpha_vee_CR;
        CM_0_w2_CR = Storing_AERO_DATA_1.Aero.CM_0_vee_CR;
        CL_max_w2_CR = Storing_AERO_DATA_1.Aero.CL_max_vee_CR;
        alpha_max_w2_CR = Storing_AERO_DATA_1.Aero.alpha_max_vee_CR;
        E_w2_max_E_prop = Storing_AERO_DATA_1.Aero.E_vee_max_E_prop;

        CL_0_w3_CR = 0;
        CL_alpha_w3_CR = 0;
        CM_0_w3_CR = 0;
        CL_max_w3_CR = 0;
        alpha_max_w3_CR = 0;
        E_w3_max_E_prop = 0;
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        CL_0_w2_CR = Storing_AERO_DATA_1.Aero.CL_0_can_CR;
        CL_alpha_w2_CR = Storing_AERO_DATA_1.Aero.CL_alpha_can_CR;
        CM_0_w2_CR = Storing_AERO_DATA_1.Aero.CM_0_can_CR;
        CL_max_w2_CR = Storing_AERO_DATA_1.Aero.CL_max_can_CR;
        alpha_max_w2_CR = Storing_AERO_DATA_1.Aero.alpha_max_can_CR;
        E_w2_max_E_prop = Storing_AERO_DATA_1.Aero.E_can_max_E_prop;

        CL_0_w3_CR = Storing_AERO_DATA_1.Aero.CL_0_vee_CR;
        CL_alpha_w3_CR = Storing_AERO_DATA_1.Aero.CL_alpha_vee_CR;
        CM_0_w3_CR = Storing_AERO_DATA_1.Aero.CM_0_vee_CR;
        CL_max_w3_CR = Storing_AERO_DATA_1.Aero.CL_max_vee_CR;
        alpha_max_w3_CR = Storing_AERO_DATA_1.Aero.alpha_max_vee_CR;
        E_w3_max_E_prop = Storing_AERO_DATA_1.Aero.E_vee_max_E_prop;
    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        CL_0_w2_CR = Storing_AERO_DATA_1.Aero.CL_0_can_CR;
        CL_alpha_w2_CR = Storing_AERO_DATA_1.Aero.CL_alpha_can_CR;
        CM_0_w2_CR = Storing_AERO_DATA_1.Aero.CM_0_can_CR;
        CL_max_w2_CR = Storing_AERO_DATA_1.Aero.CL_max_can_CR;
        alpha_max_w2_CR = Storing_AERO_DATA_1.Aero.alpha_max_can_CR;
        E_w2_max_E_prop = Storing_AERO_DATA_1.Aero.E_can_max_E_prop;
        CL_0_w3_CR = 0;
        CL_alpha_w3_CR = 0;
        CM_0_w3_CR = 0;
        CL_max_w3_CR = 0;
        alpha_max_w3_CR = 0;
        E_w3_max_E_prop = 0;
    case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        CL_0_w2_CR = Storing_AERO_DATA_1.Aero.CL_0_vee_CR;
        CL_alpha_w2_CR = Storing_AERO_DATA_1.Aero.CL_alpha_vee_CR;
        CM_0_w2_CR = Storing_AERO_DATA_1.Aero.CM_0_vee_CR;
        CL_max_w2_CR = Storing_AERO_DATA_1.Aero.CL_max_vee_CR;
        alpha_max_w2_CR = Storing_AERO_DATA_1.Aero.alpha_max_vee_CR;
        E_w2_max_E_prop = Storing_AERO_DATA_1.Aero.E_vee_max_E_prop;

        CL_0_w3_CR = 0;
        CL_alpha_w3_CR = 0;
        CM_0_w3_CR = 0;
        CL_max_w3_CR = 0;
        alpha_max_w3_CR = 0;
        E_w3_max_E_prop = 0;
    case 8 % AC_type = 8 - 2 surface: wing + V-tail + V-tail
        CL_0_w2_CR = Storing_AERO_DATA_1.Aero.CL_0_vee_CR;
        CL_alpha_w2_CR = Storing_AERO_DATA_1.Aero.CL_alpha_vee_CR;
        CM_0_w2_CR = Storing_AERO_DATA_1.Aero.CM_0_vee_CR;
        CL_max_w2_CR = Storing_AERO_DATA_1.Aero.CL_max_vee_CR;
        alpha_max_w2_CR = Storing_AERO_DATA_1.Aero.alpha_max_vee_CR;
        E_w2_max_E_prop = Storing_AERO_DATA_1.Aero.E_vee_max_E_prop;

        CL_0_w3_CR = Storing_AERO_DATA_1.Aero.CL_0_vee2_CR;
        CL_alpha_w3_CR = Storing_AERO_DATA_1.Aero.CL_alpha_vee2_CR;
        CM_0_w3_CR = Storing_AERO_DATA_1.Aero.CM_0_vee2_CR;
        CL_max_w3_CR = Storing_AERO_DATA_1.Aero.CL_max_vee2_CR;
        alpha_max_w3_CR = Storing_AERO_DATA_1.Aero.alpha_max_vee2_CR;
        E_w3_max_E_prop = Storing_AERO_DATA_1.Aero.E_vee2_max_E_prop;

end

Aero = [Storing_AERO_DATA_1.Aero.CL_0_w1_CR;...
    Storing_AERO_DATA_1.Aero.CL_alpha_w1_CR;...
    Storing_AERO_DATA_1.Aero.CM_0_w1_CR;...
    Storing_AERO_DATA_1.Aero.CL_max_w1_CR;...
    Storing_AERO_DATA_1.Aero.alpha_max_w1_CR;...

    Storing_AERO_DATA_1.Aero.alpha_w1_CR_ope;...
    Storing_AERO_DATA_1.Aero.CL_w1_CR_ope;...
    Storing_AERO_DATA_1.Aero.E_w1_max_R_jet;...
    Storing_AERO_DATA_1.Aero.E_w1_max_E_jet;...
    Storing_AERO_DATA_1.Aero.E_w1_max_E_prop;...

    CL_0_w2_CR;...
    CL_alpha_w2_CR;...
    CM_0_w2_CR;...
    CL_max_w2_CR;...
    alpha_max_w2_CR;...
    E_w2_max_E_prop;...
    CL_0_w3_CR;...
    CL_alpha_w3_CR;...
    CM_0_w3_CR;...
    CL_max_w3_CR;...
    alpha_max_w3_CR;...
    E_w3_max_E_prop;...
    Storing_AERO_DATA_1.Aero.Polar.C_D0;... % FLOW
    Storing_AERO_DATA_1.Aero.Polar.C_D1;...
    Storing_AERO_DATA_1.Aero.Polar.C_D2;...
    Storing_AERO_DATA_1.Aero_TH.CD0; % FUSE FLOW and CBM
    Storing_AERO_DATA_1.Aero_TH.CD1;
    Storing_AERO_DATA_1.Aero_TH.CD2;
    Storing_AERO_DATA_1.Aero_TH.CD0_CBM;... % CBM
    Storing_AERO_DATA_1.Aero_TH.CD1_CBM;...
    Storing_AERO_DATA_1.Aero_TH.CD2_CBM];

if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
    switch AC_type
        case 1 % AC_type = 1 - flying wing
            CL_w2 = 0;
            CL_w3 = 0;

            CL0_w1_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_w1_e_corrected;
            CL0_w2_e_corrected = 0;
            CL0_w3_e_corrected = 0;

            CLalpha_w1_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_w1_e_pw;
            CLalpha_w2_e_pw = 0;
            CLalpha_w3_e_pw = 0;

            CM_0_w1_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_w1;
            CM_0_w2_wrt_XCG = 0;
            CM_0_w3_wrt_XCG = 0;

        case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_HTP;
            CL_w3 = 0;

            CL0_w1_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_w1_e_corrected;
            CL0_w2_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_HTP_e_corrected;
            CL0_w3_e_corrected = 0;

            CLalpha_w1_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_w1_e_pw;
            CLalpha_w2_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_HTP_e_pw;
            CLalpha_w3_e_pw = 0;

            CM_0_w1_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_w1;
            CM_0_w2_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_HTP;
            CM_0_w3_wrt_XCG = 0;

        case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_can;
            CL_w3 = Storing_STABILITY_DATA_1.Trim_ITER.CL_HTP;

            CL0_w1_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_w1_e_corrected;
            CL0_w2_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_can_e_corrected;
            CL0_w3_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_HTP_e_corrected;

            CLalpha_w1_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_w1_e_pw;
            CLalpha_w2_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_can_e_pw;
            CLalpha_w3_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_HTP_e_pw;

            CM_0_w1_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_w1;
            CM_0_w2_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_can;
            CM_0_w3_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_HTP;

        case 4 % AC_type = 4 - 2 surface: wing + V-tail
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_vee;
            CL_w3 = 0;

            CL0_w1_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_w1_e_corrected;
            CL0_w2_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_vee_e_corrected;
            CL0_w3_e_corrected = 0;

            CLalpha_w1_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_w1_e_pw;
            CLalpha_w2_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_vee_e_pw;
            CLalpha_w3_e_pw = 0;

            CM_0_w1_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_w1;
            CM_0_w2_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_vee;
            CM_0_w3_wrt_XCG = 0;

        case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_can;
            CL_w3 = Storing_STABILITY_DATA_1.Trim_ITER.CL_vee;

            CL0_w1_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_w1_e_corrected;
            CL0_w2_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_can_e_corrected;
            CL0_w3_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_vee_e_corrected;

            CLalpha_w1_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_w1_e_pw;
            CLalpha_w2_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_can_e_pw;
            CLalpha_w3_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_vee_e_pw;

            CM_0_w1_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_w1;
            CM_0_w2_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_can;
            CM_0_w3_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_vee;

        case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_can;
            CL_w3 = 0;

            CL0_w1_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_w1_e_corrected;
            CL0_w2_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_can_e_corrected;
            CL0_w3_e_corrected = 0;

            CLalpha_w1_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_w1_e_pw;
            CLalpha_w2_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_can_e_pw;
            CLalpha_w3_e_pw = 0;

            CM_0_w1_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_w1;
            CM_0_w2_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_can;
            CM_0_w3_wrt_XCG = 0;

        case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP

            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_vee;
            CL_w3 = 0;

            CL0_w1_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_w1_e_corrected;
            CL0_w2_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_vee_e_corrected;
            CL0_w3_e_corrected = 0;

            CLalpha_w1_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_w1_e_pw;
            CLalpha_w2_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_vee_e_pw;
            CLalpha_w3_e_pw = 0;

            CM_0_w1_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_w1;
            CM_0_w2_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_vee;
            CM_0_w3_wrt_XCG = 0;

        case 8 % AC_type = 8 - 2 surface: wing + V-tail + V-tail
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_vee;
            CL_w3 = Storing_STABILITY_DATA_1.Trim_ITER.CL_vee2;

            CL0_w1_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_w1_e_corrected;
            CL0_w2_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_vee_e_corrected;
            CL0_w3_e_corrected = Storing_STABILITY_DATA_1.Stab_Der_parts.CL0_vee2_e_corrected;

            CLalpha_w1_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_w1_e_pw;
            CLalpha_w2_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_vee_e_pw;
            CLalpha_w3_e_pw = Storing_STABILITY_DATA_1.Stab_Der_parts.CLalpha_vee2_e_pw;

            CM_0_w1_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_w1;
            CM_0_w2_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_vee;
            CM_0_w3_wrt_XCG = Storing_STABILITY_DATA_1.Stab_Der_parts.CM_0_vee2;
    end
    % Selects the variable to be stotred depending if study is conducted
    if  OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
        data_STUDY_Long_dyn = [real(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole1);...
            imag(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole1);...
            real(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole2);...
            imag(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole2);...
            real(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole3);...
            imag(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole3);...
            real(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole4);...
            imag(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole4);...
            % Longitudinal modes
            Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.SP_wn;...
            Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.SP_damp;...
            Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.PH_wn;...
            Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.PH_damp];
    else
        data_STUDY_Long_dyn = [0;0;0;0;0;0;0;0;0;0;0;0]
    end
    % Selects the variable to be stotred depending if study is conducted
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn ==1
        % Lateral-Directional Dynamics
        data_STUDY_LatDir_dyn = [real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole1);...
            imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole1);...
            real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole2);...
            imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole2);...
            real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole3);...
            imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole3);...
            real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole4);...
            imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole4);...
            real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole5);...
            imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole5);...
            Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.DR_wn;...
            Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.DR_damp;...
            Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.ROL_T2;...
            Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.ESP_T2];
    else
        data_STUDY_LatDir_dyn = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    end

    Stab = [Storing_STABILITY_DATA_1.TRIM_RESULTS.m_TOW;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.V;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.X_NP;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.x_XCG;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.SM;...
        % Trim
        Storing_STABILITY_DATA_1.TRIM_RESULTS.CL_alpha_ac;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.CM_alpha_ac;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.CL0_ac;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.CM0_ac;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.CL_delta;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.CM_delta;...
        Storing_STABILITY_DATA_1.Trim_ITER.CL_w1;...
        CL_w2;...
        CL_w3;...
        Storing_STABILITY_DATA_1.Trim_ITER.CL;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.trim_alpha_deg;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.trim_delta_e_deg;...
        Storing_STABILITY_DATA_1.TRIM_RESULTS.delta_T;...
        % Values of corrected aerodynamic properties
        CL0_w1_e_corrected;...
        CL0_w2_e_corrected;...
        CL0_w3_e_corrected;...
        CLalpha_w1_e_pw;...
        CLalpha_w2_e_pw;...
        CLalpha_w3_e_pw;...
        CM_0_w1_wrt_XCG;...
        CM_0_w2_wrt_XCG;...
        CM_0_w3_wrt_XCG;...
        % Longitudinal Stability derivatives
        % theta derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_teta;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_teta;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_teta;...
        % alpha derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_alpha;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_alpha_ac;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_alpha_ac;...
        % u derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CDu;...
        Storing_STABILITY_DATA_1.Stab_Der.CLu;...
        Storing_STABILITY_DATA_1.Stab_Der.CMu;...
        % q derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CDq;...
        Storing_STABILITY_DATA_1.Stab_Der.CLq;...
        Storing_STABILITY_DATA_1.Stab_Der.CMq;...
        % alpha-dot derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_alphapunto;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_alphapunto;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_alphapunto;...
        % Longitudinal control derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_delta_e;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_delta_e;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_delta_e;...
        % Longitudinal control derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_delta_elevon;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_delta_elevon;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_delta_elevon;...
        % Longitudinal control derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_delta_rv;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_delta_rv;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_delta_rv;...
        % Longitudinal control derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_delta_rv2;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_delta_rv2;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_delta_rv2;...
        % Longitudinal control derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_delta_can;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_delta_can;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_delta_can;...
        % Longitudinal Elevator Trim Tab
        Storing_STABILITY_DATA_1.Stab_Der.CD_delta_e_Tab;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_delta_e_Tab;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_delta_e_Tab;...
        % Longitudinal-Directional Dynamics
        data_STUDY_Long_dyn;...
        % if  OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
        % real(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole1);...
        %     imag(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole1);...
        %     real(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole2);...
        %     imag(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole2);...
        %     real(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole3);...
        %     imag(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole3);...
        %     real(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole4);...
        %     imag(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole4);...
        %     real(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole4);...
        %     imag(Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.pole4);...
        %     % Longitudinal modes
        % Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.SP_wn;...
        %     Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.SP_damp;...
        %     Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.PH_wn;...
        %     Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.PH_damp;...
        % else
        % 0;...
        %     0;...
        %     0;...
        %     0;...
        %     0;...
        %     0;...
        %     0;...
        %     0;...
        %     0;...
        %     0;...
        %     % Longitudinal modes
        % 0;...
        %     0;...
        %     0;...
        %     0;...
        % end
        %
        % Lateral Stability Derivatives
        % beta derivatives
        Storing_STABILITY_DATA_1.Stab_Der.Cyb;...
        Storing_STABILITY_DATA_1.Stab_Der.Clb;...
        Storing_STABILITY_DATA_1.Stab_Der.Cnb;...
        % roll rate derivatives
        Storing_STABILITY_DATA_1.Stab_Der.Cyp;...
        Storing_STABILITY_DATA_1.Stab_Der.Clp;...
        Storing_STABILITY_DATA_1.Stab_Der.Cnp;...
        % yaw rate derivatives
        Storing_STABILITY_DATA_1.Stab_Der.Cyr;...
        Storing_STABILITY_DATA_1.Stab_Der.Clr;...
        Storing_STABILITY_DATA_1.Stab_Der.Cnr;...
        % y-dot derivatives
        Storing_STABILITY_DATA_1.Stab_Der.Cybpunto;...
        Storing_STABILITY_DATA_1.Stab_Der.Clbpunto;...
        Storing_STABILITY_DATA_1.Stab_Der.Cnbpunto;...
        % Lateral-Directional control derivatives
        Storing_STABILITY_DATA_1.Stab_Der.Cydeltaa;...
        Storing_STABILITY_DATA_1.Stab_Der.Cldeltaa;...
        Storing_STABILITY_DATA_1.Stab_Der.Cndeltaa;...
        % Rudder
        Storing_STABILITY_DATA_1.Stab_Der.Cydeltar;...
        Storing_STABILITY_DATA_1.Stab_Der.Cldeltar;...
        Storing_STABILITY_DATA_1.Stab_Der.Cndeltar;...
        % Rudder-vator
        Storing_STABILITY_DATA_1.Stab_Der.Cydeltarv;...
        Storing_STABILITY_DATA_1.Stab_Der.Cldeltarv;...
        Storing_STABILITY_DATA_1.Stab_Der.Cndeltarv;...
        % Rudder-vator
        Storing_STABILITY_DATA_1.Stab_Der.Cydeltarv2;...
        Storing_STABILITY_DATA_1.Stab_Der.Cldeltarv2;...
        Storing_STABILITY_DATA_1.Stab_Der.Cndeltarv2;...
        % Aileron Trim Tab
        Storing_STABILITY_DATA_1.Stab_Der.Cy_delta_a_Tab;...
        Storing_STABILITY_DATA_1.Stab_Der.Cl_delta_a_Tab;...
        Storing_STABILITY_DATA_1.Stab_Der.Cn_delta_a_Tab;...
        % Rudder Trim Tab
        Storing_STABILITY_DATA_1.Stab_Der.Cy_delta_r_Tab;...
        Storing_STABILITY_DATA_1.Stab_Der.Cl_delta_r_Tab;...
        Storing_STABILITY_DATA_1.Stab_Der.Cn_delta_r_Tab;...
        % Lateral-Directional Dynamics
        data_STUDY_LatDir_dyn];
    % if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn ==1
    %
    %     real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole1);...
    %         imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole1);...
    %         real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole2);...
    %         imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole2);...
    %         real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole3);...
    %         imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole3);...
    %         real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole4);...
    %         imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole4);...
    %         real(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole5);...
    %         imag(Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.pole5);...
    %         Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.DR_wn;...
    %         Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.DR_damp;...
    %         Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.ROL_T2;...
    %         Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.ESP_T2];
    % else
    %     % Lateral-Directional Dynamics
    %     0;...
    %         0;...
    %         0;...
    %         0;...
    %         0;...
    %         0;...
    %         0;...
    %         0;...
    %         0;...
    %         0;...
    %         0;...
    %         0;...
    %         0;...
    %         0];
    %
    % end
else
    Stab = DUMMY;
end

if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1

    % Storing_PERFORMANCE_DATA_1.Segments = Segments;
    % Storing_PERFORMANCE_DATA_1.handles = handles;
    % Storing_PERFORMANCE_DATA_1.seg = seg;
    % Storing_PERFORMANCE_DATA_1.Weights_AP = Weights_AP;
    % Storing_PERFORMANCE_DATA_1.Total_Datos = Total_Datos;
    % Storing_PERFORMANCE_DATA_1.datos = datos;

    % Storing_PERFORMANCE_DATA_1

    
    % if  OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP == 1
    %     % Segments = Storing_PERFORMANCE_DATA_1.Segments;
    %     % handles = Storing_PERFORMANCE_DATA_1.handles;
    %     % seg = Storing_PERFORMANCE_DATA_1.seg;
    %     Weights_AP = Storing_PERFORMANCE_DATA_1.Weights_AP;
    %     Total_Datos = Storing_PERFORMANCE_DATA_1.Total_Datos;
    %     datos = Storing_PERFORMANCE_DATA_1.datos;
    % 
    %     % Resultados globales:
    %     if handles.propul(1) == 4
    %         m_bat = Weights_AP.m_bat;
    %         m_f = 0;
    %         C_fuel = 0;
    %     else
    %         m_bat = 0;
    %         m_f = Weights_AP.m_f;
    %         % C_fuel = datos.TOTAL.C_fuel;
    %         C_fuel = 0;
    % 
    %     end
    % 
    %     % DOC_PLOT = Storing_PERFORMANCE_DATA_1{1}{6}.TOTAL.DOC_PLOT;
    %     % ASM_PLOT = Storing_PERFORMANCE_DATA_1{1}{6}.TOTAL.ASM_PLOT;
    %     % CAPM_PLOT = Storing_PERFORMANCE_DATA_1{1}{6}.TOTAL.CAPM_PLOT;
    % 
    %     % Storing_PERFORMANCE_DATA_1{1} = Segments;
    %     % Storing_PERFORMANCE_DATA_1{2} = handles;
    %     % Storing_PERFORMANCE_DATA_1{3} = seg;
    %     % Storing_PERFORMANCE_DATA_1{4} = Weights_AP;
    %     % Storing_PERFORMANCE_DATA_1{5} = Total_Datos;
    %     % Storing_PERFORMANCE_DATA_1{6} = datos;
    % 
    %     Perfo = [seg.datos.Mach;...
    %         Total_Datos.tiempo_total/60;...
    %         Total_Datos.distancia_total/1000;...
    %         Weights_AP.m_TOW;...
    %         Weights_AP.m_e;...
    %         Weights_AP.m_p;...
    %         Weights_AP.m_F;...
    %         m_f;...
    %         m_bat;...
    %         C_fuel];
    % else
        Perfo = DUMMY;
    % end
else
    Perfo = DUMMY;
end

if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim_lat_Trim_TAB == 1

    % Misc = [Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.rho;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.V;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.m_TOW;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.deltaa_deg;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.deltar_deg;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.deltaaT_deg;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.deltarT_deg;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.beta_deg;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.phi_deg;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.Fa;...
    %     Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C.Fr];

    Misc = [Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.rho;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.V;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.m_TOW;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.deltaa_deg;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.deltar_deg;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.deltaaT_deg;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.deltarT_deg;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.beta_deg;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.phi_deg;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.Fa;...
        Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D.Fr];
else
    Misc = DUMMY;
end



% Automatizes the range of cells for each data-write up
% Range for Geometry write up in Excel
initial = 1;
N_elements_geo = 29;
range_geo_i = initial + 1;
range_geo_f = initial + N_elements_geo;
RANGE_geo_i = num2str(range_geo_i);
RANGE_geo_f = num2str(range_geo_f);
st_geo_i = strcat(RANGE_geo_i);
st_geo_f = strcat(RANGE_geo_f);

% Range for Aerodynamics write up in Excel
N_elements_aero = 31;
range_aero_i = range_geo_f + 1;
range_aero_f = range_geo_f + N_elements_aero;
RANGE_aero_i = num2str(range_aero_i);
RANGE_aero_f = num2str(range_aero_f);
st_aero_i = strcat(RANGE_aero_i);
st_aero_f = strcat(RANGE_aero_f);

% Range for Stability write up in Excel
N_elements_stab = 116;
range_stab_i = range_aero_f + 1;
range_stab_f = range_aero_f + N_elements_stab;
RANGE_stab_i = num2str(range_stab_i);
RANGE_stab_f = num2str(range_stab_f);
st_stab_i = strcat(RANGE_stab_i);
st_stab_f = strcat(RANGE_stab_f);

% Range for Performance write up in Excel
N_elements_perfo = 10;
range_perfo_i = range_stab_f + 1;
range_perfo_f = range_stab_f + N_elements_perfo;
RANGE_perfo_i = num2str(range_perfo_i);
RANGE_perfo_f = num2str(range_perfo_f);
st_perfo_i = strcat(RANGE_perfo_i);
st_perfo_f = strcat(RANGE_perfo_f);

% Range for Performance write up in Excel
N_elements_misc = 11;
range_misc_i = range_perfo_f + 1;
range_misc_f = range_perfo_f + N_elements_misc;
RANGE_misc_i = num2str(range_misc_i);
RANGE_misc_f = num2str(range_misc_f);
st_misc_i = strcat(RANGE_misc_i);
st_misc_f = strcat(RANGE_misc_f);

% Generates the name for the range of write up cells
name_geo   = strcat(prefix,st_geo_i,':',prefix,st_geo_f);
name_aero   = strcat(prefix,st_aero_i,':',prefix,st_aero_f);
name_stab   = strcat(prefix,st_stab_i,':',prefix,st_stab_f);
name_perfo   = strcat(prefix,st_perfo_i,':',prefix,st_perfo_f);
name_misc   = strcat(prefix,st_misc_i,':',prefix,st_misc_f);

