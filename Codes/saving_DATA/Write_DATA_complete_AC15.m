function Write_DATA_complete_AC15(Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_AERO_DATA_1,Storing_PERFORMANCE_DATA_1,Storing_STABILITY_DATA_1,...
        Storing_STABILITY_DATA_2,Storing_STABILITY_DATA_2B,Storing_STABILITY_DATA_3,Storing_STABILITY_DATA_4,Storing_STABILITY_DATA_5,...
        Storing_PROPULSION_DATA,M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX)

filename = '../Results/EMERGENTIA_v2/Results_AC15.xlsx';
Sheet_AC = 'AC17';

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
        AR_w2 = Storing_GEO_DATA_1.Geo_tier.AR_w2;
        b_w2 = Storing_GEO_DATA_1.Geo_tier.b_w2;
        cR_w2 = Storing_GEO_DATA_1.Geo_tier.cR_w2;
        cT_w2 = Storing_GEO_DATA_1.Geo_tier.cT_w2;
        S_w2 = Storing_GEO_DATA_1.Geo_tier.S_w2;
        cmac_w2 = Storing_GEO_DATA_1.Geo_tier.cmac_w2;
        xbar_w2 = Storing_GEO_DATA_1.Geo_tier.xbar_w2;
        ybar_w2 = Storing_GEO_DATA_1.Geo_tier.ybar_w2;
        x_xbar_w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w2;
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w2 - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;        
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
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w1 - Storing_GEO_DATA_1.Geo_tier.x_xbar_w2;        
        AR_w3 = Storing_GEO_DATA_1.Geo_tier.AR_w2;
        b_w3 = Storing_GEO_DATA_1.Geo_tier.b_w2;
        cR_w3 = Storing_GEO_DATA_1.Geo_tier.cR_w2;
        cT_w3 = Storing_GEO_DATA_1.Geo_tier.cT_w2;
        S_w3 = Storing_GEO_DATA_1.Geo_tier.S_w2;
        cmac_w3 = Storing_GEO_DATA_1.Geo_tier.cmac_w2;
        xbar_w3 = Storing_GEO_DATA_1.Geo_tier.xbar_w2;
        ybar_w3 = Storing_GEO_DATA_1.Geo_tier.ybar_w2;
        x_xbar_w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w2;
        l_arm_w1w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w2 - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;        
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        AR_w2 = Storing_GEO_DATA_1.Geo_tier.AR_w2;
        b_w2 = Storing_GEO_DATA_1.Geo_tier.b_w2;
        cR_w2 = Storing_GEO_DATA_1.Geo_tier.cR_w2;
        cT_w2 = Storing_GEO_DATA_1.Geo_tier.cT_w2;
        S_w2 = Storing_GEO_DATA_1.Geo_tier.S_w2;
        cmac_w2 = Storing_GEO_DATA_1.Geo_tier.cmac_w2;
        xbar_w2 = Storing_GEO_DATA_1.Geo_tier.xbar_w2;
        ybar_w2 = Storing_GEO_DATA_1.Geo_tier.ybar_w2;
        x_xbar_w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w2;
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w2 - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;        
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
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w1 - Storing_GEO_DATA_1.Geo_tier.x_xbar_w2;        
        AR_w3 = Storing_GEO_DATA_1.Geo_tier.AR_w2;
        b_w3 = Storing_GEO_DATA_1.Geo_tier.b_w2;
        cR_w3 = Storing_GEO_DATA_1.Geo_tier.cR_w2;
        cT_w3 = Storing_GEO_DATA_1.Geo_tier.cT_w2;
        S_w3 = Storing_GEO_DATA_1.Geo_tier.S_w2;
        cmac_w3 = Storing_GEO_DATA_1.Geo_tier.cmac_w2;
        xbar_w3 = Storing_GEO_DATA_1.Geo_tier.xbar_w2;
        ybar_w3 = Storing_GEO_DATA_1.Geo_tier.ybar_w2;
        x_xbar_w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w2;
        l_arm_w1w3 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w2 - Storing_GEO_DATA_1.Geo_tier.x_xbar_w1;        
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
        l_arm_w1w2 = Storing_GEO_DATA_1.Geo_tier.x_xbar_w1 - Storing_GEO_DATA_1.Geo_tier.x_xbar_w2;        
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
        CL_0_w2_CR = Storing_AERO_DATA_1.Aero.CL_0_w2_CR;
        CL_alpha_w2_CR = Storing_AERO_DATA_1.Aero.CL_alpha_w2_CR;
        CM_0_w2_CR = Storing_AERO_DATA_1.Aero.CM_0_w2_CR;
        CL_max_w2_CR = Storing_AERO_DATA_1.Aero.CL_max_w2_CR;
        alpha_max_w2_CR = Storing_AERO_DATA_1.Aero.alpha_max_w2_CR;
        E_w2_max_E_prop = Storing_AERO_DATA_1.Aero.E_w2_max_E_prop;
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
        CL_0_w3_CR = Storing_AERO_DATA_1.Aero.CL_0_w2_CR;
        CL_alpha_w3_CR = Storing_AERO_DATA_1.Aero.CL_alpha_w2_CR;
        CM_0_w3_CR = Storing_AERO_DATA_1.Aero.CM_0_w2_CR;
        CL_max_w3_CR = Storing_AERO_DATA_1.Aero.CL_max_w2_CR;
        alpha_max_w3_CR = Storing_AERO_DATA_1.Aero.alpha_max_w2_CR;
        E_w3_max_E_prop = Storing_AERO_DATA_1.Aero.E_w2_max_E_prop;
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        CL_0_w2_CR = Storing_AERO_DATA_1.Aero.CL_0_w2_CR;
        CL_alpha_w2_CR = Storing_AERO_DATA_1.Aero.CL_alpha_w2_CR;
        CM_0_w2_CR = Storing_AERO_DATA_1.Aero.CM_0_w2_CR;
        CL_max_w2_CR = Storing_AERO_DATA_1.Aero.CL_max_w2_CR;
        alpha_max_w2_CR = Storing_AERO_DATA_1.Aero.alpha_max_w2_CR;
        E_w2_max_E_prop = Storing_AERO_DATA_1.Aero.E_w2_max_E_prop;
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
        CL_0_w3_CR = Storing_AERO_DATA_1.Aero.CL_0_w2_CR;
        CL_alpha_w3_CR = Storing_AERO_DATA_1.Aero.CL_alpha_w2_CR;
        CM_0_w3_CR = Storing_AERO_DATA_1.Aero.CM_0_w2_CR;
        CL_max_w3_CR = Storing_AERO_DATA_1.Aero.CL_max_w2_CR;
        alpha_max_w3_CR = Storing_AERO_DATA_1.Aero.alpha_max_w2_CR;
        E_w3_max_E_prop = Storing_AERO_DATA_1.Aero.E_w2_max_E_prop;
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
end
Aero = [Storing_AERO_DATA_1.Aero.CL_0_w1_CR;...
    Storing_AERO_DATA_1.Aero.CL_alpha_w1_CR;...
    Storing_AERO_DATA_1.Aero.CM_0_w1_CR;...
    Storing_AERO_DATA_1.Aero.CL_max_w1_CR;...
    Storing_AERO_DATA_1.Aero.alpha_max_w1_CR;...
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
    Storing_AERO_DATA_1.Aero.Polar.C_D0;...
    Storing_AERO_DATA_1.Aero.Polar.C_D1;...
    Storing_AERO_DATA_1.Aero.Polar.C_D2;...
    Storing_AERO_DATA_1.Aero.CD0_ac;...
    Storing_AERO_DATA_1.Aero.CD1_ac;...
    Storing_AERO_DATA_1.Aero.CD2_ac;...
    Storing_AERO_DATA_1.Aero_TH.CD0_CBM;...
    Storing_AERO_DATA_1.Aero_TH.CD1_CBM;...
    Storing_AERO_DATA_1.Aero_TH.CD2_CBM];

if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
    switch AC_type
        case 1 % AC_type = 1 - flying wing
            CL_w2 = 0;
            CL_w3 = 0;
            Storing_STABILITY_DATA_1.Trim_ITER.CL_HTP;
        case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_HTP;
            CL_w3 = 0;
        case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_can;
            CL_w3 = Storing_STABILITY_DATA_1.Trim_ITER.CL_HTP;
        case 4 % AC_type = 4 - 2 surface: wing + V-tail
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_Vee;
            CL_w3 = 0;
        case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_can;
            CL_w3 = Storing_STABILITY_DATA_1.Trim_ITER.CL_Vee;
        case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
            CL_w2 = Storing_STABILITY_DATA_1.Trim_ITER.CL_can;
            CL_w3 = 0;
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
        % Longitudinal Stability derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_teta;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_teta;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_teta;...
        Storing_STABILITY_DATA_1.Stab_Der.CD_alpha;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_alpha_ac;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_alpha_ac;...
        Storing_STABILITY_DATA_1.Stab_Der.CDu;...
        Storing_STABILITY_DATA_1.Stab_Der.CLu;...
        Storing_STABILITY_DATA_1.Stab_Der.CMu;...
        Storing_STABILITY_DATA_1.Stab_Der.CDq;...
        Storing_STABILITY_DATA_1.Stab_Der.CLq;...
        Storing_STABILITY_DATA_1.Stab_Der.CMq;...
        Storing_STABILITY_DATA_1.Stab_Der.CD_alphapunto;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_alphapunto;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_alphapunto;...
        % Longitudinal control derivatives
        Storing_STABILITY_DATA_1.Stab_Der.CD_delta_e;...
        Storing_STABILITY_DATA_1.Stab_Der.CL_delta_e;...
        Storing_STABILITY_DATA_1.Stab_Der.CM_delta_e;...
        % Longitudinal-Directional Dynamics
        Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.SP_wn;...
        Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.SP_damp;...
        Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.PH_wn;...
        Storing_STABILITY_DATA_1.Stab_Dyn_Long.long.PH_damp;...
        % Lateral Stability Derivatives
        Storing_STABILITY_DATA_1.Stab_Der.Cyb;...
        Storing_STABILITY_DATA_1.Stab_Der.Clb;...
        Storing_STABILITY_DATA_1.Stab_Der.Cnb;...
        Storing_STABILITY_DATA_1.Stab_Der.Cyp;...
        Storing_STABILITY_DATA_1.Stab_Der.Clp;...
        Storing_STABILITY_DATA_1.Stab_Der.Cnp;...
        Storing_STABILITY_DATA_1.Stab_Der.Cyr;...
        Storing_STABILITY_DATA_1.Stab_Der.Clr;...
        Storing_STABILITY_DATA_1.Stab_Der.Cnr;...
        Storing_STABILITY_DATA_1.Stab_Der.Cybpunto;...
        Storing_STABILITY_DATA_1.Stab_Der.Clbpunto;...
        Storing_STABILITY_DATA_1.Stab_Der.Cnbpunto;...
        % Lateral-Directional control derivatives
        Storing_STABILITY_DATA_1.Stab_Der.Cydeltaa;...
        Storing_STABILITY_DATA_1.Stab_Der.Cldeltaa;...
        Storing_STABILITY_DATA_1.Stab_Der.Cndeltaa;...
        Storing_STABILITY_DATA_1.Stab_Der.Cydeltar;...
        Storing_STABILITY_DATA_1.Stab_Der.Cldeltar;...
        Storing_STABILITY_DATA_1.Stab_Der.Cndeltar;...
        % Lateral-Directional Dynamics
        Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.DR_wn;...
        Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.DR_damp;...
        Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.ROL_T2;...
        Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.ESP_T2];
end

if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
%     Perfo = [Storing_STABILITY_DATA_1.Stab_Dyn_LatDir.lat.ESP_T2];
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
N_elements_aero = 27;
range_aero_i = range_geo_f + 1;
range_aero_f = range_geo_f + N_elements_aero;
RANGE_aero_i = num2str(range_aero_i);
RANGE_aero_f = num2str(range_aero_f);
st_aero_i = strcat(RANGE_aero_i);
st_aero_f = strcat(RANGE_aero_f);

% Range for Stability write up in Excel
N_elements_stab = 62;
range_stab_i = range_aero_f + 1;
range_stab_f = range_aero_f + N_elements_stab;
RANGE_stab_i = num2str(range_stab_i);
RANGE_stab_f = num2str(range_stab_f);
st_stab_i = strcat(RANGE_stab_i);
st_stab_f = strcat(RANGE_stab_f);

% % Range for Performance write up in Excel
% range_perfo_i = 68;
% range_perfo_f = 72;
% RANGE_perfo_i = num2str(range_perfo_i);
% RANGE_perfo_f = num2str(range_perfo_f);
% st_perfo_i = strcat(RANGE_perfo_i);
% st_perfo_f = strcat(RANGE_perfo_f);

switch mission_actual
    case 0
        prefix = 'D';
    case 1
        prefix = 'E';
    case 2
        prefix = 'F';
    case 3
        prefix = 'G';
    case 4
        prefix = 'H';
    case 5
        prefix = 'I';
    case 6
        prefix = 'J';
    case 7
        prefix = 'k';
    case 8
        prefix = 'L';
    case 9
        prefix = 'M';
end

% Generates the name for the range of write up cells
name_geo   = strcat(prefix,st_geo_i,':',prefix,st_geo_f);
name_aero   = strcat(prefix,st_aero_i,':',prefix,st_aero_f);
name_stab   = strcat(prefix,st_stab_i,':',prefix,st_stab_f);
% name_perfo   = strcat(prefix,st_perfo_i,':',prefix,st_perfo_f);

writematrix(Geo,filename,'Sheet',Sheet_AC,'Range',name_geo);
writematrix(Aero,filename,'Sheet',Sheet_AC,'Range',name_aero);
if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY == 1
    writematrix(Stab,filename,'Sheet',Sheet_AC,'Range',name_stab)
end
% if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1
%     writematrix(Perfo,filename,'Sheet',Sheet_AC,'Range',name_perfo)
% end