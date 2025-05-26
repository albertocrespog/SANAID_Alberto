function modelo = conv_modelo_v3(AC_CONFIGURATION,TRIM_RESULTS,conditions,Trim_ITER,Stab_Der_parts, Geo_tier,...
    Prop_data,Aero_TH,afe,OUTPUT_read_XLSX,Aero, Effects,Performance,conv_UNITS,Body_Geo)

%% Conversion of variables
% x_XCG = TRIM_RESULTS.x_XCG_des;
x_XCG = conditions.x_XCG;
trim_alpha = TRIM_RESULTS.trim_alpha;
CL_w1 = Trim_ITER.CL_w1;
CL_alpha_wb_w1 = Stab_Der_parts.CL_alpha_wb_w1;
z_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;
CLalpha_w1 = Aero.CL_alpha_w1_CR;
CLalpha_w1_e = CLalpha_w1;

prop_wash_effect = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;

if AC_CONFIGURATION.VTP == 1
    % Arm from Zcg to Zac_w2
    l_zcg_VTP = Geo_tier.z_zbar_VTP - z_XCG;
    l_xcg_VTP = Geo_tier.x_xbar_VTP - x_XCG;
    z_v = (l_zcg_VTP*cos(trim_alpha) - l_xcg_VTP*sin(trim_alpha));
end

if AC_CONFIGURATION.HTP == 1
    % Arm from Zcg to Zac_w2
    l_xcg_w2 = Geo_tier.x_xbar_w2 - x_XCG;
    l_zcg_w2 = Geo_tier.z_zbar_w2 - z_XCG;
    z_v = (l_zcg_w2*cos(trim_alpha) - l_xcg_w2*sin(trim_alpha));
end

if AC_CONFIGURATION.Vee == 1
    % Arm from Zcg to Zac_vee
    l_xcg_vee = Geo_tier.x_xbar_vee - x_XCG;
    l_zcg_vee = Geo_tier.z_zbar_vee - z_XCG;
    z_v = (l_zcg_vee*cos(trim_alpha) - l_xcg_vee*sin(trim_alpha));
end

if AC_CONFIGURATION.Vee2 == 1
    % Arm from Zcg to Zac_vee2
    l_xcg_vee2 = Geo_tier.x_xbar_vee2 - x_XCG;
    l_zcg_vee2 = Geo_tier.z_zbar_vee2 - z_XCG;
    z_v2 = (l_zcg_vee2*cos(trim_alpha) - l_xcg_vee2*sin(trim_alpha));
end

if AC_CONFIGURATION.Can == 1
    % Arm from Zcg to Zac_can
    l_xcg_can = x_XCG - Geo_tier.x_xbar_can;
    l_zcg_can = z_XCG - Geo_tier.z_zbar_can;
    z_v = (l_zcg_can*cos(trim_alpha) - l_xcg_can*sin(trim_alpha));
end

% Aero
% Aero_TH
% Geo_tier
% Weight_tier
% Prop_data
% conv_UNITS
% Body_Geo
% Design_criteria
% Posicion_Palanca
% Performance
% OUTPUT_read_XLSX
% AC_CONFIGURATION
% Performance_preliminar

model_conversion = 1;
if model_conversion == 1
    % oswald factor
    taper_e = 1;
    a1 = 0.0004;
    a2 = -0.0080;
    a3 = 0.0501;
    a4 =0.8642;
    lambda1 = Geo_tier.AR_w1_e*taper_e/(cos(Geo_tier.Lambda_LE_w1_e));
    R = a1*lambda1^3 + a2*lambda1^2 + a3*lambda1 + a4;
    e_e = 1.1*CLalpha_w1_e/(R*CLalpha_w1_e + (1-R)*pi*Geo_tier.AR_w1_e);
    CL = Trim_ITER.CL;

    deps_dalpha_vee = Effects.deps_dalpha_vee;
    deps_dalpha_can = Effects.deps_dalpha_can;
    deps_dalpha_h = Effects.deps_dalpha_h;

    if AC_CONFIGURATION.HTP == 1
        modelo.general.downwash = deps_dalpha_h;
    elseif AC_CONFIGURATION.Vee == 1
        modelo.general.downwash = deps_dalpha_vee;
    else
        modelo.general.downwash = 0;
    end

    modelo.general.upwash = deps_dalpha_can;

    w_R_30 = Trim_ITER.w_R_30;
    w_R_60 = Trim_ITER.w_R_60;
    w_R_90 = Trim_ITER.w_R_90;
    rho = Performance.rho;
    V = conditions.V;
    a = Performance.a;
    modelo.general.qinf = 0.5*rho*V^2;
    modelo.general.CL = Trim_ITER.CL;
    modelo.general.CL_w = Trim_ITER.CL_w1;
    %     modelo.general.CL_h = CL_HTP;
    modelo.general.mtow = conditions.m_TOW;
    % Fracción de Peso
    modelo.general.w_w0 = 1;
    modelo.general.Sref = Geo_tier.S_ref;
    modelo.general.Minf = V/a;
    modelo.general.Vinf = V;
    modelo.general.rhoinf = rho;
    modelo.general.Xcg = x_XCG;
    modelo.general.Zcg = z_XCG;
    modelo.general.CLa_WB = Stab_Der_parts.CL_alpha_w1;
    modelo.general.l_fus = Body_Geo.l_fus;

    modelo.conversion.m22ft2 = conv_UNITS.m22ft2;
    modelo.conversion.m2ft = conv_UNITS.m2ft;
    modelo.propulsion.engines = Prop_data.n_eng; %numero de motores
    modelo.propulsion.blades = 2; %numero de palas
    modelo.propulsion.beta = 0.25; %[deg]
    modelo.propulsion.X = Geo_tier.x_d_T;
    modelo.propulsion.wR30 = w_R_30;
    modelo.propulsion.wR60 = w_R_60;
    modelo.propulsion.wR90 = w_R_90;
    modelo.propulsion.D = Prop_data.D_prop;

    modelo.ala.S = Geo_tier.S_w1;
    modelo.ala.AR = Geo_tier.AR_w1;
    modelo.ala.ARwe = Geo_tier.AR_w1_e;
    modelo.ala.XLE = Geo_tier.x_w1_LE;
    modelo.ala.Xca = Geo_tier.x_xbar_w1;
    %vertical distance between the fuselage line and wing root (positive if wing above center line)
    % modelo.ala.Zca1 = z_w1_LE + ybar_w*tan(Sigma_w*D2R);
    modelo.ala.Zca1 = Geo_tier.z_zbar_w1;
    %vertical distance between the fuselage line and wing root (positive if wing below center line)
    modelo.ala.Zca = - modelo.ala.Zca1;
    modelo.ala.LAMc2 = Geo_tier.Lambda_c2_w1;
    modelo.ala.LAMc4 = Geo_tier.Lambda_c4_w1;
    modelo.ala.diedro = Geo_tier.dihedral_w1;
    modelo.ala.b = Geo_tier.b_w1;
    modelo.ala.TR = Geo_tier.lambda_w1;
    modelo.ala.TR_we = Geo_tier.lambda_w1_e;
    modelo.ala.xca = Geo_tier.xbar_w1;
    modelo.ala.le_y = (Geo_tier.b_w1/2)*tan(Geo_tier.Lambda_LE_w1);
    modelo.ala.ct = Geo_tier.cT_w1;
    modelo.ala.MAC = Geo_tier.cmac_w1;
    modelo.ala.MAC_e = Geo_tier.cmac_w1_e;

    y_1R_y1_ail = Geo_tier.y_1R_y1_ail;
    y_1R_y2_ail = Geo_tier.y_1R_y2_ail;
    cf_ail = Geo_tier.cf_ail;

    modelo.ala.y1_b2 = y_1R_y2_ail/(Geo_tier.b_w1/2); %outer aleron coordinate
    modelo.ala.y0_b2 = y_1R_y1_ail/(Geo_tier.b_w1/2); %inner aleron coordinate

    % modelo.ala.cm_c = c_ail/cmac_w1;
    modelo.ala.cm_c = cf_ail;
    modelo.ala.t_c = 0.15;
    %% OJO assumes is the same as 3D
    modelo.ala.Cla = Stab_Der_parts.CL_alpha_w1;
    modelo.ala.CLa = Stab_Der_parts.CL_alpha_w1;
    modelo.ala.CLa_we = Stab_Der_parts.CL_alpha_w1;
    modelo.ala.oswald = e_e;
    %% CORRECTION - Unifying method used to identify Polar
    modelo.ala.CD0 = Aero.Polar.C_D0;
    modelo.ala.CD1 = Aero.Polar.C_D1;
    modelo.ala.CD2 = Aero.Polar.C_D2;

    modelo.fuselaje.x = Body_Geo.x_Area_body;
    modelo.fuselaje.D_x = Body_Geo.D_x;
    modelo.fuselaje.vol = Body_Geo.Vol_TOT;
    modelo.fuselaje.CLa = Stab_Der_parts.CL_alpha_fus;
    modelo.fuselaje.D = Body_Geo.h_Area_b_max;
    modelo.fuselaje.l = Body_Geo.l_fus;
    modelo.fuselaje.Sside = Body_Geo.Area_side;
    modelo.fuselaje.W = Body_Geo.Area_b_max;

    % Real area, span and AR for the VTAIL
    S_w1_s = Geo_tier.S_w1_s;
    b_w1_s = Geo_tier.b_w1_s;
    AR_w1_s = Geo_tier.AR_w1_s;

    %     S_w2_s = Geo_tier.S_w2_s;
    %     b_w2_s = Geo_tier.b_w2_s;
    %     AR_w2_s = Geo_tier.AR_w2_s;
    %     S_w2_pv = Geo_tier.S_w2_pv;
    %     S_w2_ph = Geo_tier.S_w2_ph;

    % Defines equivalent parameters for VTAIL as projection of HTP and VTP

    switch AC_CONFIGURATION.AC_type
        case 1 % AC_type = 1 - flying wing
            S_w2_s = 0;
            b_w2_s = 0;
            AR_w2_s = 0;
            S_w2_pv = 0;
            S_w2_ph = 0;

            if AC_CONFIGURATION.VTP == 1

                if prop_wash_effect == 1
                    %                     eta_VTP_no_afe = afe.eta_VTP_no_afe;
                    % Cancel out cntribution of propwash to VTP
                    eta_VTP_no_afe = 1;
                else
                    % Assumes that the ratio of dynamic pressures is equal
                    % to 1
                    eta_VTP_no_afe = 1;
                end

                K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
                K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
                % VTP
                modelo.vertical.S = Geo_tier.S_VTP;
                modelo.vertical.eta = eta_VTP_no_afe;
                modelo.vertical.Cla = Stab_Der_parts.CLalpha_VTP;
                modelo.vertical.CLa = Stab_Der_parts.CLalpha_VTP;
                modelo.vertical.b = (Geo_tier.b_VTP_s);
                modelo.vertical.Xca = Geo_tier.x_xbar_VTP;
                modelo.vertical.Zca = l_zcg_VTP;
                modelo.vertical.AR = Geo_tier.AR_VTP;
                modelo.vertical.TR = Geo_tier.lambda_VTP;
                modelo.vertical.cm_c = 0.25;
                %                 modelo.vertical.t_c = 0.12; %naca 0012
                modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_VTP;
                modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_VTP;
                modelo.vertical.diedro = Geo_tier.dihedral_VTP;
                modelo.vertical.xca = Geo_tier.xbar_VTP;
                modelo.vertical.le_y = (Geo_tier.b_VTP)*tan(Geo_tier.Lambda_LE_VTP);
                modelo.vertical.ct = Geo_tier.cT_VTP;
                modelo.vertical.MAC = Geo_tier.cmac_VTP;
                modelo.vertical.MAC_e = Geo_tier.cmac_VTP_e;
                modelo.vertical.y0_b2 = K_y1_rudder_VTP;
                modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            end

        case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP

            S_w2_s = Geo_tier.S_w2_s;
            b_w2_s = Geo_tier.b_w2_s;
            AR_w2_s = Geo_tier.AR_w2_s;
            S_w2_pv = Geo_tier.S_w2_pv;
            S_w2_ph = Geo_tier.S_w2_ph;

            if AC_CONFIGURATION.VTP == 1
                if prop_wash_effect == 1
                    %                     eta_VTP_no_afe = afe.eta_VTP_no_afe;
                    % Cancel out cntribution of propwash to VTP
                    eta_VTP_no_afe = 1;
                else
                    % Assumes that the ratio of dynamic pressures is equal
                    % to 1
                    eta_VTP_no_afe = 1;
                end

                K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
                K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
                % VTP
                modelo.vertical.S = Geo_tier.S_VTP;
                modelo.vertical.eta = eta_VTP_no_afe;
                modelo.vertical.Cla = Stab_Der_parts.CLalpha_VTP;
                modelo.vertical.CLa = Stab_Der_parts.CLalpha_VTP;
                modelo.vertical.b = (Geo_tier.b_VTP_s);
                modelo.vertical.Xca = Geo_tier.x_xbar_VTP;
                modelo.vertical.Zca = l_zcg_VTP;
                modelo.vertical.AR = Geo_tier.AR_VTP;
                modelo.vertical.TR = Geo_tier.lambda_VTP;
                modelo.vertical.cm_c = 0.25;
                %                 modelo.vertical.t_c = 0.12; %naca 0012
                modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_VTP;
                modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_VTP;
                modelo.vertical.diedro = Geo_tier.dihedral_VTP;
                modelo.vertical.xca = Geo_tier.xbar_VTP;
                modelo.vertical.le_y = (Geo_tier.b_VTP)*tan(Geo_tier.Lambda_LE_VTP);
                modelo.vertical.ct = Geo_tier.cT_VTP;
                modelo.vertical.MAC = Geo_tier.cmac_VTP;
                modelo.vertical.MAC_e = Geo_tier.cmac_VTP_e;
                modelo.vertical.y0_b2 = K_y1_rudder_VTP;
                modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            end

            if AC_CONFIGURATION.HTP == 1
                if prop_wash_effect == 1
                    eta_w2_no_afe = afe.eta_w2_no_afe;
                else
                    % Assumes that the ratio of dynamic pressures is equal
                    % to 1
                    eta_w2_no_afe = 1;
                end

                K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
                K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;

                modelo.horizontal.S = Geo_tier.S_w2;
                modelo.horizontal.eta = eta_w2_no_afe;
                modelo.horizontal.Cla = Stab_Der_parts.CL_alpha_HTP;
                modelo.horizontal.CLa = Stab_Der_parts.CL_alpha_HTP;
                modelo.horizontal.b = (Geo_tier.b_w2_s);
                modelo.horizontal.Xca = Geo_tier.x_xbar_w2;
                modelo.horizontal.Zca = l_zcg_w2;
                modelo.horizontal.AR = Geo_tier.AR_w2;
                modelo.horizontal.TR = Geo_tier.lambda_w2;
                modelo.horizontal.cm_c = 0.25;
                %                 modelo.horizontal.t_c = 0.12; %naca 0012
                modelo.horizontal.LAMc2 = Geo_tier.Lambda_c2_w2;
                modelo.horizontal.LAMc4 = Geo_tier.Lambda_c4_w2;
                modelo.horizontal.diedro = Geo_tier.dihedral_w2;
                modelo.horizontal.xca = Geo_tier.xbar_w2;
                modelo.horizontal.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
                modelo.horizontal.ct = Geo_tier.cT_w2;
                modelo.horizontal.MAC = Geo_tier.cmac_w2;
                modelo.horizontal.MAC_e = Geo_tier.cmac_w2_e;
                modelo.horizontal.y0_b2 = K_y1_ele_w2;
                modelo.horizontal.y1_b2 = K_y2_ele_w2;
                %                 modelo.horizontal.y0_b2 = z_1R_y1_ele;
                %                 modelo.horizontal.y1_b2 = z_1R_y2_ele;

            end


        case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP

            % S_w2_s = Geo_tier.S_w2_s;
            % b_w2_s = Geo_tier.b_w2_s;
            % AR_w2_s = Geo_tier.AR_w2_s;
            % S_w2_pv = Geo_tier.S_w2_pv;
            % S_w2_ph = Geo_tier.S_w2_ph;

            if AC_CONFIGURATION.VTP == 1
                if prop_wash_effect == 1
                    eta_VTP_no_afe = afe.eta_VTP_no_afe;
                else
                    % Assumes that the ratio of dynamic pressures is equal
                    % to 1
                    eta_VTP_no_afe = 1;
                end

                K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
                K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;

                % VTP
                modelo.vertical.S = Geo_tier.S_VTP;
                modelo.vertical.eta = eta_VTP_no_afe;
                modelo.vertical.Cla = Stab_Der_parts.CLalpha_VTP;
                modelo.vertical.CLa = Stab_Der_parts.CLalpha_VTP;
                modelo.vertical.b = (Geo_tier.b_VTP_s);
                modelo.vertical.Xca = Geo_tier.x_xbar_VTP;
                modelo.vertical.Zca = l_zcg_VTP;
                modelo.vertical.AR = Geo_tier.AR_VTP;
                modelo.vertical.TR = Geo_tier.lambda_VTP;
                modelo.vertical.cm_c = 0.25;
                %                 modelo.vertical.t_c = 0.12; %naca 0012
                modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_VTP;
                modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_VTP;
                modelo.vertical.diedro = Geo_tier.dihedral_VTP;
                modelo.vertical.xca = Geo_tier.xbar_VTP;
                modelo.vertical.le_y = (Geo_tier.b_VTP)*tan(Geo_tier.Lambda_LE_VTP);
                modelo.vertical.ct = Geo_tier.cT_VTP;
                modelo.vertical.MAC = Geo_tier.cmac_VTP;
                modelo.vertical.MAC_e = Geo_tier.cmac_VTP_e;

                modelo.vertical.y0_b2 = K_y1_rudder_VTP;
                modelo.vertical.y1_b2 = K_y2_rudder_VTP;
                %                 modelo.vertical.y0_b2 = z_1R_y1_rudder;
                %                 modelo.vertical.y1_b2 = z_1R_y2_rudder;
            end

            if AC_CONFIGURATION.HTP == 1

                if prop_wash_effect == 1
                    eta_w2_no_afe = afe.eta_w2_no_afe;
                else
                    % Assumes that the ratio of dynamic pressures is equal
                    % to 1
                    eta_w2_no_afe = 1;
                end

                %                 CLalpha_wb_HTP
                K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
                K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;

                modelo.horizontal.S = Geo_tier.S_w2;
                modelo.horizontal.eta = eta_w2_no_afe;
                modelo.horizontal.Cla = Stab_Der_parts.CL_alpha_htp;
                modelo.horizontal.CLa = Stab_Der_parts.CL_alpha_htp;
                modelo.horizontal.b = (Geo_tier.b_w2_s);
                modelo.horizontal.Xca = Geo_tier.x_xbar_w2;
                modelo.horizontal.Zca = l_zcg_w2;
                modelo.horizontal.AR = Geo_tier.AR_w2;
                modelo.horizontal.TR = Geo_tier.lambda_w2;
                modelo.horizontal.cm_c = 0.25;
                %                 modelo.horizontal.t_c = 0.12; %naca 0012
                modelo.horizontal.LAMc2 = Geo_tier.Lambda_c2_w2;
                modelo.horizontal.LAMc4 = Geo_tier.Lambda_c4_w2;
                modelo.horizontal.diedro = Geo_tier.dihedral_w2;
                modelo.horizontal.xca = Geo_tier.xbar_w2;
                modelo.horizontal.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
                modelo.horizontal.ct = Geo_tier.cT_w2;
                modelo.horizontal.MAC = Geo_tier.cmac_w2;
                modelo.horizontal.MAC_e = Geo_tier.cmac_w2_e;
                %                 modelo.horizontal.y0_b2 = z_1R_y1_ele;
                %                 modelo.horizontal.y1_b2 = z_1R_y2_ele;
                modelo.horizontal.y0_b2 = K_y1_ele_w2;
                modelo.horizontal.y1_b2 = K_y2_ele_w2;

            end

        case 4 % AC_type = 4 - 2 surface: wing + V-tail

            S_vee_s = Geo_tier.S_vee_s;
            b_vee_s = Geo_tier.b_vee_s;
            AR_vee_s = Geo_tier.AR_vee_s;
            S_vee_pv = Geo_tier.S_vee_pv;
            S_vee_ph = Geo_tier.S_vee_ph;

            % Estimates properties as proyectos from virtual HTP and virtual
            % VTP
            if prop_wash_effect == 1
                eta_vee_no_afe = afe.eta_vee_no_afe;
            else
                % Assumes that the ratio of dynamic pressures is equal
                % to 1
                eta_vee_no_afe = 1;
            end

            K_y1_rudvtr_vee = Geo_tier.K_y1_rudvtr_vee;
            K_y2_rudvtr_vee = Geo_tier.K_y2_rudvtr_vee;

            %% si que utiliza la información en el modelo vertical para derterminar derivadas en beta. Es correcto???
            modelo.vertical.S = Geo_tier.S_vee_pv;
            modelo.vertical.eta = eta_vee_no_afe;
            modelo.vertical.Cla = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.vertical.CLa = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.vertical.b = (Geo_tier.b_vee_s/2)*sin(Geo_tier.dihedral_vee);
            modelo.vertical.Xca = Geo_tier.x_xbar_vee;
            modelo.vertical.Zca = l_zcg_vee;
            modelo.vertical.AR = (((Geo_tier.b_vee_s/2)*sin(Geo_tier.dihedral_vee))^2)/(Geo_tier.S_vee_pv/2);
            modelo.vertical.TR = Geo_tier.lambda_vee;
            modelo.vertical.cm_c = 0.25;

            %             modelo.vertical.t_c = 0.12; %naca 0012
            modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_vee;
            modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_vee;
            modelo.vertical.diedro = Geo_tier.dihedral_vee;
            modelo.vertical.xca = Geo_tier.xbar_vee;
            modelo.vertical.le_y = (Geo_tier.b_vee/2)*tan(Geo_tier.Lambda_LE_vee);
            modelo.vertical.ct = Geo_tier.cT_vee;
            modelo.vertical.MAC = Geo_tier.cmac_vee;
            modelo.vertical.MAC_e = Geo_tier.cmac_vee_e;
            %             modelo.vertical.y0_b2 = K_y1_rudder_VTP;
            %             modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            modelo.vertical.y0_b2 = Geo_tier.z_1R_y1_rudvtr;
            modelo.vertical.y1_b2 = Geo_tier.z_1R_y2_rudvtr;

            modelo.horizontal.S = Geo_tier.S_vee_ph;
            modelo.horizontal.eta = eta_vee_no_afe;
            modelo.horizontal.Cla = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.horizontal.CLa = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.horizontal.b = (Geo_tier.b_vee_s)*cos(Geo_tier.dihedral_vee);
            modelo.horizontal.Xca = Geo_tier.x_xbar_vee;
            modelo.horizontal.Zca = l_zcg_vee;
            modelo.horizontal.AR = (((Geo_tier.b_vee_s)*cos(Geo_tier.dihedral_vee))^2)/(Geo_tier.S_vee_ph);
            modelo.horizontal.TR = Geo_tier.lambda_vee;
            modelo.horizontal.cm_c = 0.25;

            %             modelo.horizontal.t_c = 0.12; %naca 0012
            modelo.horizontal.LAMc2 = Geo_tier.Lambda_c2_vee;
            modelo.horizontal.LAMc4 = Geo_tier.Lambda_c4_vee;
            modelo.horizontal.diedro = Geo_tier.dihedral_vee;
            modelo.horizontal.xca = Geo_tier.xbar_vee;
            modelo.horizontal.le_y = (Geo_tier.b_vee/2)*tan(Geo_tier.Lambda_LE_vee);
            modelo.horizontal.ct = Geo_tier.cT_vee;
            modelo.horizontal.MAC = Geo_tier.cmac_vee;
            modelo.horizontal.MAC_e = Geo_tier.cmac_vee_e;
            modelo.horizontal.y0_b2 = Geo_tier.z_1R_y1_rudvtr;
            modelo.horizontal.y1_b2 = Geo_tier.z_1R_y2_rudvtr;
            modelo.horizontal.y0_b2 = K_y1_rudvtr_vee;
            modelo.horizontal.y1_b2 = K_y2_rudvtr_vee;

            modelo.vee.S = Geo_tier.S_vee_s;
            modelo.vee.b = Geo_tier.b_vee_s;
            modelo.vee.Xca = Geo_tier.x_xbar_vee;
            modelo.vee.Zca = l_zcg_vee;
            modelo.vee.dihedral_vee = Geo_tier.dihedral_vee;

        case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail

            S_vee_s = Geo_tier.S_vee_s;
            b_vee_s = Geo_tier.b_vee_s;
            AR_vee_s = Geo_tier.AR_vee_s;
            S_vee_pv = Geo_tier.S_vee_pv;
            S_vee_ph = Geo_tier.S_vee_ph;

            % Estimates properties as proyectos from virtual HTP and virtual
            % VTP
            if prop_wash_effect == 1
                eta_vee_no_afe = afe.eta_vee_no_afe;
            else
                % Assumes that the ratio of dynamic pressures is equal
                % to 1
                eta_vee_no_afe = 1;
            end

            K_y1_rudvtr_vee = Geo_tier.K_y1_rudvtr_vee;
            K_y2_rudvtr_vee = Geo_tier.K_y2_rudvtr_vee;

            %% si que utiliza la información en el modelo vertical para derterminar derivadas en beta. Es correcto???
            modelo.vertical.S = Geo_tier.S_vee_pv;
            modelo.vertical.eta = eta_vee_no_afe;
            modelo.vertical.Cla = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.vertical.CLa = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.vertical.b = (Geo_tier.b_vee_s/2)*sin(Geo_tier.dihedral_vee);
            modelo.vertical.Xca = Geo_tier.x_xbar_vee;
            modelo.vertical.Zca = l_zcg_vee;
            modelo.vertical.AR = (((Geo_tier.b_vee_s/2)*sin(Geo_tier.dihedral_vee))^2)/(Geo_tier.S_vee_pv/2);
            modelo.vertical.TR = Geo_tier.lambda_vee;
            modelo.vertical.cm_c = 0.25;

            %             modelo.vertical.t_c = 0.12; %naca 0012
            modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_vee;
            modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_vee;
            modelo.vertical.diedro = Geo_tier.dihedral_vee;
            modelo.vertical.xca = Geo_tier.xbar_vee;
            modelo.vertical.le_y = (Geo_tier.b_vee/2)*tan(Geo_tier.Lambda_LE_vee);
            modelo.vertical.ct = Geo_tier.cT_vee;
            modelo.vertical.MAC = Geo_tier.cmac_vee;
            modelo.vertical.MAC_e = Geo_tier.cmac_vee_e;
            %             modelo.vertical.y0_b2 = K_y1_rudder_VTP;
            %             modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            modelo.vertical.y0_b2 = Geo_tier.z_1R_y1_rudvtr;
            modelo.vertical.y1_b2 = Geo_tier.z_1R_y2_rudvtr;

            modelo.horizontal.S = Geo_tier.S_vee_ph;
            modelo.horizontal.eta = eta_vee_no_afe;
            modelo.horizontal.Cla = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.horizontal.CLa = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.horizontal.b = (Geo_tier.b_vee_s)*cos(Geo_tier.dihedral_vee);
            modelo.horizontal.Xca = Geo_tier.x_xbar_vee;
            modelo.horizontal.Zca = l_zcg_vee;
            modelo.horizontal.AR = (((Geo_tier.b_vee_s)*cos(Geo_tier.dihedral_vee))^2)/(Geo_tier.S_vee_ph);
            modelo.horizontal.TR = Geo_tier.lambda_vee;
            modelo.horizontal.cm_c = 0.25;

            %             modelo.horizontal.t_c = 0.12; %naca 0012
            modelo.horizontal.LAMc2 = Geo_tier.Lambda_c2_vee;
            modelo.horizontal.LAMc4 = Geo_tier.Lambda_c4_vee;
            modelo.horizontal.diedro = Geo_tier.dihedral_vee;
            modelo.horizontal.xca = Geo_tier.xbar_vee;
            modelo.horizontal.le_y = (Geo_tier.b_vee/2)*tan(Geo_tier.Lambda_LE_vee);
            modelo.horizontal.ct = Geo_tier.cT_vee;
            modelo.horizontal.MAC = Geo_tier.cmac_vee;
            modelo.horizontal.MAC_e = Geo_tier.cmac_vee_e;
            modelo.horizontal.y0_b2 = Geo_tier.z_1R_y1_rudvtr;
            modelo.horizontal.y1_b2 = Geo_tier.z_1R_y2_rudvtr;
            modelo.horizontal.y0_b2 = K_y1_rudvtr_vee;
            modelo.horizontal.y1_b2 = K_y2_rudvtr_vee;

            modelo.vee.S = Geo_tier.S_vee_s;
            modelo.vee.b = Geo_tier.b_vee_s;
            modelo.vee.Xca = Geo_tier.x_xbar_vee;
            modelo.vee.Zca = l_zcg_vee;
            modelo.vee.dihedral_vee = Geo_tier.dihedral_vee;

            %             modelo.vee.Cla = CL_alpha_vee;
            %             modelo.vee.CLa = CL_alpha_vee;

        case 6 % AC_type = 5 - 3 surface: cannard + wing + VTP


            S_w2_s = 0;
            b_w2_s = 0;
            AR_w2_s = 0;
            S_w2_pv = 0;
            S_w2_ph = 0;

            if AC_CONFIGURATION.VTP == 1
                if prop_wash_effect == 1
                    %                     eta_VTP_no_afe = afe.eta_VTP_no_afe;
                    % Cancel out cntribution of propwash to VTP
                    eta_VTP_no_afe = 1;
                else
                    % Assumes that the ratio of dynamic pressures is equal
                    % to 1
                    eta_VTP_no_afe = 1;
                end

                K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
                K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;

                % VTP
                modelo.vertical.S = Geo_tier.S_VTP;
                modelo.vertical.eta = eta_VTP_no_afe;
                modelo.vertical.Cla = Stab_Der_parts.CLalpha_VTP;
                modelo.vertical.CLa = Stab_Der_parts.CLalpha_VTP;
                modelo.vertical.b = (Geo_tier.b_VTP_s);
                modelo.vertical.Xca = Geo_tier.x_xbar_VTP;
                modelo.vertical.Zca = l_zcg_VTP;
                modelo.vertical.AR = Geo_tier.AR_VTP;
                modelo.vertical.TR = Geo_tier.lambda_VTP;
                modelo.vertical.cm_c = 0.25;
                %                 modelo.vertical.t_c = 0.12; %naca 0012
                modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_VTP;
                modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_VTP;
                modelo.vertical.diedro = Geo_tier.dihedral_VTP;
                modelo.vertical.xca = Geo_tier.xbar_VTP;
                modelo.vertical.le_y = (Geo_tier.b_VTP)*tan(Geo_tier.Lambda_LE_VTP);
                modelo.vertical.ct = Geo_tier.cT_VTP;
                modelo.vertical.MAC = Geo_tier.cmac_VTP;
                modelo.vertical.MAC_e = Geo_tier.cmac_VTP_e;
                modelo.vertical.y0_b2 = K_y1_rudder_VTP;
                modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            end


            % S_w2_s = Geo_tier.S_w2_s;
            % b_w2_s = Geo_tier.b_w2_s;
            % AR_w2_s = Geo_tier.AR_w2_s;
            % S_w2_pv = Geo_tier.S_w2_pv;
            % S_w2_ph = Geo_tier.S_w2_ph;
            %
            % % Estimates properties as proyectos from virtual HTP and virtual
            % % VTP
            % if prop_wash_effect == 1
            %     eta_w2_no_afe = afe.eta_w2_no_afe;
            % else
            %     % Assumes that the ratio of dynamic pressures is equal
            %     % to 1
            %     eta_w2_no_afe = 1;
            % end
            %
            % K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
            % K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
            % K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
            % K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;
            %
            % modelo.vertical.S = Geo_tier.S_w2_pv;
            % modelo.vertical.eta = eta_w2_no_afe;
            % modelo.vertical.Cla = Stab_Der_parts.CL_alpha_wb_vee;
            % modelo.vertical.CLa = Stab_Der_parts.CL_alpha_wb_vee;
            % modelo.vertical.b = (Geo_tier.b_w2_s/2)*sin(Geo_tier.dihedral_w2);
            % modelo.vertical.Xca = Geo_tier.x_xbar_w2;
            % modelo.vertical.Zca = l_zcg_w2;
            % modelo.vertical.AR = (((Geo_tier.b_w2_s/2)*sin(Geo_tier.dihedral_w2))^2)/(Geo_tier.S_w2_pv/2);
            % modelo.vertical.TR = Geo_tier.lambda_w2;
            % modelo.vertical.cm_c = 0.25;
            % %             modelo.vertical.t_c = 0.12; %naca 0012
            % modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_w2;
            % modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_w2;
            % modelo.vertical.diedro = Geo_tier.dihedral_w2;
            % modelo.vertical.xca = Geo_tier.xbar_w2;
            % modelo.vertical.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
            % modelo.vertical.ct = Geo_tier.cT_w2;
            % modelo.vertical.MAC = Geo_tier.cmac_w2;
            % modelo.vertical.MAC_e = Geo_tier.cmac_w2_e;
            % modelo.vertical.y0_b2 = K_y1_rudder_VTP;
            % modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            % %             modelo.vertical.y0_b2 = z_1R_y1_rudvtr;
            % %             modelo.vertical.y1_b2 = z_1R_y2_rudvtr;
            %
            % modelo.horizontal.S = Geo_tier.S_w2_ph;
            % modelo.horizontal.eta = eta_w2_no_afe;
            % modelo.horizontal.Cla = Stab_Der_parts.CL_alpha_wb_vee;
            % modelo.horizontal.CLa = Stab_Der_parts.CL_alpha_wb_vee;
            % modelo.horizontal.b = (Geo_tier.b_w2_s)*cos(Geo_tier.dihedral_w2);
            % modelo.horizontal.Xca = Geo_tier.x_xbar_w2;
            % modelo.horizontal.Zca = l_zcg_w2;
            % modelo.horizontal.AR = (((Geo_tier.b_w2_s)*cos(Geo_tier.dihedral_w2))^2)/(Geo_tier.S_w2_ph);
            % modelo.horizontal.TR = Geo_tier.lambda_w2;
            % modelo.horizontal.cm_c = 0.25;
            % %             modelo.horizontal.t_c = 0.12; %naca 0012
            % modelo.horizontal.LAMc2 = Geo_tier.Lambda_c2_w2;
            % modelo.horizontal.LAMc4 = Geo_tier.Lambda_c4_w2;
            % modelo.horizontal.diedro = Geo_tier.dihedral_w2;
            % modelo.horizontal.xca = Geo_tier.xbar_w2;
            % modelo.horizontal.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
            % modelo.horizontal.ct = Geo_tier.cT_w2;
            % modelo.horizontal.MAC = Geo_tier.cmac_w2;
            % modelo.horizontal.MAC_e = Geo_tier.cmac_w2_e;
            % %             modelo.horizontal.y0_b2 = z_1R_y1_ele;
            % %             modelo.horizontal.y1_b2 = z_1R_y2_ele;
            % modelo.horizontal.y0_b2 = K_y1_ele_w2;
            % modelo.horizontal.y1_b2 = K_y2_ele_w2;

        case 8 % AC_type = 8 - 2 surface: wing + V-tail + V-tail
            % AC_type = 7 - 3 surface: wing + V-tail + VTP

            % Vee 1
            S_vee_s = Geo_tier.S_vee_s;
            b_vee_s = Geo_tier.b_vee_s;
            AR_vee_s = Geo_tier.AR_vee_s;
            S_vee_pv = Geo_tier.S_vee_pv;
            S_vee_ph = Geo_tier.S_vee_ph;

            % Vee 2
            S_vee2_s = Geo_tier.S_vee2_s;
            b_vee2_s = Geo_tier.b_vee2_s;
            AR_vee2_s = Geo_tier.AR_vee2_s;
            S_vee2_pv = Geo_tier.S_vee2_pv;
            S_vee2_ph = Geo_tier.S_vee2_ph;

            % Estimates properties as proyectos from virtual HTP and virtual
            % VTP
            if prop_wash_effect == 1
                eta_vee_no_afe = afe.eta_vee_no_afe;
                eta_vee2_no_afe = afe.eta_vee2_no_afe;
            else
                % Assumes that the ratio of dynamic pressures is equal
                % to 1
                eta_vee_no_afe = 1;
                eta_vee2_no_afe = 1;
            end

            K_y1_rudvtr_vee = Geo_tier.K_y1_rudvtr_vee;
            K_y2_rudvtr_vee = Geo_tier.K_y2_rudvtr_vee;
            K_y1_rudvtr_vee2 = Geo_tier.K_y1_rudvtr_vee2;
            K_y2_rudvtr_vee2 = Geo_tier.K_y2_rudvtr_vee2;

            %% si que utiliza la información en el modelo vertical para derterminar derivadas en beta. Es correcto???
            % Vee 1
            modelo.vertical.S = Geo_tier.S_vee_pv;
            modelo.vertical.eta = eta_vee_no_afe;
            modelo.vertical.Cla = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.vertical.CLa = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.vertical.b = (Geo_tier.b_vee_s/2)*sin(Geo_tier.dihedral_vee);
            modelo.vertical.Xca = Geo_tier.x_xbar_vee;
            modelo.vertical.Zca = l_zcg_vee;
            modelo.vertical.AR = (((Geo_tier.b_vee_s/2)*sin(Geo_tier.dihedral_vee))^2)/(Geo_tier.S_vee_pv/2);
            modelo.vertical.TR = Geo_tier.lambda_vee;
            modelo.vertical.cm_c = 0.25;

            % Vee 2
            modelo.vertical2.S = Geo_tier.S_vee2_pv;
            modelo.vertical2.eta = eta_vee2_no_afe;
            modelo.vertical2.Cla = Stab_Der_parts.CL_alpha_wb_vee2;
            modelo.vertical2.CLa = Stab_Der_parts.CL_alpha_wb_vee2;
            modelo.vertical2.b = (Geo_tier.b_vee2_s/2)*sin(Geo_tier.dihedral_vee2);
            modelo.vertical2.Xca = Geo_tier.x_xbar_vee2;
            modelo.vertical2.Zca = l_zcg_vee2;
            modelo.vertical2.AR = (((Geo_tier.b_vee2_s/2)*sin(Geo_tier.dihedral_vee2))^2)/(Geo_tier.S_vee2_pv/2);
            modelo.vertical2.TR = Geo_tier.lambda_vee2;
            modelo.vertical2.cm_c = 0.25;

            % modelo.vertical.t_c = 0.12; %naca 0012
            % Vee 1
            modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_vee;
            modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_vee;
            modelo.vertical.diedro = Geo_tier.dihedral_vee;
            modelo.vertical.xca = Geo_tier.xbar_vee;
            modelo.vertical.le_y = (Geo_tier.b_vee/2)*tan(Geo_tier.Lambda_LE_vee);
            modelo.vertical.ct = Geo_tier.cT_vee;
            modelo.vertical.MAC = Geo_tier.cmac_vee;
            modelo.vertical.MAC_e = Geo_tier.cmac_vee_e;
            modelo.vertical.y0_b2 = Geo_tier.z_1R_y1_rudvtr;
            modelo.vertical.y1_b2 = Geo_tier.z_1R_y2_rudvtr;

            % Vee 2
            modelo.vertical2.LAMc2 = Geo_tier.Lambda_c2_vee2;
            modelo.vertical2.LAMc4 = Geo_tier.Lambda_c4_vee2;
            modelo.vertical2.diedro = Geo_tier.dihedral_vee2;
            modelo.vertical2.xca = Geo_tier.xbar_vee2;
            modelo.vertical2.le_y = (Geo_tier.b_vee2/2)*tan(Geo_tier.Lambda_LE_vee2);
            modelo.vertical2.ct = Geo_tier.cT_vee2;
            modelo.vertical2.MAC = Geo_tier.cmac_vee2;
            modelo.vertical2.MAC_e = Geo_tier.cmac_vee2_e;
            modelo.vertical2.y0_b2 = Geo_tier.z_1R_y1_rudvtr2;
            modelo.vertical2.y1_b2 = Geo_tier.z_1R_y2_rudvtr2;

            % Vee 1
            modelo.horizontal.S = Geo_tier.S_vee_ph;
            modelo.horizontal.eta = eta_vee_no_afe;
            modelo.horizontal.Cla = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.horizontal.CLa = Stab_Der_parts.CL_alpha_wb_vee;
            modelo.horizontal.b = (Geo_tier.b_vee_s)*cos(Geo_tier.dihedral_vee);
            modelo.horizontal.Xca = Geo_tier.x_xbar_vee;
            modelo.horizontal.Zca = l_zcg_vee;
            modelo.horizontal.AR = (((Geo_tier.b_vee_s)*cos(Geo_tier.dihedral_vee))^2)/(Geo_tier.S_vee_ph);
            modelo.horizontal.TR = Geo_tier.lambda_vee;
            modelo.horizontal.cm_c = 0.25;

            %             modelo.horizontal.t_c = 0.12; %naca 0012
            modelo.horizontal.LAMc2 = Geo_tier.Lambda_c2_vee;
            modelo.horizontal.LAMc4 = Geo_tier.Lambda_c4_vee;
            modelo.horizontal.diedro = Geo_tier.dihedral_vee;
            modelo.horizontal.xca = Geo_tier.xbar_vee;
            modelo.horizontal.le_y = (Geo_tier.b_vee/2)*tan(Geo_tier.Lambda_LE_vee);
            modelo.horizontal.ct = Geo_tier.cT_vee;
            modelo.horizontal.MAC = Geo_tier.cmac_vee;
            modelo.horizontal.MAC_e = Geo_tier.cmac_vee_e;
            modelo.horizontal.y0_b2 = Geo_tier.z_1R_y1_rudvtr;
            modelo.horizontal.y1_b2 = Geo_tier.z_1R_y2_rudvtr;
            modelo.horizontal.y0_b2 = K_y1_rudvtr_vee;
            modelo.horizontal.y1_b2 = K_y2_rudvtr_vee;

            modelo.vertical.S = Geo_tier.S_vee_s;
            modelo.vertical.b = Geo_tier.b_vee_s;
            modelo.vertical.Xca = Geo_tier.x_xbar_vee;
            modelo.vertical.Zca = l_zcg_vee;
            modelo.vertical.dihedral_vee = Geo_tier.dihedral_vee;

            modelo.vertical2.S = Geo_tier.S_vee_s2;
            modelo.vertical2.b = Geo_tier.b_vee_s2;
            modelo.vertical2.Xca = Geo_tier.x_xbar_vee2;
            modelo.vertical2.Zca = l_zcg_vee2;
            modelo.vertical2.dihedral_vee = Geo_tier.dihedral_vee2;

            %             modelo.vee.Cla = CL_alpha_vee;
            %             modelo.vee.CLa = CL_alpha_vee;
    end
end

