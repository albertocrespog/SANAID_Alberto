function modelo = conv_modelo(AC_CONFIGURATION,TRIM_RESULTS,conditions,Trim_ITER,Stab_Der_parts, Geo_tier,...
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

if AC_CONFIGURATION.HTP == 1 | AC_CONFIGURATION.Vee == 1
    % Arm from Zcg to Zac_w2
    l_xcg_w2 = Geo_tier.x_xbar_w2 - x_XCG;
    l_zcg_w2 = Geo_tier.z_zbar_w2 - z_XCG;
    z_v = (l_zcg_w2*cos(trim_alpha) - l_xcg_w2*sin(trim_alpha));
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
    
    deps_dalpha_Vee = Effects.deps_dalpha_Vee;
    deps_dalpha_can = Effects.deps_dalpha_can;
    deps_dalpha_h = Effects.deps_dalpha_h;
    
    if AC_CONFIGURATION.HTP == 1
        modelo.general.downwash = deps_dalpha_h;
    elseif AC_CONFIGURATION.Vee == 1
        modelo.general.downwash = deps_dalpha_Vee;
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

            S_w2_s = Geo_tier.S_w2_s;
            b_w2_s = Geo_tier.b_w2_s;
            AR_w2_s = Geo_tier.AR_w2_s;
            S_w2_pv = Geo_tier.S_w2_pv;
            S_w2_ph = Geo_tier.S_w2_ph;

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

            if AC_CONFIGURATION.HTP == 1 | AC_CONFIGURATION.Vee == 1
    
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
            
            S_w2_s = Geo_tier.S_w2_s;
            b_w2_s = Geo_tier.b_w2_s;
            AR_w2_s = Geo_tier.AR_w2_s;
            S_w2_pv = Geo_tier.S_w2_pv;
            S_w2_ph = Geo_tier.S_w2_ph;
            
            % Estimates properties as proyectos from virtual HTP and virtual
            % VTP
            if prop_wash_effect == 1
                eta_w2_no_afe = afe.eta_w2_no_afe;
            else
                % Assumes that the ratio of dynamic pressures is equal
                % to 1
                eta_w2_no_afe = 1;
            end

            K_y1_rudvtr_vee = Geo_tier.K_y1_rudvtr_vee;
            K_y2_rudvtr_vee = Geo_tier.K_y2_rudvtr_vee;

            %% si que utiliza la información en el modelo vertical para derterminar derivadas en beta. Es correcto???
            modelo.vertical.S = Geo_tier.S_w2_pv;
            modelo.vertical.eta = eta_w2_no_afe;
            modelo.vertical.Cla = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.vertical.CLa = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.vertical.b = (Geo_tier.b_w2_s/2)*sin(Geo_tier.dihedral_w2);
            modelo.vertical.Xca = Geo_tier.x_xbar_w2;
            modelo.vertical.Zca = l_zcg_w2;
            modelo.vertical.AR = (((Geo_tier.b_w2_s/2)*sin(Geo_tier.dihedral_w2))^2)/(Geo_tier.S_w2_pv/2);
            modelo.vertical.TR = Geo_tier.lambda_w2;
            modelo.vertical.cm_c = 0.25;
            %             modelo.vertical.t_c = 0.12; %naca 0012
            modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_w2;
            modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_w2;
            modelo.vertical.diedro = Geo_tier.dihedral_w2;
            modelo.vertical.xca = Geo_tier.xbar_w2;
            modelo.vertical.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
            modelo.vertical.ct = Geo_tier.cT_w2;
            modelo.vertical.MAC = Geo_tier.cmac_w2;
            modelo.vertical.MAC_e = Geo_tier.cmac_w2_e;
            %             modelo.vertical.y0_b2 = K_y1_rudder_VTP;
            %             modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            modelo.vertical.y0_b2 = Geo_tier.z_1R_y1_rudvtr;
            modelo.vertical.y1_b2 = Geo_tier.z_1R_y2_rudvtr;

            modelo.horizontal.S = Geo_tier.S_w2_ph;
            modelo.horizontal.eta = eta_w2_no_afe;
            modelo.horizontal.Cla = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.horizontal.CLa = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.horizontal.b = (Geo_tier.b_w2_s)*cos(Geo_tier.dihedral_w2);
            modelo.horizontal.Xca = Geo_tier.x_xbar_w2;
            modelo.horizontal.Zca = l_zcg_w2;
            modelo.horizontal.AR = (((Geo_tier.b_w2_s)*cos(Geo_tier.dihedral_w2))^2)/(Geo_tier.S_w2_ph);
            modelo.horizontal.TR = Geo_tier.lambda_w2;
            modelo.horizontal.cm_c = 0.25;
            %             modelo.horizontal.t_c = 0.12; %naca 0012
            modelo.horizontal.LAMc2 = Geo_tier.Lambda_c2_w2;
            modelo.horizontal.LAMc4 = Geo_tier.Lambda_c4_w2;
            modelo.horizontal.diedro = Geo_tier.dihedral_w2;
            modelo.horizontal.xca = Geo_tier.xbar_w2;
            modelo.horizontal.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
            modelo.horizontal.ct = Geo_tier.cT_w2;
            modelo.horizontal.MAC = Geo_tier.cmac_w2;
            modelo.horizontal.MAC_e = Geo_tier.cmac_w2_e;
            modelo.horizontal.y0_b2 = Geo_tier.z_1R_y1_rudvtr;
            modelo.horizontal.y1_b2 = Geo_tier.z_1R_y2_rudvtr;
            modelo.horizontal.y0_b2 = K_y1_rudvtr_vee;
            modelo.horizontal.y1_b2 = K_y2_rudvtr_vee;

%             modelo.vee.S = Geo_tier.S_w2_s;
%             modelo.vee.b = Geo_tier.b_w2_s;
%             modelo.vee.Xca = Geo_tier.x_xbar_w2;
%             modelo.vee.Zca = l_zcg_w2;
%             modelo.vee.dihedral_w2 = Geo_tier.dihedral_w2;

        case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail

            S_w2_s = Geo_tier.S_w2_s;
            b_w2_s = Geo_tier.b_w2_s;
            AR_w2_s = Geo_tier.AR_w2_s;
            S_w2_pv = Geo_tier.S_w2_pv;
            S_w2_ph = Geo_tier.S_w2_ph;

            % Estimates properties as proyectos from virtual HTP and virtual
            % VTP
            if prop_wash_effect == 1
                eta_w2_no_afe = afe.eta_w2_no_afe;
            else
                % Assumes that the ratio of dynamic pressures is equal
                % to 1
                eta_w2_no_afe = 1;
            end

            K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
            K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
            K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
            K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;
            
            %% si que utiliza la información en el modelo vertical para derterminar derivadas en beta. Es correcto???
            modelo.vertical.S = Geo_tier.S_w2_pv;
            modelo.vertical.eta = eta_w2_no_afe;
            modelo.vertical.Cla = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.vertical.CLa = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.vertical.b = (Geo_tier.b_w2_s/2)*sin(Geo_tier.dihedral_w2);
            modelo.vertical.Xca = Geo_tier.x_xbar_w2;
            modelo.vertical.Zca = l_zcg_w2;
            modelo.vertical.AR = (((Geo_tier.b_w2_s/2)*sin(Geo_tier.dihedral_w2))^2)/(Geo_tier.S_w2_pv/2);
            modelo.vertical.TR = Geo_tier.lambda_w2;
            modelo.vertical.cm_c = 0.25;
            %             modelo.vertical.t_c = 0.12; %naca 0012
            modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_w2;
            modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_w2;
            modelo.vertical.diedro = Geo_tier.dihedral_w2;
            modelo.vertical.xca = Geo_tier.xbar_w2;
            modelo.vertical.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
            modelo.vertical.ct = Geo_tier.cT_w2;
            modelo.vertical.MAC = Geo_tier.cmac_w2;
            modelo.vertical.MAC_e = Geo_tier.cmac_w2_e;
            modelo.vertical.y0_b2 = K_y1_rudder_VTP;
            modelo.vertical.y1_b2 = K_y2_rudder_VTP;

            
            %             modelo.vertical.y0_b2 = z_1R_y1_rudvtr;
            %             modelo.vertical.y1_b2 = z_1R_y2_rudvtr;

            modelo.horizontal.S = Geo_tier.S_w2_ph;
            modelo.horizontal.eta = eta_w2_no_afe;
            modelo.horizontal.Cla = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.horizontal.CLa = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.horizontal.b = (Geo_tier.b_w2_s)*cos(Geo_tier.dihedral_w2);
            modelo.horizontal.Xca = Geo_tier.x_xbar_w2;
            modelo.horizontal.Zca = l_zcg_w2;
            modelo.horizontal.AR = (((Geo_tier.b_w2_s)*cos(Geo_tier.dihedral_w2))^2)/(Geo_tier.S_w2_ph);
            modelo.horizontal.TR = Geo_tier.lambda_w2;
            modelo.horizontal.cm_c = 0.25;
            %             modelo.horizontal.t_c = 0.12; %naca 0012
            modelo.horizontal.LAMc2 = Geo_tier.Lambda_c2_w2;
            modelo.horizontal.LAMc4 = Geo_tier.Lambda_c4_w2;
            modelo.horizontal.diedro = Geo_tier.dihedral_w2;
            modelo.horizontal.xca = Geo_tier.xbar_w2;
            modelo.horizontal.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
            modelo.horizontal.ct = Geo_tier.cT_w2;
            modelo.horizontal.MAC = Geo_tier.cmac_w2;
            modelo.horizontal.MAC_e = Geo_tier.cmac_w2_e;
            %             modelo.horizontal.y0_b2 = z_1R_y1_ele;
            %             modelo.horizontal.y1_b2 = z_1R_y2_ele;
            modelo.horizontal.y0_b2 = K_y1_ele_w2;
            modelo.horizontal.y1_b2 = K_y2_ele_w2;

            modelo.vee.S = Geo_tier.S_w2_s;
            modelo.vee.b = Geo_tier.b_w2_s;
            modelo.vee.Xca = Geo_tier.x_xbar_w2;
            modelo.vee.Zca = l_zcg_w2;
            modelo.vee.dihedral_w2 = Geo_tier.dihedral_w2;
            
%             modelo.vee.Cla = CL_alpha_Vee;
%             modelo.vee.CLa = CL_alpha_Vee;

        case 6 % AC_type = 5 - 3 surface: cannard + wing + VTP
            
            S_w2_s = Geo_tier.S_w2_s;
            b_w2_s = Geo_tier.b_w2_s;
            AR_w2_s = Geo_tier.AR_w2_s;
            S_w2_pv = Geo_tier.S_w2_pv;
            S_w2_ph = Geo_tier.S_w2_ph;

            % Estimates properties as proyectos from virtual HTP and virtual
            % VTP
            if prop_wash_effect == 1
                eta_w2_no_afe = afe.eta_w2_no_afe;
            else
                % Assumes that the ratio of dynamic pressures is equal
                % to 1
                eta_w2_no_afe = 1;
            end

            K_y1_rudder_VTP = Geo_tier.K_y1_rudder_VTP;
            K_y2_rudder_VTP = Geo_tier.K_y2_rudder_VTP;
            K_y1_ele_w2 = Geo_tier.K_y1_ele_w2;
            K_y2_ele_w2 = Geo_tier.K_y2_ele_w2;

            modelo.vertical.S = Geo_tier.S_w2_pv;
            modelo.vertical.eta = eta_w2_no_afe;
            modelo.vertical.Cla = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.vertical.CLa = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.vertical.b = (Geo_tier.b_w2_s/2)*sin(Geo_tier.dihedral_w2);
            modelo.vertical.Xca = Geo_tier.x_xbar_w2;
            modelo.vertical.Zca = l_zcg_w2;
            modelo.vertical.AR = (((Geo_tier.b_w2_s/2)*sin(Geo_tier.dihedral_w2))^2)/(Geo_tier.S_w2_pv/2);
            modelo.vertical.TR = Geo_tier.lambda_w2;
            modelo.vertical.cm_c = 0.25;
            %             modelo.vertical.t_c = 0.12; %naca 0012
            modelo.vertical.LAMc2 = Geo_tier.Lambda_c2_w2;
            modelo.vertical.LAMc4 = Geo_tier.Lambda_c4_w2;
            modelo.vertical.diedro = Geo_tier.dihedral_w2;
            modelo.vertical.xca = Geo_tier.xbar_w2;
            modelo.vertical.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
            modelo.vertical.ct = Geo_tier.cT_w2;
            modelo.vertical.MAC = Geo_tier.cmac_w2;
            modelo.vertical.MAC_e = Geo_tier.cmac_w2_e;
            modelo.vertical.y0_b2 = K_y1_rudder_VTP;
            modelo.vertical.y1_b2 = K_y2_rudder_VTP;
            %             modelo.vertical.y0_b2 = z_1R_y1_rudvtr;
            %             modelo.vertical.y1_b2 = z_1R_y2_rudvtr;

            modelo.horizontal.S = Geo_tier.S_w2_ph;
            modelo.horizontal.eta = eta_w2_no_afe;
            modelo.horizontal.Cla = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.horizontal.CLa = Stab_Der_parts.CL_alpha_wb_Vee;
            modelo.horizontal.b = (Geo_tier.b_w2_s)*cos(Geo_tier.dihedral_w2);
            modelo.horizontal.Xca = Geo_tier.x_xbar_w2;
            modelo.horizontal.Zca = l_zcg_w2;
            modelo.horizontal.AR = (((Geo_tier.b_w2_s)*cos(Geo_tier.dihedral_w2))^2)/(Geo_tier.S_w2_ph);
            modelo.horizontal.TR = Geo_tier.lambda_w2;
            modelo.horizontal.cm_c = 0.25;
            %             modelo.horizontal.t_c = 0.12; %naca 0012
            modelo.horizontal.LAMc2 = Geo_tier.Lambda_c2_w2;
            modelo.horizontal.LAMc4 = Geo_tier.Lambda_c4_w2;
            modelo.horizontal.diedro = Geo_tier.dihedral_w2;
            modelo.horizontal.xca = Geo_tier.xbar_w2;
            modelo.horizontal.le_y = (Geo_tier.b_w2/2)*tan(Geo_tier.Lambda_LE_w2);
            modelo.horizontal.ct = Geo_tier.cT_w2;
            modelo.horizontal.MAC = Geo_tier.cmac_w2;
            modelo.horizontal.MAC_e = Geo_tier.cmac_w2_e;
            %             modelo.horizontal.y0_b2 = z_1R_y1_ele;
            %             modelo.horizontal.y1_b2 = z_1R_y2_ele;
            modelo.horizontal.y0_b2 = K_y1_ele_w2;
            modelo.horizontal.y1_b2 = K_y2_ele_w2;

    end

end
end
