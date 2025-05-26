function [Trim_ITER,Stab_Der_parts,Stab_Der] = getCM0_v3(AC_CONFIGURATION,Geo_tier,afe,Effects,Stab_Der_parts,Body_Geo,...
    Aero,conv_UNITS,Design_criteria,conditions,Stab_Der,Trim_ITER,OUTPUT_read_XLSX)

AC_type = AC_CONFIGURATION.AC_type;
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

prop_wash_effect      = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;
flagwingspan2bodydiam = OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam; 

length_fus = Body_Geo.l_fus;
length_x_position = Body_Geo.length_x_position;
width_x_position = Body_Geo.width_x_position;
x_interp = Body_Geo.x_interp;
w_Area_b_max = Body_Geo.w_Area_b_max;
Angle_fus_interp = Body_Geo.Angle_fus_interp;
S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
S_w1_e = Geo_tier.S_w1_e;
alpha_CL_0_w1 = Aero.alpha_CL_0_w1_CR;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;
i_w1 = Design_criteria.i_w1;
CMT1 = Stab_Der.CMt1;
cmac_w1 = Geo_tier.cmac_w1;
x_XCG = conditions.x_XCG;
x_xbar_w1 = Geo_tier.x_xbar_w1;
CL0_w1_e_corrected = Stab_Der_parts.CL0_w1;
CL_alpha_w1 = Stab_Der_parts.CL_alpha_w1;
S_ref = Geo_tier.S_ref;

wingspan2bodydiam = b_w1/w_Area_b_max;
%% Modified to count for the propwash
if prop_wash_effect == 1
    if W1 == 1
        eta_w1_afe_S_w1_afe_S_ref = afe.eta_w1_afe_S_w1_afe_S_ref;
        eta_w1_no_afe_S_w1_no_afe_S_ref = afe.eta_w1_no_afe_S_w1_no_afe_S_ref;
    end
    if HTP == 1
        eta_HTP_afe_S_HTP_afe_S_ref = afe.eta_HTP_afe_S_HTP_afe_S_ref;
        eta_HTP_no_afe_S_HTP_no_afe_S_ref = afe.eta_HTP_no_afe_S_HTP_no_afe_S_ref;
    end
    if Vee == 1
        eta_vee_afe_S_vee_afe_S_ref = afe.eta_vee_afe_S_vee_afe_S_ref;
        eta_vee_no_afe_S_vee_no_afe_S_ref = afe.eta_vee_no_afe_S_vee_no_afe_S_ref;
    end
    if Vee2 == 1
        eta_vee2_afe_S_vee2_afe_S_ref = afe.eta_vee2_afe_S_vee2_afe_S_ref;
        eta_vee2_no_afe_S_vee2_no_afe_S_ref = afe.eta_vee2_no_afe_S_vee2_no_afe_S_ref;
    end
    if Can == 1
        eta_can_afe_S_can_afe_S_ref = afe.eta_can_afe_S_can_afe_S_ref;
        eta_can_no_afe_S_can_no_afe_S_ref = afe.eta_can_no_afe_S_can_no_afe_S_ref;
    end
else
end

if HTP == 1
    cmac_HTP  = Geo_tier.cmac_HTP;
    x_xbar_HTP = Geo_tier.x_xbar_HTP;
    CL0_HTP_e_corrected = Stab_Der_parts.CL0_HTP;
    CL_alpha_HTP = Stab_Der_parts.CL_alpha_HTP;
    i_HTP = Design_criteria.i_HTP;
    eps_HTP = Effects.eps_HTP;
    eps_HTP_0 = Effects.eps_HTP_0;
end

if Vee == 1
    cmac_vee  = Geo_tier.cmac_vee;
    x_xbar_vee = Geo_tier.x_xbar_vee;
    CL0_vee_e_corrected = Stab_Der_parts.CL0_vee;
    CL_alpha_vee = Stab_Der_parts.CL_alpha_vee;
    i_vee = Design_criteria.i_vee;
    eps_vee = Effects.eps_vee;
    eps_vee_0 = Effects.eps_vee_0;
end

if Vee2 == 1
    cmac_vee2  = Geo_tier.cmac_vee2;
    x_xbar_vee2 = Geo_tier.x_xbar_vee2;
    CL0_vee2_e_corrected = Stab_Der_parts.CL0_vee2;
    CL_alpha_vee2 = Stab_Der_parts.CL_alpha_vee2;
    i_vee2 = Design_criteria.i_vee2;
    eps_vee2 = Effects.eps_vee2;
    eps_vee2_0 = Effects.eps_vee2_0;
end

if Can == 1
    cmac_can = Geo_tier.cmac_can;
    x_xbar_can = Geo_tier.x_xbar_can;
    CL0_can_e_corrected = Stab_Der_parts.CL0_can;
    CL_alpha_can = Stab_Der_parts.CL_alpha_can;
    i_can = Design_criteria.i_can;
    eps_can = Effects.eps_can;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cm_0_fus%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nose section of Fuselage
% Digitalizacion Fig 3.9(b) Pamadi
%CALCULO Cmo eq 3.8 PAMADI

    % wing_body_CLalpha = 0.0785;
    int_bf = @(x) interp1(length_x_position,width_x_position,x);
    alpha_CL_0_we_fus = (alpha_CL_0_w1*D2R + i_w1); %rad
    alpha_CL_0_we_fus = alpha_CL_0_we_fus*R2D; %deg : pag. 171 Pamadi: Note that in eq 3.8
    % both alpha_0_w and i_cl_B are in degrees.
    int_low_anglef = 0;
    int_up_anglef = length_fus;
    f_anglef = @(x) interp1(x_interp,Angle_fus_interp,x);
    int_Cm0 = quad(@(x) ((int_bf(x)).^2.*(alpha_CL_0_we_fus + f_anglef(x))),int_low_anglef,int_up_anglef);
    int_Cm0 = integral(@(x) ((int_bf(x)).^2.*(alpha_CL_0_we_fus + f_anglef(x))),int_low_anglef,int_up_anglef);

    Fineness_Ratio = length_fus/w_Area_b_max;
    %digitaliazacion figura PAMADI CAP3 Fig 3.6
    x_f_k2_k1 = [4.,5.,6.,8.,10.,12.,14.,16.,18.,20.];
    y_f_k2_k1 = [.77,.825,.865,.91,.94,.955,.965,.97,.973,.975];
    f_k2_k1  = interp1(x_f_k2_k1,y_f_k2_k1,Fineness_Ratio,'spline');
    CM0_fus = (f_k2_k1/(36.5*S_w1*cmac_w1))*int_Cm0;

%% Contribution of Fuselage
wingspan2bodydiam = b_w1/w_Area_b_max;
flagwingspan2bodydiam = OUTPUT_read_XLSX.Stability_flags.flagwingspan2bodydiam; 
Munk_fuselage_constribution = OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution;
if   wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1 || Munk_fuselage_constribution == 1
    CM0_fus = (f_k2_k1/(36.5*S_w1*cmac_w1))*int_Cm0;
else
    CM0_fus = 0;
end

Stab_Der_parts.CM0_fus = CM0_fus; 
Trim_ITER.CM0_f = CM0_fus;

%%%%%%%%%%%%%%%%%%%%%%%CALCULO CM0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if W1 == 1
    CM_0_w1_ae = Aero.CM_0_w1_CR;
    if prop_wash_effect == 1
        CM0_w1_e_corrected = (eta_w1_afe_S_w1_afe_S_ref*CM_0_w1_ae + eta_w1_no_afe_S_w1_no_afe_S_ref*CM_0_w1_ae)*(cmac_w1/cmac_w1);
    else
        if wingspan2bodydiam <= 2 || flagwingspan2bodydiam==1
            CM0_w1_e_corrected = CM_0_w1_ae*(S_w1_e/S_w1);
        else
            CM0_w1_e_corrected = CM_0_w1_ae;
        end
        
    end
    CM0_w1_CL0_w1 = ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CL_alpha_w1*i_w1);
    CM_0_w1_wrt_XCG = CM_0_w1_ae + CM0_w1_CL0_w1;
    Trim_ITER.CM0_w1_e_corrected = CM0_w1_e_corrected;
    Trim_ITER.CM0_w1_CL0_w1 = CM0_w1_CL0_w1;
    Stab_Der_parts.CM0_w1_ae = CM0_w1_e_corrected;
    Stab_Der_parts.CM0_w1_CL0 = CM0_w1_CL0_w1;
    Stab_Der_parts.CM_0_w1 = CM_0_w1_wrt_XCG;
end

if HTP == 1    
    % Reference wing
    S_HTP = Geo_tier.S_HTP;
    S_HTP_e = Geo_tier.S_HTP_e;
    b_HTP = Geo_tier.b_HTP;
    w_HTP = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_HTP);
    wingspan2bodydiam_HTP = b_HTP/w_HTP;
    if wingspan2bodydiam_HTP <= 2
        CM_0_HTP_ae = Aero.CM_0_HTP_CR*(S_HTP_e/S_HTP);
    else
       CM_0_HTP_ae = Aero.CM_0_HTP_CR;
    end
    
    conv_HTP = S_HTP/S_w1; %adimensionalizados con la S_HTP_total en flow5, no con la expuesta;
     
    if prop_wash_effect == 1
        CM0_HTP_e_corrected = (eta_HTP_afe_S_HTP_afe_S_ref*CM_0_HTP_ae + eta_HTP_no_afe_S_HTP_no_afe_S_ref*CM_0_HTP_ae)*(cmac_HTP/cmac_w1);
    else
        CM0_HTP_e_corrected = CM_0_HTP_ae*(cmac_HTP/cmac_w1)*conv_HTP;
    end
    
    CM0_HTP_CL0_HTP = ((x_XCG - x_xbar_HTP)/cmac_w1)*(CL0_HTP_e_corrected + CL_alpha_HTP*(i_HTP - eps_HTP_0));
    CM_0_HTP_wrt_XCG = CM_0_HTP_ae + CM0_HTP_CL0_HTP;
    CM_0_HTP_wrt_XCG = CM0_HTP_e_corrected + CM0_HTP_CL0_HTP;
    Trim_ITER.CM0_HTP_CL0_HTP = CM0_HTP_CL0_HTP;
    Trim_ITER.CM0_HTP_e_corrected = CM0_HTP_e_corrected;
    Trim_ITER.CM0_HTP_CL0_HTP = CM0_HTP_CL0_HTP;
    Stab_Der_parts.CM0_HTP_ae = CM0_HTP_e_corrected;
    Stab_Der_parts.CM0_HTP_CL0 = CM0_HTP_CL0_HTP;
    Stab_Der_parts.CM_0_HTP = CM_0_HTP_wrt_XCG;
else
    Stab_Der_parts.CM0_HTP_ae = 0;
    Stab_Der_parts.CM0_HTP_CL0 = 0;
    Stab_Der_parts.CM_0_HTP = 0;
end

if Vee == 1
    b_vee = Geo_tier.b_vee;
    w_vee = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee);
    wingspan2bodydiam_vee = b_vee/w_vee;
    % Reference wing
    S_vee = Geo_tier.S_vee;
    S_vee_e = Geo_tier.S_vee_e;

    if wingspan2bodydiam_vee <= 2
        CM_0_vee_ae = Aero.CM_0_vee_CR*(S_vee_e/S_vee);
    else
        CM_0_vee_ae = Aero.CM_0_vee_CR;
    end

    conv_vee = S_vee/S_w1; %adimensionalizados con la S_vee_total en flow5, no con la expuesta;
     
    if prop_wash_effect == 1
        CM0_vee_e_corrected = (eta_vee_afe_S_vee_afe_S_ref*CM_0_vee_ae + eta_vee_no_afe_S_vee_no_afe_S_ref*CM_0_vee_ae)*(cmac_vee/cmac_w1);
    else
        CM0_vee_e_corrected = CM_0_vee_ae*(cmac_vee/cmac_w1)*conv_vee;
    end
    
    CM0_vee_CL0_vee = ((x_XCG - x_xbar_vee)/cmac_w1)*(CL0_vee_e_corrected + CL_alpha_vee*(i_vee - eps_vee_0));
    CM_0_vee_wrt_XCG = CM0_vee_e_corrected + CM0_vee_CL0_vee;
    Trim_ITER.CM0_vee_CL0_vee = CM0_vee_CL0_vee;
    Trim_ITER.CM0_vee_e_corrected = CM0_vee_e_corrected;
    Trim_ITER.CM0_vee_CL0_vee = CM0_vee_CL0_vee;
    Stab_Der_parts.CM0_vee_ae = CM0_vee_e_corrected;
    Stab_Der_parts.CM0_vee_CL0 = CM0_vee_CL0_vee;
    Stab_Der_parts.CM_0_vee = CM_0_vee_wrt_XCG;
else
        Stab_Der_parts.CM0_vee_ae = 0;
    Stab_Der_parts.CM0_vee_CL0 = 0;
    Stab_Der_parts.CM_0_vee = 0;
end

if Vee2 == 1
    b_vee2 = Geo_tier.b_vee2;
    w_vee2 = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_vee2);
    wingspan2bodydiam_vee2 = b_vee2/w_vee2;
    % Reference wing
    S_vee2 = Geo_tier.S_vee2;
    S_vee2_e = Geo_tier.S_vee2_e;
    
    if wingspan2bodydiam_vee2 <= 2
        CM_0_vee2_ae = Aero.CM_0_vee2_CR*(S_vee2_e/S_vee2);
    else
        CM_0_vee2_ae = Aero.CM_0_vee2_CR;
    end

    conv_vee2 = S_vee2/S_w1; %adimensionalizados con la S_vee2_total en flow5, no con la expuesta;
     
    if prop_wash_effect == 1
        CM0_vee2_e_corrected = (eta_vee2_afe_S_vee2_afe_S_ref*CM_0_vee2_ae + eta_vee2_no_afe_S_vee2_no_afe_S_ref*CM_0_vee2_ae)*(cmac_vee2/cmac_w1);
    else
        CM0_vee2_e_corrected = CM_0_vee2_ae*(cmac_vee2/cmac_w1)*conv_vee2;
    end
    
    CM0_vee2_CL0_vee2 = ((x_XCG - x_xbar_vee2)/cmac_w1)*(CL0_vee2_e_corrected + CL_alpha_vee2*(i_vee2 - eps_vee2_0));
    CM_0_vee2_wrt_XCG = CM0_vee2_e_corrected + CM0_vee2_CL0_vee2;
    Trim_ITER.CM0_vee2_CL0_vee2 = CM0_vee2_CL0_vee2;
    Trim_ITER.CM0_vee2_e_corrected = CM0_vee2_e_corrected;
    Trim_ITER.CM0_vee2_CL0_vee2 = CM0_vee2_CL0_vee2;
    Stab_Der_parts.CM0_vee2_ae = CM0_vee2_e_corrected;
    Stab_Der_parts.CM0_vee2_CL0 = CM0_vee2_CL0_vee2;
    Stab_Der_parts.CM_0_vee2 = CM_0_vee2_wrt_XCG;
else
    Stab_Der_parts.CM0_vee2_ae = 0;
    Stab_Der_parts.CM0_vee2_CL0 = 0;
    Stab_Der_parts.CM_0_vee2 = 0;
end

if Can == 1
    % Reference wing
    S_can = Geo_tier.S_can;
    S_can_e = Geo_tier.S_can_e;
    b_can = Geo_tier.b_can;
    w_can = interp1(Body_Geo.length_x_position, Body_Geo.width_x_position, Geo_tier.x_xbar_can);
    wingspan2bodydiam_can = b_can/w_can;    conv_can = S_can/S_w1; %adimensionalizados con la S_can_total en flow5, no con la expuesta;
    
    if wingspan2bodydiam_can <= 2
        CM_0_can_ae = Aero.CM_0_can_CR*(S_can_e/S_can);
    else
        CM_0_can_ae = Aero.CM_0_can_CR;
    end
    if prop_wash_effect == 1
        CM0_can_e_corrected = (eta_can_afe_S_can_afe_S_ref*CM_0_can_ae + eta_can_no_afe_S_can_no_afe_S_ref*CM_0_can_ae)*(cmac_can/cmac_w1);
    else
        CM0_can_e_corrected = CM_0_can_ae*(cmac_can/cmac_w1)*conv_can;
    end
    CM0_can_CL0_can = ((x_XCG - x_xbar_can)/cmac_w1)*(CL0_can_e_corrected + CL_alpha_can*(i_can + eps_can));
    CM_0_can_wrt_XCG = CM0_can_e_corrected + CM0_can_CL0_can;
    Trim_ITER.CM0_can_CL0_can = CM0_can_CL0_can;
    Trim_ITER.CM0_can_e_corrected = CM0_can_e_corrected;
    Stab_Der_parts.CM0_can_ae = CM0_can_e_corrected;
    Stab_Der_parts.CM0_can_CL0 = CM0_can_CL0_can;
    Stab_Der_parts.CM_0_can = CM_0_can_wrt_XCG;
else
    Stab_Der_parts.CM0_can_ae = 0;
    Stab_Der_parts.CM0_can_CL0 = 0;
    Stab_Der_parts.CM_0_can = 0;
end

if Nac == 1
    Stab_Der_parts.CM0_nac_ae = 0;
    Stab_Der_parts.CM0_nac_CL0 = 0;
    Stab_Der_parts.CM_0_nac = 0;
else
    Stab_Der_parts.CM0_nac_ae = 0;
    Stab_Der_parts.CM0_nac_CL0 = 0;
    Stab_Der_parts.CM_0_nac = 0;
end

CMT1 = Trim_ITER.CMT1;

% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
% AC_type = 6 - 3 surface: cannard + wing + VTP
if AC_type == 1
%     CM0_ac = CM0_fus + CM_0_w1_wrt_XCG;
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CMT1;
elseif AC_type == 2
%     CM0_fus
%     CM_0_w1_wrt_XCG
%     CM_0_HTP_wrt_XCG
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_HTP_wrt_XCG + CMT1;    
%     CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_HTP_wrt_XCG;    
elseif AC_type == 3
    %     CM0_w1_wrt_XCG = CM0_w1_e_corrected + ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CL_alpha_w1_corrected*i_w1);
    %     CM0_can_wrt_XCG = CM0_can_e_corrected*(cmac_can/cmac_w1) + ((x_XCG - x_xbar_can)/cmac_w1)*(CL0_can_e_corrected + CLalpha_can*(i_can + eps_can));
    %     CM0_HTP_wrt_XCG = CM0_HTP_e_corrected*(cmac_HTP/cmac_w1) + ((x_XCG - x_xbar_HTP)/cmac_w1)*(CL0_HTP_e_corrected + CL_alpha_wb_HTP*(i_HTP - eps_HTP));
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_can_wrt_XCG + CM_0_HTP_wrt_XCG + CMT1;  
elseif AC_type == 4
    %     CM0_w1_wrt_XCG = CM0_w1_e_corrected + ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CL_alpha_w1_corrected*i_w1);
    %     CM0_HTP_wrt_XCG = CM0_HTP_e_corrected*(cmac_HTP/cmac_w1) + ((x_XCG - x_xbar_HTP)/cmac_w1)*(CL0_vee_e_corrected + CL_alpha_wb_vee*(i_HTP - eps_HTP));
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_vee_wrt_XCG + CMT1;
%     CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_HTP_wrt_XCG;
elseif AC_type == 5
    %     CM0_w1_wrt_XCG = CM0_w1_e_corrected + ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CL_alpha_w1_corrected*i_w1);
    %     CM0_can_wrt_XCG = CM0_can_e_corrected*(cmac_can/cmac_w1) + ((x_XCG - x_xbar_can)/cmac_w1)*(CL0_can_e_corrected + CLalpha_can*(i_can + eps_can));
    %     CM0_HTP_wrt_XCG = CM0_HTP_e_corrected*(cmac_HTP/cmac_w1) + ((x_XCG - x_xbar_HTP)/cmac_w1)*(CL0_vee_e_corrected + CL_alpha_wb_vee*(i_HTP - eps_HTP));
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_can_wrt_XCG + CM_0_vee_wrt_XCG + CMT1;
%     CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_can_wrt_XCG + CM_0_HTP_wrt_XCG;
elseif AC_type == 6
    %     CM0_w1_wrt_XCG = CM0_w1_e_corrected + ((x_XCG - x_xbar_w1)/cmac_w1)*(CL0_w1_e_corrected + CL_alpha_w1_corrected*i_w1);
    %     CM0_can_wrt_XCG = CM0_can_e_corrected*(cmac_can/cmac_w1) + ((x_XCG - x_xbar_can)/cmac_w1)*(CL0_can_e_corrected + CLalpha_can*(i_can + eps_can));
    %     CM0_HTP_wrt_XCG = CM0_HTP_e_corrected*(cmac_HTP/cmac_w1) + ((x_XCG - x_xbar_HTP)/cmac_w1)*(CL0_vee_e_corrected + CL_alpha_wb_vee*(i_HTP - eps_HTP));
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_can_wrt_XCG + CMT1;
%     CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_can_wrt_XCG;
elseif AC_type == 7
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_vee_wrt_XCG + CMT1;
elseif AC_type == 8
    CM0_ac = CM0_fus + CM_0_w1_wrt_XCG + CM_0_vee_wrt_XCG + CM_0_vee2_wrt_XCG + CMT1;
end
Trim_ITER.CM0_ac = CM0_ac;
Stab_Der.CM0_ac = CM0_ac;
Stab_Der_parts.CM0_ac = CM0_ac;
end
