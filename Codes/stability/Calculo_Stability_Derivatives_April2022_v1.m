function [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts,Propulsion, Effects] = ...
    Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
    Body_Geo,Design_criteria,Posicion_Palanca,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar)
SM_des = OUTPUT_read_XLSX.Stability_flags.SM_des;

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

%Calculo de parametros propulsivos
[Propulsion,Stab_Der] = get_propulsion(AC_CONFIGURATION,Aero_TH,Aero,Prop_data,Geo_tier, conditions, Performance, conv_UNITS, OUTPUT_read_XLSX);

%Calculo de down-wash y upwash
% Effects = effects(AC_CONFIGURATION,Geo_tier,Design_criteria,conditions,Aero_TH,Aero,conv_UNITS);
% Effects = effects_v2(AC_CONFIGURATION,Geo_tier,Design_criteria,conditions,Aero_TH,Aero,conv_UNITS);
Effects = effects_v3(AC_CONFIGURATION,Geo_tier,Design_criteria,conditions,Aero_TH,Aero,conv_UNITS);

% Upwash influencing in prop
% [Stab_Der_parts, afe] = propwash_influence(AC_CONFIGURATION,Geo_tier,Posicion_Palanca,conditions,Propulsion,Performance,Aero,Effects,OUTPUT_read_XLSX,Body_Geo);
% [Stab_Der_parts, afe] = propwash_influence_v2(AC_CONFIGURATION,Geo_tier,Posicion_Palanca,conditions,Propulsion,Performance,Aero,Effects,OUTPUT_read_XLSX,Body_Geo);
[Stab_Der_parts, afe] = propwash_influence_v3(AC_CONFIGURATION,Geo_tier,Posicion_Palanca,conditions,Propulsion,Performance,Aero,Effects,OUTPUT_read_XLSX,Body_Geo);

%NOTA IMPORTANTE: En Aero estan los CL0, CM0 y CLalpha originales, tal y como vienen de FLOW
%                 En Stab_Der_parts estan CL0 y CLalpha adimensionalizados con la Sref.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Lift Curve Slope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
body_interference = 0;
% [Stab_Der_parts,Trim_ITER] = getCLalpha(AC_CONFIGURATION,Body_Geo,Geo_tier,body_interference,Stab_Der_parts, Effects,Design_criteria, OUTPUT_read_XLSX);
% [Stab_Der_parts,Trim_ITER] = getCLalpha_v2(AC_CONFIGURATION,Body_Geo,Geo_tier,body_interference,Stab_Der_parts, Effects,Design_criteria, OUTPUT_read_XLSX);
[Stab_Der_parts,Trim_ITER] = getCLalpha_v3(AC_CONFIGURATION,Body_Geo,Geo_tier,body_interference,Stab_Der_parts, Effects,Design_criteria, OUTPUT_read_XLSX);

%Con body_interference=0, lo que hace es llamar a las
%obtenidas en propwash_influence, es decir, a las adimensionalizadas con Sref y calcula CLalpha y CL_0 de la aeronave completa.

%% Asumes symetrical fuselage and nacelle that should be upgraded with experimental/CFD results
CL0_fus = 0;
CL0_nac = 0;
Stab_Der_parts.CL0_fus = CL0_fus;
Stab_Der_parts.CL0_nac = CL0_nac;
% Angle_fus_x_fus = -Angle_fus_x_fus; %DUDA, esto para que?
% Angle_fus_interp = -Angle_fus_interp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Propulsive Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der,Trim_ITER] = get_propulsive_derivatives(AC_CONFIGURATION, Propulsion, Aero_TH,Aero,Prop_data,Geo_tier,Stab_Der_parts,afe,conv_UNITS,conditions,Performance,Trim_ITER,Stab_Der);
%Se calculan las derivadas propulsivas, no he tocado nada.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%CALCULO CM0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Trim_ITER,Stab_Der_parts,Stab_Der] = getCM0(AC_CONFIGURATION,Geo_tier,afe,Effects,Stab_Der_parts,Body_Geo,....
%     Aero,conv_UNITS,Design_criteria,conditions,Stab_Der,Trim_ITER,OUTPUT_read_XLSX);
% [Trim_ITER,Stab_Der_parts,Stab_Der] = getCM0_v2(AC_CONFIGURATION,Geo_tier,afe,Effects,Stab_Der_parts,Body_Geo,....
%     Aero,conv_UNITS,Design_criteria,conditions,Stab_Der,Trim_ITER,OUTPUT_read_XLSX);
[Trim_ITER,Stab_Der_parts,Stab_Der] = getCM0_v3(AC_CONFIGURATION,Geo_tier,afe,Effects,Stab_Der_parts,Body_Geo,....
    Aero,conv_UNITS,Design_criteria,conditions,Stab_Der,Trim_ITER,OUTPUT_read_XLSX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Pitch Moment Alpha%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Cmalpha_wf, Cmalpha_WB,Cmalpha_wb_nelson] = get_cmalpha_wf_comparacion(AC_type,AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,Effects,afe)
% [Trim_ITER,Stab_Der_parts,TRIM_RESULTS,Stab_Der] = get_cmalpha_wf_def(AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,Effects,afe,SM_des,...
%     OUTPUT_read_XLSX,Stab_Der,Trim_ITER);
% [Trim_ITER,Stab_Der_parts,TRIM_RESULTS,Stab_Der] = get_cmalpha_wf_def_v2(AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,Effects,afe,SM_des,...
%     OUTPUT_read_XLSX,Stab_Der,Trim_ITER);
[Trim_ITER,Stab_Der_parts,TRIM_RESULTS,Stab_Der] = get_cmalpha_wf_def_v3(AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,Effects,afe,SM_des,...
    OUTPUT_read_XLSX,Stab_Der,Trim_ITER);
%Se calcula el punto neutro, el SM deseado y no deseado, el CMalpha del fus
%y de la aeronave completa, todos los cálculos con y sin tener en cuenta el fuselaje.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%LONGITUDINAL CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Stab_Der,Trim_ITER] = get_long_control(AC_CONFIGURATION,Geo_tier,Stab_Der_parts,Aero_TH,conditions,afe,Performance,Stab_Der,Trim_ITER,OUTPUT_read_XLSX);
% [Stab_Der,Trim_ITER] = get_long_control_v2(AC_CONFIGURATION,Geo_tier,Stab_Der_parts,Aero_TH,conditions,afe,Performance,Stab_Der,Trim_ITER,OUTPUT_read_XLSX);
[Stab_Der,Trim_ITER] = get_long_control_v3(AC_CONFIGURATION,Geo_tier,Stab_Der_parts,Aero_TH,conditions,afe,Performance,Stab_Der,Trim_ITER,OUTPUT_read_XLSX);

% if only_trim == 0
Trim_ITER.SM_des = SM_des;
%Trim_ITER.x_XCG = conditions.x_XCG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Pitch Derivatives%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Stab_Der, Stab_Der_parts] = get_pitch_derivatives(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,Body_Geo,Aero,OUTPUT_read_XLSX);
% [Stab_Der, Stab_Der_parts] = get_pitch_derivatives_v2(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,Body_Geo,Aero,OUTPUT_read_XLSX);
[Stab_Der, Stab_Der_parts] = get_pitch_derivatives_v2(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,Body_Geo,Aero,OUTPUT_read_XLSX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%CD_Alpha%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stab_Der = get_CDalpha(Aero, conditions, Stab_Der, conv_UNITS, Performance, Geo_tier);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Alpha punto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Stab_Der,Stab_Der_parts] = get_alphadot_derivatives(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,Effects,Aero_TH,OUTPUT_read_XLSX);
% [Stab_Der,Stab_Der_parts] = get_alphadot_derivatives_v2(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,Effects,Aero_TH,OUTPUT_read_XLSX);
[Stab_Der,Stab_Der_parts] = get_alphadot_derivatives_v3(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,Effects,Aero_TH,OUTPUT_read_XLSX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPEED DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der] = get_speed_derivatives(Stab_Der,conditions,Performance,Aero, conv_UNITS,Geo_tier);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%LONGITUDINAL TRIM CONDITIONS%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TRIM_RESULTS,Trim_ITER,Stab_Der] = resolution_trim_conditions(AC_CONFIGURATION,Stab_Der,Stab_Der_parts,Design_criteria,Effects,...
%     Geo_tier,conv_UNITS,conditions,Performance,Aero_TH,Aero,TRIM_RESULTS,Propulsion,Trim_ITER);
% [TRIM_RESULTS,Trim_ITER,Stab_Der] = resolution_trim_conditions_v2(AC_CONFIGURATION,Stab_Der,Stab_Der_parts,Design_criteria,Effects,...
%     Geo_tier,conv_UNITS,conditions,Performance,Aero_TH,Aero,TRIM_RESULTS,Propulsion,Trim_ITER);
[TRIM_RESULTS,Trim_ITER,Stab_Der] = resolution_trim_conditions_v3(OUTPUT_read_XLSX,AC_CONFIGURATION,Stab_Der,Stab_Der_parts,Design_criteria,Effects,...
    Geo_tier,conv_UNITS,conditions,Performance,Aero_TH,Aero,TRIM_RESULTS,Propulsion,Trim_ITER,Prop_data,Weight_tier);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PITCH ANGLE DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stab_Der = get_pitch_angle_derivatives(Stab_Der, TRIM_RESULTS);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%MODEL CONVERSION FOR LATERAL DIRECTIONAL STUDIES%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modelo = conv_modelo(AC_CONFIGURATION,TRIM_RESULTS,conditions,Trim_ITER,Stab_Der_parts, Geo_tier,...
%     Prop_data,Aero_TH,afe,OUTPUT_read_XLSX,Aero, Effects,Performance,conv_UNITS,Body_Geo);
% Actualizado Junio 2022
% modelo = conv_modelo_v2(AC_CONFIGURATION,TRIM_RESULTS,conditions,Trim_ITER,Stab_Der_parts, Geo_tier,...
%     Prop_data,Aero_TH,afe,OUTPUT_read_XLSX,Aero, Effects,Performance,conv_UNITS,Body_Geo);
% modelo = conv_modelo_v3(AC_CONFIGURATION,TRIM_RESULTS,conditions,Trim_ITER,Stab_Der_parts, Geo_tier,...
%    Prop_data,Aero_TH,afe,OUTPUT_read_XLSX,Aero, Effects,Performance,conv_UNITS,Body_Geo);
modelo = conv_modelo_v4(AC_CONFIGURATION,TRIM_RESULTS,conditions,Trim_ITER,Stab_Der_parts, Geo_tier,...
    Prop_data,Aero_TH,afe,OUTPUT_read_XLSX,Aero, Effects,Performance,conv_UNITS,Body_Geo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%COEFICIENTES Stab_Der ESTABILIDAD ESTáTICA Stab_DerLATERAL.DIRECCIONAL %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%derivadas en funcion de beta%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Stab_Der_parts,Stab_Der] = get_beta_derivatives(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,conditions);
[Stab_Der_parts,Stab_Der] = get_beta_derivatives_v2(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,Geo_tier,TRIM_RESULTS,Trim_ITER,Body_Geo,Aero,conditions,OUTPUT_read_XLSX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%control alerones%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der_parts Stab_Der] = get_aileron_latdir_der(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,OUTPUT_read_XLSX);
Cydeltaa = Stab_Der.Cydeltaa;
Cldeltaa = Stab_Der.Cldeltaa;
Cndeltaa = Stab_Der.Cndeltaa;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%control rudder%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%derivadas en deltar
if AC_CONFIGURATION.d_rudder == 1
    % Stab_Der = get_deltar_latdir_deriv(AC_CONFIGURATION,modelo,Stab_Der,Geo_tier,afe,TRIM_RESULTS.trim_alpha,OUTPUT_read_XLSX, Aero);
    [Stab_Der_parts Stab_Der] = get_deltar_latdir_deriv_v2(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,Geo_tier,afe,TRIM_RESULTS.trim_alpha,OUTPUT_read_XLSX, Aero);
    Cydeltar = Stab_Der.Cydeltar;
    Cldeltar = Stab_Der.Cldeltar;
    Cndeltar = Stab_Der.Cndeltar;
else
    % No rudder
    Stab_Der.Cydeltar = 0;
    Stab_Der.Cldeltar = 0;
    Stab_Der.Cndeltar = 0;
    
    Cydeltar = Stab_Der.Cydeltar;
    Cldeltar = Stab_Der.Cldeltar;
    Cndeltar = Stab_Der.Cndeltar;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%control ruddervator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if AC_CONFIGURATION.d_rudvtr == 1
    if Vee ==1
        [Stab_Der_parts Stab_Der] = get_delta_rv_latdir_deriv(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,afe,TRIM_RESULTS.trim_alpha,Aero,Geo_tier,OUTPUT_read_XLSX);
        Cydeltarv = Stab_Der.Cydeltarv;
        Cldeltarv = Stab_Der.Cldeltarv;
        Cndeltarv = Stab_Der.Cndeltarv;
  
        Stab_Der.Cydeltar = Cydeltarv;
        Stab_Der.Cldeltar = Cldeltarv;       
        Stab_Der.Cndeltar = Cndeltarv;
    else
        Cydeltarv = 0;
        Cldeltarv = 0;
        Cndeltarv = 0;
   
        Stab_Der.Cydeltarv = Cydeltarv;
        Stab_Der.Cldeltarv = Cldeltarv;       
        Stab_Der.Cndeltarv = Cndeltarv;
    end

    if Vee2 ==1
        [Stab_Der_parts Stab_Der] = get_delta_rv_latdir_deriv_v2(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,afe,TRIM_RESULTS.trim_alpha,Aero,Geo_tier,OUTPUT_read_XLSX);
        Cydeltarv2 = Stab_Der.Cydeltarv2;
        Cldeltarv2 = Stab_Der.Cldeltarv2;
        Cndeltarv2 = Stab_Der.Cndeltarv2;

        Stab_Der.Cydeltar = Cydeltarv2;
        Stab_Der.Cldeltar = Cldeltarv2;       
        Stab_Der.Cndeltar = Cndeltarv2;
    else
        Cydeltarv2 = 0;
        Cldeltarv2 = 0;
        Cndeltarv2 = 0;

        Stab_Der.Cydeltarv2 = Cydeltarv2;
        Stab_Der.Cldeltarv2 = Cldeltarv2;       
        Stab_Der.Cndeltarv2 = Cndeltarv2;
    end

    if Vee ==1 && Vee2 ==1
        Cydeltar = Cydeltarv + Cydeltarv2;    
        Cldeltar = Cldeltarv + Cldeltarv2;
        Cndeltar = Cndeltarv + Cndeltarv2;
        
        Stab_Der.Cydeltar = Cydeltar;
        Stab_Der.Cldeltar = Cldeltar;       
        Stab_Der.Cndeltar = Cndeltar;
    end

else
    % No Vee-tail
    % Stab_Der.Cydeltar = 0;
    % Stab_Der.Cldeltar = 0;
    % Stab_Der.Cndeltar = 0;

    Stab_Der.Cydeltarv = 0;
    Stab_Der.Cldeltarv = 0;
    Stab_Der.Cndeltarv = 0;

    Stab_Der.Cydeltarv2 = 0;
    Stab_Der.Cldeltarv2 = 0;
    Stab_Der.Cndeltarv2 = 0;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% hinge derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der_parts Stab_Der] = get_hinge_derivatives(AC_CONFIGURATION,modelo,TRIM_RESULTS,Stab_Der,Stab_Der_parts,Trim_ITER,OUTPUT_read_XLSX,Geo_tier,conditions,Performance);
Ch_delta_a = Stab_Der.Ch_delta_a;
Ch_delta_r = Stab_Der.Ch_delta_r;
Ch_delta_at = Stab_Der.Ch_delta_at;
Ch_delta_rt = Stab_Der.Ch_delta_rt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Longitudinal Forces due to Vtail incidence Derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der_parts Stab_Der] = get_tail_incidence_derivatives(AC_CONFIGURATION,modelo,TRIM_RESULTS,Stab_Der,Stab_Der_parts,Trim_ITER,OUTPUT_read_XLSX);
% CD_Vtail_i = Stab_Der.CD_Vtail_i;
% CL_Vtail_i = Stab_Der.CL_Vtail_i;
% CM_Vtail_i = Stab_Der.CM_Vtail_i;
% Cy_Vtail_i = Stab_Der.Cy_Vtail_i;
% Cl_Vtail_i = Stab_Der.Cl_Vtail_i;
% Cn_Vtail_i = Stab_Der.Cn_Vtail_i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Tab Derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der_parts Stab_Der] = get_tab_derivatives(modelo,OUTPUT_read_XLSX, AC_CONFIGURATION,Geo_tier,afe,TRIM_RESULTS.trim_alpha,Aero,Stab_Der,Stab_Der_parts);
% Elevator Trim Tab
CL_delta_e_Tab = Stab_Der.CL_delta_e_Tab;
CD_delta_e_Tab = Stab_Der.CD_delta_e_Tab;
CM_delta_e_Tab = Stab_Der.CM_delta_e_Tab;
% Aileron Trim Tab
% Stab_Der.Ch_deltaa_Tab = Ch_deltaa_Tab;
Cy_delta_a_Tab = Stab_Der.Cy_delta_a_Tab;
Cl_delta_a_Tab = Stab_Der.Cl_delta_a_Tab;
Cn_delta_a_Tab = Stab_Der.Cn_delta_a_Tab;
% Rudder Trim Tab
% Stab_Der.Ch_deltaa_Tab = Ch_deltaa_Tab;
Cy_delta_r_Tab = Stab_Der.Cy_delta_r_Tab;
Cl_delta_r_Tab = Stab_Der.Cl_delta_r_Tab;
Cn_delta_r_Tab = Stab_Der.Cn_delta_r_Tab;
% Spoiler
CL_delta_sp = Stab_Der.CL_delta_sp;
CD_delta_sp = Stab_Der.CD_delta_sp;
CM_delta_sp = Stab_Der.CM_delta_sp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Derivadas en funcion de p%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ASpro Cyp
% Stab_Der = getCyp(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts,Trim_ITER, Geo_tier, OUTPUT_read_XLSX);
[Stab_Der_parts Stab_Der] = getCyp_v2(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts,Trim_ITER, Geo_tier, OUTPUT_read_XLSX);
% Cyp_w_diedro = Stab_Der.Cyp_w_diedro;
% Cyp_w = Stab_Der.Cyp_w;
% Cyp_v = Stab_Der.Cyp_v;
% Cyp_v2 = Stab_Der.Cyp_v2;
% Cyp = Stab_Der.Cyp;

%ASpro Clp
% Stab_Der = getClp(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts, Geo_tier, OUTPUT_read_XLSX);
[Stab_Der_parts Stab_Der] = getClp_v2(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts, Geo_tier, OUTPUT_read_XLSX);
% Clp_w = Stab_Der.Clp_w;
% Clp_v = Stab_Der.Clp_v;
% Clp = Stab_Der.Clp;

%ASpro Cnp
% Stab_Der = getCnp(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts,Trim_ITER);
[Stab_Der_parts Stab_Der] = getCnp_v2(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts,Trim_ITER, Geo_tier, OUTPUT_read_XLSX);
% Cnp_w = Stab_Der.Cnp_w;
% Cnp_v = Stab_Der.Cnp_v;
% Cnp_v2 = Stab_Der.Cnp_v2;
% Cnp = Stab_Der.Cnp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Derivadas en respecto a r%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ASpro Cyr
% Stab_Der = getCyr(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts);
[Stab_Der_parts Stab_Der] = getCyr_v2(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts);
% Cyr_w = Stab_Der.Cyr_w;
% Cyr_v = Stab_Der.Cyr_v;
% Cyr_v2 = Stab_Der.Cyr_v2;
% Cyr = Stab_Der.Cyr;

%ASpro Clr
% Stab_Der = getClr(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER);
[Stab_Der_parts Stab_Der] = getClr_v2(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER,Geo_tier, OUTPUT_read_XLSX);
% Clr_w = Stab_Der.Clr_w;
% Clr_v = Stab_Der.Clr_v;
% Clr_v2 = Stab_Der.Clr_v2;
% Clr = Stab_Der.Clr;

%ASpro Cnr
% Stab_Der = getCnr(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER);
[Stab_Der_parts Stab_Der] = getCnr_v2(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER,Geo_tier, OUTPUT_read_XLSX);
% Cnr_w = Stab_Der.Cnr_w;
% Cnr_v = Stab_Der.Cnr_v;
% Cnr_v2 = Stab_Der.Cnr_v2;
% Cnr = Stab_Der.Cnr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%derivadas respecto a beta punto%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Roskam page 149: Except for airplanes in the high subsonic speed
%range, the beta_dot derivatives are frequently considered
%negligible.

% ASpro Cybpunto
% Stab_Der = getCybdot(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,AC_CONFIGURATION, Geo_tier, Aero);
[Stab_Der_parts Stab_Der] = getCybdot_v2(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Geo_tier, Aero,Trim_ITER,Stab_Der_parts);
% Cybpunto = Stab_Der.Cybpunto;
% Cybpunto_v = Stab_Der.Cybpunto_v;
% Cybpunto_v2 = Stab_Der.Cybpunto_v2;


%ASpro Clbpunto
% Stab_Der = getClbdot(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der);
[Stab_Der_parts Stab_Der] = getClbdot_v2(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Geo_tier, Aero,Trim_ITER,Stab_Der_parts);
% Clbpunto = Stab_Der.Clbpunto;
% Clbpunto_v = Stab_Der.Clbpunto_v;
% Clbpunto_v2 = Stab_Der.Clbpunto_v2;

%ASpro nbpunto
% Stab_Der = getCnbdot(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der);
[Stab_Der_parts Stab_Der] = getCnbdot_v2(AC_CONFIGURATION,modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Geo_tier, Aero,Trim_ITER,Stab_Der_parts);
% Cnbpunto = Stab_Der.Cnbpunto;
% Cnbpunto_v = Stab_Der.Cnbpunto_v;
% Cnbpunto_v2 = Stab_Der.Cnbpunto_v2;

    % Conversor Stability Derivatives
%     Conversion_Der = Conversor_Derivatives_SANAID2AAA(Performance_preliminar,Stab_Der,Weight_tier,Geo_tier,conv_UNITS,TRIM_RESULTS,OUTPUT_read_XLSX,Trim_ITER,AC_CONFIGURATION);
%     Stab_Der_parts.Conversion_Der = Conversion_Der;      
% end

