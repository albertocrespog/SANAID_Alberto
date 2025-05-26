function [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts,Propulsion] = ...
    Calculo_Stability_Derivatives_Noviembre2021(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
    Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar)

%Calculo de parametros propulsivos
[Propulsion,Stab_Der] = get_propulsion(AC_CONFIGURATION,Aero_TH,Aero,Prop_data,Geo_tier, conditions, Performance, conv_UNITS, OUTPUT_read_XLSX);

%Calculo de down-wash y upwash
Effects = effects(AC_CONFIGURATION,Geo_tier,Design_criteria,conditions,Aero_TH,Aero,conv_UNITS);

% Upwash influencing in prop
[Stab_Der_parts, afe] = propwash_influence(AC_CONFIGURATION,Geo_tier,Posicion_Palanca,conditions,Propulsion,Performance,Aero,Effects,OUTPUT_read_XLSX);
%NOTA IMPORTANTE: En Aero estan los CL0, CM0 y CLalpha originales, tal y como vienen de FLOW
%                 En Stab_Der_parts estan CL0 y CLalpha adimensionalizados con la Sref.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Lift Curve Slope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
body_interference = 0;
[Stab_Der_parts,  Trim_ITER] = getCLalpha(AC_CONFIGURATION,Body_Geo,Geo_tier,body_interference,Stab_Der_parts, Effects,Design_criteria);
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
%%%%%%%%%%%%%%%%%%%%%%%%CALCULO CM0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Trim_ITER,Stab_Der_parts,Stab_Der] = getCM0(AC_CONFIGURATION,Geo_tier,afe,Effects,Stab_Der_parts,Body_Geo,....
    Aero,conv_UNITS,Design_criteria,conditions,Stab_Der,Trim_ITER,OUTPUT_read_XLSX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Pitch Moment Alpha%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Cmalpha_wf, Cmalpha_WB,Cmalpha_wb_nelson] = get_cmalpha_wf_comparacion(AC_type,AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,Effects,afe)
[Trim_ITER,Stab_Der_parts,TRIM_RESULTS,Stab_Der] = get_cmalpha_wf_def(AC_CONFIGURATION,Stab_Der_parts,Body_Geo,Geo_tier,conditions,Performance,conv_UNITS,Effects,afe,SM_des,...
    OUTPUT_read_XLSX,Stab_Der,Trim_ITER);
%Se calcula el punto neutro, el SM deseado y no deseado, el CMalpha del fus
%y de la aeronave completa, todos los cálculos con y sin tener en cuenta el fuselaje.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%LONGITUDINAL CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Stab_Der,Trim_ITER] = get_long_control(AC_CONFIGURATION,Geo_tier,Stab_Der_parts,Aero_TH,conditions,afe,Performance,Stab_Der,Trim_ITER,OUTPUT_read_XLSX);


% if only_trim == 0
    
    Trim_ITER.SM_des = SM_des;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%Pitch Derivatives%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    [Stab_Der, Stab_Der_parts] = get_pitch_derivatives(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,Body_Geo,Aero,OUTPUT_read_XLSX);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Alpha punto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Stab_Der,Stab_Der_parts] = get_alphadot_derivatives(AC_CONFIGURATION, Stab_Der, Stab_Der_parts,afe,conditions, Performance,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,Effects,Aero_TH,OUTPUT_read_XLSX);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PITCH ANGLE DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    Stab_Der = get_pitch_angle_derivatives(Stab_Der, TRIM_RESULTS);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPEED DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Stab_Der] = get_speed_derivatives(Stab_Der,conditions,Performance);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%LONGITUDINAL TRIM CONDITIONS%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[TRIM_RESULTS,Trim_ITER,Stab_Der] = resolution_trim_conditions(AC_CONFIGURATION,Stab_Der,Stab_Der_parts,Design_criteria,Effects,...
    Geo_tier,conv_UNITS,conditions,Performance,Aero_TH,Aero,TRIM_RESULTS,Propulsion,Trim_ITER);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%MODEL CONVERSION FOR LATERAL DIRECTIONAL STUDIES%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelo = conv_modelo(AC_CONFIGURATION,TRIM_RESULTS,conditions,Trim_ITER,Stab_Der_parts, Geo_tier,...
    Prop_data,Aero_TH,afe,OUTPUT_read_XLSX,Aero, Effects,Performance,conv_UNITS,Body_Geo);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%COEFICIENTES Stab_Der ESTABILIDAD ESTÁTICA Stab_DerLATERAL.DIRECCIONAL %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%derivadas en funcion de beta%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Stab_Der_parts,Stab_Der] = get_beta_derivatives(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,Geo_tier,TRIM_RESULTS,Body_Geo,Aero,conditions);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%control alerones%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Stab_Der = get_aileron_latdir_der(modelo,Stab_Der,OUTPUT_read_XLSX);
    Cydeltaa = Stab_Der.Cydeltaa;
    Cldeltaa = Stab_Der.Cldeltaa;
    Cndeltaa = Stab_Der.Cndeltaa;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%control rudder%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %derivadas en deltar
    if AC_CONFIGURATION.d_rudder ==1
        Stab_Der = get_deltar_latdir_deriv(modelo,Stab_Der,Geo_tier,afe,TRIM_RESULTS.trim_alpha,OUTPUT_read_XLSX, Aero);
        Cydeltar = Stab_Der.Cydeltar;
        Cldeltar = Stab_Der.Cldeltar;
        Cndeltar = Stab_Der.Cndeltar;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%control ruddervator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    if AC_CONFIGURATION.d_rudvtr == 1
        Stab_Der = get_delta_rv_latdir_deriv(modelo,Stab_Der,afe,TRIM_RESULTS.trim_alpha,Aero,Geo_tier,OUTPUT_read_XLSX);
        
        Cydeltarv = Stab_Der.Cydeltarv;
        Cydeltar = Cydeltarv;
        Stab_Der.Cydeltar = Cydeltar;
        
        Cldeltarv = Stab_Der.Cldeltarv;
        Cldeltar = Cldeltarv;
        Stab_Der.Cldeltar = Cldeltar;
        
        Cndeltarv = Stab_Der.Cndeltarv;
        Cndeltar = Cndeltarv;
        Stab_Der.Cndeltar = Cndeltar;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% hinge derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Stab_Der = get_hinge_derivatives(modelo,TRIM_RESULTS,Stab_Der,Stab_Der_parts,Trim_ITER,OUTPUT_read_XLSX);
    Ch_alpha_deltaa= Stab_Der.Ch_alpha_deltaa;
    Ch_delta_deltaa = Stab_Der.Ch_delta_deltaa;
    Ch_beta_deltarv = Stab_Der.Ch_beta_deltarv;
    Ch_delta_deltarv = Stab_Der.Ch_delta_deltarv;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Longitudinal Forces due to Vtail incidence Derivative
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Stab_Der = get_tail_incidence_derivatives(modelo,TRIM_RESULTS,Stab_Der,Stab_Der_parts,Trim_ITER,OUTPUT_read_XLSX);
    CD_Vtail_i = Stab_Der.CD_Vtail_i;
    CL_Vtail_i = Stab_Der.CL_Vtail_i;
    CM_Vtail_i = Stab_Der.CM_Vtail_i;
    Cy_Vtail_i = Stab_Der.Cy_Vtail_i;
    Cl_Vtail_i = Stab_Der.Cl_Vtail_i;
    Cn_Vtail_i = Stab_Der.Cn_Vtail_i;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Tab Derivative
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Stab_Der = get_tab_derivatives(modelo,TRIM_RESULTS,Stab_Der,Stab_Der_parts,Trim_ITER,OUTPUT_read_XLSX);
    Cy_deltaa_Tab = Stab_Der.Cy_deltaa_Tab;
    Cl_deltaa_Tab = Stab_Der.Cl_deltaa_Tab;
    Cn_deltaa_Tab = Stab_Der.Cn_deltaa_Tab;
    Ch_deltaa_Tab = Stab_Der.Ch_deltaa_Tab;
    CL_deltarv_Tab = Stab_Der.CL_deltarv_Tab;
    Cm_deltarv_Tab = Stab_Der.Cm_deltarv_Tab;
    Ch_deltarv_Tab = Stab_Der.Ch_deltarv_Tab;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%Derivadas en funcion de p%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %ASpro Cyp
    Stab_Der = getCyp(modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts,Trim_ITER, Geo_tier, OUTPUT_read_XLSX);
    Cyp_w_diedro = Stab_Der.Cyp_w_diedro;
    Cyp_w = Stab_Der.Cyp_w;
    Cyp_v = Stab_Der.Cyp_v;
    Cyp = Stab_Der.Cyp;
    
    %ASpro Clp
    Stab_Der = getClp(modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts, Geo_tier, OUTPUT_read_XLSX);
    Clp_w = Stab_Der.Clp_w;
    Clp_v = Stab_Der.Clp_v;
    Clp = Stab_Der.Clp;
    
    %ASpro Cnp
    Stab_Der = getCnp(modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts,Trim_ITER);
    Cnp_w = Stab_Der.Cnp_w;
    Cnp_v = Stab_Der.Cnp_v;
    Cnp = Stab_Der.Cnp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%Derivadas en respecto a r%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %ASpro Cyr
    Stab_Der = getCyr(modelo,TRIM_RESULTS.trim_alpha,Stab_Der,Stab_Der_parts);
    Cyr_w = Stab_Der.Cyr_w;
    Cyr_v = Stab_Der.Cyr_v;
    Cyr = Stab_Der.Cyr;
    
    %ASpro Clr
    Stab_Der = getClr(modelo,TRIM_RESULTS.trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER);
    Clr_w = Stab_Der.Clr_w;
    Clr_v = Stab_Der.Clr_v;
    Clr = Stab_Der.Clr;
    
    %ASpro Cnr
    Stab_Der = getCnr(modelo,TRIM_RESULTS.trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER);
    Cnr_w = Stab_Der.Cnr_w;
    Cnr_v = Stab_Der.Cnr_v;
    Cnr = Stab_Der.Cnr;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%derivadas respecto a beta punto%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Roskam page 149: Except for airplanes in the high subsonic speed
    %range, the beta_dot derivatives are frequently considered
    %negligible.
    
    % ASpro Cybpunto
    Stab_Der = getCybdot(modelo,TRIM_RESULTS.trim_alpha,Stab_Der,AC_CONFIGURATION, Geo_tier, Aero);
    Cybpunto = Stab_Der.Cybpunto;
    
    
    %ASpro Clbpunto
    Stab_Der = getClbdot(modelo,TRIM_RESULTS.trim_alpha,Stab_Der);
    Clbpunto = Stab_Der.Clbpunto;
    
    %ASpro nbpunto
    Stab_Der = getCnbdot(modelo,TRIM_RESULTS.trim_alpha,Stab_Der);
    Cnbpunto = Stab_Der.Cnbpunto;
    
    % Conversor Stability Derivatives
%     Conversion_Der = Conversor_Derivatives_SANAID2AAA(Performance_preliminar,Stab_Der,Weight_tier,Geo_tier,conv_UNITS,TRIM_RESULTS,OUTPUT_read_XLSX,Trim_ITER,AC_CONFIGURATION);
%     Stab_Der_parts.Conversion_Der = Conversion_Der;      
% end

