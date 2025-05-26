function Conversion_Der = Conversor_Derivatives_SANAID2AAA(Performance_preliminar,Stab_Der,Weight_tier,Geo_tier,conv_UNITS...
    ,TRIM_RESULTS,OUTPUT_read_XLSX,Trim_ITER,AC_CONFIGURATION)

% get_add_path
% filename = '../AIRCRAFT/AAA/AAA_Comparisson_AC.xlsx';
% excel_file = 'AAA_Comparisson_AC.xlsx';
% writecell(C,excel_file,'Sheet',2,'Stability_Derivatives','A3:C5')

filename_read = '../AIRCRAFT/AAA/AAA_320-200.xlsx';
excel_file_read = 'AAA_320-200.xlsx';
Sheet1 = readcell(filename_read,'Sheet','1MATLAB','Range','A1:B76');

filename_write = '../AIRCRAFT/AAA/AAA_Comparisson_AC.xlsx';
cabecera_A320_200 = 'AAA A320-200';
writematrix(cabecera_A320_200,filename_write,'Sheet',1,'Range','B1')

filename_write = '../AIRCRAFT/AAA/AAA_Comparisson_AC.xlsx';
writecell(Sheet1,filename_write,'Sheet',1,'Range','A2:B77')

cabecera_SANAID = 'SANAID A320-200';
writematrix(cabecera_SANAID,filename_write,'Sheet',1,'Range','C1')

% write_cabeceras_comparisson_Excel(Sheet1)
   
Conversion_Der.OpFlightAlt	= Performance_preliminar.h;
Conversion_Der.OpDeltaTemperature	= Performance_preliminar.Temp;
Conversion_Der.OpFlightSpeedKTS = Performance_preliminar.V;
Conversion_Der.OpCurrentWeight	= Weight_tier.m_TOW;

rho_init = Performance_preliminar.rho;
S_ref = Geo_tier.S_ref;
g = conv_UNITS.g;

Conversion_Der.OpAlpha = TRIM_RESULTS.trim_alpha;
Conversion_Der.OpCLAirplane = Stab_Der.CL;

q_dynamic = 0.5*rho_init*(Performance_preliminar.V);
Lift = q_dynamic*Conversion_Der.OpCLAirplane*S_ref;

Conversion_Der.OpLoadFactor = Lift/Conversion_Der.OpCurrentWeight*g;
Conversion_Der.OpXcg = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
Conversion_Der.OpZcg = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;

Initial_Conditions = [Conversion_Der.OpFlightAlt;Conversion_Der.OpDeltaTemperature; Conversion_Der.OpFlightSpeedKTS;Conversion_Der.OpCurrentWeight;...
    Conversion_Der.OpAlpha; Conversion_Der.OpCLAirplane; Conversion_Der.OpLoadFactor; Conversion_Der.OpXcg; Conversion_Der.OpZcg];

writematrix(Initial_Conditions,filename_write,'Sheet',1,'Range','C2:C10')

if AC_CONFIGURATION.HTP == 1
    Conversion_Der.OpHorTailDownwashAngle = 0;
    Conversion_Der.OpHorTailPressRatioPowerOff	 = 1;   
    Conversion_Der.OpHorTailPressRatio	= 1;
    Tail_Conditions = [Conversion_Der.OpHorTailDownwashAngle; Conversion_Der.OpHorTailPressRatioPowerOff;Conversion_Der.OpHorTailPressRatio];
    writematrix(Tail_Conditions,filename_write,'Sheet',1,'Range','C11:C13')
end


if AC_CONFIGURATION.Vee == 1
    Conversion_Der.OpVeeTailDownwashAngle = 0;
    Conversion_Der.OpVeeTailPressRatioPowerOff	 = 1;   
    Conversion_Der.OpVeeTailPressRatio	= 1;
    Tail_Conditions = [Conversion_Der.OpVeeTailDownwashAngle; Conversion_Der.OpVeeTailPressRatioPowerOff;Conversion_Der.OpVeeTailPressRatio];
    writematrix(Tail_Conditions,filename_write,'Sheet',1,'Range','C11:C13')
end

% Propulsive Derivatives
Conversion_Der.OpCTx1	 = Stab_Der.CTx1;
Conversion_Der.OpCmT1	 = Stab_Der.CMt1;
Propulsive_Conditions = [Conversion_Der.OpCTx1; Conversion_Der.OpCmT1];

% Forward Speed Derivatives
Conversion_Der.OpCDu	 = Stab_Der.CDu;
Conversion_Der.OpCLu	 = Stab_Der.CLu;
Conversion_Der.OpCmu	 = Stab_Der.CMu;
Conversion_Der.OpCTxu	 = Stab_Der.CTxu;
Conversion_Der.OpCmTu	 = Stab_Der.CMtu;
ForwardSpeed_Conditions = [Conversion_Der.OpCDu;Conversion_Der.OpCLu;Conversion_Der.OpCmu;Conversion_Der.OpCTxu;Conversion_Der.OpCmTu];

% Angle of Attack Derivatives
Conversion_Der.OpCDAlpha  = Stab_Der.CD_alpha;
Conversion_Der.OpCLAlpha	 = Stab_Der.CL_alpha_ac;
Conversion_Der.OpCmAlpha	 = Stab_Der.CM_alpha_ac;
Conversion_Der.OpCmTAlpha	 = Stab_Der.CMtalpha;
AngleAttack_Conditions = [Conversion_Der.OpCDAlpha;Conversion_Der.OpCLAlpha;Conversion_Der.OpCmAlpha;Conversion_Der.OpCmTAlpha];

% Angle of Attack derivative Derivatives
Conversion_Der.OpCDAlphaDot  = Stab_Der.CD_alphapunto;
Conversion_Der.OpCLAlphaDot  = Stab_Der.CL_alphapunto;
Conversion_Der.OpCmAlphaDot  = Stab_Der.CM_alphapunto;
AngleAttackDot_Conditions = [Conversion_Der.OpCDAlphaDot;Conversion_Der.OpCLAlphaDot;Conversion_Der.OpCmAlphaDot];

% Pitch Rate Derivatives
Conversion_Der.OpCDq	 = Stab_Der.CDq;
Conversion_Der.OpCLq	 = Stab_Der.CLq;
Conversion_Der.OpCmq	 = Stab_Der.CMq;
PitchRate_Conditions = [Conversion_Der.OpCDq;Conversion_Der.OpCLq;Conversion_Der.OpCmq];

% Beta Derivatives
Conversion_Der.OpCyBeta  = Stab_Der.Cyb;
Conversion_Der.OpClBeta  = Stab_Der.Clb;
Conversion_Der.OpCnBeta  = Stab_Der.Cnb;
Conversion_Der.OpCnTBeta	 = Stab_Der.CNTb;
Beta_Conditions = [Conversion_Der.OpCyBeta;Conversion_Der.OpClBeta;Conversion_Der.OpCnBeta;Conversion_Der.OpCnTBeta];

% Beta dot Derivatives
Conversion_Der.OpCyBetadot	 = Stab_Der.Cybpunto;
Conversion_Der.OpClBetadot	 = Stab_Der.Clbpunto;
Conversion_Der.OpCnBetadot	 = Stab_Der.Cnbpunto;
BetaDot_Conditions = [Conversion_Der.OpCyBetadot;Conversion_Der.OpClBetadot;Conversion_Der.OpCnBetadot];

% Roll Rate Derivatives
Conversion_Der.OpCyp	 = Stab_Der.Cyp;
Conversion_Der.OpClp	 = Stab_Der.Clp;
Conversion_Der.OpCnp	 = Stab_Der.Cnp;
RollRate_Conditions = [Conversion_Der.OpCyp;Conversion_Der.OpClp;Conversion_Der.OpCnp];

% Yaw Rate Derivatives
Conversion_Der.OpCyr	 = Stab_Der.Cyr;
Conversion_Der.OpClr	 = Stab_Der.Clr;
Conversion_Der.OpCnr	 = Stab_Der.Cnr;
YawRate_Conditions = [Conversion_Der.OpCyr;Conversion_Der.OpClr;Conversion_Der.OpCnr];

Derivatives_Conditions = [Propulsive_Conditions;ForwardSpeed_Conditions;AngleAttack_Conditions;AngleAttackDot_Conditions;PitchRate_Conditions;Beta_Conditions;...
    BetaDot_Conditions;RollRate_Conditions;YawRate_Conditions];
writematrix(Derivatives_Conditions,filename_write,'Sheet',1,'Range','C14:C43')

if AC_CONFIGURATION.Vee == 1
    % Vtail incidence Derivative
    Conversion_Der.OpCDVeeTailIncidence = Stab_Der.CD_Vtail_i;
    Conversion_Der.OpCLVeeTailIncidence = Stab_Der.CL_Vtail_i;
    Conversion_Der.OpCmVeeTailIncidence = Stab_Der.CM_Vtail_i;
    
    Conversion_Der.OpCyVeeTailIncidence = Stab_Der.Cy_Vtail_i;
    Conversion_Der.OpCrollVeeTailIncidence = Stab_Der.Cl_Vtail_i;
    Conversion_Der.OpCnVeeTailIncidence = Stab_Der.Cn_Vtail_i;

    Conversion_Der.OpCDdeltaRuddervator = 0;
    Conversion_Der.OpCLdeltaRuddervatorZero = 0;
    Conversion_Der.OpCLdeltaRuddervator = 0;
    Conversion_Der.OpCmdeltaRuddervatorZero = 0;
    Conversion_Der.OpCmdeltaRuddervator = 0;

    Conversion_Der.OpCydeltaRuddervatorZero = Stab_Der.Cydeltarv;
    Conversion_Der.OpCydeltaRuddervator = Stab_Der.Cydeltarv;
    Conversion_Der.OpCndeltaRuddervatorZero = Stab_Der.Cndeltarv;
    Conversion_Der.OpCndeltaRuddervator = Stab_Der.Cndeltarv;
    Conversion_Der.OpCrolldeltaRuddervatorZero = Stab_Der.Cldeltarv;
    Conversion_Der.OpCrolldeltaRuddervator = Stab_Der.Cldeltarv;
    Conversion_Der.OpChBetaRuddervator	= Stab_Der.Ch_beta_deltarv;
    Conversion_Der.OpChdeltaRuddervator = Stab_Der.Ch_delta_deltarv;
    Conversion_Der.OpCLdeltaRuddervatorTab = Stab_Der.CL_deltarv_Tab;
    Conversion_Der.OpCmdeltaRuddervatorTab = Stab_Der.Cm_deltarv_Tab;
    Conversion_Der.OpChdeltaRuddervatorTab = Stab_Der.Ch_deltarv_Tab;
end

% if AC_CONFIGURATION.HTP == 1
%     Conversion_Der.OpCDHorTailIncidence = Stab_Der.CD_HTP_i;
%     Conversion_Der.OpCLHorTailIncidence = Stab_Der.CL_HTP_i;
%     Conversion_Der.OpCmHorTailIncidence = Stab_Der.CM_HTP_i;
%     % Longitudinal Rudder vator derivatives
%     Conversion_Der.OpCDdeltaElevator = Stab_Der.CD_delta_e;
%     Conversion_Der.OpCLdeltaElevatorZero = Stab_Der.CL_delta_e;
%     Conversion_Der.OpCLdeltaElevator = Stab_Der.CL_delta_e;
%     Conversion_Der.OpCmdeltaElevatorZero = Stab_Der.CM_delta_e;
%     Conversion_Der.OpCmdeltaElevator = Stab_Der.CM_delta_e;
%     Conversion_Der.OpChAlphaElevator = Stab_Der.Ch_alpha_e;
%     Conversion_Der.OpChdeltaElevator = Stab_Der.Ch_delta_e;
% end

% if AC_CONFIGURATION.VTP == 1
%     Conversion_Der.OpCyVertTailIncidence = Cy_VTP_i;
%     Conversion_Der.OpClVertTailIncidence = Cl_VTP_i;
%     Conversion_Der.OpCnVertTailIncidence = Cn_VTP_i;
%     Conversion_Der.OpCydeltaRudderZero = Cy_delta_r;
%     Conversion_Der.OpCydeltaRudder = Cy_delta_r;
%     Conversion_Der.OpCldeltaRudderZero = Cl_delta_r;
%     Conversion_Der.OpCldeltaRudder = Cl_delta_r;
%     Conversion_Der.OpCndeltaRudderZero = Cn_delta_r;
%     Conversion_Der.OpCndeltaRudder = Cn_delta_r;
%     Conversion_Der.OpChBetaRudder = Ch_betta_r;
%     Conversion_Der.OpChdeltaRudder = Ch_delta_r;
% end

% % Lateral Directional Rudder vator derivatives
% if AC_CONFIGURATION.d_rudder == 1
%     Conversion_Der.OpCydeltaRudderZero = Stab_Der.Cydeltarv;
%     Conversion_Der.OpCydeltaRudder = Stab_Der.Cydeltarv;
%     Conversion_Der.OpCndeltaRudderZero = Stab_Der.Cndeltarv;
%     Conversion_Der.OpCndeltaRudder = Stab_Der.Cndeltarv;    
%     Conversion_Der.OpChBetaRudder = Stab_Der.Ch_beta_deltarv;
%     Conversion_Der.OpChdeltaRuddervator = Stab_Der.Ch_delta_deltarv;
%     Conversion_Der.OpCLdeltaRuddervatorTab = Stab_Der.CL_deltarv_Tab;
%     Conversion_Der.OpCmdeltaRuddervatorTab = Stab_Der.Cm_deltarv_Tab;
%     Conversion_Der.OpChdeltaRudderv= Stab_Der.Ch_deltarv_Tab;
%     Conversion_Der.OpChBetaRudder = Stab_Der.Ch_deltarv_Tab;
% end

% Aileron Derivatives
Conversion_Der.OpCydeltaAileron  = Stab_Der.Cydeltaa;
Conversion_Der.OpCldeltaAileron  = Stab_Der.Cldeltaa;
Conversion_Der.OpCndeltaAileron  = Stab_Der.Cndeltaa;
% Hinge Derivatives
Conversion_Der.OpChAlphaAileron  = Stab_Der.Ch_alpha_deltaa;
Conversion_Der.OpChdeltaAileron  = Stab_Der.Ch_delta_deltaa;

AC_CONFIGURATION.d_spoiler = 0;
if AC_CONFIGURATION.d_spoiler == 1
    % Tab Derivatives
    Conversion_Der.OpCydeltaSpoiler = Stab_Der.Cy_delta_spo;
    Conversion_Der.OpCldeltaSpoiler = Stab_Der.Cl_delta_spo;
    Conversion_Der.OpCndeltaSpoiler = Stab_Der.Cn_delta_spo;
end

AC_CONFIGURATION.d_ail_tab = 0;
if AC_CONFIGURATION.d_ail_tab == 1
    % Tab Derivatives
    Conversion_Der.OpCydeltaAileronTab = Stab_Der.Cy_deltaa_Tab;
    Conversion_Der.OpCldeltaAileronTab= Stab_Der.Cl_deltaa_Tab;
    Conversion_Der.OpCndeltaAileronTab = Stab_Der.Cn_deltaa_Tab;
    Conversion_Der.OpChdeltaAileronTab = Stab_Der.Ch_deltaa_Tab;
end

Conversion_Der.OpCLzeroAirplaneClean = TRIM_RESULTS.CL0_ac;
Conversion_Der.OpCLzeroAirplane = TRIM_RESULTS.CL0_ac;
Conversion_Der.OpCmzerobarWingFus = Trim_ITER.CM0_w1_e_corrected;
Conversion_Der.OpCmzerobarNoEmpennage = Trim_ITER.CM0_w1_e_corrected;
Conversion_Der.OpCmzeroAirplane = TRIM_RESULTS.CM0_ac;
Trim_Conditions = [Conversion_Der.OpCLzeroAirplaneClean;Conversion_Der.OpCLzeroAirplane;Conversion_Der.OpCmzerobarWingFus;Conversion_Der.OpCmzerobarNoEmpennage;Conversion_Der.OpCmzeroAirplane];
writematrix(Trim_Conditions,filename_write,'Sheet',1,'Range','C73:C77')