function Aero_True = get_Aero_MILVUS(Aero,Aero_TH,Geo_tier)

% Cruise configuration inventado
Aero_True.CR.C_D0 = 0.03802;
Aero_True.CR.C_D1 = -0.04747;
Aero_True.CR.C_D2 = 0.06970;

C_D0_TO = Aero_approx.TO.C_D0 = Aero_True.CR.C_D0;
C_D1_TO = Aero_approx.TO.C_D1 = Aero_True.CR.C_D1;
C_D2_TO = Aero_approx.TO.C_D2;
%         elseif Flight_climb
C_D0_LN = Aero_approx.LN.C_D0 = Aero_True.CR.C_D0;
C_D1_LN = Aero_approx.LN.C_D1 = Aero_True.CR.C_D1;
C_D2_LN = Aero_approx.LN.C_D2 = Aero_True.CR.C_D2;
