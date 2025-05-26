function Aero_True = get_Aero_MILVUS_POLAR(Aero,Aero_TH,Geo_tier)

% Corrección CD0 CFD
Aero_True.CR.C_D0 = 0.099691137712358;   
Aero_True.CR.C_D1 = -0.062222749483185 ;
Aero_True.CR.C_D2 = 0.064419413009501;


% Take Off MILVUS 
Aero_True.TO.C_D0 = Aero_True.CR.C_D0;
Aero_True.TO.C_D1 = Aero_True.CR.C_D1;
Aero_True.TO.C_D2 = Aero_True.CR.C_D2;

% Take Off MILVUS 
Aero_True.LN.C_D0 = Aero_True.CR.C_D0;
Aero_True.LN.C_D1 = Aero_True.CR.C_D1;
Aero_True.LN.C_D2 = Aero_True.CR.C_D0;

CD0_ac = Aero_True.CR.C_D0;
CD1_ac = Aero_True.CR.C_D1;
CD2_ac = Aero_True.CR.C_D2;

Aero_True.CD0_ac = CD0_ac;
Aero_True.CD1_ac = CD1_ac;
Aero_True.CD2_ac = CD2_ac;
