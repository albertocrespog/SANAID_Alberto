function Aero_True = get_Aero_T75(Aero,Aero_TH,Geo_tier)

% Cruise configuration Tarsis 75 - AICIA
Aero_True.CR.C_D0 = 0.03802;
Aero_True.CR.C_D1 = -0.04747;
Aero_True.CR.C_D2 = 0.06970;

% Corrección CD0 experimentos
Aero_True.CR.C_D0 =  0.072916427615719;
Aero_True.CR.C_D1 = -0.04747;
Aero_True.CR.C_D2 = 0.06970;

% Take Off configuration Tarsis 75
Aero_True.TO.C_D0 = 0.1324;
Aero_True.TO.C_D1 = -0.1545;
Aero_True.TO.C_D2 = 0.1131;

% LANDING configuration Tarsis 75
Aero_True.LN.C_D0 = 0.1324;
Aero_True.LN.C_D1 = -0.1545;
Aero_True.LN.C_D2 = 0.1131;

% % Climb configuration Tarsis 75
% Aero_True.CLIMB.C_D0 = 0.03802;
% Aero_True.CLIMB.C_D1 = -0.04747;
% Aero_True.CLIMB.C_D2 = 0.06970;

% Coeficiente de sustentación máximo en configuración aerodinámica de despegue/aterrizaje: 1.675
% Coeficientes de resistencia en configuración aerodinámica de despegue/aterrizaje (en configuración con tren de aterrizaje):

correction_polar = 1;
if correction_polar == 1
    CD0_w1 = Aero.CD0_w1;
    CD1_w1 = Aero.CD1_w1;
    CD2_w1 = Aero.CD2_w1;

    CD0_w1T75 = 0.015385165212135;
    CD1_w1T75 = -0.025935913273452;
    CD2_w1T75 = 0.047001230383594;

    S_ref_T75 = 2.197620000000000;

    CD0_ac_NoWing = Aero_True.CR.C_D0 - CD0_w1T75;
    CD0_ac = CD0_ac_NoWing*S_ref_T75/Geo_tier.S_w1 + CD0_w1;

    CD1_ac_NoWing = Aero_True.CR.C_D1 - CD1_w1T75;
    CD1_ac = CD1_ac_NoWing*S_ref_T75/Geo_tier.S_w1 + CD1_w1;

    CD2_ac_NoWing = Aero_True.CR.C_D2 - CD2_w1T75;
    CD2_ac = CD2_ac_NoWing*S_ref_T75/Geo_tier.S_w1 + CD2_w1;

    Aero_True.CD0_ac = CD0_ac;
    Aero_True.CD1_ac = CD1_ac;
    Aero_True.CD2_ac = CD2_ac;
end