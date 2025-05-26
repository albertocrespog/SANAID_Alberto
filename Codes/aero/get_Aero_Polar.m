function Polar = get_Aero_Polar(Aero,Aero_TH,OUTPUT_read_XLSX,Performance,Geo_tier)

Flight_cruise = OUTPUT_read_XLSX.Performance_pre_flags.Flight_cruise; %
Flight_takeoff = OUTPUT_read_XLSX.Performance_pre_flags.Flight_takeoff; %
Flight_climb = OUTPUT_read_XLSX.Performance_pre_flags.Flight_climb; %

polar_model = OUTPUT_read_XLSX.Aerodynamic_Data_flags.polar_model;
switch polar_model
    case 1 % Composite Build Up Method
%         if Flight_cruise
            C_D0 = Aero_TH.CR.CD0;
            C_D1 = Aero_TH.CR.CD1;
            C_D2 = Aero_TH.CR.CD2;
%         elseif Flight_takeoff
            C_D0_TO = Aero_TH.TO.CD0;
            C_D1_TO = Aero_TH.TO.CD1;
            C_D2_TO = Aero_TH.TO.CD2;
%         elseif Flight_landing
            C_D0_LN = Aero_TH.LND.CD0;
            C_D1_LN = Aero_TH.LND.CD1;
            C_D2_LN = Aero_TH.LND.CD2;
%         elseif Flight_climb
%             C_D0_CLIMB = Aero_TH.CLIMB.CD0;
%             C_D1_CLIMB = Aero_TH.CLIMB.CD1;
%             C_D2_CLIMB = Aero_TH.CLIMB.CD2;
%         end
    case 2 % CFD or similar
%         if Flight_cruise
            C_D0 = Aero_TH.CR.CD0;
            C_D1 = Aero_TH.CR.CD1;
            C_D2 = Aero_TH.CR.CD2;
%         elseif Flight_takeoff
            C_D0_TO = Aero_TH.TO.CD0;
            C_D1_TO = Aero_TH.TO.CD1;
            C_D2_TO = Aero_TH.TO.CD2;
%         elseif Flight_landing
            C_D0_LN = Aero_TH.LND.CD0;
            C_D1_LN = Aero_TH.LND.CD1;
            C_D2_LN = Aero_TH.LND.CD2;
%         elseif Flight_climb
%             C_D0_CLIMB = Aero_TH.CLIMB.CD0;
%             C_D1_CLIMB = Aero_TH.CLIMB.CD1;
%             C_D2_CLIMB = Aero_TH.CLIMB.CD2;
%         end
    case 3 % Approximated values
        case_AC = OUTPUT_read_XLSX.AC_Data_flags.case_AC;
        Aero_approx = get_approximate_polar_values(case_AC,Aero,Aero_TH,Geo_tier);
        Polar.Aero_approx = Aero_approx;
%         if Flight_cruise
            C_D0 = Aero_approx.CR.C_D0;
            C_D1 = Aero_approx.CR.C_D1;
            C_D2 = Aero_approx.CR.C_D2;
%         elseif Flight_takeoff
            C_D0_TO = Aero_approx.TO.C_D0;
            C_D1_TO = Aero_approx.TO.C_D1;
            C_D2_TO = Aero_approx.TO.C_D2;
%         elseif Flight_climb
            C_D0_LN = Aero_approx.LN.C_D0;
            C_D1_LN = Aero_approx.LN.C_D1;
            C_D2_LN = Aero_approx.LN.C_D2;
%         end
    case 4 % Experimental results
%         if Flight_cruise
            C_D0 = Aero_TH.CR.CD0;
            C_D1 = Aero_TH.CR.CD1;
            C_D2 = Aero_TH.CR.CD2;
%         elseif Flight_takeoff
            C_D0_TO = Aero_TH.TO.CD0;
            C_D1_TO = Aero_TH.TO.CD1;
            C_D2_TO = Aero_TH.TO.CD2;
%         elseif Flight_landing
            C_D0_TO = Aero_TH.LND.CD0;
            C_D1_TO = Aero_TH.LND.CD1;
            C_D2_TO = Aero_TH.LND.CD2;
%         elseif Flight_climb
%             C_D0_CLIMB = Aero_TH.CLIMB.CD0;
%             C_D1_CLIMB = Aero_TH.CLIMB.CD1;
%             C_D2_CLIMB = Aero_TH.CLIMB.CD2;
%         end
    case 5 % Fusion CBM and Numerical results
%         if Flight_cruise
            C_D0 = Aero_TH.CR.CD0;
            C_D1 = Aero_TH.CR.CD1;
            C_D2 = Aero_TH.CR.CD2;
%         elseif Flight_takeoff
            C_D0_TO = Aero_TH.TO.CD0;
            C_D1_TO = Aero_TH.TO.CD1;
            C_D2_TO = Aero_TH.TO.CD2;
%         elseif Flight_landing
            C_D0_TO = Aero_TH.LND.CD0;
            C_D1_TO = Aero_TH.LND.CD1;
            C_D2_TO = Aero_TH.LND.CD2;
%         elseif Flight_climb
%             C_D0_CLIMB = Aero_TH.CLIMB.CD0;
%             C_D1_CLIMB = Aero_TH.CLIMB.CD1;
%             C_D2_CLIMB = Aero_TH.CLIMB.CD2;
%         end
end

% Stores Values for CRUISE
Polar.C_D0 = C_D0;
Polar.C_D1 = C_D1;
Polar.C_D2 = C_D2;
% Stores Values for TAKE OFF
Polar.C_D0_TO = C_D0_TO;
Polar.C_D1_TO = C_D1_TO;
Polar.C_D2_TO = C_D2_TO;
% Stores Values for LANDING
Polar.C_D0_LN = C_D0_LN;
Polar.C_D1_LN = C_D1_LN;
Polar.C_D2_LN = C_D2_LN;