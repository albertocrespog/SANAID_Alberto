% geometry King AIR
% NOse
x = 0;
y = 0 ;
z = 0 ;

%% Wing root
% wing root - LE
x_w1_CR_LE = 3854.796;
y_w1_CR_LE = 766.873;
z_w1_CR_LE = -331.779;

% wing root - TE
x_w1_CR_TE = 5738.067;
y_w1_CR_TE = 767.403;
z_w1_CR_TE = -460.726;

c_w1_CR = x_w1_CR_TE - x_w1_CR_LE;
i_w1_CR = atan((z_w1_CR_LE-z_w1_CR_TE)/(abs(x_w1_CR_LE-x_w1_CR_TE)));
y_w1_CR = (y_w1_CR_TE + y_w1_CR_LE)/2;


%% Wing root
% wing root - LE
x_w1_CR_LE = 3928.909;
y_w1_CR_LE = 0;
z_w1_CR_LE = -399.075;

% wing root - TE
x_w1_CR_TE = 6298.607;
y_w1_CR_TE = 0;
z_w1_CR_TE = -399.075;

c_w1_CR = x_w1_CR_TE - x_w1_CR_LE;
i_w1_CR = atan((z_w1_CR_LE-z_w1_CR_TE)/(abs(x_w1_CR_LE-x_w1_CR_TE)));
y_w1_CR = (y_w1_CR_TE + y_w1_CR_LE)/2;


%% Wing root kink 1
% wing kink 1 LE
x_w1_s1_LE = 3997.179;
y_w1_s1_LE = 980.608;
z_w1_s1_LE = -315.322;

% wing kink 1 TE
x_w1_s1_TE = 6298.607;
y_w1_s1_TE = 979.241;
z_w1_s1_TE = -395.107;

% Chord
c_w1_s1 = x_w1_s1_TE - x_w1_s1_LE;
% Incidence angle
i_w1_s1 = atan((z_w1_s1_LE-z_w1_s1_TE)/(abs(x_w1_s1_LE-x_w1_s1_TE)));
% Normalized wingspan
y_w1_s1 = (y_w1_s1_TE + y_w1_s1_LE)/2;
% Dihedral beteen this section and previous
dihedral_s1 = atan((z_w1_s1_LE-z_w1_CR_LE)/(abs(y_w1_s1-y_w1_CR)));
% Sweep beteen this section and previous
sweep_s1_LE = atan((x_w1_s1_LE-x_w1_CR_LE)/(abs(y_w1_s1-y_w1_CR)));

%% Wing root engine  side interior
% wing kink 2 LE
x_w1_s2_LE = 4096.921;
y_w1_s2_LE = 2179.343;
z_w1_s2_LE = -196.927;

% wing kink 2 TE
x_w1_s2_TE = x_w1_s1_TE; % Assume no TE sweep, same value to avid errors
y_w1_s2_TE = 2182.026;
z_w1_s2_TE = -258.082;

% Chord
c_w1_s2 = x_w1_s2_TE - x_w1_s2_LE;
% Incidence angle
i_w1_s2 = atan((z_w1_s2_LE-z_w1_s2_TE)/(abs(x_w1_s2_LE-x_w1_s2_TE)));
% Normalized wingspan
y_w1_s2 = (y_w1_s2_TE + y_w1_s2_LE)/2;
% Dihedral beteen this section and previous
dihedral_s2 = atan((z_w1_s2_LE-z_w1_s1_LE)/(abs(y_w1_s2-y_w1_s1)));
% Sweep beteen this section and previous
sweep_s2_LE = atan((x_w1_s2_LE-x_w1_s1_LE)/(abs(y_w1_s2-y_w1_s1)));

%% Wing root engine  side exterior
% wing kink 2 LE
x_w1_s3_LE = 4044.216;
y_w1_s3_LE = 2975.632;
z_w1_s3_LE = -101.903;

% wing kink 2 TE
x_w1_s3_TE = x_w1_s2_TE;
y_w1_s3_TE = 2936.587;
z_w1_s3_TE = -173.684;

% Chord
c_w1_s3 = x_w1_s3_TE - x_w1_s3_LE;
% Incidence angle
i_w1_s3 = atan((z_w1_s3_LE-z_w1_s3_TE)/(abs(x_w1_s3_LE-x_w1_s3_TE)));
% Normalized wingspan
y_w1_s3 = (y_w1_s3_TE + y_w1_s3_LE)/2;
% Dihedral beteen this section and previous
dihedral_s3 = atan((z_w1_s3_LE-z_w1_s2_LE)/(abs(y_w1_s3-y_w1_s2)));
% Sweep beteen this section and previous
sweep_s3_LE = atan((x_w1_s3_LE-x_w1_s2_LE)/(abs(y_w1_s3-y_w1_s2)));

%% Wing beginning of flap
% wing kink 2 LE
x_w1_s4_LE = 4118.90;
y_w1_s4_LE = 3130.92;
z_w1_s4_LE = -64.668;

% wing kink 2 TE
x_w1_s4_TE = x_w1_s3_TE;
y_w1_s4_TE = 3124.36;
z_w1_s4_TE = -150.562;

% Chord
c_w1_s4 = x_w1_s4_TE - x_w1_s4_LE;
% Incidence angle
i_w1_s4 = atan((z_w1_s4_LE-z_w1_s4_TE)/(abs(x_w1_s4_LE-x_w1_s4_TE)));
% Normalized wingspan
y_w1_s4 = (y_w1_s4_TE + y_w1_s4_LE)/2;
% Dihedral beteen this section and previous
dihedral_s4 = atan((z_w1_s4_LE-z_w1_s3_LE)/(abs(y_w1_s4-y_w1_s3)));
% Sweep beteen this section and previous
sweep_s4_LE = atan((x_w1_s4_LE-x_w1_s3_LE)/(abs(y_w1_s4-y_w1_s3)));


%% Wing Leading Edge Kink 
% wing kink LE
x_w1_s5_LE = 4357.415;
y_w1_s5_LE = 3685.877;
z_w1_s5_LE = -27.496;

% wing kink TE
x_w1_s5_TE = 6203.854;
y_w1_s5_TE = 3691.991;
z_w1_s5_TE = -87.319;

% Chord
c_w1_s5 = x_w1_s5_TE - x_w1_s5_LE;
% Incidence angle
i_w1_s5 = atan((z_w1_s5_LE-z_w1_s5_TE)/(abs(x_w1_s5_LE-x_w1_s5_TE)));
% Normalized wingspan
y_w1_s5 = (y_w1_s5_TE + y_w1_s5_LE)/2;
% Dihedral beteen this section and previous
dihedral_s5 = atan((z_w1_s5_LE-z_w1_s4_LE)/(abs(y_w1_s5-y_w1_s4)));
% Sweep beteen this section and previous
sweep_s5_LE = atan((x_w1_s5_LE-x_w1_s4_LE)/(abs(y_w1_s5-y_w1_s4)));


%% Wing rLeading Edge Kink 
%% flap/aileron
x_w1_s6_LE = 4413.917;
y_w1_s6_LE = 4748.28;
z_w1_s6_LE = 95.311;

% wing kink TE
x_w1_s6_TE = 6036.32;
y_w1_s6_TE = y_w1_s6_LE;
z_w1_s6_TE = 34.24;

% Chord
c_w1_s6 = x_w1_s6_TE - x_w1_s6_LE;
% Incidence angle
i_w1_s6 = atan((z_w1_s6_LE-z_w1_s6_TE)/(abs(x_w1_s6_LE-x_w1_s6_TE)));
% Normalized wingspan
y_w1_s6 = (y_w1_s6_TE + y_w1_s6_LE)/2;
% Dihedral beteen this section and previous
dihedral_s6 = atan((z_w1_s6_LE-z_w1_s5_LE)/(abs(y_w1_s6-y_w1_s5)));
% Sweep beteen this section and previous
sweep_s6_LE = atan((x_w1_s6_LE-x_w1_s5_LE)/(abs(y_w1_s6-y_w1_s5)));


%% Wing rLeading Edge Kink 
%% aileron Trim Tab - Outer
x_w1_s6B_LE = 4469.615;
y_w1_s6B_LE = 5753.624;
z_w1_s6B_LE = 206.205;

% wing kink TE
x_w1_s6B_TE = 5876.946;
y_w1_s6B_TE = 5751.979;
z_w1_s6B_TE = 150.206;

% Chord
c_w1_s6B = x_w1_s6B_TE - x_w1_s6B_LE;
% Incidence angle
i_w1_s6B = atan((z_w1_s6B_LE-z_w1_s6B_TE)/(abs(x_w1_s6B_LE-x_w1_s6B_TE)));
% Normalized wingspan
y_w1_s6B = (y_w1_s6B_TE + y_w1_s6B_LE)/2;
% Dihedral beteen this section and previous
dihedral_s6B = atan((z_w1_s6B_LE-z_w1_s6_LE)/(abs(y_w1_s6B-y_w1_s6)));
% Sweep beteen this section and previous
sweep_s6B_LE = atan((x_w1_s6B_LE-x_w1_s6_LE)/(abs(y_w1_s6B-y_w1_s6)));


%% Wing rLeading Edge Kink 
%% aileron out
x_w1_s7_LE = 4572.013;
y_w1_s7_LE = 7761.13;
z_w1_s7_LE = 95.311;

% wing kink TE
x_w1_s7_TE = 5560.25;
y_w1_s7_TE = y_w1_s7_LE;
z_w1_s7_TE = 34.24;

% Chord
c_w1_s7 = x_w1_s7_TE - x_w1_s7_LE;
% Incidence angle
i_w1_s7 = atan((z_w1_s7_LE-z_w1_s7_TE)/(abs(x_w1_s7_LE-x_w1_s7_TE)));
% Normalized wingspan
y_w1_s7 = (y_w1_s7_TE + y_w1_s7_LE)/2;
% Dihedral beteen this section and previous
dihedral_s7 = atan((z_w1_s7_LE-z_w1_s6_LE)/(abs(y_w1_s7-y_w1_s6)));
% Sweep beteen this section and previous
sweep_s7_LE = atan((x_w1_s7_LE-x_w1_s6_LE)/(abs(y_w1_s7-y_w1_s6)));

%% Wing tip chord
% Wing tip chord LE
x_w1_CT_LE = 4587.58;
y_w1_CT_LE = 8097.97;
z_w1_CT_LE = 442.16;

% wing kink TE
x_w1_CT_TE = 5507.64;
y_w1_CT_TE = 8098.66;
z_w1_CT_TE = 417.51;

% Chord
c_w1_CT = x_w1_CT_TE - x_w1_CT_LE;
% Incidence angle
i_w1_CT = atan((z_w1_CT_LE-z_w1_CT_TE)/(abs(x_w1_CT_LE-x_w1_CT_TE)));
% Normalized wingspan
y_w1_CT = (y_w1_CT_TE + y_w1_CT_LE)/2;
% Dihedral beteen this section and previous
dihedral_CT = atan((z_w1_CT_LE-z_w1_s7_LE)/(abs(y_w1_CT-y_w1_s7)));
% Sweep beteen this section and previous
sweep_CT_LE = atan((x_w1_CT_LE-x_w1_s7_LE)/(abs(y_w1_CT-y_w1_s7)));

%% For Flow 5
%% Winglet - 1
x_w1_winglet1_LE = 4590.434;
y_w1_winglet1_LE = 8177.11;
z_w1_winglet1_LE = 452.514;

% wing kink TE
x_w1_winglet1_TE = 5496.706;
y_w1_winglet1_TE = 8181.17;
z_w1_winglet1_TE = 428.335;

% Chord
c_w1_winglet1 = x_w1_winglet1_TE - x_w1_winglet1_LE;
% Incidence angle
i_w1_winglet1 = atan((z_w1_winglet1_LE-z_w1_winglet1_TE)/(abs(x_w1_winglet1_LE-x_w1_winglet1_TE)));
% Normalized wingspan
y_w1_winglet1 = (y_w1_winglet1_TE + y_w1_winglet1_LE)/2;
% Dihedral beteen this section and previous
dihedral_winglet1 = atan((z_w1_winglet1_LE-z_w1_CT_LE)/(abs(y_w1_winglet1-y_w1_CT)));
% Sweep beteen this section and previous
sweep_winglet1_LE = atan((x_w1_winglet1_LE-x_w1_CT_LE)/(abs(y_w1_winglet1-y_w1_CT)));

%% Winglet - 2
x_w1_winglet2_LE = 4598.316;
y_w1_winglet2_LE = 8268.498;
z_w1_winglet2_LE = 463.825;

% wing kink TE
x_w1_winglet2_TE = 5482.905;
y_w1_winglet2_TE = 8268.926;
z_w1_winglet2_TE = 440.428;

% Chord
c_w1_winglet2 = x_w1_winglet2_TE - x_w1_winglet2_LE;
% Incidence angle
i_w1_winglet2 = atan((z_w1_winglet2_LE-z_w1_winglet2_TE)/(abs(x_w1_winglet2_LE-x_w1_winglet2_TE)));
% Normalized wingspan
y_w1_winglet2 = (y_w1_winglet2_TE + y_w1_winglet2_LE)/2;
% Dihedral beteen this section and previous
dihedral_winglet2 = atan((z_w1_winglet2_LE-z_w1_winglet1_LE)/(abs(y_w1_winglet2-y_w1_winglet1)));
% Sweep beteen this section and previous
sweep_winglet2_LE = atan((x_w1_winglet2_LE-x_w1_winglet1_LE)/(abs(y_w1_winglet2-y_w1_winglet1)));

%% Winglet - 3
x_w1_winglet3_LE = 4690.448;
y_w1_winglet3_LE = 8510.529;
z_w1_winglet3_LE = 503.181;

% wing kink TE
x_w1_winglet3_TE = 5436.155;
y_w1_winglet3_TE = 8521.997;
z_w1_winglet3_TE = 481.578;

% Chord
c_w1_winglet3 = x_w1_winglet3_TE - x_w1_winglet3_LE;
% Incidence angle
i_w1_winglet3 = atan((z_w1_winglet3_LE-z_w1_winglet3_TE)/(abs(x_w1_winglet3_LE-x_w1_winglet3_TE)));
% Normalized wingspan
y_w1_winglet3 = (y_w1_winglet3_TE + y_w1_winglet3_LE)/2;
% Dihedral beteen this section and previous
dihedral_winglet3 = atan((z_w1_winglet3_LE-z_w1_winglet2_LE)/(abs(y_w1_winglet3-y_w1_winglet2)));
% Sweep beteen this section and previous
sweep_winglet3_LE = atan((x_w1_winglet3_LE-x_w1_winglet2_LE)/(abs(y_w1_winglet3-y_w1_winglet2)));

%% Winglet - 3B
x_w1_winglet3B_LE = 4792.62;
y_w1_winglet3B_LE = 8617.57;
z_w1_winglet3B_LE = 525.83;

% wing kink TE
x_w1_winglet3B_TE = 5404.157;
y_w1_winglet3B_TE = 8617.987;
z_w1_winglet3B_TE = 504.733;

% Chord
c_w1_winglet3B = x_w1_winglet3B_TE - x_w1_winglet3B_LE;
% Incidence angle
i_w1_winglet3B = atan((z_w1_winglet3B_LE-z_w1_winglet3B_TE)/(abs(x_w1_winglet3B_LE-x_w1_winglet3B_TE)));
% Normalized wingspan
y_w1_winglet3B = (y_w1_winglet3B_TE + y_w1_winglet3B_LE)/2;
% Dihedral beteen this section and previous
dihedral_winglet3B = atan((z_w1_winglet3B_LE-z_w1_winglet3_LE)/(abs(y_w1_winglet3B-y_w1_winglet3)));
% Sweep beteen this section and previous
sweep_winglet3B_LE = atan((x_w1_winglet3B_LE-x_w1_winglet3_LE)/(abs(y_w1_winglet3B-y_w1_winglet3)));


%% Winglet - 4
x_w1_winglet4_LE = 5044.581;
y_w1_winglet4_LE = 8707.037;
z_w1_winglet4_LE = 653.499;

% wing kink TE
x_w1_winglet4_TE = 5422.12;
y_w1_winglet4_TE = 8711.363;
z_w1_winglet4_TE = 653.594;

chord = 296.199;

% Chord
c_w1_winglet4 = x_w1_winglet4_TE - x_w1_winglet4_LE;
% Incidence angle
i_w1_winglet4 = atan((z_w1_winglet4_LE-z_w1_winglet4_TE)/(abs(x_w1_winglet4_LE-x_w1_winglet4_TE)));
% Normalized wingspan
y_w1_winglet4 = (y_w1_winglet4_TE + y_w1_winglet4_LE)/2;
% Dihedral beteen this section and previous
dihedral_winglet4 = atan((z_w1_winglet4_LE-z_w1_winglet3_LE)/(abs(y_w1_winglet4-y_w1_winglet3)));
% Sweep beteen this section and previous
sweep_winglet4_LE = atan((x_w1_winglet4_LE-x_w1_winglet3_LE)/(abs(y_w1_winglet4-y_w1_winglet3)));

%% Winglet - 5
x_w1_winglet5_LE = 5355.466;
y_w1_winglet5_LE = 8748.292;
z_w1_winglet5_LE = 1086.116;

% wing kink TE
x_w1_winglet5_TE = 5532.249;
y_w1_winglet5_TE = 8747.652;
z_w1_winglet5_TE = 1086.314;

% Chord
c_w1_winglet5 = x_w1_winglet5_TE - x_w1_winglet5_LE;
% Incidence angle
i_w1_winglet5 = atan((z_w1_winglet5_LE-z_w1_winglet5_TE)/(abs(x_w1_winglet5_LE-x_w1_winglet5_TE)));
% Normalized wingspan
y_w1_winglet5 = (y_w1_winglet5_TE + y_w1_winglet5_LE)/2;
% Dihedral beteen this section and previous
dihedral_winglet5 = atan((z_w1_winglet5_LE-z_w1_winglet4_LE)/(abs(y_w1_winglet5-y_w1_winglet4)));
% Sweep beteen this section and previous
sweep_winglet5_LE = atan((x_w1_winglet5_LE-x_w1_winglet4_LE)/(abs(y_w1_winglet5-y_w1_winglet4)));

% Recolection of geometry
% Wing for the Beechcraft King Air 350
x_w1 = [x_w1_CR_LE x_w1_s1_LE x_w1_s2_LE x_w1_s3_LE x_w1_s4_LE x_w1_s5_LE x_w1_s6_LE x_w1_s6B_LE ...
    x_w1_s7_LE x_w1_CT_LE x_w1_winglet1_LE x_w1_winglet2_LE x_w1_winglet3_LE x_w1_winglet3B_LE x_w1_winglet4_LE x_w1_winglet5_LE]/1000
y_w1 = [y_w1_CR y_w1_s1 y_w1_s2 y_w1_s3 y_w1_s4 y_w1_s5 y_w1_s6 y_w1_s6B y_w1_s7 y_w1_CT...
    y_w1_winglet1 y_w1_winglet2 y_w1_winglet3 y_w1_winglet3B y_w1_winglet4 y_w1_winglet5]/1000
z_w1 = [z_w1_CR_LE z_w1_s1_LE z_w1_s2_LE z_w1_s3_LE z_w1_s4_LE z_w1_s5_LE z_w1_s6_LE z_w1_s6B_LE ...
    z_w1_s7_LE z_w1_CT_LE z_w1_winglet1_LE z_w1_winglet2_LE z_w1_winglet3_LE z_w1_winglet3B_LE z_w1_winglet4_LE z_w1_winglet5_LE]/1000

%% Correction for the entry of data on FLOW5 according to the dihedral
b1 = sqrt((y_w1_s1-y_w1_CR)^2 + (z_w1_s1_LE-z_w1_CR_LE)^2)/1000;
b2 = sqrt((y_w1_s2-y_w1_s1)^2 + (z_w1_s2_LE-z_w1_s1_LE)^2)/1000;
b3 = sqrt((y_w1_s3-y_w1_s2)^2 + (z_w1_s3_LE-z_w1_s2_LE)^2)/1000;
b4 = sqrt((y_w1_s4-y_w1_s3)^2 + (z_w1_s4_LE-z_w1_s3_LE)^2)/1000;
b5 = sqrt((y_w1_s5-y_w1_s4)^2 + (z_w1_s5_LE-z_w1_s4_LE)^2)/1000;
b6 = sqrt((y_w1_s6-y_w1_s5)^2 + (z_w1_s6_LE-z_w1_s5_LE)^2)/1000;
b6B = sqrt((y_w1_s6B-y_w1_s6)^2 + (z_w1_s6B_LE-z_w1_s6_LE)^2)/1000;
b7 = sqrt((y_w1_s7-y_w1_s6B)^2 + (z_w1_s7_LE-z_w1_s6B_LE)^2)/1000;
bCT = sqrt((y_w1_CT-y_w1_s7)^2 + (z_w1_CT_LE-z_w1_s7_LE)^2)/1000;
bwing1 = sqrt((y_w1_winglet1-y_w1_CT)^2 + (z_w1_winglet1_LE-z_w1_CT_LE)^2)/1000;
bwing2 = sqrt((y_w1_winglet2-y_w1_winglet1)^2 + (z_w1_winglet2_LE-z_w1_winglet1_LE)^2)/1000;
bwing3 = sqrt((y_w1_winglet3-y_w1_winglet2)^2 + (z_w1_winglet3_LE-z_w1_winglet2_LE)^2)/1000;
bwing3B = sqrt((y_w1_winglet3B-y_w1_winglet3)^2 + (z_w1_winglet3B_LE-z_w1_winglet3_LE)^2)/1000;
bwing4 = sqrt((y_w1_winglet4-y_w1_winglet3B)^2 + (z_w1_winglet4_LE-z_w1_winglet3B_LE)^2)/1000;
bwing5 = sqrt((y_w1_winglet5-y_w1_winglet4)^2 + (z_w1_winglet5_LE-z_w1_winglet4_LE)^2)/1000;

span_FLOW_w1 = [ 0, b1, b2, b3, b4, b5, b6, b6B, b7, bCT, bwing1, bwing2, bwing3, bwing3B, bwing4, bwing5]

b_WING(1) = 0;
N = length(span_FLOW_w1);
for i=1:N-1
    b_WING(i+1) = b_WING(i) + span_FLOW_w1(i+1);
end

% New wingspan modified as the wingspan along the surface
b_WING

chord_w1 = [c_w1_CR c_w1_s1 c_w1_s2 c_w1_s3 c_w1_s4 c_w1_s5 c_w1_s6 c_w1_s6B ...
    c_w1_s7 c_w1_CT c_w1_winglet1 c_w1_winglet2 c_w1_winglet3 c_w1_winglet3B c_w1_winglet4 c_w1_winglet5]/1000
offset_w1 = ([x_w1_CR_LE x_w1_s1_LE x_w1_s2_LE x_w1_s3_LE x_w1_s4_LE x_w1_s5_LE x_w1_s6_LE x_w1_s6B_LE ...
    x_w1_s7_LE x_w1_CT_LE x_w1_winglet1_LE x_w1_winglet2_LE x_w1_winglet3_LE x_w1_winglet3B_LE x_w1_winglet4_LE x_w1_winglet5_LE] - x_w1_CR_LE)/1000
dihedral_w1 =[0 dihedral_s1 dihedral_s2 dihedral_s3  dihedral_s4 dihedral_s5 dihedral_s6 dihedral_s6B dihedral_s7 dihedral_CT...
    dihedral_winglet1 dihedral_winglet2 dihedral_winglet3 dihedral_winglet3B dihedral_winglet4 dihedral_winglet5]*180/pi
sweep_w1 = [0 sweep_s1_LE sweep_s2_LE sweep_s3_LE sweep_s4_LE sweep_s5_LE sweep_s6_LE sweep_s6B_LE sweep_s7_LE sweep_CT_LE...
    sweep_winglet1_LE sweep_winglet2_LE sweep_winglet3_LE sweep_winglet3B_LE sweep_winglet4_LE sweep_winglet5_LE]*180/pi

% Recolection of geometry
% Wing for the Beechcraft King Air 200B
x_w1 = [x_w1_CR_LE x_w1_s1_LE x_w1_s2_LE x_w1_s3_LE x_w1_s4_LE x_w1_s5_LE x_w1_s6_LE x_w1_s6B_LE ...
    x_w1_s7_LE x_w1_CT_LE x_w1_winglet1_LE x_w1_winglet2_LE]/1000
y_w1 = [y_w1_CR y_w1_s1 y_w1_s2 y_w1_s3 y_w1_s4 y_w1_s5 y_w1_s6 y_w1_s6B y_w1_s7 y_w1_CT...
    y_w1_winglet1 y_w1_winglet2 ]/1000
z_w1 = [z_w1_CR_LE z_w1_s1_LE z_w1_s2_LE z_w1_s3_LE z_w1_s4_LE z_w1_s5_LE z_w1_s6_LE z_w1_s6B_LE ...
    z_w1_s7_LE z_w1_CT_LE z_w1_winglet1_LE z_w1_winglet2_LE]/1000

chord_w1 = [c_w1_CR c_w1_s1 c_w1_s2 c_w1_s3 c_w1_s4 c_w1_s5 c_w1_s6 c_w1_s6B ...
    c_w1_s7 c_w1_CT c_w1_winglet1 c_w1_winglet2]/1000
offset_w1 = ([x_w1_CR_LE x_w1_s1_LE x_w1_s2_LE x_w1_s3_LE x_w1_s4_LE x_w1_s5_LE x_w1_s6_LE x_w1_s6B_LE ...
    x_w1_s7_LE x_w1_CT_LE x_w1_winglet1_LE x_w1_winglet2_LE ] - x_w1_CR_LE)/1000
dihedral_w1 =[0 dihedral_s1 dihedral_s2 dihedral_s3  dihedral_s4 dihedral_s5 dihedral_s6 dihedral_s6B dihedral_s7 dihedral_CT...
    dihedral_winglet1 dihedral_winglet2]*180/pi
sweep_w1 = [0 sweep_s1_LE sweep_s2_LE sweep_s3_LE sweep_s4_LE sweep_s5_LE sweep_s6_LE sweep_s6B_LE sweep_s7_LE sweep_CT_LE...
    sweep_winglet1_LE sweep_winglet2_LE]*180/pi

Geo_tier.cR_w1 = c_w1_CR/1000;
Geo_tier.cB_k1_w1 = c_w1_s3/1000;
Geo_tier.cB_k2_w1 = c_w1_s4/1000;
Geo_tier.cT_w1 = c_w1_winglet2/1000;


%% Controsl surfaces WING
%% FLAP
% Flap LE inboard
x_w1_flap_y1_LE = 5866.933;
y_w1_flap_y1_LE = y_w1_s4;

% Flap TE inboard
x_w1_flap_y1_TE = x_w1_s4_TE;
y_w1_flap_y1_TE = y_w1_s4;
z_w1_flap_y1_TE = z_w1_s4_TE;

% Chord of control surface
c_w1_flap_y1 = x_w1_flap_y1_TE - x_w1_flap_y1_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_flap_y1_c_s4 = c_w1_flap_y1/c_w1_s4;
% Delta z location of wing
Delta_z_w1_y1_LE = (z_w1_s4_LE - z_w1_s4_TE)*c_flap_y1_c_s4;
z_w1_flap_y1_LE = z_w1_s4_TE + Delta_z_w1_y1_LE;

% Flap LE outboard
x_w1_flap_y2_LE = 5733.936;
y_w1_flap_y2_LE = y_w1_s6;

% Flap TE outboard
x_w1_flap_y2_TE = x_w1_s6_TE;
y_w1_flap_y2_TE = y_w1_s6;
z_w1_flap_y2_TE = z_w1_s6_TE;

% Chord of control surface
c_w1_flap_y2 = x_w1_flap_y2_TE - x_w1_flap_y2_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_flap_y2_c_s6 = c_w1_flap_y2/c_w1_s6;
% Delta z location of wing
Delta_z_w1_y2_LE = (z_w1_s6_LE - z_w1_s6_TE)*c_flap_y2_c_s6;
z_w1_flap_y2_LE = z_w1_s6_TE + Delta_z_w1_y2_LE;


b_flap = y_w1_flap_y2_TE - y_w1_flap_y1_TE;

b_w1_e = (y_w1_winglet2 - y_w1_CR);
Geo_tier.K_y1_flap_w1 = (y_w1_flap_y1_TE - y_w1_CR)/b_w1_e;
Geo_tier.K_y2_flap_w1 = (y_w1_flap_y2_TE - y_w1_CR)/b_w1_e;

% % Flap LE outboard
% 5733.936 mm
% 4741.283 mm
% 93.748 mm
% % Flap TE outboard
% 6036.318 mm
% 4747.384 mm
% 33.982 mm

%% Aileron
% Aileron LE inboard
x_w1_ail_y1_LE = 5511.184;
y_w1_ail_y1_LE = y_w1_s6;

% Aileron TE inboard
x_w1_ail_y1_TE = x_w1_s6_TE;
y_w1_ail_y1_TE = y_w1_s6;
z_w1_ail_y1_TE = z_w1_s6_TE;

% Chord of control surface
c_w1_ail_y1 = x_w1_ail_y1_TE - x_w1_ail_y1_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_ail_y1_c_s6 = c_w1_ail_y1/c_w1_s6;
% Delta z location of wing
Delta_z_w1_y1_LE = (z_w1_s6_LE - z_w1_s6_TE)*c_ail_y1_c_s6;
z_w1_ail_y1_LE = z_w1_s6_TE + Delta_z_w1_y1_LE;

% Aileron LE outboard
x_w1_ail_y2_LE = 5238.682;
y_w1_ail_y2_LE = y_w1_s7;

% Aileron TE outboard
x_w1_ail_y2_TE = x_w1_s7_TE;
y_w1_ail_y2_TE = y_w1_s7;
z_w1_ail_y2_TE = z_w1_s7_TE;

% Chord of control surface
c_w1_ail_y2 = x_w1_ail_y2_TE - x_w1_ail_y2_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_ail_y2_c_s7 = c_w1_ail_y2/c_w1_s7;
% Delta z location of wing
Delta_z_w1_y2_LE = (z_w1_s7_LE - z_w1_s7_TE)*c_ail_y2_c_s7;
z_w1_ail_y2_LE = z_w1_s7_TE + Delta_z_w1_y2_LE;

b_ail = y_w1_ail_y2_TE - y_w1_ail_y1_TE;

b_w1_e = (y_w1_winglet2 - y_w1_CR);
Geo_tier.K_y1_ail_w1 = (y_w1_ail_y1_TE - y_w1_CR)/b_w1_e;
Geo_tier.K_y2_ail_w1 = (y_w1_ail_y2_TE - y_w1_CR)/b_w1_e;

cf_ail_tmp = [c_ail_y1_c_s6, c_ail_y2_c_s7];
Geo_tier.cf_ail = mean(cf_ail_tmp) ;
 
%% Aileron - TRIM-TAB
% LE inboard
x_w1_ail_TT_y1_LE = 5927.904;
y_w1_ail_TT_y1_LE = y_w1_s6;

% Aileron TE inboard
x_w1_ail_TT_y1_TE = x_w1_s6_TE;
y_w1_ail_TT_y1_TE = y_w1_s6;
z_w1_ail_TT_y1_TE = z_w1_s6_TE;

% Chord of control surface
c_w1_ail_TT_y1 = x_w1_ail_TT_y1_TE - x_w1_ail_TT_y1_LE;
cf_TT_ail_y1 = c_w1_ail_TT_y1/c_w1_ail_y1;

% Ratio con constrol surface chord with respecto to surface chord
c_ail_TT_y1_c_s6 = c_w1_ail_TT_y1/c_w1_s6;
% Delta z location of wing
Delta_z_w1_y1_LE = (z_w1_s6_LE - z_w1_s6_TE)*c_ail_TT_y1_c_s6;
z_w1_ail_TT_y1_LE = z_w1_s6_TE + Delta_z_w1_y1_LE;

% LE outboard
x_w1_ail_TT_y2_LE = 5783.454;
y_w1_ail_TT_y2_LE = y_w1_s6B;

% TE outboard
x_w1_ail_TT_y2_TE = x_w1_s6B_TE;
y_w1_ail_TT_y2_TE = y_w1_s6B_TE;
z_w1_ail_TT_y2_TE = z_w1_s6B_TE;

% Chord of control surface
c_w1_ail_TT_y2 = x_w1_ail_TT_y2_TE - x_w1_ail_TT_y2_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_ail_TT_y2_c_s6B = c_w1_ail_TT_y2/c_w1_s6B;
cf_TT_ail_y2 = c_w1_ail_TT_y2/c_w1_ail_y2;

cf_ail_TT_tmp = [cf_TT_ail_y1, cf_TT_ail_y2];
Geo_tier.cf_TT_ail = mean(cf_ail_TT_tmp) ;

% Delta z location of wing
Delta_z_w1_y2_LE = (z_w1_s6B_LE - z_w1_s6B_TE)*c_ail_TT_y2_c_s6B;
z_w1_ail_TT_y2_LE = z_w1_s6B_TE + Delta_z_w1_y2_LE;

b_ail_TT = y_w1_ail_TT_y2_TE - y_w1_ail_TT_y1_TE;

b_w1_e = (y_w1_winglet2 - y_w1_CR);
Geo_tier.K_y1_ail_TT_w1 = (y_w1_ail_TT_y1_TE - y_w1_CR)/b_w1_e;
Geo_tier.K_y2_ail_TT_w1 = (y_w1_ail_TT_y2_TE - y_w1_CR)/b_w1_e;

cf_ail_TT_tmp = [cf_TT_ail_y2, cf_TT_ail_y2];
Geo_tier.cf_ail = mean(cf_ail_TT_tmp) ;


%% HTP
%% Wing root
% wing root - LE
x_HTP_CR_LE = 12223.347;
y_HTP_CR_LE = 0;
z_HTP_CR_LE = 3096.483;

% wing root - TE
x_HTP_CR_TE = 13742.864;
y_HTP_CR_TE = 0;
z_HTP_CR_TE = 3096.483;

c_HTP_CR = x_HTP_CR_TE - x_HTP_CR_LE;
i_HTP_CR = atan((z_HTP_CR_LE-z_HTP_CR_TE)/(abs(x_HTP_CR_LE-x_HTP_CR_TE)));
y_HTP_CR = (y_HTP_CR_TE + y_HTP_CR_LE)/2;

%% HTP knk 1. inner elevator
% wing root - LE
x_HTP_s1_LE = 12279.335;
y_HTP_s1_LE = 140.18;
z_HTP_s1_LE = z_HTP_CR_LE;

% wing root - TE
x_HTP_s1_TE = 13758.15;
y_HTP_s1_TE = 140.18;
z_HTP_s1_TE = z_HTP_CR_LE;

% Chord
c_HTP_s1 = x_HTP_s1_TE - x_HTP_s1_LE;
% Incidence angle
i_HTP_s1 = atan((z_HTP_s1_LE-z_HTP_s1_TE)/(abs(x_HTP_s1_LE-x_HTP_s1_TE)));
% Normalized wingspan
y_HTP_s1 = (y_HTP_s1_TE + y_HTP_s1_LE)/2;
% Dihedral beteen this section and previous
dihedral_s1 = atan((z_HTP_s1_LE-z_HTP_CR_LE)/(abs(y_HTP_s1-y_HTP_CR)));
% Sweep beteen this section and previous
sweep_s1_LE = atan((x_HTP_s1_LE-x_HTP_CR_LE)/(abs(y_HTP_s1-y_HTP_CR)));

%% HTP kink 2 inner trim tab
% wing root - LE
x_HTP_s2_LE = 12661.092;
y_HTP_s2_LE = 1134.87;
z_HTP_s2_LE = z_HTP_CR_LE;

% wing root - TE
x_HTP_s2_TE = 13876.81;
y_HTP_s2_TE = 1134.87;
z_HTP_s2_TE = z_HTP_CR_LE;

% Chord
c_HTP_s2 = x_HTP_s2_TE - x_HTP_s2_LE;
% Incidence angle
i_HTP_s2 = atan((z_HTP_s2_LE-z_HTP_s2_TE)/(abs(x_HTP_s2_LE-x_HTP_s2_TE)));
% Normalized wingspan
y_HTP_s2 = (y_HTP_s2_TE + y_HTP_s2_LE)/2;
% Dihedral beteen this section and previous
dihedral_s2 = atan((z_HTP_s2_LE-z_HTP_s1_LE)/(abs(y_HTP_s2-y_HTP_s1)));
% Sweep beteen this section and previous
sweep_s2_LE = atan((x_HTP_s2_LE-x_HTP_s1_LE)/(abs(y_HTP_s2-y_HTP_s1)));


%% HTP tip
% wing root - LE
x_HTP_CT_LE = 13421.60;
y_HTP_CT_LE = 2865.253;
z_HTP_CT_LE = z_HTP_CR_LE;

% wing root - TE
x_HTP_CT_TE = 14071.9875;
y_HTP_CT_TE = 2885.558;
z_HTP_CT_TE = z_HTP_CR_LE;

% Chord
c_HTP_CT = x_HTP_CT_TE - x_HTP_CT_LE;
% Incidence angle
i_HTP_CT = atan((z_HTP_CT_LE-z_HTP_CT_TE)/(abs(x_HTP_CT_LE-x_HTP_CT_TE)));
% Normalized wingspan
y_HTP_CT = (y_HTP_CT_TE + y_HTP_CT_LE)/2;
% Dihedral beteen this section and previous
dihedral_CT = atan((z_HTP_CT_LE-z_HTP_CR_LE)/(abs(y_HTP_CT-y_HTP_s2)));
% Sweep beteen this section and previous
sweep_CT_LE = atan((x_HTP_CT_LE-x_HTP_CR_LE)/(abs(y_HTP_CT-y_HTP_s2)));

x_HTP = [x_HTP_CR_LE x_HTP_s1_LE x_HTP_s2_LE y_HTP_CT_LE]/1000
y_HTP = [y_HTP_CR y_HTP_s1 y_HTP_s2 y_HTP_CT]/1000
z_HTP = [z_HTP_CR_LE z_HTP_s1_LE z_HTP_s2_LE z_HTP_CT_LE]/1000
chord_HTP =[c_HTP_CR c_HTP_s1 c_HTP_s2 c_HTP_CT]/1000
offset_HTP = ([x_HTP_CR_LE x_HTP_s1_LE x_HTP_s2_LE x_HTP_CT_LE] - x_HTP_CR_LE)/1000
dihedral_HTP =[0 dihedral_s1 dihedral_s2 dihedral_CT]*180/pi
sweep_HTP =[0 sweep_s1_LE sweep_s2_LE sweep_CT_LE]*180/pi

%% Correction for the entry of data on FLOW5 according to the dihedral
b1 = sqrt((y_HTP_s1-y_HTP_CR)^2 + (z_HTP_s1_LE-z_HTP_CR_LE)^2)/1000;
b2 = sqrt((y_HTP_s2-y_HTP_s1)^2 + (z_HTP_s2_LE-z_HTP_s1_LE)^2)/1000;
bCT = sqrt((y_w1_CT-y_HTP_s2)^2 + (z_HTP_CT_LE-z_HTP_s2_LE)^2)/1000;

span_FLOW5_HTP = [ 0, b1, b2, bCT]

b_HTP(1) = 0;
N = length(span_FLOW5_HTP);
for i=1:N-1
    b_HTP(i+1) = b_HTP(i) + span_FLOW5_HTP(i+1);
end

% New wingspan modified as the wingspan along the surface
b_HTP

Geo_tier.cR_HTP = c_HTP_CR/1000;
Geo_tier.cB_k1_HTP = c_HTP_s1/1000;
Geo_tier.cB_k2_HTP = c_HTP_s2/1000;
Geo_tier.cT_HTP = c_HTP_CT/1000;

%% Elevator
% Elevator LE inboard
x_HTP_ele_y1_LE = 13163.774;
y_HTP_ele_y1_LE = y_HTP_s1;

% Elevator TE inboard
x_HTP_ele_y1_TE = x_HTP_s1_TE;
y_HTP_ele_y1_TE = y_HTP_s1;
z_HTP_ele_y1_TE = z_HTP_s1_TE;

% Chord of control surface
c_HTP_ele_y1 = x_HTP_ele_y1_TE - x_HTP_ele_y1_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_ele_y1_c_s1 = c_HTP_ele_y1/c_HTP_s1;
% Delta z location of wing
Delta_z_HTP_y1_LE = (z_HTP_s1_LE - z_HTP_s1_TE)*c_ele_y1_c_s1;
z_HTP_ele_y1_LE = z_HTP_s1_TE + Delta_z_HTP_y1_LE;

% Elevator LE outboard
x_HTP_ele_y2_LE = 13697.402;
y_HTP_ele_y2_LE = y_HTP_CT;

% Elevator TE outboard
x_HTP_ele_y2_TE = x_HTP_CT_TE;
y_HTP_ele_y2_TE = y_HTP_CT;
z_HTP_ele_y2_TE = z_HTP_CT_TE;

% Chord of control surface
c_HTP_ele_y2 = x_HTP_ele_y2_TE - x_HTP_ele_y2_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_ele_y2_c_CT = c_HTP_ele_y2/c_HTP_CT;
% Delta z location of wing
Delta_z_HTP_y2_LE = (z_HTP_CT_LE - z_HTP_CT_TE)*c_ele_y2_c_CT;
z_HTP_ele_y2_LE = z_HTP_CT_TE + Delta_z_HTP_y2_LE;

b_ele = y_HTP_ele_y2_TE - y_HTP_ele_y1_TE;

b_HTP_e = (y_HTP_CT_LE - y_HTP_CR);
Geo_tier.K_y1_ele_HTP = (y_HTP_ele_y1_TE - y_HTP_CR)/b_HTP_e;
Geo_tier.K_y2_ele_HTP = (y_HTP_ele_y2_TE - y_HTP_CR)/b_HTP_e;

cf_ele_tmp = [c_ele_y1_c_s1, c_ele_y2_c_CT];
Geo_tier.cf_ele = mean(cf_ele_tmp) ;


%% Elevator Trim Tab
% Elevator LE inboard
x_HTP_ele_TT_y1_LE = 13602.50;
y_HTP_ele_TT_y1_LE = y_HTP_s1;

% Elevator TE inboard
x_HTP_ele_TT_y1_TE = x_HTP_s1_TE;
y_HTP_ele_TT_y1_TE = y_HTP_s1;
z_HTP_ele_TT_y1_TE = z_HTP_s1_TE;

% Chord of control surface
c_HTP_ele_TT_y1 = x_HTP_ele_TT_y1_TE - x_HTP_ele_TT_y1_LE;
cf_TT_ele_y1 = c_HTP_ele_TT_y1/c_HTP_ele_y1;

% Ratio con constrol surface chord with respecto to surface chord
c_ele_TT_y1_c_s1 = c_HTP_ele_TT_y1/c_HTP_s1;
% Delta z location of wing
Delta_z_HTP_y1_LE = (z_HTP_s1_LE - z_HTP_s1_TE)*c_ele_TT_y1_c_s1;
z_HTP_ele_TT_y1_LE = z_HTP_s1_TE + Delta_z_HTP_y1_LE;

% Elevator LE outboard
x_HTP_ele_TT_y2_LE = 13769.2292;
y_HTP_ele_TT_y2_LE = y_HTP_s2;

% Elevator TE outboard
x_HTP_ele_TT_y2_TE = x_HTP_s2_TE;
y_HTP_ele_TT_y2_TE = y_HTP_s2;
z_HTP_ele_TT_y2_TE = z_HTP_s2_TE;

% Chord of control surface
c_HTP_ele_TT_y2 = x_HTP_ele_TT_y2_TE - x_HTP_ele_TT_y2_LE;
cf_TT_ele_y2 = c_HTP_ele_TT_y2/c_HTP_ele_y2;

% Ratio con constrol surface chord with respecto to surface chord
c_ele_TT_y2_c_s2 = c_HTP_ele_TT_y2/c_HTP_s2;
% Delta z location of wing
Delta_z_HTP_y2_LE = (z_HTP_s2_LE - z_HTP_s2_TE)*c_ele_TT_y2_c_s2;
z_HTP_ele_TT_y2_LE = z_HTP_s2_TE + Delta_z_HTP_y2_LE;

b_ele = y_HTP_ele_y2_TE - y_HTP_ele_y1_TE;

b_HTP_e = (y_HTP_CT_LE - y_HTP_CR);
Geo_tier.K_y1_TT_ele_HTP = (y_HTP_ele_TT_y1_TE - y_HTP_CR)/b_HTP_e;
Geo_tier.K_y2_TT_ele_HTP = (y_HTP_ele_TT_y2_TE - y_HTP_CR)/b_HTP_e;

cf_TT_ele_tmp = [cf_TT_ele_y1, cf_TT_ele_y2];
Geo_tier.cf_TT_ele = mean(cf_TT_ele_tmp) ;

%% VTP
% VTP Root points
%% HTP tip
% wing root - LE
x_fin_CR_LE = 8025.981;
y_fin_CR_LE = 0;
z_fin_CR_LE = 1337.503;

% wing root - TE
x_fin_CR_TE = 10793.235;
y_fin_CR_TE = 0;
z_fin_CR_TE = 1122.535;

% Fin LE TIP - only since triangle
% wing root - TE
x_fin_CR1_TE = 10951.322;
y_fin_CR1_TE = 0;
z_fin_CR1_TE = 1646.189;

%% VTP
%% Wing root
% wing root - LE
x_VTP_CR_LE = 10562.961;
y_VTP_CR_LE = 0;
z_VTP_CR_LE = 1145.417;

% wing root - TE
x_VTP_CR_TE = 12777.834;
y_VTP_CR_TE = 0;
z_VTP_CR_TE = 938.215;
y_VTP_CR = 0;

c_VTP_CR = x_VTP_CR_TE - x_VTP_CR_LE;
i_VTP_CR = atan((y_VTP_CR_LE-y_VTP_CR_TE)/(abs(x_VTP_CR_LE-x_VTP_CR_TE)));
z_VTP_CR = (z_VTP_CR_TE + z_VTP_CR_LE)/2;

%% VTP knk 1 outter Trim Tab
% wing root - LE
x_VTP_s1_LE = 11258.558;
y_VTP_s1_LE = 0;
z_VTP_s1_LE = 2012.937;

% wing root - TE
x_VTP_s1_TE = 13248.979;
y_VTP_s1_TE = 0;
z_VTP_s1_TE = 1983.97;
y_VTP_s1 = 0;

% Chord
c_VTP_s1 = x_VTP_s1_TE - x_VTP_s1_LE;
% Incidence angle
i_VTP_s1 = atan((y_VTP_s1_LE-y_VTP_s1_TE)/(abs(x_VTP_s1_LE-x_VTP_s1_TE)));
% Normalized wingspan
z_VTP_s1 = (z_VTP_s1_TE + z_VTP_s1_LE)/2;
% Dihedral beteen this section and previous
dihedral_s1 = atan((y_VTP_s1_LE-y_VTP_CR_LE)/(abs(y_VTP_s1-y_VTP_CR)));
% Sweep beteen this section and previous
sweep_s1_LE = atan((x_VTP_s1_LE-x_VTP_CR_LE)/(abs(z_VTP_s1-z_VTP_CR)));

%% VTP tip
% wing root - LE
x_VTP_CT_LE = 11940.869;
y_VTP_CT_LE = 2864.849;
z_VTP_CT_LE = 2925.23;

% wing root - TE
x_VTP_CT_TE = 13675.00;
y_VTP_CT_TE = 2864.849;
z_VTP_CT_TE = z_VTP_CT_LE;

y_VTP_CT = (y_VTP_CT_LE + y_VTP_CT_TE)/2;

% Chord
c_VTP_CT = x_VTP_CT_TE - x_VTP_CT_LE;
% Incidence angle
i_VTP_CT = atan((y_VTP_CT_LE-y_VTP_CT_TE)/(abs(x_VTP_CT_LE-x_VTP_CT_TE)));
% Normalized wingspan
z_VTP_CT = (z_VTP_CT_TE + z_VTP_CT_LE)/2;
% Dihedral beteen this section and previous
dihedral_CT = atan((y_VTP_CT_LE-y_VTP_s1_LE)/(abs(y_VTP_CT-y_VTP_s1)));
% Sweep beteen this section and previous
sweep_CT_LE = atan((x_VTP_CT_LE-x_VTP_s1_LE)/(abs(z_VTP_CT-z_VTP_s1)));

y_VTP = [y_VTP_CR y_VTP_s1 y_VTP_CT]/1000
x_VTP = [x_VTP_CR_LE x_VTP_s1_LE y_VTP_CT_LE]/1000
z_VTP = [z_VTP_CR_LE z_VTP_s1_LE z_VTP_CT_LE]/1000
z_VTP = ([z_VTP_CR_LE z_VTP_s1_LE z_VTP_CT_LE] - z_VTP_CR_LE) /1000
chord_VTP =[c_VTP_CR c_VTP_s1 c_VTP_CT]/1000
offset_VTP = ([x_VTP_CR_LE x_VTP_s1_LE x_VTP_CT_LE] - x_VTP_CR_LE)/1000
dihedral_VTP =[0 dihedral_s1 dihedral_CT]*180/pi
sweep_VTP =[0 sweep_s1_LE sweep_CT_LE]*180/pi

%% Correction for the entry of data on FLOW5 according to the dihedral
b1 = sqrt((y_VTP_s1-y_VTP_CR)^2 + (z_VTP_s1_LE-z_VTP_CR_LE)^2)/1000;
bCT = sqrt((y_w1_CT-y_VTP_s1)^2 + (z_VTP_CT_LE-z_VTP_s1_LE)^2)/1000;

span_FLOW5_VTP = [ 0, b1, bCT]

b_VTP(1) = 0;
N = length(span_FLOW5_VTP);
for i=1:N-1
    b_VTP(i+1) = b_VTP(i) + span_FLOW5_VTP(i+1);
end

% New wingspan modified as the wingspan along the surface
b_VTP

Geo_tier.cR_VTP = c_VTP_CR/1000;
Geo_tier.cB_k1_VTP = c_VTP_s1/1000;
Geo_tier.cT_VTP = c_VTP_CT/1000;

%% Rudder
% Rudder LE inboard
x_VTP_rud_y1_LE = 11761.199;
y_VTP_rud_y1_LE = y_VTP_CR_LE;
z_VTP_rud_y1_LE = 1026.972;

% Rudder TE inboard
x_VTP_rud_y1_TE = x_VTP_CR_TE;
y_VTP_rud_y1_TE = y_VTP_CR_TE;
z_VTP_rud_y1_TE = z_VTP_CR_TE;

% Chord of control surface
c_VTP_rud_y1 = x_VTP_rud_y1_TE - x_VTP_rud_y1_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_rud_y1_c_CR = c_VTP_rud_y1/c_VTP_CR;
% Delta z location of wing

% Rudder LE outboard
x_VTP_rud_y2_LE = 12935.788;
y_VTP_rud_y2_LE = y_VTP_CT_LE;
z_VTP_rud_y2_LE = 1026.972;

% Rudder TE inboard
x_VTP_rud_y2_TE = x_VTP_CT_TE;
y_VTP_rud_y2_TE = y_VTP_CT_LE;
z_VTP_rud_y2_TE = z_VTP_CT_TE;

% Chord of control surface
c_VTP_rud_y2 = x_VTP_rud_y2_TE - x_VTP_rud_y2_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_rud_y2_c_CT = c_VTP_rud_y2/c_VTP_CT;
% Delta z location of wing

b_rudder = y_VTP_rud_y2_TE - y_VTP_rud_y1_TE;

b_VTP_e = (y_VTP_CT_LE - y_VTP_CR);
Geo_tier.K_y1_rudder_VTP = (z_VTP_rud_y1_TE - z_VTP_CR)/b_VTP_e;
Geo_tier.K_y2_rudder_VTP = (z_VTP_rud_y2_TE - z_VTP_CR)/b_VTP_e;

cf_rudder_tmp = [c_rud_y1_c_CR, c_rud_y2_c_CT];
Geo_tier.cf_rudder = mean(cf_rudder_tmp) ;

%% Rudder Trim Tab
% Rudder LE inboard
x_VTP_rud_TT_y1_LE = 12511.19;
y_VTP_rud_TT_y1_LE = 0;
z_VTP_rud_TT_y1_LE = 1026.972;

% Rudder TE inboard
x_VTP_rud_TT_y1_TE = x_VTP_CR_TE;
y_VTP_rud_TT_y1_TE = 0;
z_VTP_rud_TT_y1_TE = z_VTP_CR_TE;

% Chord of control surface
c_VTP_rud_TT_y1 = x_VTP_rud_TT_y1_TE - x_VTP_rud_TT_y1_LE;
% Ratio con constrol surface chord with respecto to surface chord
c_rud_TT_y1_c_CR = c_VTP_rud_TT_y1/c_VTP_CR;
cf_TT_rudder_y1 = c_VTP_rud_TT_y1/c_VTP_rud_y1;

% Delta z location of wing

% Rudder LE outboard
x_VTP_rud_TT_y2_LE = 13123.508;
y_VTP_rud_TT_y2_LE = 0;
z_VTP_rud_TT_y2_LE = z_VTP_s1_LE;

% Rudder TE inboard
x_VTP_rud_TT_y2_TE = x_VTP_s1_TE;
y_VTP_rud_TT_y2_TE = 0;
z_VTP_rud_TT_y2_TE = z_VTP_s1_TE;

% Chord of control surface
c_VTP_rud_TT_y2 = x_VTP_rud_TT_y2_TE - x_VTP_rud_TT_y2_LE;
cf_TT_rudder_y2 = c_VTP_rud_TT_y2/c_VTP_rud_y2;

% Ratio con constrol surface chord with respecto to surface chord
c_rud_TT_y2_c_CR = c_VTP_rud_TT_y2/c_VTP_s1;
% Delta z location of wing

b_rudder_TT = y_VTP_rud_TT_y2_TE - y_VTP_rud_TT_y1_TE;

b_VTP_e = (y_VTP_CT_LE - y_VTP_CR);
Geo_tier.K_y1_TT_rudder_VTP = (z_VTP_rud_TT_y1_TE - z_VTP_CR)/b_VTP_e;
Geo_tier.K_y2_TT_rudder_VTP = (z_VTP_rud_TT_y2_TE - z_VTP_CR)/b_VTP_e;

cf_rudder_tmp = [cf_TT_rudder_y1, cf_TT_rudder_y2];
Geo_tier.cf_TT_rudder = mean(cf_rudder_tmp) ;

Geo_tier