%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conversion Units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function conv_UNITS = conversion_UNITS
% Radians 2 degrees
R2D = 180/pi;
% Degrees to radians
D2R = pi/180; 
% feet 2 m
ft2m = 0.3048;
% meter 2 feet
m2ft = 3.28084;
% Watt to hp
W2hp = 0.001341022;
% Watt to hp
hp2W = 1/W2hp;
% m/s 2 ft/s
mps2ftps = 3.28084;
% ft/s 3 m/s
ftps2mps = 1/3.28084;
% kg 2 lb
kg2lb = 2.20462;
% lb 2 kg
lb2kg= 1/2.20462;
% ft/m 2 m/s
ftpm2mps = 0.00508;
% m/s 2 ft/m
mps2ftpm = 1/0.00508;
% Gravity
g = 9.80665;
g_imp  = 32.174; %[ft/s^2]
% in 2 m
in2m = 0.0254;
% in 2 m
m2in = 1/0.0254;
% m/s 2 knot
mps2knot = 1.943844;
% knot 2 m/s
knot2mps = 1/1.943844;

% Wats to pft/sec
W2pftsec = 0.7375621;
% pft/sec 2 Wats
pftsec2W = 1/0.7375621;
% m^2 to ft^2
m22ft2 = 10.76391;
% ft^2 2 m^2
ft22m2 = 1/10.76391;
% density SI to density IMP
rho_SI2rho_IMP = (0.0023772/1.2250);
% density IMP to density SI
rho_IMP2rho_SI = 1/(0.0023772/1.2250);
% dynamic pressure in SI to IMP
qmet2qimp = m22ft2*rho_SI2rho_IMP;
% dynamic pressure in IMP to SI
qimp2qmet = 1/(m22ft2*rho_SI2rho_IMP);
% Newton to lb force
N2lbf = 0.2248089;
% lb force 2 Newton
lbf2N = 1/0.2248089;

% kg 2 Tm
kg2Tm = 1/1000;
% Tm 2 kg 
Tm2kg = 1000;
% m 2 km
m2km = 1/1000;
% km 2 m
km2m = 1000;
% seg 2 hours
seg2hrs = 1/3600;
% hours 2 seg
hrs2seg = 3600;
% km 2 nm
km2nm = 0.5399568;
% nm 2 km
nm2km = 1/0.5399568;
% nm 2 m
nm2m = 1852; %1/0.5399568;
% litter 2 gal
l2gal = 0.2641721 ;
% gal 2 litter
gal2l = 1/0.2641721 ;

% Slugs ft^2 to kg m^2 (Inertias)
slft2_2_kgm2 = 1.355817962; 
density_fuel = 6.70; %lb/gallon

conv_UNITS.g = g;
conv_UNITS.g_imp = g_imp;
% Rads & degrees
conv_UNITS.R2D = R2D;
conv_UNITS.D2R = D2R;

% Distance
conv_UNITS.ft2m = ft2m;
conv_UNITS.m2ft = m2ft;
conv_UNITS.in2m = in2m;
conv_UNITS.m2in = m2in;
conv_UNITS.km2nm = km2nm;
conv_UNITS.nm2km = nm2km;
conv_UNITS.nm2m = nm2m;
conv_UNITS.m2km = m2km;
conv_UNITS.km2m = km2m;

% Power
conv_UNITS.W2hp = W2hp;
conv_UNITS.hp2W = hp2W;
conv_UNITS.W2pftsec = W2pftsec;
conv_UNITS.pftsec2W = pftsec2W;

% Speed
conv_UNITS.mps2ftps = mps2ftps;
conv_UNITS.ftps2mps = ftps2mps;
conv_UNITS.mps2knot = mps2knot;
conv_UNITS.knot2mps = knot2mps;
conv_UNITS.ftpm2mps = ftpm2mps;
conv_UNITS.mps2ftpm = mps2ftpm;

% WEight
conv_UNITS.kg2lb = kg2lb;
conv_UNITS.lb2kg = lb2kg;
conv_UNITS.kg2Tm = kg2Tm;
conv_UNITS.Tm2kg = Tm2kg;

% Area
conv_UNITS.m22ft2 = m22ft2;
conv_UNITS.ft22m2 = ft22m2;
conv_UNITS.rho_SI2rho_IMP = rho_SI2rho_IMP;
conv_UNITS.rho_IMP2rho_SI = rho_IMP2rho_SI;

% Volumen
conv_UNITS.l2gal = l2gal;
conv_UNITS.gal2l = gal2l;

conv_UNITS.qmet2qimp = qmet2qimp;
conv_UNITS.qimp2qmet = qimp2qmet;
conv_UNITS.N2lbf = N2lbf;
conv_UNITS.lbf2N = lbf2N;
conv_UNITS.seg2hrs = seg2hrs;
conv_UNITS.hrs2seg = hrs2seg;

conv_UNITS.slft2_2_kgm2 = slft2_2_kgm2;

% Fuel - Jet-A 1
conv_UNITS.density_fuel = density_fuel;

% save('data/Data_conv_UNITS.mat', 'conv_UNITS')

