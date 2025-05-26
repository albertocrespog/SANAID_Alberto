function Geo_tier = Geometric_Data_tier_UPDATE(Geo_tier,MAC)

% Geometric_Data
% Loads the mat file in order to actualize with the new variables
load Geo_tier.mat

% Correction of location of mean aerodynamic center (MAC) using XFLR5 
xbar_w1 = MAC.xbar_w1;
xbar_w2 = MAC.xbar_w2;
xbar_v = MAC.xbar_v;

xbar_w1_e = MAC.xbar_w1;
xbar_w2_e = MAC.xbar_w2;

% Distances
x_w1_LE = Geo_tier.x_w1_LE;
x_w1_xbar = x_w1_LE + xbar_w1;
x_w1_xbar_e = x_w1_LE + xbar_w1_e;

x_w2_LE = Geo_tier.x_w2_LE;
x_w2_xbar = x_w2_LE + xbar_w2;
x_w2_xbar_e = x_w2_LE + xbar_w2_e;

x_v_LE = Geo_tier.x_v_LE;
x_v_xbar = x_v_LE + xbar_v;

% Storing DATA
Geo_tier.x_w1_xbar = x_w1_xbar;
Geo_tier.x_w1_xbar_e = x_w1_xbar_e;
Geo_tier.x_w2_xbar = x_w2_xbar;
Geo_tier.x_w2_xbar_e = x_w2_xbar_e;
Geo_tier.x_v_xbar = x_v_xbar;

save Geo_tier.mat