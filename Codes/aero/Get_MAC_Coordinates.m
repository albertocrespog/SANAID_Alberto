function XYZ_MAC = Get_MAC_Coordinates(b_w,lambda,cR,dihedral,Lambda_c4)

% From ETKIN
% Y-locatrion of MAC
ybar_w = (b_w/2)*((1 + 2*lambda)/(3*(1+lambda)));
XYZ_MAC.ybar_w = ybar_w;

% X-locatrion of MAC
n = 1/4; % assume is at 25% chord
xbar_w = n*cR + ybar_w*tan(Lambda_c4);
XYZ_MAC.xbar_w = xbar_w;

% Z-locatrion of MAC
zbar_w   = ybar_w*tan(dihedral);
XYZ_MAC.zbar_w = zbar_w; 

cbar = cR*(2/3)*((1 + lambda + lambda^2)/(1+lambda));
XYZ_MAC.cbar = cbar;