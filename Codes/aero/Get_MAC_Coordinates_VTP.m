function XYZ_MAC = Get_MAC_Coordinates_VTP(b_w,lambda,cR,dihedral,Lambda_c4,b_VTP)

% From ETKIN
% Y-locatrion of MAC
zbar_w = (b_w)*((1 + 2*lambda)/(3*(1+lambda)));
XYZ_MAC.zbar_w = zbar_w;

% X-locatrion of MAC
n = 1/4; % assume is at 25% chord
xbar_w = n*cR + zbar_w*tan(Lambda_c4);
XYZ_MAC.xbar_w = xbar_w;

% Z-locatrion of MAC
ybar_w   = zbar_w*sin(dihedral);
XYZ_MAC.ybar_w = ybar_w; 

cbar = cR*(2/3)*((1 + lambda + lambda^2)/(1+lambda));
XYZ_MAC.cbar = cbar;