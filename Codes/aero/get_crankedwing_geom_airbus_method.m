function [S_e, S, S_s, S_fus, S_ph, S_pv, AR, AR_e, AR_s, XYZ_MAC, XYZ_MAC_e] = get_crankedwing_geom_airbus_method(cr, cB1, cB2, ct,...
    yB0, yB1, yB2, yt, Lambda_LE, Lambda_LE_1,Lambda_LE_2,dihedral, dihedral_1, dihedral_2) 

% Wing geometry and MAC Airbus Method
% Commercial Airplane Design Principles by Pasquale M. Sforza,
% Elsevier Science & Technology, 2014-04-05
% Wing area and mac calculation.
% Only valid for airbus method!!

%Inputs:

% cr ;  %root chord
% cB1; %kink1 chord; %If 0, no kinks
% cB2; %kink2 chord; %If 0, just one kink
% ct;   %tip chord;
% 
% yB0;  %y loc of wing root chord LE position
% yB1;   %y loc of wing kink1 chord LE position*
% yB2;   %y loc of wing kink2 chord LE position*
% yt;  %y loc of wing tip chord LE position
% 
% Lambda_LE; %sweep of wing (LE)
% Lambda_LE_1; %sweep of kink1 (LE)
% Lambda_LE_2; %sweep of kink2 (LE)
% 
% dihedral;    %dihedral of wing
% dihedral_1;    %dihedral of kink1
% dihedral_2;    %dihedral of kink2


cB0           = cr;   %root chord (Airbus Wing Definition)

if  cB1 == 0
    cB1 = [];
    yB1 = [];
    Lambda_LE_1 = [];
    dihedral_1= [];
    cB2 = [];
    yB2 = [];
    Lambda_LE_2 = [];
    dihedral_2 = [];
    
elseif cB2 == 0
    cB2 = [];
    yB2 = [];
    Lambda_LE_2 = [];
    dihedral_2 = [];
end
    

c_vec         = [cr, cB0, cB1, cB2, ct];
yB_vec        = [0, yB0, yB1, yB2, yt];
Lambda_LE_vec = [0,Lambda_LE, Lambda_LE_1, Lambda_LE_2];
dihedral_vec  = [0, dihedral, dihedral_1, dihedral_2];
ll            = length(yB_vec);


%Outputs:
b   = yB_vec(end)*2;   %span
b_e = (yB_vec(end)-yB_vec(2))*2; %exposed wing span

b_s_vec = ((yB_vec(2:end)-yB_vec(1:end-1))./cos(dihedral_vec(1:end)));
b_s = sum(b_s_vec)*2;

S_vec(1:ll-1) = 0.5*(yB_vec(2:ll)-yB_vec(1:ll-1)).*(c_vec(1:ll-1) + c_vec(2:ll));
S_vec_s(1:ll-1) = 0.5*(yB_vec(2:ll)-yB_vec(1:ll-1))./cos(dihedral_vec(1:ll-1)).*(c_vec(1:ll-1) + c_vec(2:ll));
S_vec_s(1:ll-1) = 0.5*(yB_vec(2:ll)-yB_vec(1:ll-1))./cos(dihedral_vec(1:ll-1)).*(c_vec(1:ll-1) + c_vec(2:ll));
S_ph(1:ll-1)  = S_vec_s.*cos(dihedral_vec(1:ll-1));
S_pv(1:ll-1)  = S_vec_s.*sin(dihedral_vec(1:ll-1));

Shalf = sum(S_vec);
Shalf_e = sum(S_vec(2:end));
Shalf_s = sum(S_vec_s(2:end));
Shalf_ph = sum(S_ph(2:end));
Shalf_pv = sum(S_pv(2:end));

% Solves error of negative area
S_vec(1:ll-1) = abs(0.5*(yB_vec(2:ll)-yB_vec(1:ll-1)).*(c_vec(1:ll-1) + c_vec(2:ll)));
S_vec_s(1:ll-1) = abs(0.5*(yB_vec(2:ll)-yB_vec(1:ll-1))./cos(dihedral_vec(1:ll-1)).*(c_vec(1:ll-1) + c_vec(2:ll)));
S_vec_s(1:ll-1) = abs(0.5*(yB_vec(2:ll)-yB_vec(1:ll-1))./cos(dihedral_vec(1:ll-1)).*(c_vec(1:ll-1) + c_vec(2:ll)));
S_ph(1:ll-1)  = abs(S_vec_s.*cos(dihedral_vec(1:ll-1)));
S_pv(1:ll-1)  = abs(S_vec_s.*sin(dihedral_vec(1:ll-1)));

Shalf = sum(S_vec);
Shalf_e = sum(S_vec(2:end));
Shalf_s = sum(S_vec_s(2:end));
Shalf_ph = sum(S_ph(2:end));
Shalf_pv = sum(S_pv(2:end));

S   = 2*Shalf;
S_e = 2*Shalf_e;
S_s = 2*Shalf_s;
S_fus = c_vec(1)*2*yB0;
S_ph  = 2*Shalf_ph;
S_pv  = 2*Shalf_pv;


AR   = b^2/S;
AR_e = b_e^2/S_e;
AR_s = b_s^2/S_s;

lambda_vec(1:ll-1) = c_vec(2:ll)./c_vec(1:ll-1);


cmac_vec(1:ll-1) = 2/3*c_vec(1:ll-1).*((1 + lambda_vec(1:ll-1) + lambda_vec(1:ll-1).^2)./(1 + lambda_vec(1:ll-1)));

cmac   = sum(cmac_vec.*S_vec/Shalf);
cmac_e = sum(cmac_vec(2:end).*S_vec(2:end)/Shalf_e);


ymac_vec(1:ll-1) = (yB_vec(2:ll) - yB_vec(1:ll-1))./3.*((1 + 2*lambda_vec(1:ll-1))./(1 + lambda_vec(1:ll-1)));

ymac   = sum((ymac_vec(1:ll-1) + yB_vec(1:ll-1)).*S_vec(1:ll-1)/Shalf);
ymac_e = sum((ymac_vec(2:ll-1) + yB_vec(2:ll-1)).*S_vec(2:ll-1)/Shalf_e);


xmac     = sum((yB_vec(2:ll-1).*tan(Lambda_LE_vec(1:ll-2)) + ymac_vec(2:ll-1).*tan(Lambda_LE_vec(2:ll-1))).*S_vec(2:ll-1)/Shalf);
xmac_e   = sum((yB_vec(2:ll-1).*tan(Lambda_LE_vec(1:ll-2)) + ymac_vec(2:ll-1).*tan(Lambda_LE_vec(2:ll-1))).*S_vec(2:ll-1)/Shalf_e);

zmac     = sum((yB_vec(2:ll-1).*tan(dihedral_vec(1:ll-2)) + ymac_vec(2:ll-1).*tan(dihedral_vec(2:ll-1))).*S_vec(2:ll-1)/Shalf);
zmac_e   = sum((yB_vec(2:ll-1).*tan(dihedral_vec(1:ll-2)) + ymac_vec(2:ll-1).*tan(dihedral_vec(2:ll-1))).*S_vec(2:ll-1)/Shalf_e);


x25mac   = xmac + 0.25*cmac; % assume is at 25% chord
x25mac_e = xmac_e + 0.25*cmac_e; % assume is at 25% chord

XYZ_MAC.xbar_w = x25mac;
XYZ_MAC_e.xbar_w = x25mac_e;
XYZ_MAC.ybar_w = ymac;
XYZ_MAC_e.ybar_w = ymac_e;
XYZ_MAC.zbar_w = zmac;
XYZ_MAC_e.zbar_w = zmac_e;
XYZ_MAC.cbar_w = cmac;
XYZ_MAC_e.cbar_w = cmac_e;
end