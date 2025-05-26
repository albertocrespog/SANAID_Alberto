function Stab_Der = getCnp(AC_CONFIGURATION,modelo,alpha,Stab_Der,Stab_Der_parts,Trim_ITER)

AR_we      = modelo.ala.ARwe;
LAMc4_w    = modelo.ala.LAMc4;
CL_w1         = Trim_ITER.CL_w1;
CLa_we      = Stab_Der_parts.CLalpha_w1_e_pw;
b_w        = modelo.ala.b;
eOswald_w  = modelo.ala.oswald;
Xca_w      = modelo.ala.Xca;
cMAC_w     = modelo.ala.MAC;

Zca_v      = modelo.vertical.Zca;
Xca_v      = modelo.vertical.Xca;

Xcg        = modelo.general.Xcg;
Minf       = modelo.general.Minf;

Cy_beta_vert = Stab_Der_parts.Cy_beta_vert;
Cl_p_w       = Stab_Der.Clp_w;

beta = sqrt(1 - Minf^2);


%% Cn_p
% Ala
aw1         = CLa_we/(pi*AR_we*eOswald_w);
K           = (1-aw1)/(1-eOswald_w*aw1);    % Pamadi Ecuacion 4.549
xi         = (Xca_w - Xcg)/cMAC_w;

Cnp_CL_CL0_M0   = -(AR_we + 6*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)*(xi/AR_we + tan(LAMc4_w)/12))/(6*(AR_we + 4*cos(LAMc4_w)));
Cnp_CL_CL0      = ((AR_we + 4*cos(LAMc4_w))/((AR_we*beta + 4*cos(LAMc4_w))))*...
                ((AR_we*beta + 0.5*((AR_we*beta + cos(LAMc4_w)))*(tan(LAMc4_w)^2))/...
                (AR_we + 0.5*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)^2))*Cnp_CL_CL0_M0;

Cn_p_w          = Cl_p_w*tan(alpha)*(K-1) + K*Cnp_CL_CL0*CL_w1; %Pamadi 4.594

% Vertical
Cn_p_vert       = -2/b_w*(Zca_v*sin(alpha) + (Xca_v - Xcg)*cos(alpha))*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)/b_w*Cy_beta_vert;

if isnan(Cn_p_vert)
    Cn_p_vert = 0;
end

% Twin Vertical Tail configuration
if AC_CONFIGURATION.twin_VTP == 1
    Cn_p_vert = 2*Cn_p_vert;
end


% DERIVADA TOTAL
Cn_p            = Cn_p_vert + Cn_p_w;

Stab_Der.Cnp_w = Cn_p_w;
Stab_Der.Cnp_v = Cn_p_vert;
Stab_Der.Cnp = Cn_p;

end
