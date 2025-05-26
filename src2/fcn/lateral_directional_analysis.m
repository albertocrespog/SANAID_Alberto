function Stab_Dyn_LatDir = lateral_directional_analysis(V,rho,S_w,b_w,cmac_w,m_TO_1,Stab_Der,CG,...
    q_inf,theta1,conv_UNITS,Stability_Pamadi,CL)

g = conv_UNITS.g;

Ixx = CG.Ixx;
Iyy = CG.Iyy;
Izz = CG.Izz;
Ixz = CG.Ixz;

c_w = cmac_w;
cmed = c_w/(2*V);                  %adimensionalización de la cuerda
m1 = 2*m_TO_1/(rho*S_w*V);               %adimensionalización de la masa

% Store Stability Derivatives
Cyb = Stab_Der.Cyb;
Clb = Stab_Der.Clb;
Cnb = Stab_Der.Cnb;

Cyp = Stab_Der.Cyp;
Clp = Stab_Der.Clp;
Cnp = Stab_Der.Cnp;

Cyr = Stab_Der.Cyr;
Clr = Stab_Der.Clr;
Cnr = Stab_Der.Cnr;

Cybpunto = Stab_Der.Cybpunto;
Clbpunto = Stab_Der.Clbpunto;
Cnbpunto = Stab_Der.Cnbpunto;

CyTb = Stab_Der.CyTb;
CNTb = Stab_Der.CNTb;

Cydeltaa = Stab_Der.Cydeltaa;
Cldeltaa = Stab_Der.Cldeltaa;
Cndeltaa = Stab_Der.Cndeltaa;

Cydeltar = Stab_Der.Cydeltar;
Cldeltar = Stab_Der.Cldeltar;
Cndeltar = Stab_Der.Cndeltar;

%inercias de catia
%inercias adimensionales
Ix1=Ixx/(0.5*rho*V^2*S_w*b_w); 
Iz1=Izz/(0.5*rho*V^2*S_w*b_w); 
Ixz1=Ixz/(0.5*rho*V^2*S_w*b_w); 
Iy1=Iyy/(0.5*rho*V^2*S_w*c_w);       %adimensionalización de la inercia

%inercias de trabajo
Ix1prima=Ix1/(Ix1*Iz1-Ixz1^2);
Iz1prima=Iz1/(Ix1*Iz1-Ixz1^2);
Ixz1prima=Ixz1/(Ix1*Iz1-Ixz1^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Lateral Dimensional Stability Derivatives      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    span =b_w;
    Ybeta = (q_inf*S_w*Cyb)/m_TO_1;
    Yp = (q_inf*S_w*span*Cyp)/(2*m_TO_1*V);
    Yr = (q_inf*S_w*span*Cyr)/(2*m_TO_1*V);
    Yda = (q_inf*S_w*Cydeltaa)/(m_TO_1);
    Ydr = (q_inf*S_w*Cydeltar)/(m_TO_1);
    
    Lbeta = (q_inf*S_w*span*Clb)/Ixx;
    Lp = (q_inf*S_w*span*span*Clp)/(2*Ixx*V);
    Lr = (q_inf*S_w*span*span*Clr)/(2*Ixx*V);
    Lda = (q_inf*S_w*span*Cldeltaa)/(Ixx);
    Ldr = (q_inf*S_w*span*Cldeltar)/(Ixx);
    
    CNTb =0; 
    Nbeta = (q_inf*S_w*span*Cnb)/Izz;
    Ntbeta = (q_inf*S_w*span*CNTb)/Izz;
    Np = (q_inf*S_w*span*span*Cnp)/(2*Izz*V);
    Nr = (q_inf*S_w*span*span*Cnr)/(2*Izz*V);
    Nda = (q_inf*S_w*span*Cndeltaa)/(Izz);
    Ndr = (q_inf*S_w*span*Cndeltar)/(Izz);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % terms of the Matrix that finds the poles of the system response %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A1=Ixz/Ixx;
    B1=Ixz/Izz;
    
    aal = (Ybeta)/V;
    bbl = Yp;
    ccl = Yr-V;
    ddl = g*cos(theta1);
    ddl2= 0;
    
    aal_r = Ybeta;
    bbl_r = Yp;
    ccl_r = Yr-V;
    ddl_r = g*cos(theta1);
    ddl2_r= 0;
        
    eel_r = Lbeta;
    ffl_r = Lp;
    ggl_r = Lr;
    hhl_r = 0;
    hhl2_r = 0;

    iil_r = (Nbeta + Ntbeta);
    jjl_r = Np;
    kkl_r = Nr;
    lll_r = 0;
    lll2_r = 0;
    
    mml = 0;
    nnl = 1;
    ool = tan(theta1);
    ppl = 0;
    ppl2= 0;
    
    fifth1 = 0;
    fifth2 = 0;
    fifth3 = 1/cos(theta1);
    fifth4 = 0;
    fifth5 = 0;
    
    M_lat = [V 0 0 0 0; 0 1 -A1 0 0; 0 -B1 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
    Ra_lat = [aal_r bbl_r ccl_r ddl_r ddl2_r; eel_r ffl_r ggl_r hhl_r hhl2_r;...
        iil_r jjl_r kkl_r lll_r lll2_r; mml nnl ool ppl ppl2; fifth1 fifth2 fifth3 fifth4 fifth5];
    
    A_lat = inv(M_lat)*Ra_lat;

    qql_r = Yda;
    rrl_r = Ydr;
    ssl_r = Lda;
    ttl_r = Ldr;
    uul_r = Nda;
    vvl_r = Ndr;
    wwl_r = 0;
    xxl_r = 0;
    yyl_r = 0;
    zzl_r = 0;
    
    Rb_lat = [qql_r rrl_r; ssl_r ttl_r; uul_r vvl_r; wwl_r xxl_r; yyl_r zzl_r];
    B_lat = inv(M_lat)*Rb_lat;

% poles_lat = eig(A_lat_simp)
poles_lat = eig(A_lat);
k=length(poles_lat);

%checks which set of poles are the shorp period and which one the Phugoid
pole1=poles_lat(1);
pole2=poles_lat(2);
pole3=poles_lat(3);
pole4=poles_lat(4);
pole5=poles_lat(5);

Stab_Dyn_LatDir.pole1 = pole1;
Stab_Dyn_LatDir.pole2 = pole2;
Stab_Dyn_LatDir.pole3 = pole3;
Stab_Dyn_LatDir.pole4 = pole4;
Stab_Dyn_LatDir.pole5 = pole5;

for i=1:length(poles_lat)
    real_poles(i) = real(poles_lat(i));
end

pole_check = 0;
for i=1:length(real_poles)
    for j=1:length(real_poles)
        if real_poles(i) == real_poles(j)
            check_ij(i,j)=1;
        else
            check_ij(i,j)=0;
        end
    end
end

for i=1:length(real_poles)
    same_poles(i) = sum(check_ij(i,:));
end

i=0;
j=0;
k=0;
for i=1:length(poles_lat)
    if same_poles(i) == 2
        j=j+1;
        DR_poles(j) = poles_lat(i);
    else 
        k=k+1;
        Rest_poles(k)= poles_lat(i);
    end
end

for i=1:length(Rest_poles)
    Rest_poles_vec(i) = real(Rest_poles(i));
end

[B,I] = sort(Rest_poles_vec);

Yaw_pole = Rest_poles(I(3));
Spiral_pole = Rest_poles(I(2));
Roll_pole = Rest_poles(I(1));

Stab_Dyn_LatDir.DR_poles = DR_poles(1);
Stab_Dyn_LatDir.Spiral_pole = Spiral_pole;
Stab_Dyn_LatDir.Roll_pole = Roll_pole;
Stab_Dyn_LatDir.Yaw_pole = Yaw_pole;

% real and imaginary parts
re_dr=real(DR_poles(1));
im_dr=imag(DR_poles(1));

% damping ratio
ZETA_DR=abs(re_dr/sqrt(re_dr^2 + im_dr^2));
% natural frequency
WN_DR=sqrt(re_dr^2 + im_dr^2);

Stab_Dyn_LatDir.ZETA_DR = ZETA_DR;
Stab_Dyn_LatDir.WN_DR = WN_DR;

% period
T_DR=(2*pi)/(WN_DR*sqrt(1-ZETA_DR^2));
%time to half or double
ta_DR=0.693/abs(re_dr);

Stab_Dyn_LatDir.T_DR = T_DR;
Stab_Dyn_LatDir.ta_DR = ta_DR;

if Roll_pole < 0
    ta_half_roll = 0.693/abs(Roll_pole);
    Stab_Dyn_LatDir.ta_half_roll = ta_half_roll;
else
    ta_double_roll = 0.693/abs(Roll_pole);
    Stab_Dyn_LatDir.ta_double_roll = ta_double_roll;
end

if Spiral_pole < 0
    ta_half_spiral = 0.693/abs(Spiral_pole);
    Stab_Dyn_LatDir.ta_half_spiral = ta_half_spiral;
else
    ta_double_spiral = 0.693/abs(Spiral_pole);
    Stab_Dyn_LatDir.ta_double_spiral = ta_double_spiral;
end

if Stability_Pamadi == 0
    % Aproximations
    WN_DR_approx = sqrt(Nbeta + (1/V)*(Ybeta*Nr - Nbeta*Yr));
    ZETA_DR_approx = - (Nr + Ybeta/V)/(2*WN_DR);
    Spiral_approx = (Lbeta*Nr-Nbeta*Lr)/(Lbeta + Nbeta*A1);
end

% Side Speed Stability
crit2 = Cyb - CyTb;
if crit2 < 0
    Sideslip_Speed_Stability =1;
else
    Sideslip_Speed_Stability =0;
end

% Angle of Sideslip Stability 
crit5 = Cnb + CNTb;
if crit5 > 0
    Beta_Stability =1;
else
    Beta_Stability =0;
end

% Roll Rate Stablity
crit6 = Clp;
if crit6 < 0
    P_Stability =1;
else
    P_Stability =0;
end

% Yaw Rate Stablity
crit8 = Cnr;
if crit8 < 0
    R_Stability =1;
else
    R_Stability =0;
end

% Sideslip on Rolling Moment Stablity
crit10 = Clb;
if crit10 < 0
    beta_L_Stability =1;
else
    beta_L_Stability =0;
end

% Sideslip on Yawing Momwnt Stablity
crit11 = Cnb;
if crit11 > 0
    beta_N_Stability =1;
else
    beta_N_Stability =0;
end

% Spiral Root Stablity
crit11 = (Lbeta*Nr-Nbeta*Lr);
if crit11 > 0
    spiral_Stability =1;
else
    spiral_Stability =0;
end

Stab_Dyn_LatDir.Sideslip_Speed_Stability = Sideslip_Speed_Stability;
Stab_Dyn_LatDir.Beta_Stability = Beta_Stability;
Stab_Dyn_LatDir.P_Stability = P_Stability;
Stab_Dyn_LatDir.R_Stability = R_Stability;
Stab_Dyn_LatDir.beta_L_Stability = beta_L_Stability;
Stab_Dyn_LatDir.beta_N_Stability = beta_N_Stability;
Stab_Dyn_LatDir.spiral_Stability = spiral_Stability;
