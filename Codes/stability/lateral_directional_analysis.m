function Stab_Dyn_LatDir = lateral_directional_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,Weight_tier,TRIM_RESULTS,...
    Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen,filenameS)

% Obtains only if variable study
if Variable_Study == 1
    V = conditions.V;
    rho = Performance.rho;
    q_inf = 0.5*rho*V^2;
    m_TOW = conditions.m_TOW;
else
    V = Performance.V;
    q_inf = Performance.q_inf;
    rho = Performance.rho;
    m_TOW = Weight_tier.m_TOW;
end

% theta1 = TRIM_RESULTS.trim_alpha;
theta1 = 0;

in2m = conv_UNITS.in2m;
g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

% Performance
% V = Performance.V;
% rho = Performance.rho;
% q_inf = Performance.q_inf;

% Geometric data
S_ref = Geo_tier.S_ref;
S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
cmac_w1 = Geo_tier.cmac_w1;
% m_TOW = Weight_tier.m_TOW;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%COEFICIENTES Stab_DerLIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ixx = Weight_tier.Ixx;
Iyy = Weight_tier.Iyy;
Izz = Weight_tier.Izz;
Ixz = Weight_tier.Ixz;

cmed = cmac_w1/(2*V);                  %adimensionalización de la cuerda
m1 = 2*m_TOW/(rho*S_w1*V);               %adimensionalización de la masa

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

% --------------------------------------
% Cyb = -0.6031;
% Clb = -0.0567;
% Cnb = 0.1282;

% Cyp = -0.1348;
% Clp = -0.4682;
% Cnp = -0.2786;
% 
% Cyr = 0.2317;
% Clr = 0.3218;
% Cnr = -0.1373;

% --------------------------------------
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

Cydeltarv = Stab_Der.Cydeltarv;
Cldeltarv = Stab_Der.Cldeltarv;
Cndeltarv = Stab_Der.Cndeltarv;

Cydeltarv2 = Stab_Der.Cydeltarv2;
Cldeltarv2 = Stab_Der.Cldeltarv2;
Cndeltarv2 = Stab_Der.Cndeltarv2;

% Aileron Trim Tab
% Stab_Der.Ch_deltaa_Tab = Ch_deltaa_Tab;
Cy_delta_a_Tab = Stab_Der.Cy_delta_a_Tab;
Cl_delta_a_Tab = Stab_Der.Cl_delta_a_Tab;
Cn_delta_a_Tab = Stab_Der.Cn_delta_a_Tab;
% Rudder Trim Tab
% Stab_Der.Ch_deltaa_Tab = Ch_deltaa_Tab;
Cy_delta_r_Tab = Stab_Der.Cy_delta_r_Tab;
Cl_delta_r_Tab = Stab_Der.Cl_delta_r_Tab;
Cn_delta_r_Tab = Stab_Der.Cn_delta_r_Tab;

%inercias de catia
%inercias adimensionales
Ix1=Ixx/(0.5*rho*V^2*S_w1*b_w1); 
Iz1=Izz/(0.5*rho*V^2*S_w1*b_w1); 
Ixz1=Ixz/(0.5*rho*V^2*S_w1*b_w1); 
Iy1=Iyy/(0.5*rho*V^2*S_w1*cmac_w1);       %adimensionalización de la inercia
%inercias de trabajo
Ix1prima=Ix1/(Ix1*Iz1-Ixz1^2);
Iz1prima=Iz1/(Ix1*Iz1-Ixz1^2);
Ixz1prima=Ixz1/(Ix1*Iz1-Ixz1^2);

switch StabilityModel
    
    case 1    
    span =b_w1;
    
    chi1=Iz1prima*Clbpunto+Ixz1prima*Cnbpunto;
    chi2=Ix1prima*Cnbpunto+Ixz1prima*Clbpunto;
    
    bmed = b_w1/(2*V);
    
    %los coeficientes son por orden:
    a11=Cyb/(m1-bmed*Cybpunto);
    Cyphi=0.2;      %0.003125                      %esto es una hipótesis
    a12=Cyphi/(m1-bmed*Cybpunto);
    a13=(Cyp*bmed)/(m1-bmed*Cybpunto);
    a14=0;
    a15=-(m1-bmed*Cyr)/(m1-bmed*Cybpunto);
    a21=0;
    a22=0;
    a23=1;
    a24=0;
    a25=0;
    a31=Clb*Iz1prima+Cnb*Ixz1prima+chi1*bmed*a11;
    a32=chi1*bmed*a12;
    a33=Clp*bmed*Iz1prima+Cnp*Ixz1prima*bmed+chi1*bmed*a13;
    a34=0;
    a35=Clr*bmed*Iz1prima+Cnr*Ixz1prima*bmed+chi1*bmed*a15;
    a41=0;
    a42=0;
    a43=0;
    a44=0;
    a45=1;
    a51=Ix1prima*Cnb+Ixz1prima*Clb+bmed*chi2*a11;
    a52=chi2*bmed*a22;
    a53=bmed*(Cnp*Ix1prima+Clp*Ixz1prima+chi2*a13);
    a54=0;
    a55=bmed*(Cnr*Ix1prima+Clr*Ixz1prima+chi2*a15);
% 
%     Ybeta = (q_inf*S_w1*Cyb)/m_TOW;
%     Yp = (q_inf*S_w1*span*Cyp)/(2*m_TOW*V);
%     Yr = (q_inf*S_w1*span*Cyr)/(2*m_TOW*V);
%     Yda = (q_inf*S_w1*Cydeltaa)/(m_TOW);
%     Ydr = (q_inf*S_w1*Cydeltar)/(m_TOW);
% 
%     Lbeta = (q_inf*S_w1*span*Clb)/Ixx;
%     Lp = (q_inf*S_w1*span*span*Clp)/(2*Ixx*V);
%     Lr = (q_inf*S_w1*span*span*Clr)/(2*Ixx*V);
%     Lda = (q_inf*S_w1*span*Cldeltaa)/(Ixx);
%     Ldr = (q_inf*S_w1*span*Cldeltar)/(Ixx);
% 
% %     CNTb =0; 
% %     CyTb =0;
% 
%     Nbeta = (q_inf*S_w1*span*Cnb)/Izz;
%     Ntbeta = (q_inf*S_w1*span*CNTb)/Izz;
%     Np = (q_inf*S_w1*span*span*Cnp)/(2*Izz*V);
%     Nr = (q_inf*S_w1*span*span*Cnr)/(2*Izz*V);
%     Nda = (q_inf*S_w1*span*Cndeltaa)/(Izz);
%     Ndr = (q_inf*S_w1*span*Cndeltar)/(Izz);
% 


    A1=Ixz/Ixx;
    B1=Ixz/Izz;

    Matrix_lat=[a11 a12 a13 a14 a15;...
        a21 a22 a23 a24 a25;...
        a31 a32 a33 a34 a35;...
        a41 a42 a43 a44 a45;...
        a51 a52 a53 a54 a55];
   
    %para la matriz B los coeficientes son:
    b11=Cydeltaa/(m1-bmed*Cybpunto);
    b12=Cydeltar/(m1-bmed*Cybpunto);
    b21=0;
    b22=0;
    b31=Cldeltaa*Iz1prima+Cndeltaa*Ixz1prima+chi1*bmed*b11;
    b32=Cldeltar*Iz1prima+Cndeltar*Ixz1prima+chi1*bmed*b12;
    b41=0;
    b42=0;
    b51=Cndeltaa*Ix1prima+Cldeltaa*Ixz1prima+chi2*bmed*b11;
    b52=Cndeltar*Ix1prima+Cldeltar*Ixz1prima+chi2*bmed*b12;
    
    Deflection_lat=[b11 b12;b21 b22;b31 b32;b41 b42;b51 b52];
    
    
    if Variable_Study == 1
        % no saving for variable study
    else
%         save Lateral.mat  Matrix_lat Deflection_lat m_TOW V
%         save('data/Lateral.mat', 'Matrix_lat','Deflection_lat','m_TOW','V')
        Saving_LateralDirectional_Dynamics(Matrix_lat,Deflection_lat,m_TOW,V,OUTPUT_read_XLSX);
    end
        
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Lateral Dimensional Stability Derivatives      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        span = b_w1;
        Ybeta = (q_inf*S_w1*Cyb)/m_TOW;
        Yp = (q_inf*S_w1*span*Cyp)/(2*m_TOW*V);
        Yr = (q_inf*S_w1*span*Cyr)/(2*m_TOW*V);
        
        
        Lbeta = (q_inf*S_w1*span*Clb)/Ixx;
        Lp = (q_inf*S_w1*span*span*Clp)/(2*Ixx*V);
        Lr = (q_inf*S_w1*span*span*Clr)/(2*Ixx*V);
        
        
        % CNTb = 0;
        Nbeta = (q_inf*S_w1*span*Cnb)/Izz;
        Ntbeta = (q_inf*S_w1*span*CNTb)/Izz;
        Np = (q_inf*S_w1*span*span*Cnp)/(2*Izz*V);
        Nr = (q_inf*S_w1*span*span*Cnr)/(2*Izz*V);
        
        Yda = (q_inf*S_w1*Cydeltaa)/(m_TOW);
        Lda = (q_inf*S_w1*span*Cldeltaa)/(Ixx);
        Nda = (q_inf*S_w1*span*Cndeltaa)/(Izz);

        Ydr = (q_inf*S_w1*Cydeltar)/(m_TOW);
        Ldr = (q_inf*S_w1*span*Cldeltar)/(Ixx);
        Ndr = (q_inf*S_w1*span*Cndeltar)/(Izz);

        Ydrv = (q_inf*S_w1*Cydeltarv)/(m_TOW);
        Ldrv = (q_inf*S_w1*span*Cldeltarv)/(Ixx);
        Ndrv = (q_inf*S_w1*span*Cndeltarv)/(Izz);

        Ydrv2 = (q_inf*S_w1*Cydeltarv2)/(m_TOW);
        Ldrv2 = (q_inf*S_w1*span*Cldeltarv2)/(Ixx);
        Ndrv2 = (q_inf*S_w1*span*Cndeltarv2)/(Izz);

        Ydr_aTab = (q_inf*S_w1*Cy_delta_a_Tab)/(m_TOW);
        Ldr_aTab = (q_inf*S_w1*span*Cl_delta_a_Tab)/(Ixx);
        Ndr_aTab = (q_inf*S_w1*span*Cn_delta_a_Tab)/(Izz);

        Ydr_rTab = (q_inf*S_w1*Cy_delta_r_Tab)/(m_TOW);
        Ldr_rTab = (q_inf*S_w1*span*Cl_delta_r_Tab)/(Ixx);
        Ndr_rTab = (q_inf*S_w1*span*Cn_delta_r_Tab)/(Izz);
      

        %     Yda = (q1*Sref*Cy_da)/(m);
        %     Ydr = (q1*Sref*Cy_dr)/(m);
        %     Lda = (q1*Sref*b_w*Cl_da)/(Ixx);
        %     Ldr = (q1*Sref*b_w*Cl_dr)/(Ixx);
        %     Nda = (q1*Sref*b_w*Cn_da)/(Izz);
        %     Ndr = (q1*Sref*b_w*Cn_dr)/(Izz);

        Stab_Dyn_LatDir.Ybeta = Ybeta;
        Stab_Dyn_LatDir.Yp = Yp;
        Stab_Dyn_LatDir.Yr = Yr;

        Stab_Dyn_LatDir.Lbeta = Lbeta;
        Stab_Dyn_LatDir.Ntbeta = Ntbeta;
        Stab_Dyn_LatDir.Np = Np;
        Stab_Dyn_LatDir.Nr = Nr;

        Stab_Dyn_LatDir.Nbeta = Nbeta;
        Stab_Dyn_LatDir.Lp = Lp;
        Stab_Dyn_LatDir.Lr = Lr;

        Stab_Dyn_LatDir.Yda = Yda;
        Stab_Dyn_LatDir.Lda = Lda;
        Stab_Dyn_LatDir.Nda = Nda;

        Stab_Dyn_LatDir.Ydr = Ydr;
        Stab_Dyn_LatDir.Ldr = Ldr;
        Stab_Dyn_LatDir.Ndr = Ndr;

        Stab_Dyn_LatDir.Ydrv = Ydrv;
        Stab_Dyn_LatDir.Ldrv = Ldrv;
        Stab_Dyn_LatDir.Ndrv = Ndrv;

        Stab_Dyn_LatDir.Ydrv2 = Ydrv2;
        Stab_Dyn_LatDir.Ldrv2 = Ldrv2;
        Stab_Dyn_LatDir.Ndrv2 = Ndrv2;

        Stab_Dyn_LatDir.Ydr_aTab = Ydr_aTab;
        Stab_Dyn_LatDir.Ldr_aTab = Ldr_aTab;
        Stab_Dyn_LatDir.Ndr_aTab = Ndr_aTab;

        Stab_Dyn_LatDir.Ydr_rTab = Ydr_rTab;
        Stab_Dyn_LatDir.Ldr_rTab = Ldr_rTab;
        Stab_Dyn_LatDir.Ndr_rTab = Ndr_rTab;
        
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
        
        eel = ((Lbeta + A1*(Nbeta+Ntbeta))/(1-A1*B1))/V;
        ffl = (Lp + A1*Np)/(1-A1*B1);
        ggl = (Lr + A1*Nr)/(1-A1*B1);
        hhl = 0;
        hhl2 = 0;
        
        eel_r = Lbeta;
        ffl_r = Lp;
        ggl_r = Lr;
        hhl_r = 0;
        hhl2_r = 0;
        
        iil = ((Nbeta + Ntbeta + B1*Lbeta)/(1-A1*B1))/V;
        jjl = (Np + B1*Lp)/(1-A1*B1);
        kkl = (Nr + B1*Lr)/(1-A1*B1);
        lll = 0;
        lll2 = 0;
        
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
        
        
        
        Ma_lat = [1 0 0 0 0; 0 1 -A1 0 0; 0 -B1 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
        Ra_lat = [aal_r/V bbl_r ccl_r ddl_r ddl2_r; eel_r/V ffl_r ggl_r hhl_r hhl2_r;...
            iil_r/V jjl_r kkl_r lll_r lll2_r; mml nnl ool ppl ppl2; fifth1 fifth2 fifth3 fifth4 fifth5];
        
        Matrix_lat = inv(Ma_lat)*Ra_lat;
        %     A_lat=[aal bbl ccl ddl ddl2; eel ffl ggl hhl hhl2; ...
        %         iil jjl kkl lll lll2; mml nnl ool ppl ppl2; fifth1 fifth2 fifth3 fifth4 fifth5]
        Matrix_lat       =   Ma_lat\Ra_lat;
        
        qql = Yda;
        rrl = Ydr;
        ssl = (Lda+A1*Nda)/(1-A1*B1);
        ttl = (Ldr+A1*Ndr)/(1-A1*B1);
        uul = (B1*Lda+Nda)/(1-A1*B1);
        vvl = (B1*Ldr+Ndr)/(1-A1*B1);
        wwl = 0;
        xxl = 0;
        yyl = 0;
        zzl = 0;
        
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
        %
        Rb_lat = [qql_r rrl_r; ssl_r ttl_r; uul_r vvl_r; wwl_r xxl_r; yyl_r zzl_r];
        Deflection_lat = inv(Ma_lat)*Rb_lat;
        
        %     B_lat=[qql rrl; ssl ttl; uul vvl; wwl xxl; yyl zzl]
        %     pause
            
        if Variable_Study == 1
            % no saving for variable study
        else
%             save Lateral.mat  Matrix_lat Deflection_lat m_TOW V
%             save('data/Lateral.mat', 'Matrix_lat','Deflection_lat','m_TOW','V')
            Saving_LateralDirectional_Dynamics(Matrix_lat,Deflection_lat,m_TOW,V,OUTPUT_read_XLSX,filenameS);
        end
end


Stab_Dyn_LatDir.Matrix_lat = Matrix_lat;
Stab_Dyn_LatDir.Deflection_lat = Deflection_lat;

lat_poles   = eig(Matrix_lat);

%checks which set of poles are the dutch roll, spiral or roll
pole1=lat_poles(1);
pole2=lat_poles(2);
pole3=lat_poles(3);
pole4=lat_poles(4);
pole5=lat_poles(5);

Stab_Dyn_LatDir.lat.pole1 = pole1;
Stab_Dyn_LatDir.lat.pole2 = pole2;
Stab_Dyn_LatDir.lat.pole3 = pole3;
Stab_Dyn_LatDir.lat.pole4 = pole4;
Stab_Dyn_LatDir.lat.pole5 = pole5;

real_poles = real(lat_poles);

pole_check = 0;

for ii=1:length(real_poles)
    for jj=1:length(real_poles)
        if real_poles(ii) == real_poles(jj)
            check_ij(ii,jj)=1;
        else
            check_ij(ii,jj)=0;
        end
    end
end

same_poles = sum(check_ij');
%modification SER
[GC,GR] = groupcounts(same_poles');
ii=0;
jj=0;
k=0;
DR_poles = [];
Rest_poles = [];
% for ii=1:length(lat_poles)
%     if same_poles(ii) == 2
%         jj=jj+1;
%         DR_poles(jj) = lat_poles(ii);
%     else
%         k=k+1;
%         Rest_poles(k)= lat_poles(ii);
%     end
% end

for ii=1:length(lat_poles)
    if GC(1) == 5
        k=k+1;
        Rest_poles(k)= lat_poles(ii);
    else
        if GR(2)== 2
            if same_poles(ii) == 2
                jj=jj+1;
                DR_poles(jj) = lat_poles(ii);
            else
                k=k+1;
                Rest_poles(k)= lat_poles(ii);
            end
        elseif GR(2)== 4
            if same_poles(ii) == 2
                jj=jj+1;
                DR_poles(jj) = lat_poles(ii);
            else
                k=k+1;
                Rest_poles(k)= lat_poles(ii);
            end
        end
    end
end

for ii=1:length(Rest_poles)
    Rest_poles_vec(ii) = real(Rest_poles(ii));
end

if ~isempty(DR_poles) && ~isempty(Rest_poles)
    % If 2 pair of poles encountered
    if GC(2)== 2
        [B,I] = sort(abs(Rest_poles_vec));
        
        Spiral_approx = (Lbeta*Nr-Nbeta*Lr)/(Lbeta + Nbeta*A1);
        Rolling_approx = Lp;
        
        Yaw_pole = Rest_poles(I(1));
        
        if abs(Spiral_approx) < abs(Rolling_approx)
            Spiral_pole = Rest_poles(I(2));
            Roll_pole = Rest_poles(I(3));
        else
            Spiral_pole = Rest_poles(I(3));
            Roll_pole = Rest_poles(I(2));
        end
        
        Stab_Dyn_LatDir.lat.DR_poles1 = DR_poles(1);
        Stab_Dyn_LatDir.lat.DR_poles2 = DR_poles(2);
        Stab_Dyn_LatDir.lat.Spiral_pole = Spiral_pole;
        Stab_Dyn_LatDir.lat.Roll_pole = Roll_pole;
        Stab_Dyn_LatDir.lat.Yaw_pole = Yaw_pole;
        
        % real and imaginary parts
        re_dr=real(DR_poles(1));
        im_dr=imag(DR_poles(1));
        
        % damping ratio
        ZETA_DR=abs(re_dr/sqrt(re_dr^2 + im_dr^2));
        % natural frequency
        WN_DR=sqrt(re_dr^2 + im_dr^2);
        
        Stab_Dyn_LatDir.lat.DR_damp     = ZETA_DR;
        Stab_Dyn_LatDir.lat.DR_wn       = WN_DR;
        
        % period
        T_DR=(2*pi)/(WN_DR*sqrt(1-ZETA_DR^2));
        %time to half or double
        ta_DR=0.693/abs(re_dr);
        
        % cycles to double or half
        cycles_DR = 0.110*(sqrt(1-ZETA_DR^2))/abs(ZETA_DR);
        % Logarithmic decrement
        decrement_DR = 0.693/cycles_DR;
        
        Stab_Dyn_LatDir.lat.DR_T    = T_DR;
        Stab_Dyn_LatDir.lat.DR_T2   = ta_DR;
        
        Stab_Dyn_LatDir.lat.cycles_DR    = cycles_DR;
        Stab_Dyn_LatDir.lat.decrement_DR   = decrement_DR;
        
        Stab_Dyn_LatDir.lat.ROL_T2  = 0.693/abs(Roll_pole);
        Stab_Dyn_LatDir.lat.ESP_T2  = 0.693/abs(Spiral_pole);
    elseif GC(2)== 4
        % logic to determine which of the poles is closest to the Dutch
        % Roll
        Spiral_approx = (Lbeta*Nr-Nbeta*Lr)/(Lbeta + Nbeta*A1);
        Rolling_approx = Lp;
        
        Yaw_pole = Rest_poles;
        Spiral_pole = [];
        Roll_pole = [];
        
        Stab_Dyn_LatDir.lat.DR_poles1 = DR_poles(1);
        Stab_Dyn_LatDir.lat.DR_poles2 = DR_poles(2);
        Stab_Dyn_LatDir.lat.DR_poles1 = DR_poles(3);
        Stab_Dyn_LatDir.lat.DR_poles2 = DR_poles(4);
        Stab_Dyn_LatDir.lat.Spiral_pole = Spiral_pole;
        Stab_Dyn_LatDir.lat.Roll_pole = Roll_pole;
        Stab_Dyn_LatDir.lat.Yaw_pole = Yaw_pole;
        
        % real and imaginary parts
        re_dr=real(DR_poles(1));
        im_dr=imag(DR_poles(1));
        % real and imaginary parts
        re_dr2=real(DR_poles(3));
        im_dr2=imag(DR_poles(3));
        
        WN_DR_approx = sqrt(Nbeta + (1/V)*(Ybeta*Nr - Nbeta*Yr));
        % ZETA_DR;
        ZETA_DR_approx = - (Nr + Ybeta/V)/(2*WN_DR_approx);
        re_dr_approx = -WN_DR_approx*ZETA_DR_approx;
        
        % damping ratio
        ZETA_DR=abs(re_dr/sqrt(re_dr^2 + im_dr^2));
        % natural frequency
        WN_DR=sqrt(re_dr^2 + im_dr^2);
        
        % damping ratio
        ZETA_DR2=abs(re_dr2/sqrt(re_dr2^2 + im_dr2^2));
        % natural frequency
        WN_DR2=sqrt(re_dr2^2 + im_dr2^2);
        
        Stab_Dyn_LatDir.lat.DR_damp     = ZETA_DR;
        Stab_Dyn_LatDir.lat.DR_wn       = WN_DR;
        
        Stab_Dyn_LatDir.lat.DR_damp2     = ZETA_DR2;
        Stab_Dyn_LatDir.lat.DR_wn2       = WN_DR2;
        
        % period
        T_DR=(2*pi)/(WN_DR*sqrt(1-ZETA_DR^2));
        %time to half or double
        ta_DR=0.693/abs(re_dr);

        % period
        T_DR2=(2*pi)/(WN_DR2*sqrt(1-ZETA_DR2^2));
        %time to half or double
        ta_DR2=0.693/abs(re_dr2);
        
        % cycles to double or half
        cycles_DR = 0.110*(sqrt(1-ZETA_DR^2))/abs(ZETA_DR);
        % Logarithmic decrement
        decrement_DR = 0.693/cycles_DR;
        % cycles to double or half
        cycles_DR2 = 0.110*(sqrt(1-ZETA_DR2^2))/abs(ZETA_DR2);
        % Logarithmic decrement
        decrement_DR2 = 0.693/cycles_DR2;
        
        Stab_Dyn_LatDir.lat.DR_T    = T_DR;
        Stab_Dyn_LatDir.lat.DR_T2   = ta_DR;
        
        Stab_Dyn_LatDir.lat.cycles_DR    = cycles_DR;
        Stab_Dyn_LatDir.lat.decrement_DR   = decrement_DR;
        
        Stab_Dyn_LatDir.lat.DR_T2    = T_DR2;
        Stab_Dyn_LatDir.lat.DR_T22   = ta_DR2;
        
        Stab_Dyn_LatDir.lat.cycles_DR2    = cycles_DR2;
        Stab_Dyn_LatDir.lat.decrement_DR2   = decrement_DR2;
        
        %         [B,I] = sort(abs(DR_poles));
        %         if I(1)
        %             a = [abs(DR_poles(1)), abs(DR_poles(3))];
        %             b = abs(re_dr_approx);
        %             % this will work even if 'a' is  in random order
        %             d=sort(abs(b-a));
        %             lowest=find(abs(b-a)==d(1));
        %             sec_lowest=find(abs(b-a)==d(2));
        %         end
        %         if lowest == 1
        %             DR_polesDR(1) = DR_poles()
        %         end
        
    end

else
    re_dr = [];
    im_dr = [];

    ZETA_DR = [];
    WN_DR = [];

    T_DR = [];
    ta_DR = [];
    
    cycles_DR = [];
    decrement_DR = [];

    Spiral_pole = [];
    Roll_pole = [];
    

    Stab_Dyn_LatDir.lat.DR_damp     = [];
    Stab_Dyn_LatDir.lat.DR_wn       = [];
    
    Stab_Dyn_LatDir.lat.DR_T        = [];
    Stab_Dyn_LatDir.lat.DR_T2       = [];
    
    Stab_Dyn_LatDir.lat.ROL_T2      = [];
    Stab_Dyn_LatDir.lat.ESP_T2      = [];
    
    Stab_Dyn_LatDir.lat.Spiral_pole = [];
    Stab_Dyn_LatDir.lat.Roll_pole   = [];
    Stab_Dyn_LatDir.lat.Yaw_pole    = [];
    
end

Poles_msg = strcat('The poles of the lateral-directional open-loop system are');
Poles_msg2 = strcat('-------------------------------------');

% Aproximations
% WN_DR;
WN_DR_approx = sqrt(Nbeta + (1/V)*(Ybeta*Nr - Nbeta*Yr));
% ZETA_DR;
ZETA_DR_approx = - (Nr + Ybeta/V)/(2*WN_DR_approx);
re_dr_approx = -WN_DR_approx*ZETA_DR_approx;
im_dr_approx = WN_DR_approx*sqrt(1 - ZETA_DR_approx^2);
Spiral_approx = (Lbeta*Nr-Nbeta*Lr)/(Lbeta + Nbeta*A1);
T_spiral_approx = 1/abs(Spiral_approx);
Rolling_approx = Lp;
T_rolling_approx = abs(1/Lp);

% Eigenvalues
DR_msg = strcat('The Dutch Roll (DR) Analysis');
DR_msg2 = strcat('------------------------------');

% if all poles are different and no Dutch Roles are encountered
if ~isempty(DR_poles) && ~isempty(Rest_poles)
    if GC(2)== 2
        poles_dr_msg = strcat('The DR Poles are: real part: ',num2str(re_dr),' & imaginary part: ',num2str(im_dr),'i');
    elseif GC(2)== 4
        poles_dr_msg1 = strcat('The DR Poles are: real part: ',num2str(re_dr),' & imaginary part: ',num2str(im_dr),'i');
        poles_dr_msg2 = strcat('The DR Poles are: real part: ',num2str(re_dr2),' & imaginary part: ',num2str(im_dr2),'i');
    end
else
    poles_dr_msg = strcat('There are No Dutch Roll Poles');
end

% if all poles are different and no Dutch Roles are encountered
if ~isempty(DR_poles) && ~isempty(Rest_poles)
    poles_dr_approx_msg = strcat('The DR Poles aproximation are: real part: ',num2str(re_dr_approx),' & imaginary part: ',num2str(im_dr_approx),'i');
    WN_dr_msg = strcat('The DR natural frequency is: ',num2str(WN_DR),' and the approximation: ',num2str(WN_DR_approx));
    DAMP_dr_msg = strcat('The DR damping is: ',num2str(ZETA_DR),' and the approximation: ',num2str(ZETA_DR_approx));
    Period_dr_msg = strcat('The DR period is: ',num2str(T_DR),' s');
    ta_dr_msg = strcat('The DR time to half/double is: ',num2str(ta_DR),' s');
    cycles_dr_msg = strcat('The DR Cycles to double or half is: ',num2str(cycles_DR));
    decrement_dr_msg = strcat('The DR Logarithmic decrement is: ',num2str(decrement_DR));
    if GC(2)== 2
        poles_rolling_msg = strcat('The Roll pole is: ',num2str(Roll_pole), 'and the Roll pole approximation is: ',num2str(Rolling_approx));
        Period_rolling_msg = strcat('The Rolling period is: ',num2str(T_rolling_approx),' s');
        poles_spiral_msg = strcat('The Spiral pole is: ',num2str(Spiral_pole), 'and the Spiral pole approximation is: ',num2str(Spiral_approx));
        Period_spiral_msg = strcat('The Spiral period is: ',num2str(T_spiral_approx),' s');
    elseif GC(2)== 4
        poles_rolling_msg = strcat('No Roll pole');
        Period_rolling_msg = strcat('No Rolling period');
        poles_spiral_msg = strcat('No Spiral pole');
        Period_spiral_msg = strcat('No Spiral period');
    end
else
   
    poles_dr_approx_msg = strcat('The DR Poles aproximation are: real part: ',num2str(re_dr_approx),' & imaginary part: ',num2str(im_dr_approx),'i');
    WN_dr_msg = strcat('The DR natural frequency approximation: ',num2str(WN_DR_approx));
    DAMP_dr_msg = strcat('The DR damping approximation: ',num2str(ZETA_DR_approx));
    
    poles_rolling_msg = strcat('There are no Roll poles for the system, and the Theoretical Roll pole approximation is: ',num2str(Rolling_approx));
    Period_rolling_msg = strcat('The Rolling period is: ',num2str(T_rolling_approx),' s');
    poles_spiral_msg = strcat('The Spiral pole approximation is: ',num2str(Spiral_approx));
    Period_spiral_msg = strcat('The Spiral period is: ',num2str(T_spiral_approx),' s');
    LAT_POLES_msg = strcat('The lateral-directional Poles are: ',num2str(T_spiral_approx),' s');
end

% Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
% disp(Warning)

% Side Speed Stability
crit2 = Cyb - CyTb;
if crit2 < 0
    Sideslip_Speed_Stability =1;
    Sideslip_Speed_Stability_msg = 'Satisfies with Side Speed Stability';
else
    Sideslip_Speed_Stability =0;
    Sideslip_Speed_Stability_msg = 'Does NOT Satisfy with Side Speed Stability';
end

% Angle of Sideslip Stability
crit5 = Cnb + CNTb;
if crit5 > 0
    Beta_Stability =1;
    Beta_Stability_msg = 'Satisfies with Angle of Sideslip Stability';
else
    Beta_Stability =0;
    Beta_Stability_msg = 'Does NOT Satisfy with Angle of Sideslip Stability';
end

% Roll Rate Stablity
crit6 = Clp;
if crit6 < 0
    P_Stability =1;
    P_Stability_msg = 'Satisfies with Roll Rate Stability';
else
    P_Stability =0;
    P_Stability_msg = 'Does NOT Satisfy with Roll Rate Stability';
end

% Yaw Rate Stablity
crit8 = Cnr;
if crit8 < 0
    R_Stability =1;
    R_Stability_msg = 'Satisfies with Yaw Rate Stability';
else
    R_Stability =0;
    R_Stability_msg = 'Does NOT Satisfy with Yaw Rate Stability';
end

% Sideslip on Rolling Moment Stablity
crit10 = Clb;
if crit10 < 0
    beta_L_Stability =1;
    beta_L_Stability_msg = 'Satisfies with Sideslip on Rolling Moment Stability';
else
    beta_L_Stability =0;
    beta_L_Stability_msg = 'Does NOT Satisfy with Sideslip on Rolling Moment Stability';
end

% Sideslip on Yawing Momwnt Stablity
crit11 = Cnb;
if crit11 > 0
    beta_N_Stability =1;
    beta_N_Stability_msg = 'Satisfies with Sideslip on Yawing Moment Stability';
else
    beta_N_Stability =0;
    beta_N_Stability_msg = 'Does NOT Satisfy with Sideslip on Yawing Moment Stability';
end

% Spiral Root Stablity
crit11 = (Lbeta*Nr-Nbeta*Lr);
if crit11 > 0
    spiral_Stability =1;
    spiral_Stability_msg = 'Satisfies with Spiral Root Stablity';
else
    spiral_Stability =0;
    spiral_Stability_msg = 'Does NOT Satisfy with Spiral Root Stablity';
end

% Stability Analysis
Lat_criteria_msg = strcat('Lateral-Directiona Criteria for Longitudinal Stability');
lat_msg2         = strcat('-----------------------------------------------------');

% Does not display message for variable study
if Variable_Study == 1
    % No message displayed
else
    if show_messages_screen == 1
        % Shos the poles
        disp(Poles_msg)
        disp(lat_poles)

        % Lateral Stability
        disp(DR_msg)
        disp(DR_msg2)
        disp(poles_dr_msg)
        disp(poles_dr_approx_msg)
        disp(WN_dr_msg)
        disp(DAMP_dr_msg)
        disp(Period_dr_msg)
        disp(ta_dr_msg)
        disp(cycles_dr_msg)
        disp(decrement_dr_msg)
        disp(poles_rolling_msg)
        disp(Period_rolling_msg)
        disp(poles_spiral_msg)
        disp(Period_spiral_msg)
        disp(DR_msg2)
        %     Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
        %     disp(Warning)

        % Stability Analysis
        disp(Lat_criteria_msg)
        disp(lat_msg2)
        disp(Sideslip_Speed_Stability_msg)
        disp(Beta_Stability_msg)
        disp(P_Stability_msg)
        disp(R_Stability_msg)
        disp(beta_L_Stability_msg)
        disp(beta_N_Stability_msg)
        disp(spiral_Stability_msg)

        %     Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
        %     disp(Warning)
    end
end
