% Modelo propulsivo SP210
function [T,c,Vnmax,delta_act,n] = propulsive_model_SP210(V,h,delta_p,prop)
% V = velocidad, h = altura, delta_p = posicion de la palanca de gases, 
% prop = hélice (3612,2812,3412,3112)

% T = tracción, c = consumo

% El motor considerado es un SP210, con un modelo de potencia cuadrático,
% con parámetros estimados mediante ajuste polinómico de datos de operación
% real en banco.

% El modelo de consumo específico depende de la posición de palanca de gases
% y de la velocidad de giro de forma cuadrática. El modelo se ha ajustado a
% partir de datos de operación en banco.

% Determinación de la penalización por altitud: Modelo de Raymer.
delta_rho = (1-6.5e-3*h/288.15)^4.25193; %Troposfera ISA.
rho       = 1.225*delta_rho;
delta_h   = delta_rho-(1-delta_rho)/7.55;

% Hélices
% cT = cT0 + cT1*J + cT2*J^2;
% cP = cP0 + cP1*J + cP2*J^2;
if prop==3612 % Hélice BIPALA 36x12, modelada en PropCalc con p/D = 0.33, NACA 6412 con "angle
              % adjust" 1º, chord ratio 1.2, y ensayada a 6000 rpm
    cT0 = 0.07897;   cT1 = -0.1111;  cT2 = -0.03181;
    cP0 = 0.0189;   cP1 = 0.01604;   cP2 = -0.06985;

    D   = 36*0.0254; % Diámetro de la hélice en [m]
    
elseif prop==2812 % Hélice TRIPALA 28x12, modelada en PropCalc con p/D = 0.4286, NACA 6412 con "angle
                  % adjust" 1º, chord ratio 1.2, y ensayada a 6000 rpm
    cT0 = 0.1062;   cT1 = -0.1024;  cT2 = -0.07823;
    cP0 = 0.03036;   cP1 = 0.03298;   cP2 = -0.1007; 

    D   = 28*0.0254; % Diámetro de la hélice en [m]

elseif prop==3412 % Hélice BIPALA 34x12, modelada en PropCalc con p/D = 0.3529, NACA 6412 con "angle
                  % adjust" 1º, chord ratio 1.2, y ensayada a 6000 rpm
    cT0 = 0.07102;   cT1 = -0.09107;  cT2 = -0.03765;
    cP0 = 0.0161;   cP1 = 0.02106;   cP2 = -0.06772;

    D   = 34*0.0254; % Diámetro de la hélice en [m]

elseif prop==3112 % Hélice TRIPALA 31x12, modelada en PropCalc con p/D = 0.3871, NACA 6412 con "angle
                  % adjust" 1º, chord ratio 1.2, y ensayada a 6000 rpm
    cT0 = 0.0807;   cT1 = -0.08869;  cT2 = -0.05565;
    cP0 = 0.02014;   cP1 = 0.0256;   cP2 = -0.07632;

    D   = 31*0.0254; % Diámetro de la hélice en [m]
else
    warning(['La hélice no se ha introducido correctamente. En "prop", introducir 3612 para la hélice bipala 36x12, ...' ...
        '3412 para la bipala 34x12, 2812 para la tripala 28x12 o 3112 para la tripala 31x12.  '])
end

%Corrección de John T. Lowry (Slow Down Efficiency Factor) por interferencia de la hélice con el fuselaje:
Z    = 2*sqrt(0.122/pi)/D; %Diámetro fuselaje / diámetro hélice.
SDEF = 1.05263-0.04185*Z-0.01481*Z^2-0.62001*Z^3;

% DATOS DEL MOTOR SP210
n_Pzero=600/60; % rev/s / potencia entregada nula --> SE HA TOMADO LA MISMA QUE EL MOTOR DA170
n_ref=105; % vel. de referencia = 6300 rpm --> SE HA TOMADO LA MISMA QUE EL MOTOR DA170 (Es válida ya que simplemente adimensionaliza la ec.)
P_ref1=10000; % potencia de referencia 10kW --> SE HA TOMADO LA MISMA QUE EL MOTOR DA170 (Es válida ya que simplemente adimensionaliza la ec.)
deltae_vec=[10,15,20,30,40,50,55,60,65,70,75,80,85,90,100]'; % valores de la posición de la palanca de gases para los que se tienen datos
AB_de=[0.0000    0.4177;0.0000    0.5158;0.0000    0.5986;0.0000    0.8221;0.0000    0.9901;0.0000    1.1322;...
    0.0000    1.1693;0.0000    1.2157;0.0000    1.2561;0.0000    1.2739;0.0005    1.2869;0.0023    1.2893;...
    0.0012    1.2936;0.0043    1.3071;0.0212    1.3182];
% AB_de=[ajuste_AB(10,n_Pzero,n_ref,P_ref1);ajuste_AB(15,n_Pzero,n_ref,P_ref1);ajuste_AB(20,n_Pzero,n_ref,P_ref1);ajuste_AB(30,n_Pzero,n_ref,P_ref1);...
%     ajuste_AB(40,n_Pzero,n_ref,P_ref1);ajuste_AB(50,n_Pzero,n_ref,P_ref1);ajuste_AB(55,n_Pzero,n_ref,P_ref1);ajuste_AB(60,n_Pzero,n_ref,P_ref1);...
%     ajuste_AB(65,n_Pzero,n_ref,P_ref1);ajuste_AB(70,n_Pzero,n_ref,P_ref1);ajuste_AB(75,n_Pzero,n_ref,P_ref1);ajuste_AB(80,n_Pzero,n_ref,P_ref1);...
%     ajuste_AB(85,n_Pzero,n_ref,P_ref1);ajuste_AB(90,n_Pzero,n_ref,P_ref1);ajuste_AB(100,n_Pzero,n_ref,P_ref1)];
P_ref=P_ref1*delta_h;
ABC=[1.6061,   -2.1198,    1.2058]; % obtenido de function [ABC]=ajuste_cE_ABC(deltae,n_ref) para deltae=70 (es indiferente)
F_vec=[3.2677,2.4076,1.5729,1.1027,0.9712,0.9243,0.9289,0.9599,0.9915,1.0000,0.9964,1.0214,1.0364,1.0460,1.0425]';
% F_vec= [ajuste_cE_F(10,n_ref,ABC);ajuste_cE_F(15,n_ref,ABC);ajuste_cE_F(20,n_ref,ABC);ajuste_cE_F(30,n_ref,ABC);ajuste_cE_F(40,n_ref,ABC);...
%     ajuste_cE_F(50,n_ref,ABC);ajuste_cE_F(55,n_ref,ABC);ajuste_cE_F(60,n_ref,ABC);ajuste_cE_F(65,n_ref,ABC);ajuste_cE_F(70,n_ref,ABC);...
%     ajuste_cE_F(75,n_ref,ABC);ajuste_cE_F(80,n_ref,ABC);ajuste_cE_F(85,n_ref,ABC);ajuste_cE_F(90,n_ref,ABC);ajuste_cE_F(100,n_ref,ABC)];

AB_de_1 = interp1(deltae_vec,AB_de(:,1),delta_p);
AB_de_2 = interp1(deltae_vec,AB_de(:,2),delta_p);
F       = interp1(deltae_vec,F_vec,delta_p);

nmax   = 7000/60;    % pág. 58 System Description. Comprobar.
nmin   = 2500/60;    % pág. 42 System Description. Comprobar (¿ 2500, 2600 ?).
Pelec  = 320;        % Estimación de consumos (Teams 30/04/2020). SUPONGO QUE SE MANTIENE
Pcomax = Pelec/0.75; % Se incluye rendimiento de transformaciones. SE SUPONE RENDIMIENTO TOTAL DE GENERACIÓN = 0.75
%%% FIN CAMBIOS PARA EL MOTOR SP 210

%Determinación de Vnmax:
Pnmax_mot = P_ref*(AB_de_2-AB_de_1*nmax/n_ref)*(nmax/n_ref-n_Pzero/n_ref)-Pcomax;

if cP0 - cP1^2/(4*cP2)<Pnmax_mot/(rho*nmax^3*D^5)
    %A cualquier velocidad se tiene el motor a n=nmax.
    Vnmax = 0;
else
    %-cP2*J^2 - cP1*J + (Pnmax_mot/(rho*nmax^3*D^5)-cP0) = 0;
    JVnmax = (cP1+sqrt(cP1^2+4*cP2*(Pnmax_mot/(rho*nmax^3*D^5)-cP0)))/(-2*cP2);
    Vnmax  = JVnmax*nmax*D;
end

%Determinación del modo de operación
J_nmax    = V/(nmax*D);
cPnmax    = cP0 + cP1*J_nmax + cP2*J_nmax^2;
Pnmax_hel = cPnmax*rho*nmax^3*D^5;

if Pnmax_mot>Pnmax_hel
    %El equilibrio motor-hélice se daría con n>n_max si delta_p no se
    %redujese -->  Operación a n=nmax;
    n     = nmax;
    J     = J_nmax;
    cT    = cT0 + cT1*J + cT2*J^2;
    T     = SDEF*cT*rho*n^2*D^4;
    cP    = cP0 + cP1*J + cP2*J^2;
    Ptot  = rho*n^3*D^5*cP;
    
    %Determinación de delta_act:
    Pnmax_v   = P_ref*(AB_de(:,2)-AB_de(:,1)*nmax/n_ref)*(nmax/n_ref-n_Pzero/n_ref)-Pcomax;
    delta_act = interp1(Pnmax_v,deltae_vec,Ptot);
    F_act     = interp1(deltae_vec,F_vec,delta_act);
    c         = (Ptot+Pcomax)/1000*F_act*(ABC(1)*(n/n_ref).^2+ABC(2)*(n/n_ref)+ABC(3));
else
    %El equilibrio motor-hélice se da a n<n_max con delta_p
    % -->  Operación a delta_p;
    delta_act = delta_p;
    
    Jref = V/(n_ref*D);
    a_x  = cP0;
    b_x  = cP1*Jref + P_ref*AB_de_1/(rho*n_ref^3*D^5);
    c_x  = cP2*Jref^2 - P_ref*(AB_de_2+AB_de_1*n_Pzero/n_ref)/(rho*n_ref^3*D^5);
    d_x  = (P_ref*(AB_de_2*n_Pzero/n_ref)+Pcomax)/(rho*n_ref^3*D^5);
    %Ecuación de equilibrio: a_x*x^3 + b_x*x^2 + c_x*x + d_x = 0; Se hace
    %el cambio x = z - b_x/(3*a_x); 
    p_z = (3*a_x*c_x-b_x^2)/(3*a_x^2);
    q_z = (2*b_x^3-9*a_x*b_x*c_x+27*a_x^2*d_x)/(27*a_x^3);
    %Ahora la ecuación es: z^3 + p_z*z + q_z = 0;
    Delta_z = -4*p_z^3-27*q_z^2;
    if Delta_z>0
        %Hay 3 soluciones reales: Una negativa (imposible!) otra inestable
        %(n cercano a nPzero) y otra estable (la mayor).
        theta = acos(sqrt(27*q_z^2/(-4*p_z^3)));
        z_1   = -2*sign(q_z)*sqrt(-p_z/3)*cos(theta/3);
        z_2   = -2*sign(q_z)*sqrt(-p_z/3)*cos(theta/3+2/3*pi);
        z_3   = -2*sign(q_z)*sqrt(-p_z/3)*cos(theta/3+4/3*pi);
        z     = max([z_1 z_2 z_3]);
        x     = z - b_x/(3*a_x);   
    else
        %Sólo hay una solución real, que tiene sentido sólo si conlleva una
        %potencia positiva:
        if p_z<0
            z = -2*sign(q_z)*sqrt(-p_z/3)*cosh(1/3*acosh(-3*abs(q_z)/(2*p_z)*sqrt(-3/p_z)));
        else
            z = -2*sqrt(p_z/3)*sinh(1/3*asinh(3*q_z/(2*p_z)*sqrt(3/p_z)));
        end
        x   = z - b_x/(3*a_x);
    end
    n    = x*n_ref;
    if n<nmin
        %El motor se cala:
        T         = NaN;
        c         = NaN;
    else
        %La operación es posible:
        J     = V/(n*D);
        cT    = cT0 + cT1*J + cT2*J^2;
        T     = SDEF*cT*rho*n^2*D^4;
        cP    = cP0 + cP1*J + cP2*J^2;
        Ptot  = rho*n^3*D^5*cP;
        %P_mot = P_ref*(AB_de_2-AB_de_1*x)*(x-n_Pzero/n_ref)-Pcomax;
        c     = (Ptot+Pcomax)/1000*F*(ABC(1)*(n/n_ref).^2+ABC(2)*(n/n_ref)+ABC(3));
    end
end
end


