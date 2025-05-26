%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%CONDICIÓN DE VUELO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Data_Der = get_derivatives_Cefiro_1214(Data_Trim,h,V)

%los parámetros fundamentales para el vuelo son:
%altura
%velocidad de vuelo
%temperatura del día

alpha = Data_Trim.alpha_1;
delta_e = Data_Trim.delta_e_1;

v=V;                                 %velocidad de crucero
g=9.8065;
% Data around the equilibrium point
Data_ATM = get_Atmospheric_Cefiro(h);
rho = Data_ATM.rho;
a_speed = Data_ATM.a_speed;
M = v/a_speed;
Q = 0.5*rho*v^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GEOMETRÍA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%longitudes en metros
%pesos en kilogramos
%%%%%%%%%%%%%%%%%%%%%%%%%DATOS INICIALES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

anchot=0.6;               %ancho de la cola
b=2.8124;                 %b es el wing span
bfmax=0.26;               %máximo ancho del fuselaje 
bv=0.8;                   %envergadura cola vertical va multiplicada por 2. 0,4
c=0.39299;                %cuerda media del ala
ca=0.14;                  %cuerda aleron
ce=0.36675;               %cuerda media efectiva
ch=0.28;                  %cuerda del estabilizador horizontal
ctip=0.295;               %cuerda en la punta del ala
Df=0.2877;                %diámetro del fuselaje en la raíz del ala
flecha=0;                 %ángulo de flecha del ala
diedro=2;                          %diedro del ala
lf=1.52;                  %longitud del fuselaje
ln=1.05;                  %distancia del morro a la raiz del ala 
lt=1.51;                  %distancia desde el cg del avión al cgt hori
lv=1.556;                 %distancia del cg del avión al cav
ltubo=1.1;                %longitud del tailboom puramente circular
S=1.088;                  %superficie alar
Sbmax=0.066;              %area máxima del fuselaje seccion transversal
Sbs=0.362;                %área lateral del fuselaje
Se=0.988;                 %superficie expuesta
Sh=ch*anchot;             %superficie estabilizador horizontal
S0=0.04908;               %sección circular del fuselaje 0.04908
Sr=0.039;                 %superficie del rudder
Sv=0.0763;                %superficie de los dos timones verticales
Vb=0.073;                 %volumen del fuselaje 0.073
zw=0.090;                 %diferencia de cotas entre cg y  c/4 del ala
zv=0.255;                 %diferencia de cotas entre cg y ca del vertical
zt=0;                     %diferencia de cotas entre cg y cg helice

Data_Der.S = S;
Data_Der.c = c;
Data_Der.g = g;
Data_Der.b = b;
%________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL DOWN-WASH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%
env=2.8124;                %envergadura del avión
hb=0;                      %altura de la cola
A=8.2569;                  %aspect ratio
Ae=7.64;                   %aspect ratio expuesto
lambda=ctip/c;             %parámetro lambda
lb=ltubo+3*(c/4)+ch/4;     %longitud corregida coje 3/4 de la cuerda ala
kh=(1-hb/env)/(2*lb/env)^(1/3);
kl=(10-3*lambda)/7;
ka=1/A-1/(1+A^(1.7));
downwash=4.44*(kh*kl*ka)^1.19;
d=1-downwash;

% Vertical
av=3.0477;              %pendiente sustentación estab vertical
etav=0.9;                          %eficiencia de las superficies de control
eficiencia=etav*(0.724+(3.06/(1+cos(flecha)))*(Sv/S)+0.4*(zw/Df)+0.0009*A);

%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL CENTRO DE GRAVEDAD%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xcgwb=1.25;                %centro de gravedad del ala 0.4
Xcgwb1=Xcgwb/c;            %centro de gravedad del ala adimensional
Xcale=ce/4;                %distancia del borde ataque al cg del alae
Xcgle=Xcale-0.0815;        %distancia del borde ataque al cg del alae

Xcawb=Xcgwb-0.0815;        %centro aerodinámico del ala
Xcawb1=Xcawb/c;            %centro aeodinámico del ala adimensional  

%OJO ese 0,0815 es la distancia entre centro de
%gravedad ala y centro aerodinamico ala recalcular si varia la cuerda del
%ala
    
Xcgt=Xcgwb+0.0806+0.0651+ltubo; %0.091037
Xcgt1=Xcgt/c; %0.1293

%OJO el 0,1293 es la distancia entre el
%borde de ataque y el cg de la cola en el estabilizador horizontal
%recalcular si se varia la forma de la cola            
Xcat=Xcgwb+0.0651+ltubo+ch/4;   
%OJO  ese 0,091037 es la long. entre el 
%cg del ala y el enganche del tailboom,
%recalcularla si secambia la cuerda del ala
Xcat1=Xcat/c;                        

%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL CENTRO DE GRAVEDAD%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%pero con el caso real%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                      
Xtrendelantero=0.535;         %brazo del tren delantero
Xtrentrasero=1.32;            %brazo del tren trasero
Xcarga1=0.4323;               %posición carga de pago delantera

m_trendelantero=3.046;          %peso del tren delantero
m_trentrasero1=9.818;           %peso del tren trasero 1
m_trentrasero2=9.622;           %peso del tren trasero 2
m_carga1=0;                     %peso de la carga delantera 3.606                 

Xavion=(Xtrendelantero*m_trendelantero + Xtrentrasero*(m_trentrasero1 + m_trentrasero2))/...
    (m_trendelantero + m_trentrasero1 + m_trentrasero2);                

Xdeposito=1.43;               %posición del depósito 0.7346  Ó 0.58
XbateriaNiMh=0.95;            %posición bateria NiMh
XbateriaLiPo=0.95;            %posición bateria LiPo

m_avion = m_trendelantero + m_trentrasero1 + m_trentrasero2;

m_deposito=0.7;                  %peso del depósito 1.2  
m_bateriaNiMh=0;               %peso bateria NiMh 0.116
m_bateriaLiPo=0;               %peso bateria LiPo 0.116

m = m_avion + m_carga1 + m_deposito + m_bateriaNiMh + m_bateriaLiPo;
W1=m*g;

Xcg=(Xavion*m_avion + Xcarga1*m_carga1 + Xdeposito*m_deposito + XbateriaNiMh*m_bateriaNiMh +...
    XbateriaLiPo*m_bateriaLiPo)/m;
Xcg1=Xcg/c;
%Xcg1=3.2422

Xcgsincarga = (Xcg*m - Xcarga1*m_carga1)/(m - m_carga1);
m_cgsincarga = m - m_carga1;


m_concarga = m;

Data_Der.W1=W1;
Data_Der.m=m;

%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%MODELO DE LA POLAR VS. VELOCIDAD%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                           
Cd0=0.02866;              %Cd0 inicial
k1=0.0012;                %coeficiente lineal de la polar
k2=0.04279;               %coeficiente cuadrático de la polar
T=64*0.25;                %17.46

%T=(q*S)*(Cd0+k1*Cl+k2*Cl^2)empuje que varía con la velocidad de vuelo
Data_Der.Cd0 = Cd0;
Data_Der.k1 = k1;
Data_Der.k2 = k2;

%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CONFIGURACIÓN PUSHER. MODELO DE HELICÓPTERO%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi=0.4572;              %diámetro de la hélice en metros, son 18 in

%Se usa el modelo helicóptero movimiento axial ascendente

Sheli=(pi*(phi/2)^2);              %superficie de la hélice
Safe=ch*phi/sqrt(2);               %superficie de la cola afectada 
Snoafe=(anchot-phi/sqrt(2))*ch;    %superficie de la cola no afectada 
bafe=Safe/S;                       %adimensionalización superficie afect
bnoafe=Snoafe/S;                   %adimensionalización superficie no afect
vio=sqrt(T/(2*rho*Sheli));          %velocidad inicial inducida por hélice
vi=-0.5*v+sqrt(0.25*v^2+vio^2);    %velocidad real inducida por la hélice
Vc=v+vi;                           %velocidad que afecta a la cola
q=0.5*rho*v^2;                      %presión dinámica inicial
qhafe=0.5*rho*Vc^2;                 %presion dinamica afectada 
qhnoafe =0.9*q;                    %presion dinamica no afectada 
nhafe=qhafe/q;
nhnoafe=qhnoafe/q;
% Aproximations
% nhafe=0.95;
% nhnoafe=0.95;

%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%CÁLCULO DE COEFICIENTES LOS COEFICIENTES DE Data_DerLIDAD%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%COEFICIENTES Data_DerLIDAD ESTÁTICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                %Cl0wb
                                %%%%%%                             
Cl0wb=0.29927;                  %sust. ala a ángulo de ataque nulo
                                %se supone que sustenta cero el fuselaje
                                
                                %Clawb
                                %%%%%%
                                
K1=0.98;                                        %tabla del finnennes ratio                                                               
Kwb=0.1714*(bfmax/b)^2+0.8326*(bfmax/b)+0.9974;
Kbw=0.781*(bfmax/b)^2+1.1976*(bfmax/b)+0.0088;
ClaN=2*K1*Sbmax/S;
Clae=3.8658*(3.8804/4.8083);                   %varia área expuesta
KN=(ClaN/Clae)*(S/Se);
Clawb=(KN+Kwb+Kbw)*Clae*(Se/S);
Clabody=2*K1*S0/Vb^(2/3);
Clawbcomprobacion=Clae+Clabody;

                                %Clha
                                %%%%%
                                
Clha=3.6814+2.27075*(anchot-1);    %sustentación de la cola
%se ha realizado una interpolación para ver cuanto disminuye si varia el
%ancho


                                %Clelev
                                %%%%%%%
                                
Clelev=2.4614+1.312*(anchot-1);    %sustentación del elevador
%se ha realizado una interpolación para ver cuanto disminuye si varia el
%ancho 

                                %Clh0
                                %%%%%
                                
Clh0=0;                           %perfil simétrico a ángulo cero no sust.

                                %Cmawb
                                %%%%%%

%beta*Ae= 7.6194  por tanto cogemos el primer tipo de interpolacion pg22
kappa=0.265;                                   %según pg28
flecha=0;
Xaccrebw=1/4+((b-bfmax)/(2*ce))*kappa*tan(flecha);
Xaccrewb=0.38;                                 %segun tabla de pg27 
%int('Sbmax*(ln-x)','x','0','X0')=.21027072000000000000000000000000e-1
XaccreN=-1/(ce*Sbmax)*.21027072000000000000000000000000e-1;                                
Clabdew=Clae*Kbw*(Se/S);
Clawdeb=Clae*Kwb*(Se/S);
Xacwbcre=(XaccreN*ClaN+Xaccrewb*Clawdeb+Xaccrebw*Clabdew)/Clawb;                                
Xcawb=Xacwbcre*(ce/c);
Cmawb=(Xcg1-Xcawb)*Clawb;
Cmawb=0;
                                %Cm0wb
                                %%%%%%

%Cmof
Cmof=(K1/(36.6*S*c))*bfmax^2*lf;      %momento de cabeceo del fuselaje
Cmow=-0.046;                          %momento de cabeceo del ala -0.046
Cm0wb=Cmow+Cmof;                      %del conjunto ala-fuselaje
                                      
%%%%%%%%%%%%%%%%Cálculo del punto neutro y margen estático%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%
                                           
N01=(Clawb*Xcawb1+Clha*d*(nhafe*bafe+nhnoafe*bnoafe)*Xcat1)/(Clawb+Clha*d*(nhafe*bafe+nhnoafe*bnoafe));
N0=N01*c;
Margenestatico=N01-Xcg1;

Data_Der.Margenestatico=Margenestatico;
Data_Der.Xcg=Xcg;
Data_Der.Xcgsincarga=Xcgsincarga;
Data_Der.m_sincarga=m_cgsincarga;
Data_Der.m_concarga=m_concarga;
Data_Der.Xcg1=Xcg1;
Data_Der.N01=N01;

%Datos de incidencias para el trimado
iwgrados=2;                     %incidencia del ala en grados
itgrados=-2;                    %incidencia de la cola en grados
iw=iwgrados*(2*pi)/360;
it=itgrados*(2*pi)/360;
epsilon0=0;     

%%%%%%%%%%%%%Cálculo de los coeficientes para el trimado%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
CL0 = Cl0wb+Clha*(it-epsilon0)*(nhafe*bafe+nhnoafe*bnoafe)+Clawb*iw;
CLalpha = Clawb+Clha*d*(nhafe*bafe+nhnoafe*bnoafe);
CLdelta_e = Clelev*(nhafe*bafe+nhnoafe*bnoafe);

CM0=Cm0wb+(Xcg1-Xcawb1)*(Cl0wb+Clawb*iw)+(nhafe*bafe+nhnoafe*bnoafe)*(Xcg1-Xcat1)*(Clh0+Clha*(it-epsilon0));
CMalpha=CLalpha*(Xcg1-N01);
CMdelta_e=-Clelev*(nhafe*bafe+nhnoafe*bnoafe)*(Xcat1-Xcg1);

Data_Der.CL0=CL0;
Data_Der.CLalpha=CLalpha;
Data_Der.CLdelta_e=CLdelta_e;

Data_Der.CM0=CM0;
Data_Der.CMalpha=CMalpha;
Data_Der.CMdelta_e=CMdelta_e;

%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%MODELO DE LA POLAR VS. VELOCIDAD%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation of Oswald Efficiency
a1=0.0004;
a2=-0.008;
a3=0.0501;
a4=0.8642;
lambda1=A*lambda;
R=a1*lambda1^3+a2*lambda1^2+a3*lambda1+a4;                            
e=(1.1*CLalpha)/(R*CLalpha+(1-R)*pi*A);
k=1/(pi*A*e);
Data_Der.k = k;

%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL TRIMADO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%

X=[-CM0;W1/(q*S)-CL0];
XX=[CMalpha CMdelta_e;CLalpha CLdelta_e];
trim=inv(XX)*X;         %el primer valor es alpha y el segundo delta

alpha=trim(1,1);
delta=trim(2,1);
alphagrados=alpha*360/(2*pi);
deltagrados=delta*360/(2*pi);

Data_Der.alpha1 = alpha;
Data_Der.delta1 = delta;

CL = CL0 + CLalpha*alpha + CLdelta_e*delta;

Data_Der.CL=CL;
Data_Der.alpha=alpha;
Data_Der.delta=delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%COEFICIENTES Data_DerLIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
Vh=(Sh*lt)/(S*c);              %volumen ratio estabilizador horizontal
Cdalpha=0.13023;               %resistencia del ala a ángulo trim
teta0=0;

cmed=c/(2*v);                  %adimensionalización de la cuerda
m1=2*m/(rho*S*v);               %adimensionalización de la masa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%CONSTRUCCIÓN DE LA MATRIZ A y B%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inercias de catia
Ix=2.727;                       %inercia en eje x
Iz=9.95;                        %inercia eje z
Ixz=-0.002;                     %inercia eje xz
Iy=7.447;                      %inercia sin adimensionalizar Iy=7.447

%inercias adimensionales
Ix1=Ix/(0.5*rho*v^2*S*b); 
Iz1=Iz/(0.5*rho*v^2*S*b); 
Ixz1=Ixz/(0.5*rho*v^2*S*b); 
Iy1=Iy/(0.5*rho*v^2*S*c);       %adimensionalización de la inercia

%inercias de trabajo
Ix1prima=Ix1/(Ix1*Iz1-Ixz1^2);
Iz1prima=Iz1/(Ix1*Iz1-Ixz1^2);
Ixz1prima=Ixz1/(Ix1*Iz1-Ixz1^2);

Data_Der.Ix=Ix;
Data_Der.Iz=Iz;
Data_Der.Ixz=Ixz;
Data_Der.Iy=Iy;
 
Data_Der.Ix1=Ix1;
Data_Der.Iz1=Iz1;
Data_Der.Ixz1=Ixz1;
Data_Der.Iy1=Iy1;


% mu=2*W/(rho*S*v);               %adimensionalización de la masa
bmed=b/(2*v);                  %adimensionalización de la envergadura

                            %Cxalfa
                            %%%%%%%

% Cl=CLalpha*alpha;                            
% Cdalfa=2*k*Cl*CLalpha;
% Cxalfa=Cl-Cdalfa;
% Cl=CLalpha*alpha;                            
Cdalfa=2*k*CL*CLalpha;
Cxalfa=CL-Cdalfa;

                            %Czalfa
                            %%%%%%%
Cd=Cd0+k*CL^2;
Czalfa=-CLalpha-Cd;

                            %Cmalfa
                            %%%%%%%
                            
Cmalfa=CMalpha;

                            %Cxu
                            %%%%
                            
Cdu=0;                      %para vuelo subsóico bajo  
Cxu=-2*Cd-Cdu; 
%Cxu=0;

                            %Czu
                            %%%%

Clu=0;                      %para vuelo subsónico  
q0=0;                       %OJO que no lo se
Czu=-2*CL -Clu - 2*m1*q0;

                            %Cmu
                            %%%%

Cmu=0;                      %0 para vuelo subsónico

                            %Cxq
                            %%%%

Cxq=0;

                            %Czq
                            %%%%


Clab=2*K1*(Sbmax/Vb^(2/3));                            
Claprima=Clab*(Vb^(2/3)/Sbmax);                            
Clqb=2*Claprima*(1-Xcg/lf);                            
chi=(0.07054)/c;                   
Clqe=(0.5+2*chi)*Clae;                            
Clqwb=(Kwb+Kbw)*((Se*ce)/(S*c))*Clqe+Clqb*((Sbmax*lf)/(S*c));                            
Clqt=2*Clha*Vh*(nhafe*bafe+nhnoafe*bnoafe);                           
Czq=-(Clqwb+Clqt);


                            %Cmq
                            %%%%

%Cmqe
B=sqrt(1-M^2*(cos(flecha))^2);
c1=A^3*(tan(flecha))^2;
c2=3/B;
c3=A*B+6*cos(flecha);
c4=A+6*cos(flecha);
c5=A+2*cos(flecha);
Cmqe02=-0.7*CLalpha*cos(flecha)*((A*(0.5*chi+2*chi^2))/(c5)+c1/(24*c4)+1/8);
Cmqe=((c1/c3+c2)/(c1/c4+3))*Cmqe02;

%Cmqb
%int('dSbmax/dx*(xcg-x)','x','0','X0')
%int('0.066/0.304*(
%1.2959-x)','x','0','0.304')=.75497399999999999999999999999998e-1
Cmab=2*K1/Vb*.75497399999999999999999999999998e-1;
Cmaprimab=Cmab*(Vb/Sbmax);
Xm1=Xcg/lf;
Xc=1/Vb*Sbmax*lf^2/2;
Xc1=Xc/lf;
Vb1=Vb/(Sbmax*lf);
Cmqb=2*Cmaprimab*(((1-Xm1)^2-Vb1*(Xc1-Xm1))/(1-Xm1-Vb1));


Cmqwb=(Kwb+Kbw)*(Se/S)*(ce/c)^2*Cmqe+Cmqb*(Sbmax/S)*(lf/c)^2;
Cmqt=-2*Clha*Vh*(nhafe*Safe+nhnoafe*Snoafe)/(Safe+Snoafe)*(lt/c);

Cmq=(Cmqwb+Cmqt);

                            %Cxalfapunto
                            %%%%%%%%%%%%
                            
                            
Cxalfapunto=0;              %vale 0 para vuelo subsónico


                            %Czalfapunto
                            %%%%%%%%%%%%
                            
beta=sqrt(1-M^2);
tau=beta*Ae;
Clg=((-pi*Ae)/(2*beta^2))*(0.0013*tau^4-0.0122*tau^3+0.0317*tau^2+0.0186*tau-0.0004);
Clalfapuntoe=1.5*(Xcawb/ce)*Clae+3*Clg;
 
%Clalfapuntob

%Clab=2*K1*(Sbmax/Vb^(2/3));              %anteriormente definido              
%Claprima=Clab*(Vb^(2/3)/Sbmax);          %anteriormente definido                  
Clalfapuntob=2*Claprima*(Vb/(Sbmax*lf)) ;                           
                            
                            
                            
Clalfapuntowb=(Kwb+Kbw)*((Se*ce)/(S*c))*Clalfapuntoe+Clalfapuntob*((Sbmax*lf)/(S*c)) ;                           
Clalfapuntot=2*Clha*Vh*(nhafe*Safe+nhnoafe*Snoafe)/(Safe+Snoafe)*downwash;         
Clalfapunto=Clalfapuntot;%+Clalfapuntowb
Czalfapunto=-Clalfapunto; %OJO 
                            
                            %Cmalfapunto
                            %%%%%%%%%%%%
                            
%Cmalfapuntoe
Cm0g=((-pi*Ae)/(2*beta^2))*(0.0008*tau^4-0.0075*tau^3+0.0185*tau^2+0.0128*tau-0.0003);
Cmalfasegundae=-(81/32)*(Xcale/ce)^2*Clae+9/2*Cm0g;
Cmalfapuntoe=Cmalfasegundae+(Xcgle/c)*Clalfapuntoe;
                            
%Cmalfapuntob
Cmalfapuntob=2*Cmaprimab*((Xc1-Xm1)/(1-Xm1-Vb1))*(Vb/(Sbmax*lf));
                                                      
Cmalfapuntowb=(Kwb+Kbw)*((Se*ce^2)/(S*c^2))*Cmalfapuntoe+Cmalfapuntob*((Sbmax*(lf^2))/(S*(c^2)));                            
Cmalfapuntot=-2*Clha*Vh*((nhafe*Safe+nhnoafe*Snoafe)/(Safe+Snoafe))*downwash*(lt/c);                            
Cmalfapunto=Cmalfapuntot;%+Cmalfapuntowb

                            %Cxdeltae
                            %%%%%%%%%                            
Cxdeltae=0;

                            %Czdeltae
                            %%%%%%%%%
                            
% Correcció Trimado
% Czdeltae=-Clha*Sh/S;
Czdeltae=-CLdelta_e;
                            %Cmdeltae
                            %%%%%%%%%                           
% Corrección Trimado
% Cmdeltae=-Clha*Vh;
Cmdeltae=CMdelta_e;

                            %Cxteta
                            %%%%%%%
                            
Cxteta=-CL;

                            %Czteta
                            %%%%%%%
                            
Czteta=0;     %0.3562

                            %Cmteta
                            %%%%%%%
                            
Cmteta=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%CONSTRUCCIÓN DE LA MATRIZ A y B%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%se procede a colocar cada coeficiente por su sitio dentro de la matriz
%además de mu, Iy1 y cmed se necesitan otros coeficientes adim

chi1=Cxalfapunto*cmed/(m1-Czalfapunto*cmed);
chi2=Cmalfapunto*cmed/(m1-Czalfapunto*cmed);

%los coeficientes son por orden:
a11=(Cxu+chi1*Czu)/(m1);
a12=(Cxalfa+chi1*Czalfa)/(m1);
a13=(Cxq*cmed+chi1*(m1+Czq*cmed))/(m1);
a14=(Cxteta+chi1*Czteta)/(m1);
a21=(Czu)/(m1-Czalfapunto*cmed);
a22=(Czalfa)/(m1-Czalfapunto*cmed);
a23=(m1+Czq*cmed)/(m1-Czalfapunto*cmed);
a24=(Czteta)/(m1-Czalfapunto*cmed);
a31=(Cmu+chi2*Czu)/(Iy1);
a32=(Cmalfa+chi2*Czalfa)/(Iy1);
a33=(Cmq*cmed+chi2*(m1+Czq*cmed))/(Iy1);
a34=(chi2*Czteta)/(Iy1);
a41=0;
a42=0;
a43=1;
a44=0;

A_long=[a11 a12 a13 a14;a21 a22 a23 a24;a31 a32 a33 a34;a41 a42 a43 a44];
autovalores_long=eig(A_long);

Data_Der.autovalores_long = autovalores_long;

frecuenciacorto=sqrt(imag(autovalores_long(1,1))^2+real(autovalores_long(1,1))^2);
frecuenciafugoide=sqrt(imag(autovalores_long(3,1))^2+real(autovalores_long(3,1))^2);

amortiguamientocorto=-real(autovalores_long(1,1))/frecuenciacorto;
amortiguamientofugoide=-real(autovalores_long(3,1))/frecuenciafugoide;

periodocorto=2*pi/(frecuenciacorto*sqrt(1-amortiguamientocorto^2));
periodofugoide=2*pi/(frecuenciafugoide*sqrt(1-amortiguamientofugoide^2));

%para la matriz B los coeficientes son:
b1=(Cxdeltae+chi1*Czdeltae)/(m1);
b2=(Czdeltae)/(m1-cmed*Czalfapunto);
b3=(Cmdeltae+chi2*Czdeltae)/(Iy1);
b4=0;

B_long=[b1;b2;b3;b4];

Data_Der.Cxalfa=Cxalfa;
Data_Der.Czalfa=Czalfa;
Data_Der.Cmalfa=Cmalfa;
Data_Der.Cxu=Cxu;
Data_Der.Czu=Czu;
Data_Der.Cmu=Cmu;
Data_Der.Cxq=Cxq;
Data_Der.Czq=Czq;
Data_Der.Cmq=Cmq;
Data_Der.Cxalfapunto=Cxalfapunto;
Data_Der.Czalfapunto=Czalfapunto;                           
Data_Der.Cmalfapunto=Cmalfapunto;
Data_Der.Cxdeltae=Cxdeltae;
Data_Der.Czdeltae=Czdeltae;
Data_Der.Cmdeltae=Cmdeltae;
Data_Der.Cxteta=Cxteta;
Data_Der.Czteta=Czteta;
Data_Der.Cmteta=Cmteta;

Data_Der.A_long=A_long;
Data_Der.B_long=B_long;

Data_Der.frecuenciacorto=frecuenciacorto;
Data_Der.frecuenciafugoide=frecuenciafugoide;
Data_Der.amortiguamientocorto=amortiguamientocorto;
Data_Der.amortiguamientofugoide=amortiguamientofugoide;
Data_Der.periodocorto=periodocorto;
Data_Der.periodofugoide=periodofugoide;
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%COEFICIENTES Data_DerLIDAD ESTÁTICA Data_DerERAL %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                        %Cyb
                        %%%%

%Cybw
Cybw_grados=-0.0001*diedro;  
Cybw=Cybw_grados*(360/(2*pi));
%radianes

%Cybv
Cybv=-2*av*eficiencia*(Sv/S);       
%radianes

%Cybf
Ki=1.32;                             
Claf=2*K1*S0/Vb^(2/3);
Cybf=-Ki*Claf*Vb^(2/3)/S;
%radianes

Cyb=Cybw+Cybf+Cybv;                 
%radianes

                        %Clb
                        %%%%

%Clbv
Clbv=-2*av*eficiencia*(Sv/S)*((zv*cos(alpha)-lv*sin(alpha))/b); 
%radianes

%Clbwb
Kf=0.98;
Kma=1;
ClbClflecha=-0.000;
ClbClaspect=-0.0005;
Clbdiedro=-0.00024;  
Kmdiedro=1;
incremento_Clbdiedro=-0.0005*sqrt(A)*(Df/b)^2;
Clbwb_grados=CL*(ClbClflecha*Kma*Kf+ClbClaspect)+diedro*(Clbdiedro*Kmdiedro+incremento_Clbdiedro)+1.2*sqrt(A)/57.3*(zw/b)*(2*Df/b);
%grados

Clbwb=Clbwb_grados*(360/(2*pi));

Clb=Clbwb+Clbv;
%radianes

                        %Cnb
                        %%%%
                                         
%Cnbdiedro
diedro_radianes=diedro*((2*pi)/360);
Cnbdiedro=-0.075*diedro_radianes*CL;                           
%radianes  

%Cnbflecha
Cnbflecha=CL^2/(4*pi*Ae);                             
%radianes

%Cnbv
Cnbv=av*(lv/b)*(Sv/S)*eficiencia ;      
%radianes                                                                                                                             "3,0477 es la pendiente de sustentacion de la cola"
                                                                                                                                                                                                            
%Cnbw
Knn=0.0015;                       %de grafica de altura/bfmax
Kri=1.35;                         %sale del numero de reynolds

Cnbwb=-Knn*Kri*(Sbs/S)*(lf/b)*(360/(2*pi));            
%radianes                                                 

Cnbw=Cnbdiedro+Cnbflecha+Cnbwb;
Cnb=(Cnbw+Cnbv);                         
%radianes
                             

                                
                           %Cydeltaa
                           %%%%%%%%%
                           

Cydeltaa=0;                          
%radianes


                           %Cldeltaa
                           %%%%%%%%%
taudeltaa=0.55;
vardeltaatau=0.4;
Cldeltaa=taudeltaa*vardeltaatau;        
%radianes    

                           %Cndeltaa
                           %%%%%%%%%


K11=-0.13;
Cndeltaa=2*K11*CL*Cldeltaa; 
%radianes           


                            %Cydeltar
                            %%%%%%%%%

tau1=0.6;
Cydeltar=tau1*2*av*(Sv/S);          
%radianes


                            %Cldeltar
                            %%%%%%%%%

Cldeltar=Cydeltar*(zv/b);           
%radianes


                            %Cndeltar
                            %%%%%%%%%

Cndeltar=-Cydeltar*(lv/b)*etav;
%radianes

%%%%%%%%%%%%%%%%%%%%%CÁLCULO DEL TRIMADO Data_DerERAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         %%%%%%%%%%%%%%%%%%%%%


%beta_grados=12;
%beta=beta_grados*((2*pi)/360);
%X1=[-Cyb*beta;-Clb*beta;-Cnb*beta];
%XX1=[(W*9.81/(q*S)) Cydeltaa Cydeltar;0 Cldeltaa Cndeltar; 0 Cndeltaa Cndeltar];
%trim=inv(XX1)*X1;         %el primer valor es alpha y el segundo delta

%alabeo=trim(1,1);
%deltaa=trim(2,1);
%deltar=trim(3,1);
%alabeo_grados=alabeo*360/(2*pi);
%deltaa_grados=deltaa*360/(2*pi);
%deltar_grados=deltar*360/(2*pi);                                                                                                                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%COEFICIENTES Data_DerLIDAD DINAMICA Data_DerERAL%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                        %Cyp
                        %%%%
                        
Cyp=0;
%grados


                        %Clp
                        %%%%
                        
%Clpv
z=zv*cos(alpha)-lv*sin(alpha);
Clpv=-2*(z/b)*((z-zv)/b)*Cybv;     
%radianes

%Clpw
Cd=Cd0+k1*CL+k2*CL^2;
Clpw=-1/6*(CLalpha+Cd);      
%radianes

Clp=Clpv+Clpw;                      
%radianes


                        %Cnp
                        %%%%

%Cnpv
Cnpv=-(2/b)*(lv*cos(alpha)+zv*sin(alpha))*((z-zv)/b)*Cybv;
%radianes

%Cnpw
a1=0.0004;
a2=-0.008;
a3=0.0501;
a4=0.8642;
lambda1=A*lambda;
R=a1*lambda1^3+a2*lambda1^2+a3*lambda1+a4;                            
e=(1.1*CLalpha)/(R*CLalpha+(1-R)*pi*A);                            
k=1/(pi*A*e);                                                        
Cdalfa=2*k*CL*CLalpha;
Cnpw=-1/6*(CL-Cdalfa);       
%radianes

Cnp=Cnpw+Cnpv;                          
%radianes


                        %Cyr
                        %%%%
                        

Cyr=0;
%grados

                        %Clr
                        %%%%
                        
%Clrv
Clrv=-(2/b^2)*(zv*sin(alpha)+lv*cos(alpha))*z*Cybv;
%radianes

%Clrw
Clrw=CL/3;
%radianes

Clr=Clrv+Clrw;
%radianes

                        %Cnr
                        %%%%
                        
                        
%Cnrv
Cnrv=(2/b^2)*(zv*sin(alpha)+lv*cos(alpha))^2*Cybv;
%radianes

%Cnrw
Cnrw=-(1/3)*Cdalfa;
%radianes

Cnr=Cnrv+Cnrw;
%radianes

                      %Cybpunto
                      %%%%%%%%%
                      
Cybpunto=0;
%grados

                      %Clbpunto
                      %%%%%%%%%
                      
Clbpunto=0;
%grados

                      %Cnbpunto
                      %%%%%%%%%
                      
Cnbpunto=0;
%grados
             
                         
%¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿BMED????????????????????????????????????????????????

chi1=Iz1prima*Clbpunto+Ixz1prima*Cnbpunto;
chi2=Ix1prima*Cnbpunto+Ixz1prima*Clbpunto;

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

A_lat=[a11 a12 a13 a14 a15;a21 a22 a23 a24 a25;a31 a32 a33 a34 a35;a41 a42 a43 a44 a45;a51 a52 a53 a54 a55];
autovalores_lat=eig(A_lat);

Data_Der.autovalores_lat = autovalores_lat;

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

B_lat=[b11 b12;b21 b22;b31 b32;b41 b42;b51 b52];

frecuencia_balanceo=sqrt(imag(autovalores_lat(3,1))^2+real(autovalores_lat(3,1))^2);
amortiguamiento_balanceo=-real(autovalores_lat(3,1))/frecuencia_balanceo;

frecuencia_convergencia=sqrt(imag(autovalores_lat(5,1))^2+real(autovalores_lat(5,1))^2);
amortiguamiento_convergencia=-real(autovalores_lat(5,1))/frecuencia_convergencia;

frecuencia_espiral=sqrt(imag(autovalores_lat(2,1))^2+real(autovalores_lat(2,1))^2);
amortiguamiento_espiral=-real(autovalores_lat(2,1))/frecuencia_espiral;

%se procede a colocar cada coeficiente por su sitio dentro de la matriz
%además de mu, Iy1 y cmed se necesitan otros coeficientes adim

Data_Der.Cyb=Cyb;
Data_Der.Cnb=Cnb;
Data_Der.Clb=Clb;
Data_Der.Cydeltaa=Cydeltaa;
Data_Der.Cndeltaa=Cndeltaa;
Data_Der.Cldeltaa=Cldeltaa;
Data_Der.Cydeltar=Cydeltar;
Data_Der.Cndeltar=Cndeltar;
Data_Der.Cldeltar=Cldeltar;
Data_Der.Cyp=Cyp;
Data_Der.Clp=Clp;
Data_Der.Cnp=Cnp;
Data_Der.Cyr=Cyr;
Data_Der.Clr=Clr;
Data_Der.Cnr=Cnr;
Data_Der.Cybpunto=Cybpunto;
Data_Der.Clbpunto=Clbpunto;
Data_Der.Cnbpunto=Cnbpunto;

Data_Der.A_lat=A_lat;
Data_Der.B_lat=B_lat;

Data_Der.frecuencia_balanceo=frecuencia_balanceo;
Data_Der.amortiguamiento_balanceo=amortiguamiento_balanceo;
Data_Der.frecuencia_convergencia=frecuencia_convergencia;
Data_Der.amortiguamiento_convergencia=amortiguamiento_convergencia;
Data_Der.frecuencia_espiral=frecuencia_espiral;
Data_Der.amortiguamiento_espiral=amortiguamiento_espiral;