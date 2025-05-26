%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%CONDICIÓN DE VUELO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Data_Trim = Trim_Cefiro_1214(h,V)

%los parámetros fundamentales para el vuelo son:
%altura
%velocidad de vuelo
%temperatura del día

v=V;                                 %velocidad de crucero
g=9.8065;
% Data around the equilibrium point
Data_ATM = get_Atmospheric_Cefiro(h);
rho = Data_ATM.rho;
a_speed = Data_ATM.a_speed;
M = v/a_speed;
Q = 0.5*rho*v^2;

Data_Trim.rho = rho;
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

Xcgsincarga = (Xcg*m - Xcarga1*m_carga1)/(m - m_carga1);
m_cgsincarga = m - m_carga1;

m_concarga = m;

Data_Trim.W1=W1;
Data_Trim.m=m;

%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%MODELO DE LA POLAR VS. VELOCIDAD%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                           
Cd0=0.02866;              %Cd0 inicial
k1=0.0012;                %coeficiente lineal de la polar
k2=0.04279;               %coeficiente cuadrático de la polar
T=64*0.25;                %17.46

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
%%%%%%%CÁLCULO DE COEFICIENTES LOS COEFICIENTES DE Data_TrimLIDAD%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%COEFICIENTES Data_TrimLIDAD ESTÁTICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%
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

Data_Trim.Margenestatico=Margenestatico;
Data_Trim.Xcg=Xcg;
Data_Trim.Xcgsincarga=Xcgsincarga;
Data_Trim.m_sincarga=m_cgsincarga;
Data_Trim.m_concarga=m_concarga;
Data_Trim.Xcg1=Xcg1;
Data_Trim.N01=N01;

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

Data_Trim.CL0=CL0;
Data_Trim.CLalpha=CLalpha;
Data_Trim.CLdelta_e=CLdelta_e;

Data_Trim.CM0=CM0;
Data_Trim.CMalpha=CMalpha;
Data_Trim.CMdelta_e=CMdelta_e;

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
Data_Trim.k = k;

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

Data_Trim.alpha_1 = alpha;
Data_Trim.delta_e_1 = delta;

CL = CL0 + CLalpha*alpha + CLdelta_e*delta;
Data_Trim.CL=CL;