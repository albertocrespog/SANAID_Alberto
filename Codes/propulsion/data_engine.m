function Engine = data_engine(k)

load('Data.mat')
load('Data_Prop.mat')

%e= Engine number concerning to excel master table. From 1 to 5.
e=6;

%---------------------------ENGINE DATA----------------------------------%
W_eng=Data.engine(e).Weight;     % Engine Weight.
n_eng=2;                         % Engine number.
n_fan=2;
n_enga=2;
D_prop=0.6043;        % Propeller Diameter [m].
D_propin=D_prop*100/2.54; % Propeller Diameter [in].
A_prop=pi*(D_prop/2)^2;          % Swept area by propeller.
RPM_max = 145000/D_propin; % Max RPM by engine builder.

eta_p=0.7;                       % Propulsive efficiency.
eta_gear=0.96;                   % Gear box efficiency.
eta_m=0.88;                      % Engine efficiency (output/input).
eta_esc=0.98;                    % Speed controller efficiency.
eta_dist=0.96;                   % Shaft efficiency.

% Prop geometry
% Datos genéricos Hélice
b_p = (22*2.54/100);
c_p = 3/100;
S_p = b_p*c_p;
AR_p = (b_p^2)/S_p;

% Polyfit Coefficients for APC 22x10
% CT_Polyfit{i} = CT{i}(3) + CT{i}(2).*XData_CT{i} + CT{i}(1).*XData_CT{i}.^2;
CT2 = CT{k}(1);
CT1 = CT{k}(2);
CT0 = CT{k}(3);

% CP_Polyfit{i} = CP{i}(4) + CP{i}(3).*XData_CP{i} + CP{i}(2).*XData_CP{i}.^2 + CP{i}(1).*XData_CP{i}.^3;  
       
CP3 = CP{k}(1);
CP2 = CP{k}(2);
CP1 = CP{k}(3);
CP0 = CP{k}(4);

% CQ_Polyfit{i} = CQ{i}(4) + CQ{i}(3).*XData_CQ{i} + CQ{i}(2).*XData_CQ{i}.^2 + CQ{i}(1).*XData_CQ{i}.^3;
CQ3 = CQ{k}(1);
CQ2 = CQ{k}(2);
CQ1 = CQ{k}(3);
CQ0 = CQ{k}(4);

% ethamp_Polyfit{i} = ethamp{i}(4) + ethamp{i}(3).*XData_ethamp{i} + ethamp{i}(2).*XData_ethamp{i}.^2 + ethamp{i}(1).*XData_ethamp{i}.^3;
etha_mp3 = ethamp{k}(1);
etha_mp2 = ethamp{k}(2);
etha_mp1 = ethamp{k}(3);
etha_mp0 = ethamp{k}(4);

Engine.W_eng = W_eng;
Engine.n_eng = n_eng;
Engine.n_fan = n_fan;
Engine.n_enga = n_enga;
Engine.D_propin = D_propin;
Engine.D_prop = D_prop;
Engine.A_prop = A_prop;
Engine.RPM_max = RPM_max;

Engine.eta_p = eta_p;
Engine.eta_gear = eta_gear;
Engine.eta_m = eta_m;
Engine.eta_esc = eta_esc;
Engine.eta_dist = eta_dist;

Engine.b_p = b_p;
Engine.c_p = c_p;
Engine.S_p = S_p;
Engine.AR_p = AR_p;

Engine.CT0 = CT0;
Engine.CT1 = CT1;
Engine.CT2 = CT2;
Engine.CP0 = CP0;
Engine.CP1 = CP1;
Engine.CP2 = CP2;
Engine.CP3 = CP3;
Engine.CQ0 = CQ0;
Engine.CQ1 = CQ1;
Engine.CQ2 = CQ2;
Engine.CQ3 = CQ3;
Engine.etha_mp0 = etha_mp0;
Engine.etha_mp1 = etha_mp1;
Engine.etha_mp2 = etha_mp2;
Engine.etha_mp3 = etha_mp3;
Engine.RPM_max = RPM_max;

% Assume AXI 2826/12 GOLD LINE
Engine.dia = 35/1000;
Engine.length = 54/1000;