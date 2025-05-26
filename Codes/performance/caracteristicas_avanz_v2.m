% function handles = caracteristicas_avanz_v1(varargin)
function handles = caracteristicas_avanz_v2(Aero_TH,Prop_data,Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Aero,Prop_sel)

%% assigns the variable input argument
% switch nargin
%     case 1 
%         Aero_TH = varargin{1};        
%     case 2
%         Prop_data = varargin{2};
%     case 3
%         Weight_tier = varargin{3};
%     case 4
%         conv_UNITS = varargin{4};
%     case 5
%         Geo_tier = varargin{5};
%     case 6
%         AC_CONFIGURATION = varargin{6};
%     case 7
%         Aero = varargin{7};
% end

%% PROPUL
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR
% 4: CONSUMO ESPECIFICO
% 5: AVION CIVIL =1/MILITAR = 2
% 6: EFICIENCIA DE LA HELICE (ETA_P)
% 7: DERIVACION(TURBOFANES)
propul(1) = 2; % Type of engine
propul(2) = AC_CONFIGURATION.n_eng; % Number of engines
propul(3) = 8500; % Thrust (lbf) or Power (shp) per engine
propul(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
propul(5) = 1; % Normativa
propul(6) = 0.82; % Prop efficiency
propul(7) = 1; % By-pass

% propul(1) = Prop_sel(1); % Type of engine
% propul(2) = Prop_sel(2); % Number of engines
% propul(3) = Prop_sel(3); % Thrust (lbf) or Power (shp) per engine
% propul(4) = Prop_sel(4); % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
% propul(5) = Prop_sel(5); % Normativa
% propul(6) = Prop_sel(6); % Prop efficiency
% propul(7) = Prop_sel(7); % By-pass

empuje = propul(3);
consumo = propul(4);

if propul(1) == 1,
    propul(3) = empuje * 4.448221615255; %[N]
    propul(4) = consumo * 2.832546065 * 10^-5; %[kg/(N*s)]
else
    propul(3) = empuje * 745.699872; %[W]
    propul(4) = consumo * 1.68965941 * 10^-7; %[kg/(W*s)]
end

handles.propul(1) = propul(1); % Type of engine
handles.propul(2) = propul(2); % Number of engines
handles.propul(3) = propul(3); % Thrust (lbf) or Power (shp) per engine
handles.propul(4) = propul(4); % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
handles.propul(5) = propul(5); % Normativa
handles.propul(6) = propul(6); % Prop efficiency
handles.propul(7) = propul(7); % By-pass

%% AERODINAMICA
% 1: SUPERFICIE
% 2: CD0
% 3: K1 = K
% 4: K2
% 5: CLMAX LIMPIO
aerodinamica(1) = Geo_tier.S_ref; % Superficie alar
aerodinamica(2) = Aero_TH.CD0; % Coeficiente de resistencia parasitaria CD0: CD = CD0 + K1*CL^2 - K2*CL
aerodinamica(3) = Aero_TH.CD2; % Coeficiente de resistencia inducida K1: CD = CD0 + K1*CL^2 - K2*CL
aerodinamica(4) = Aero_TH.CD1; % Coeficiente de resistencia inducida K2: CD = CD0 + K1*CL^2 - K2*CL
aerodinamica(5) = Aero.CL_max_w1_CR; % Coeficiente de sustentación máxima en limpio

handles.aerodinamica(1) = aerodinamica(1); % Superficie alar
handles.aerodinamica(2) = aerodinamica(2); % Coeficiente de resistencia parasitaria CD0: CD = CD0 + K1*CL^2 - K2*CL
handles.aerodinamica(3) = aerodinamica(3); % Coeficiente de resistencia inducida K1: CD = CD0 + K1*CL^2 - K2*CL
handles.aerodinamica(4) = aerodinamica(4); % Coeficiente de resistencia inducida K2: CD = CD0 + K1*CL^2 - K2*CL
handles.aerodinamica(5) = aerodinamica(5); % Coeficiente de sustentación máxima en limpio

%% AERODINAMICA_DESPEGUE
%       % 1: CD0 EN DESPEGUE (CONSIDERANDO TREN + FLAPS)
%       % 2: CL EN DESPEGUE
%       % 3: K CON EFECTO SUELO (O ALTURA DEL ALA SOBRE EL SUELO Y
%       ENVERGADURA DEL AVION)
%       % 4: CLMAX SUCIO DESPEGUE (FLAPS = 20º)
aero_despegue(1) = Aero_TH.CD0; % Coeficiente de resistencia parasitaria CD0: CD = CD0 + K1*CL^2 - K2*CL
aero_despegue(2) = Aero_TH.CD2; % Coeficiente de resistencia inducida K1: CD = CD0 + K1*CL^2 - K2*CL
aero_despegue(3) = Aero.CL_w1_CR_iw; % Coeficiente de sustentación para ángulo de ataque nulo
aero_despegue(4) = Aero.CL_max_w1_CR; % Coeficiente de sustentación máxima en sucio Take Off

handles.aero_despegue(1) = aero_despegue(1); % Coeficiente de resistencia parasitaria CD0: CD = CD0 + K1*CL^2 - K2*CL
handles.aero_despegue(2) = aero_despegue(2); % Coeficiente de resistencia inducida K1: CD = CD0 + K1*CL^2 - K2*CL
handles.aero_despegue(3) = aero_despegue(3); % Coeficiente de sustentación para ángulo de ataque nulo
handles.aero_despegue(4) = aero_despegue(4); % Coeficiente de sustentación máxima en sucio Take Off

%% AERODINAMICA_ATERRIZAJE
%       % 1: CD0 EN ATERRIZAJE (CONSIDERANDO TREN + FLAPS + SPOILERS)
%       % 2: CL EN ATERRIZAJE
%       % 3: K CON EFECTO SUELO (O ALTURA DEL ALA SOBRE EL SUELO Y
%       ENVERGADURA DEL AVION)
%       % 4: CLMAX SUCIO ATERRIZAJE (FLAPS = 30º
aero_aterrizaje(1) = Aero_TH.CD0; % Coeficiente de resistencia parasitaria CD0: CD = CD0 + K1*CL^2 - K2*CL
aero_aterrizaje(2) = Aero_TH.CD2; % Coeficiente de resistencia inducida K1: CD = CD0 + K1*CL^2 - K2*CL
aero_aterrizaje(3) = Aero.CL_w1_CR_iw; % Coeficiente de sustentación para ángulo de ataque nulo
aero_aterrizaje(4) = Aero.CL_max_w1_CR;; % Coeficiente de sustentación máxima en sucio Landing

handles.aero_aterrizaje(1) = aero_aterrizaje(1); % Coeficiente de resistencia parasitaria CD0: CD = CD0 + K1*CL^2 - K2*CL
handles.aero_aterrizaje(2) = aero_aterrizaje(2); % Coeficiente de resistencia inducida K1: CD = CD0 + K1*CL^2 - K2*CL
handles.aero_aterrizaje(3) = aero_aterrizaje(3); % Coeficiente de sustentación para ángulo de ataque nulo
handles.aero_aterrizaje(4) = aero_aterrizaje(4); % Coeficiente de sustentación máxima en sucio Landing

%% PESOS
% 1: PESO EN VACIO
% 2: CARGA DE PAGO INICIAL
% 3: PESO TRIPULACION
% 4: %FUEL RESTANTE AL ACABAR
pesos(1) = Weight_tier.m_empty; % Peso en vacío
pesos(2) = Weight_tier.m_payload; % Carga de pago
pesos(3) = 0; % Peso de tripulación
pesos(4) = 6; % % Fuel restante al final

handles.pesos(1) = pesos(1); % Peso en vacío
handles.pesos(2) = pesos(2); % Carga de pago
handles.pesos(3) = pesos(3); % Peso de tripulación
handles.pesos(4) = pesos(4); % % Fuel restante al final
