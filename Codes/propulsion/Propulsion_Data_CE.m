%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propulsion Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Propulsion = Propulsion_Data_CE(AC_CONFIGURATION)

n_eng = AC_CONFIGURATION.n_eng;

eta_prop = 0.7;
eta_ex = 0.85; % 15% asociada a perdidas tubo de escape
eta_gen = 0.85; % 15% asociada a generador de energía
eta_mp = eta_prop*eta_ex*eta_gen;

Propulsion.eta_prop = eta_prop;
Propulsion.eta_ex = eta_ex;
Propulsion.eta_gen = eta_gen; 
Propulsion.eta_mp = eta_mp;

% propeller especific fuel consumptio (lb/hr/bhp)
cbhp_1 = 0.8;
cbhp_2 = 0.9;
cbhp_3 = 1.0;
cbhp_4 = 1.1;
cbhp_5 = 1.2;
cbhp_6 = 1.36;

Propulsion.cbhp_1 = cbhp_1;
Propulsion.cbhp_2 = cbhp_2;
Propulsion.cbhp_3 = cbhp_3;
Propulsion.cbhp_4 = cbhp_4;
Propulsion.cbhp_5 = cbhp_5;
Propulsion.cbhp_6 = cbhp_6;

% Estimación de posición de palanca
delta_TO = 1.00;
delta_R = 0.85;
delta_E = 0.75;

Propulsion.delta_TO = delta_TO;
Propulsion.delta_R = delta_R;
Propulsion.delta_E = delta_E;