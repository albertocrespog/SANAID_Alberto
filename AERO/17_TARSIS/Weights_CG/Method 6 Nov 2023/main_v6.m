clear; close all; clc

caso = 2;    % [0, 1, 2]. 0 = caso control; 1 = No BAT; 2 = Sí BAT
W_PL = 5.5;  % Payload (kg)

% Lectura datos T120
T120 = xlsread('Tabla cdg e inercias T120 KSA 2.xlsx', 'D4:M13');

%Lectura datos CG depósito
FUEL = xlsread('Tabla cdg  e inercias depósito v2.xlsx', 'B3:N17');


[percentage_FUEL OUTPUT] = CG_TARSIS120_KSA_V6(W_PL, T120, FUEL, caso)

porcentaje_actual = 0.5;
[percentage_FUEL OUTPUT] = CG_TARSIS120_KSA_V6B(W_PL, T120, FUEL, caso,porcentaje_actual)