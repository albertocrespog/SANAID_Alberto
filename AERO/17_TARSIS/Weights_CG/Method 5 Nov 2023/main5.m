% Weight Configuration from AERTEC
clear; close all; clc

MOTOR    = 1; % SP210
RACK     = 0; % sin rack
DEPOSITO = 1; % deposito original
FUEL     = 0.5; % en tanto por 1
W_PL     = 0; %carga de pago a considerar 
n_MSL    = 0; % se han sustituido los datos de la configuración de 3 micros por los de la configuración KSA SÓLO BOLA, sin racks

% Lectura datos T120
CATIAT120 = xlsread('Tabla cdg e inercias T120 KSA.xlsx','D4:M25');

%Lectura datos CG depósito
CG_DEP = xlsread('Tabla cdg  e inercias depósito.xlsx');

OUTPUT = CG_TARSIS120_KSA_V5(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120,CG_DEP)