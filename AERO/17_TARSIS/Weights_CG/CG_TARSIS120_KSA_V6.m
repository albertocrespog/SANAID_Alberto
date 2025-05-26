function [percentage_FUEL OUTPUT] =  CG_TARSIS120_KSA_V6(W_PL, T120, excel_FUEL, caso)

%% Extracción por componentes [Masa, Xcg, Ycg, Zcg, Ix, Iy, Iz, Ixy, Ixz, Iyz]
FUS     = T120(1,:);
WNG_LFT = T120(2,:);
WNG_RGT = T120(3,:);
VTP_LFT = T120(4,:);
VTP_RGT = T120(5,:);
HTP     = T120(6,:);

KSA_CPR = T120(7,:);
KSA_BAT = T120(8,:);

PLD     = T120(9,:);


%% Ajuste masa PAYLOAD (En V6, PLD(1)=3.75kg)
PLD(1) = W_PL;


%% Ajuste masa FUEL (En V6, FUEL(1)=8.11kg)
if caso == 0
    % Ni BAT ni rack+contrapeso
    fuel_max_real = min(20.327, 120 - (FUS(1) + WNG_LFT(1) + WNG_RGT(1) + VTP_LFT(1) + VTP_RGT(1) + HTP(1) + PLD(1)));
elseif caso == 1
    % BAT no presente
    fuel_max_real = min(20.327, 120 - (FUS(1) + WNG_LFT(1) + WNG_RGT(1) + VTP_LFT(1) + VTP_RGT(1) + HTP(1) + KSA_CPR(1) + PLD(1)));

elseif caso == 2
    % BAT sí presente
    fuel_max_real = min(20.327, 120 - (FUS(1) + WNG_LFT(1) + WNG_RGT(1) + VTP_LFT(1) + VTP_RGT(1) + HTP(1) + KSA_CPR(1) + KSA_BAT(1) + PLD(1)));

else
    disp('Error en la configuración')
end

MTOW_fuel = fuel_max_real;


%% Ajuste inercias CATIA-EXPERIMENTALES (excepto FUEL) 

FUS(5:10)     =     FUS(5:10) * (FUS(1)/63.987);
WNG_LFT(5:10) = WNG_LFT(5:10) * (WNG_LFT(1)/7.122);
WNG_RGT(5:10) = WNG_RGT(5:10) * (WNG_RGT(1)/6.893);
VTP_LFT(5:10) = VTP_LFT(5:10) * (VTP_LFT(1)/3.327);
VTP_RGT(5:10) = VTP_RGT(5:10) * (VTP_RGT(1)/3.327);
HTP(5:10)     =     HTP(5:10) * (HTP(1)/2.604);

PLD(5:10)     =     PLD(5:10) * (PLD(1)/3.75);


%% Ajuste inercias y cg FUEL

percentage_FUEL = interp1(excel_FUEL(:,12), excel_FUEL(:,13), MTOW_fuel, 'spline')

cg_FUEL   = interp1(excel_FUEL(:,12), excel_FUEL(:,2:4), MTOW_fuel, 'spline');
iner_FUEL = interp1(excel_FUEL(:,12), excel_FUEL(:,5:10), MTOW_fuel, 'spline');


%% CENTRO DE GRAVEDAD (desde MAMPARO)

if caso == 0
    a = 0; b = 0;
elseif caso == 1
    a = 1; b = 0;
elseif caso == 2
    a = 1; b = 1;
end

m=1; x=2; y=3; z=4;

sum_masas = FUS(1) + WNG_LFT(1) + WNG_RGT(1) + VTP_LFT(1) + VTP_RGT(1) + HTP(1) + a*KSA_CPR(1) + b*KSA_BAT(1) + PLD(1) + MTOW_fuel; %BAT: [0,1]

x_cg_mamp = (FUS(m)*FUS(x) + WNG_LFT(m)*WNG_LFT(x) + WNG_RGT(m)*WNG_RGT(x) + VTP_LFT(m)*VTP_LFT(x) + VTP_RGT(m)*VTP_RGT(x) + HTP(m)*HTP(x) + ...
             a*KSA_CPR(m)*KSA_CPR(x) + b*KSA_BAT(m)*KSA_BAT(x) + PLD(m)*PLD(x) + MTOW_fuel*cg_FUEL(1)) / sum_masas; % [mm]

y_cg_mamp = (FUS(m)*FUS(y) + WNG_LFT(m)*WNG_LFT(y) + WNG_RGT(m)*WNG_RGT(y) + VTP_LFT(m)*VTP_LFT(y) + VTP_RGT(m)*VTP_RGT(y) + HTP(m)*HTP(y) + ...
             a*KSA_CPR(m)*KSA_CPR(y) + b*KSA_BAT(m)*KSA_BAT(y) + PLD(m)*PLD(y) + MTOW_fuel*cg_FUEL(2)) / sum_masas; % [mm]

z_cg_mamp = (FUS(m)*FUS(z) + WNG_LFT(m)*WNG_LFT(z) + WNG_RGT(m)*WNG_RGT(z) + VTP_LFT(m)*VTP_LFT(z) + VTP_RGT(m)*VTP_RGT(z) + HTP(m)*HTP(z) + ...
             a*KSA_CPR(m)*KSA_CPR(z) + b*KSA_BAT(m)*KSA_BAT(z) + PLD(m)*PLD(z) + MTOW_fuel*cg_FUEL(3)) / sum_masas; % [mm]

% En m.
x_cg_mamp = x_cg_mamp/1000; 
y_cg_mamp = y_cg_mamp/1000;
z_cg_mamp = z_cg_mamp/1000;


%% CENTRO DE GRAVEDAD (desde MORRO)

x_morro  = 490.224/1000; % Distancia morro-mamparo [m]
CG_morro = [x_cg_mamp + x_morro, y_cg_mamp, z_cg_mamp]; % [m]


%% Inercias en MAMPARO

I_total_mamparo = FUS(5:10) + WNG_LFT(5:10) + WNG_RGT(5:10) + VTP_LFT(5:10) + VTP_RGT(5:10) + HTP(5:10) + a*KSA_CPR(5:10) + ...
                  b*KSA_BAT(5:10) + PLD(5:10) + iner_FUEL; % [kg*m2]


%% Inercias en CG (T. Steiner)

I_total_cg(1) = I_total_mamparo(1) - sum_masas*(y_cg_mamp^2 + z_cg_mamp^2);
I_total_cg(2) = I_total_mamparo(2) - sum_masas*(x_cg_mamp^2 + z_cg_mamp^2);
I_total_cg(3) = I_total_mamparo(3) - sum_masas*(x_cg_mamp^2 + y_cg_mamp^2);

I_total_cg(4) = I_total_mamparo(4) + sum_masas*x_cg_mamp^2*y_cg_mamp^2;
I_total_cg(5) = I_total_mamparo(5) + sum_masas*x_cg_mamp^2*z_cg_mamp^2;
I_total_cg(6) = I_total_mamparo(6) + sum_masas*y_cg_mamp^2*z_cg_mamp^2;

OUTPUT=[sum_masas, CG_morro(1), CG_morro(2), CG_morro(3), I_total_cg]'; % [kg, m, m, m, kg*m2, kg*m2, kg*m2]
format long g

end