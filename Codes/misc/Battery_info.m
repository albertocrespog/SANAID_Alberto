% 3DR IRIS+ 3S 5100MAH 11.1V RC LIPO DRONE QUADCOPTER BATTERY W/XT60 PLUG BY VENOM
% https://www.venompower.com/3dr-iris-3s-5100mah-11-1v-rc-lipo-drone-quadcopter-battery-w-xt60-plug-by-venom-35028
% % Producto: Lipo Batería
% Battery Type: Lithium Polymer (LiPo Battery)
% C Rate: 8C
% Volts: 11.1
% Capacity: 5100mAh
% Cell Count: 3S
% Cell Configuration: 3S1P
% Continuous Discharge: 8C (40.8A)
% Max Burst Rate: 10C (51A)
% Max Volts per Pack: 12.6V
% Min Volts per Pack: 9V
% Charge Rate: 1C (5.1A)
% Wire Gauge: 12 AWG Soft and Flexible Low Resistance Silicone Wire
% Plug Type: XT60
% Dimensions: 131 x 44 x 24.5 mm / 5.2 x 1.7 x 1 in
% Watt Hours: 56.61
% Weight: 10.8 oz (306 g)

Capacidad= 8000; %mAh
Cells = 6; % number of cells
Discharge = 60;
% Dimensions
w_bat = 46; % wide
h_bat = 90; % height
l_bat = 158; % length
% Volumen 
Vol_bat = w_bat*h_bat*l_bat;
Weight_batt = 1242/1000;
Weight_cell = Weight_batt/Cells;
Lipo_V_cell = 3.7;
V_bat = Lipo_V_cell*Cells; 