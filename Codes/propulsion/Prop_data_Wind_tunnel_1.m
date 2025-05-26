function Prop_data_WT1 = Prop_data_Wind_tunnel_1(Prop_type)

% Loads data for Wind Tunneñ Model 1
load('Data.mat')
load('Data_Prop.mat')

%---------------------------ENGINE DATA----------------------------------%
%% Propulsive Model
% Number of Prop used
k = Prop_type; % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/

% Thrust coefficients from wind tunnel
% Polyfit Coefficients for APC 22x10
% CT_Polyfit{k} = CT{k}(3) + CT{k}(2).*XData_CT{k} + CT{k}(1).*XData_CT{k}.^2;
CT2 = CT{k}(1);
CT1 = CT{k}(2);
CT0 = CT{k}(3);
Prop_data_WT1.CT0 = CT0;
Prop_data_WT1.CT1 = CT1;
Prop_data_WT1.CT2 = CT2;

% CP_Polyfit{k} = CP{k}(4) + CP{k}(3).*XData_CP{k} + CP{k}(2).*XData_CP{k}.^2 + CP{k}(1).*XData_CP{k}.^3;         
CP3 = CP{k}(1);
CP2 = CP{k}(2);
CP1 = CP{k}(3);
CP0 = CP{k}(4);
Prop_data_WT1.CP0 = CP0;
Prop_data_WT1.CP1 = CP1;
Prop_data_WT1.CP2 = CP2;
Prop_data_WT1.CP3 = CP3;

% CQ_Polyfit{k} = CQ{k}(4) + CQ{k}(3).*XData_CQ{k} + CQ{k}(2).*XData_CQ{k}.^2 + CQ{k}(1).*XData_CQ{k}.^3;
CQ3 = CQ{k}(1);
CQ2 = CQ{k}(2);
CQ1 = CQ{k}(3);
CQ0 = CQ{k}(4);
Prop_data_WT1.CQ0 = CQ0;
Prop_data_WT1.CQ1 = CQ1;
Prop_data_WT1.CQ2 = CQ2;
Prop_data_WT1.CQ3 = CQ3;

% ethamp_Polyfit{k} = ethamp{k}(4) + ethamp{k}(3).*XData_ethamp{k} + ethamp{k}(2).*XData_ethamp{k}.^2 + ethamp{k}(1).*XData_ethamp{k}.^3;
etha_mp3 = ethamp{k}(1);
etha_mp2 = ethamp{k}(2);
etha_mp1 = ethamp{k}(3);
etha_mp0 = ethamp{k}(4);
Prop_data_WT1.etha_mp0 = etha_mp0;
Prop_data_WT1.etha_mp1 = etha_mp1;
Prop_data_WT1.etha_mp2 = etha_mp2;
Prop_data_WT1.etha_mp3 = etha_mp3;