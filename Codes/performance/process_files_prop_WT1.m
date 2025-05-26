function PROPDATA_WT1 = process_files_prop_WT1(casos_prop_WT1,Prop_data_WT1,N_props_WT1)

% Stores data
N_order_CT = 2;
N_order_CP = 3;
N_order_CQ = 3;
N_order_etamp = 3;

PROPDATA_WT1.N_order_CT = N_order_CT;
PROPDATA_WT1.N_order_CP = N_order_CP;
PROPDATA_WT1.N_order_CQ = N_order_CQ;
PROPDATA_WT1.N_order_etamp = N_order_etamp;

% Thrust coefficients from wind tunnel
% Polyfit Coefficients for APC 22x10
CT_Polyfit_WT1 = [Prop_data_WT1.CT2, Prop_data_WT1.CT1, Prop_data_WT1.CT0];
CP_Polyfit_WT1 = [Prop_data_WT1.CP3, Prop_data_WT1.CP2, Prop_data_WT1.CP1, Prop_data_WT1.CP0];
CQ_Polyfit_WT1 = [Prop_data_WT1.CQ3, Prop_data_WT1.CQ2, Prop_data_WT1.CQ1, Prop_data_WT1.CQ0];
% etap_Polyfit_WT1 = [Prop_data_WT1.etha_mp3, Prop_data_WT1.etha_mp2, Prop_data_WT1.etha_mp1, Prop_data_WT1.etha_mp0];
etamp_Polyfit_WT1 = [Prop_data_WT1.etha_mp3, Prop_data_WT1.etha_mp2, Prop_data_WT1.etha_mp1, Prop_data_WT1.etha_mp0];

% Stores data
PROPDATA_WT1.CT_Polyfit_WT1 = CT_Polyfit_WT1;
PROPDATA_WT1.CP_Polyfit_WT1 = CP_Polyfit_WT1;
PROPDATA_WT1.CQ_Polyfit_WT1 = CQ_Polyfit_WT1;
PROPDATA_WT1.etamp_Polyfit_WT1 = etamp_Polyfit_WT1;
