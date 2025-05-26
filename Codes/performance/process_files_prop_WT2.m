function PROPDATA_WT2 = process_files_prop_WT2(Prop_data_WT2)

% %% Inputs
% rpms = [3000,2500];
% angles = {{'00';'05';'10';'15';'20'};{'30';'45';'60';'75';'90'}};
% for i=1:length(rpms)
%     data.(['rpm',num2str(rpms(i))]).angles = str2num(cell2mat(angles{i}));
%     for j=1:length(angles{i})        
%         filename = ['22x12W_RPM_',num2str(rpms(i)),'_',angles{i}{j}];
%         A = load([dir.data,filename]);
%         data.(['rpm',num2str(rpms(i))]).J{j} = A.J;
%         data.(['rpm',num2str(rpms(i))]).Cp{j} = A.Cp;
%         data.(['rpm',num2str(rpms(i))]).Cq{j} = A.Cq;
%         data.(['rpm',num2str(rpms(i))]).Ct{j} = A.Ct;
%     end
% end

% For 3000 RPM and alpha = 0;
data_J_WT2 = Prop_data_WT2.rpm3000.J{1};
data_Ct_WT2 = Prop_data_WT2.rpm3000.Ct{1};
data_Cp_WT2 = Prop_data_WT2.rpm3000.Cp{1};
data_etap_WT2 = Prop_data_WT2.rpm3000.J{1}.*(Prop_data_WT2.rpm3000.Ct{1}./Prop_data_WT2.rpm3000.Cp{1});

J_limit_max = max(data_J_WT2);
n_order_ct = 3;
Ct_WT2 = polyfit(data_J_WT2(data_J_WT2<=J_limit_max),data_Ct_WT2(data_Ct_WT2<=J_limit_max),n_order_ct);
Ct_total =  Ct_WT2(n_order_ct+1);
for j=1:n_order_ct
    Ct_intermediate = Ct_WT2(j).*data_J_WT2.^(n_order_ct+1-j);
    Ct_total = Ct_total + Ct_intermediate;
end
CT_Polyfit_WT2 = Ct_total;
N_order_CT = n_order_ct;

% CP estimation
n_order_cp = 3;
Cp_WT2 = polyfit(data_J_WT2(data_J_WT2<=J_limit_max),data_Cp_WT2(data_Cp_WT2<=J_limit_max),n_order_cp);
Cp_total =  Cp_WT2(n_order_cp+1);
for j=1:n_order_cp
    Cp_intermediate = Cp_WT2(j).*data_J_WT2.^(n_order_cp+1-j);
    Cp_total = Cp_total + Cp_intermediate;
end
CP_Polyfit_WT2 = Cp_total;
N_order_CP = n_order_cp;

% eta estimation
n_order_eta = 9;
% etap_WT2 = polyfit(data_J_WT2(data_J_WT2<=J_limit_max),data_etap_WT2(data_etap_WT2<=J_limit_max),n_order_eta);
etap_WT2 = polyfit(data_J_WT2,data_etap_WT2,n_order_eta);
eta_total =  etap_WT2(n_order_eta+1);
for j=1:n_order_eta
    eta_intermediate = etap_WT2(j).*data_J_WT2.^(n_order_eta+1-j);
    eta_total = eta_total + eta_intermediate;
end
etap_Polyfit_WT2 = eta_total;
N_order_etap = n_order_eta;

PROPDATA_WT2.CT_Polyfit_WT2 = CT_Polyfit_WT2;
PROPDATA_WT2.N_order_CT = N_order_CT;
PROPDATA_WT2.CP_Polyfit_WT2 = CP_Polyfit_WT2;
PROPDATA_WT2.N_order_CP = N_order_CP;
PROPDATA_WT2.etap_Polyfit_WT2 = etap_Polyfit_WT2;
PROPDATA_WT2.N_order_etap = N_order_etap;
