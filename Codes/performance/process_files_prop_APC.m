function PROPDATA_APC = process_files_prop_APC(casos_prop_APC,Prop_data_APC)

for i=1:length(casos_prop_APC)
    cas = i;
    % Identifies the max value of J to ensure that polifit is for those
    % values of J max
    J_limit_max = max(Prop_data_APC(cas).J);
    
    % CT EStimation
    n_order_ct = 2;
    Ct_APC{i} = polyfit(Prop_data_APC(cas).J(Prop_data_APC(cas).J<=J_limit_max),Prop_data_APC(cas).Ct(Prop_data_APC(cas).Ct<=J_limit_max),n_order_ct);
    Ct_total =  Ct_APC{1}(n_order_ct+1);
    for j=1:n_order_ct
        Ct_intermediate = Ct_APC{i}(j).*Prop_data_APC(cas).J.^(n_order_ct+1-j);
        Ct_total = Ct_total + Ct_intermediate;
    end
    CT_Polyfit_APC{i} = Ct_total;
    N_order_CT{i} = n_order_ct;
    
    % CP estimation
    n_order_cp = 3;
    Cp_APC{i} = polyfit(Prop_data_APC(cas).J(Prop_data_APC(cas).J<=J_limit_max),Prop_data_APC(cas).Cp(Prop_data_APC(cas).Cp<=J_limit_max),n_order_cp);
    Cp_total =  Cp_APC{1}(n_order_cp+1);
    for j=1:n_order_cp
        Cp_intermediate = Cp_APC{i}(j).*Prop_data_APC(cas).J.^(n_order_cp+1-j);
        Cp_total = Cp_total + Cp_intermediate;
    end
    CP_Polyfit_APC{i} = Cp_total;
    N_order_CP{i} = n_order_cp;
    
    % eta estimation
    n_order_eta = 3;
%     size(Prop_data_APC(cas).J(Prop_data_APC(cas).J<=J_limit_max))
%     size(Prop_data_APC(cas).Pe(Prop_data_APC(cas).Pe<=J_limit_max))
%     size(Prop_data_APC(cas).J);
%     size(Prop_data_APC(cas).Pe);
%     etap_APC{i} = polyfit(Prop_data_APC(cas).J(Prop_data_APC(cas).J<=J_limit_max),Prop_data_APC(cas).Pe(Prop_data_APC(cas).Pe<=J_limit_max),n_order_eta);
    etap_APC{i} = polyfit(Prop_data_APC(cas).J,Prop_data_APC(cas).Pe,n_order_eta);
    eta_total =  etap_APC{1}(n_order_eta+1);
    for j=1:n_order_eta
        eta_intermediate = etap_APC{i}(j).*Prop_data_APC(cas).J.^(n_order_eta+1-j);
        eta_total = eta_total + eta_intermediate;
    end
    etap_Polyfit_APC{i} = eta_total;
    N_order_etap{i} = n_order_eta;
end

% Stores data
PROPDATA_APC.CT_Polyfit_APC = Ct_APC;
PROPDATA_APC.N_order_CT = N_order_CT;

PROPDATA_APC.CP_Polyfit_APC = Cp_APC;
PROPDATA_APC.N_order_CP = N_order_CP;
% Same coefficients as CP
PROPDATA_APC.CQ_Polyfit_APC = Cp_APC;
PROPDATA_APC.N_order_CQ = N_order_CP;

PROPDATA_APC.etap_Polyfit_APC = etap_APC;
PROPDATA_APC.N_order_etap = N_order_etap;