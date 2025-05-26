close all
clear all
casE = 1;
casos{1} = 'PER3_22x12WE_RPM1000.dat'; 
casos{2} = 'PER3_22x12WE_RPM2000.dat'; 
casos{3} = 'PER3_22x12WE_RPM3000.dat'; 
casos{4} = 'PER3_22x12WE_RPM4000.dat'; 
casos{5} = 'PER3_22x12WE_RPM5000.dat'; 
casos{6} = 'PER3_22x12WE_RPM6000.dat'; 
casos{7} = 'PER3_22x12WE_RPM7000.dat'; 

dir.current = pwd;
dir.data = [dir.current,'\data\'];

%% Units conversions

D2R = pi/180;
R2D = 180/pi;

%% Inputs
rpms = [3000,2500];
angles = {{'00';'05';'10';'15';'20'};{'30';'45';'60';'75';'90'}};
for i=1:length(rpms)
    data.(['rpm',num2str(rpms(i))]).angles = str2num(cell2mat(angles{i}));
    for j=1:length(angles{i})        
        filename = ['22x12W_RPM_',num2str(rpms(i)),'_',angles{i}{j}];
        A = load([dir.data,filename]);
        data.(['rpm',num2str(rpms(i))]).J{j} = A.J;
        data.(['rpm',num2str(rpms(i))]).Cp{j} = A.Cp;
        data.(['rpm',num2str(rpms(i))]).Cq{j} = A.Cq;
        data.(['rpm',num2str(rpms(i))]).Ct{j} = A.Ct;
    end
end

% % i=3;
% J_limit_max=0.7;
% for i=1:length(XData_CT)
%     CT{i} = polyfit(XData_CT{i}(XData_CT{i}<=J_limit_max),YData_CT{i}(XData_CT{i}<=J_limit_max),2);
%     CT_Polyfit{i} = CT{i}(3) + CT{i}(2).*XData_CT{i} + CT{i}(1).*XData_CT{i}.^2;
% end

%% Propulsive Model
% Number of Prop used
% Number of Prop used
% 1 - APC 20x8
% 2 - APC 22x10 
% 3 - APC 22x12 
% 4 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
% 5 - APC 21x13W 
% 6 - APC 21x14  
Prop_type = 4; 
k = Prop_type; % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/

load('Data.mat')
load('Data_Prop.mat')

for i =1:length(casos)
    % Extracción de datos XFLR5
    file_id = fopen(casos{1});
        aux = textscan(file_id,'%n%n%n%n%n%n%n%n','Headerlines',0);
        % New results 13 entries
        Data_P(i).raw_data= aux;
        Data_P(i).label = casos{i};
        Data_P(i).V = aux{1}*0.44704;      % m/s
        Data_P(i).J = aux{2};
        Data_P(i).Pe    = aux{3};
        Data_P(i).Ct   = aux{4};
        Data_P(i).Cp   = aux{5};
        Data_P(i).Hp    = aux{6}*745.6999; % Watts
        Data_P(i).Q    = aux{7};
        Data_P(i).T    = aux{8}*4.448222;  % Newtons
        fclose(file_id);
end


% i=3;
cas = 7;
J_limit_max= max(Data_P(cas).J);

for i=1:length(Data_P(cas))
    % CT EStimation
    CT_APC{i} = polyfit(Data_P(cas).J(Data_P(cas).J<=J_limit_max),Data_P(cas).Ct(Data_P(cas).Ct<=J_limit_max),2);
    CT_Polyfit_APC{i} = CT_APC{i}(3) + CT_APC{i}(2).*Data_P(cas).J + CT_APC{i}(1).*Data_P(cas).J.^2;
    
    % CP estimation
    CP_APC{i} = polyfit(Data_P(cas).J(Data_P(cas).J<=J_limit_max),Data_P(cas).Cp(Data_P(cas).Cp<=J_limit_max),3);
    CP_Polyfit_APC{i} = CP_APC{i}(4) + CP_APC{i}(3).*Data_P(cas).J + CP_APC{i}(2).*Data_P(cas).J.^2 + CP_APC{i}(1).*Data_P(cas).J.^3;

    % CP estimation
    n_order = 9;
    etap_APC{i} = polyfit(Data_P(cas).J(Data_P(cas).J<=J_limit_max),Data_P(cas).Pe(Data_P(cas).Pe<=J_limit_max),n_order);
    eta_total =  etap_APC{1}(n_order+1);
    for j=1:n_order
        eta_intermediate = etap_APC{i}(j).*Data_P(cas).J.^(n_order+1-j);
        eta_total = eta_total + eta_intermediate;
    end
        etap_Polyfit_APC{i} = eta_total;
%         etap_Polyfit_APC{i} = etap_APC{i}(6+1) + etap_APC{i}(5).*Data_P(cas).J + etap_APC{i}(4).*Data_P(cas).J.^2 + etap_APC{i}(3).*Data_P(cas).J.^3 +...
%         etap_APC{i}(2).*Data_P(cas).J.^4 + etap_APC{i}(1).*Data_P(cas).J.^5;
end

% for i=1:length(Data_P(1).J)
    % Thrust coefficients from wind tunnel
    % Polyfit Coefficients for APC 22x10
    CT_Polyfit{k} = CT{k}(3) + CT{k}(2).*Data_P(1).J + CT{k}(1).*Data_P(1).J.^2;
    % CT2 = CT{k}(1);
    % CT1 = CT{k}(2);
    % CT0 = CT{k}(3);
    
    CP_Polyfit{k} = CP{k}(4) + CP{k}(3).*Data_P(1).J + CP{k}(2).*Data_P(1).J.^2 + CP{k}(1).*Data_P(1).J.^3;
    % CP3 = CP{k}(1);
    % CP2 = CP{k}(2);
    % CP1 = CP{k}(3);
    % CP0 = CP{k}(4);
    
    CQ_Polyfit{k} = CQ{k}(4) + CQ{k}(3).*Data_P(1).J + CQ{k}(2).*Data_P(1).J.^2 + CQ{k}(1).*Data_P(1).J.^3;
    % CQ3 = CQ{k}(1);
    % CQ2 = CQ{k}(2);
    % CQ1 = CQ{k}(3);
    % CQ0 = CQ{k}(4);
    
    ethamp_Polyfit{k} = ethamp{k}(4) + ethamp{k}(3).*Data_P(1).J + ethamp{k}(2).*Data_P(1).J.^2 + ethamp{k}(1).*Data_P(1).J.^3;
    
    ethap_Polyfit{k} = CT_Polyfit{k}.*Data_P(1).J./CP_Polyfit{k};
    % etha_mp3 = ethamp{k}(1);
    % etha_mp2 = ethamp{k}(2);
    % etha_mp1 = ethamp{k}(3);
    % etha_mp0 = ethamp{k}(4);
% end

Fig = 1;

figure(Fig)
plot(Data_P(7).J,Data_P(7).Ct,'b')
hold on
% plot(Data_P(7).J,Data_P(7).Ct,'b:')
plot(Data_P(7).J,CT_Polyfit{k},'g--')
plot(data.rpm3000.J{1},data.rpm3000.Ct{1},'k-.')
plot(Data_P(7).J,CT_Polyfit_APC{1},'m--')
hold off
legend('APC RPM = 7000','Model-1','Model-2','APC RPM = 7000 - Polyfit')
title('C_T vs J')
xlabel('J')
xlabel('C_T')
grid on
Fig = Fig + 1;

figure(Fig)
plot(Data_P(7).J,Data_P(7).Cp,'b')
hold on
plot(Data_P(7).J,CP_Polyfit{k},'g--')
plot(data.rpm3000.J{1},data.rpm3000.Cp{1},'k-.')
plot(Data_P(7).J,CP_Polyfit_APC{1},'m--')
hold off
legend('APC RPM = 7000','Model-1','Model-2','APC RPM = 7000 - Polyfit')
title('C_P vs J')
xlabel('J')
xlabel('C_P')
grid on
Fig = Fig + 1;


eta_mod2 = (data.rpm3000.Ct{1}.*data.rpm3000.J{1})./data.rpm3000.Cp{1};
figure(Fig)
plot(Data_P(7).J,Data_P(7).Pe,'b')
hold on
plot(Data_P(7).J,ethamp_Polyfit{k},'g--')
plot(Data_P(7).J,ethap_Polyfit{k},'g*')
plot(data.rpm3000.J{1},eta_mod2,'k-.')
plot(Data_P(7).J,etap_Polyfit_APC{1},'m--')
hold off
legend('APC RPM = 7000','Model-1 - eta_{mp}','Model-1 - eta_{p}','Model-2','APC RPM = 7000 - Polyfit')
% hold off
% legend('RPM = 7000','Wind Tunnel - \eta_{p}','Wind Tunnel - \eta_{mp}')
title('\eta vs J')
xlabel('J')
xlabel('\eta')
grid on
Fig = Fig + 1;
