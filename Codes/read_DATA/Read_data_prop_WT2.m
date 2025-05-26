function [Prop_data_WT2] = Read_data_prop_WT2

dir.current = pwd;
dir.data = [dir.current,'..\Prop_data\Data_WT2\'];

%% Units conversions
D2R = pi/180;
R2D = 180/pi;
%% Inputs
rpms = [3000,2500];
angles = {{'00';'05';'10';'15';'20'};{'30';'45';'60';'75';'90'}};
for i=1:length(rpms)
    Prop_data_WT2.(['rpm',num2str(rpms(i))]).angles = str2num(cell2mat(angles{i}));
    for j=1:length(angles{i})        
        filename = ['22x12W_RPM_',num2str(rpms(i)),'_',angles{i}{j}];
%         A = load([dir.data,filename]);
        A = load(filename);
        Prop_data_WT2.(['rpm',num2str(rpms(i))]).J{j} = A.J;
        Prop_data_WT2.(['rpm',num2str(rpms(i))]).Cp{j} = A.Cp;
        Prop_data_WT2.(['rpm',num2str(rpms(i))]).Cq{j} = A.Cq;
        Prop_data_WT2.(['rpm',num2str(rpms(i))]).Ct{j} = A.Ct;
    end
end
