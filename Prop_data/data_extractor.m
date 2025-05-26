clear all, close all, clc

% This script loads the experimental data of the propeller at different
% angles of attack 

%% Directories

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

%% Plot
Cp_plot = figure;
Cq_plot = figure;
Ct_plot = figure;

% % As lines
% 
% for i=1:length(rpms)
%     for j=1:length(data.(['rpm',num2str(rpms(i))]).angles)
%         figure(Cp_plot);
%         hold all
%         plot3(data.(['rpm',num2str(rpms(i))]).J{j}, data.(['rpm',num2str(rpms(i))]).angles(j)*ones(1,10), data.(['rpm',num2str(rpms(i))]).Cp{j})
%         figure(Cq_plot);
%         hold all
%         plot3(data.(['rpm',num2str(rpms(i))]).J{j}, data.(['rpm',num2str(rpms(i))]).angles(j)*ones(1,10), data.(['rpm',num2str(rpms(i))]).Cp{j})
%         figure(Ct_plot);
%         hold all
%         plot3(data.(['rpm',num2str(rpms(i))]).J{j}, data.(['rpm',num2str(rpms(i))]).angles(j)*ones(1,10), data.(['rpm',num2str(rpms(i))]).Cp{j})
%     end
% end

% As surface
% J_vec       = [];
% angle_vec   = [];
% Cp_vec      = [];
% Cq_vec      = [];
% Ct_vec      = [];

for i=1:length(rpms)
    for j=1:length(data.(['rpm',num2str(rpms(i))]).angles)
        angle_mat(:,j+(i-1)*5)  = data.(['rpm',num2str(rpms(i))]).angles(j)*ones(10,1);
        J_mat(:,j+(i-1)*5)      = data.(['rpm',num2str(rpms(i))]).J{j}';
        Cp_mat(:,j+(i-1)*5)     = data.(['rpm',num2str(rpms(i))]).Cp{j}';
        Cq_mat(:,j+(i-1)*5)     = data.(['rpm',num2str(rpms(i))]).Cq{j}';
        Ct_mat(:,j+(i-1)*5)     = data.(['rpm',num2str(rpms(i))]).Ct{j}';
        
%         J_vec       = [J_vec; data.(['rpm',num2str(rpms(i))]).J{j}'];
%         angle_vec   = [angle_vec; data.(['rpm',num2str(rpms(i))]).angles(j)*ones(10,1)];
%         Cp_vec      = [Cp_vec; data.(['rpm',num2str(rpms(i))]).Cp{j}'];
%         Cq_vec      = [Cq_vec; data.(['rpm',num2str(rpms(i))]).Cq{j}'];
%         Ct_vec      = [Ct_vec; data.(['rpm',num2str(rpms(i))]).Ct{j}'];
    end
end

% Angles converted into rads
angle_mat = angle_mat*D2R;
%angle_vec = angle_vec*D2R;

figure(Cp_plot);
surf(J_mat,angle_mat*R2D,Cp_mat)
shading interp
colormap(jet)

ylabel('\alpha (deg)')
xlabel('J (-)')
zlabel('Cp (-)')
grid on
figure(Cq_plot);
surf(J_mat,angle_mat*R2D,Cq_mat)
ylabel('\alpha (deg)')
xlabel('J (-)')
zlabel('Cq (-)')
grid on 
figure(Ct_plot);
surf(J_mat,angle_mat*R2D,Ct_mat)
shading interp
colormap(jet)
colorbar
ylabel('\alpha (deg)')
xlabel('J (-)')
zlabel('Ct (-)')
grid on

