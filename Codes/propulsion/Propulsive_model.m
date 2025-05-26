function [ENGINE PROP]=Propulsive_model_v1(k)

load Data_Prop.mat
% APC 20x8
% APC 22x10
% APC 22x12
% APC 22x12W
% APC 21x13W
% APC 21x14

% Polyfit Coefficients for APC 22x10
k = 4;
% CT_Polyfit{i} = CT{i}(3) + CT{i}(2).*XData_CT{i} + CT{i}(1).*XData_CT{i}.^2;
CT2 = CT{k}(1);
CT1 = CT{k}(2);
CT0 = CT{k}(3);

CP3 = CP{k}(1);
CP2 = CP{k}(2);
CP1 = CP{k}(3);
CP0 = CP{k}(4);

etha_mp3 = ethamp{k}(1);
etha_mp2 = ethamp{k}(2);
etha_mp1 = ethamp{k}(3);
etha_mp0 = ethamp{k}(4);

Engine.CT0 = CT0;
Engine.CT1 = CT1;
Engine.CT2 = CT2;
Engine.CP0 = CP0;
Engine.CP1 = CP1;
Engine.CP2 = CP2;
Engine.CP3 = CP3;
Engine.etha_mp0 = etha_mp0;
Engine.etha_mp1 = etha_mp1;
Engine.etha_mp2 = etha_mp2;
Engine.etha_mp3 = etha_mp3;
Engine.RPM_max = RPM_max;

% Checks first for the range of valid CT looking for J @ CT=0
J_vec = linspace(0,1,100);
CT_vec = CT0 + CT1*J_vec + CT2*(J_vec.^2);
CP_vec = CP0 + CP1*J_vec + CP2*(J_vec.^2) + CP3*(J_vec.^3);
etha_vec = etha_mp0 + etha_mp1*J_vec + etha_mp2*(J_vec.^2) + etha_mp3*(J_vec.^3);

J_max = interp1(CT_vec,J_vec,0,'spline');

% Estimates new J vector
J_vec = linspace(0,J_max,100);
% Recalculates Vectors
CT_vec = CT0 + CT1*J_vec + CT2*(J_vec.^2);
CP_vec = CP0 + CP1*J_vec + CP2*(J_vec.^2) + CP3*(J_vec.^3);
etha_vec = etha_mp0 + etha_mp1*J_vec + etha_mp2*(J_vec.^2) + etha_mp3*(J_vec.^3);

etha_max = max(etha_vec);
J_eta_max = interp1(etha_vec,J_vec,etha_max,'spline');

PROP.CT = CT;
PROP.CP = CP;
PROP.ethamp = ethamp;

PROP.J_vec = J_vec;
PROP.CT_vec = CT_vec;
PROP.CP_vec = CP_vec;
PROP.etha_vec = etha_vec;

PROP.CT0 = CT0;
PROP.CT1 = CT1;
PROP.CT2 = CT2;

PROP.CP0 = CP0;
PROP.CP1 = CP1;
PROP.CP2 = CP2;
PROP.CP3 = CP3;

PROP.etha_mp0 = etha_mp0;
PROP.etha_mp1 = etha_mp1;
PROP.etha_mp2 = etha_mp2;
PROP.etha_mp3 = etha_mp3; 

PROP.etha_max = etha_max; 
PROP.J_eta_max = J_eta_max; 

figures_CT_model = 0;
fig = 0;
if figures_CT_model ==1
    fig = fig + 1
    figure(fig)
    plot(J_vec,CT_vec)
    title('C_T vs. J')
    xlabel('J')
    ylabel('C_T')
    grid on
    
    fig = fig + 1
    figure(fig)
    plot(J_vec,CP_vec)
    title('C_P vs. J')
    xlabel('J')
    ylabel('C_P')
    grid on
    
    fig = fig + 1
    figure(fig)
    plot(J_vec,etha_vec)
    title('\eta vs. J')
    xlabel('J')
    ylabel('\eta')
    grid on
end 