function [AB]=ajuste_AB(deltae,n_Pzero,n_ref,P_ref)

% Valores de n para los que disponemos de datos (Tabla 41, pág 54, Documento System Description SP210)
n_aux=[1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000];
n=n_aux./60; % pasamos de rpm a rev/s
% i=0;

if deltae==100
P_aux=[1,1.9,2.8,3.6,4.9,6.0,7.1,8.2,9.1,10.4,10.9];  % valores de 'y'
elseif deltae==90
    P_aux=[1,1.9,2.8,3.5,4.8,6,7,8.1,9.2,10.4,11.1];
elseif deltae==85
    P_aux=[1,1.9,2.7,3.4,4.7,5.9,6.9,8,9.2,10.4,11];
elseif deltae==80
    P_aux=[1,1.8,2.8,3.5,4.7,5.8,6.8,8,9.3,10.3,10.9];
elseif deltae==75
    P_aux=[1,1.8,2.7,3.5,4.7,5.7,6.7,8,9.2,10.4,11];
elseif deltae==70
    P_aux=[1,1.6,2.7,3.4,4.6,5.6,6.4,7.9,9.3,10.4,10.9];
elseif deltae==65
    P_aux=[0.9,1.6,2.6,3.3,4.4,5.5,6.3,7.8,9.2,10.4,10.7];
elseif deltae==60
    P_aux=[0.9,1.5,2.6,3.3,4.3,5.3,5.9,7.5,8.8,10,10.6];
elseif deltae==55
    P_aux=[0.9,1.4,2.6,3.3,4.1,5.1,5.9,7.3,8.4,9.5,10.1];
elseif deltae==50
    P_aux=[0.7,1.3,2.5,3,4,4.8,5.5,7.1,8.2,9.3,9.9];
elseif deltae==40
    P_aux=[0.6,1.3,2.3,2.7,3.4,4.1,4.8,6.2,7.2,8.2,8.6];
elseif deltae==30
    P_aux=[0.4,0.9,1.5,2,2.7,3.1,3.6,5.1,6.1,7.1,7.5];
elseif deltae==20
    P_aux=[0.3,0.5,0.7,0.9,1.4,1.8,2.2,3.9,4.7,5.5,5.9];
elseif deltae==15
    P_aux=[0.1,0.2,0.4,0.6,0.9,1.2,1.8,3.2,4.2,5.1,5.3];
elseif deltae==10
    P_aux=[0.1,0.2,0.2,0.4,0.6,0.8,1.3,2.6,3.4,3.9,4.8];
end

% if i==0 % no se ha dado deltae=10
P=P_aux.*1000; % pasamos de kW a W

% Definos funcion y matrices para resolución por minimos cuadrados
f=@(ab) P_ref*(ab(2).^2-(ab(1).^2)*(n/n_ref)).*(n./n_ref-n_Pzero/n_ref);

M1=@(ab) sum(((P-(f(ab))).*((-ab(1).*n/n_ref).*(n/n_ref-n_Pzero/n_ref))));
M2=@(ab) sum(((P-(f(ab))).*ab(2).*(n-n_Pzero)));
M=@(ab) [M1(ab),M2(ab)];
% M([1,1])

ab1=fsolve(M,[10,10]);
AB=ab1.^2;

%AB=[-0.5,0]
f=@(ab,n)  P_ref*(ab(2).^2-(ab(1).^2)*(n/n_ref)).*(n./n_ref-n_Pzero/n_ref);
plot(n,P,'or')
hold on
xx=16:150;
yy=f(ab1,xx);
plot(xx,yy,'blue')
hold off
end


