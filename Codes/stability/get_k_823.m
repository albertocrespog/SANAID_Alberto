function K_nacareport = get_k_823(AR,lambda)
% 
A = readmatrix('Fig2 - NACA Report No_823.txt');
A = A(:,1:4);
[i,k] = find(isnan(A));
A(i,:)=[];
x = A(:,1);
lambda1 = A(:,2);
lambda0_5 = A(:,3);
lambda0_25= A(:,4);
if AR < x(1)
    AR = x(1);
elseif x > x(end)
    AR = x(end);
end
if lambda < 0.25
    lambda = 0.25;
end

X = x;
Y = [1, 0.5, 0.25];
Z=A(:,2:4)';
K_nacareport = interp2(X,Y,Z,AR,lambda);

% xx = linspace(x(1),x(end),1000);
% lambda1int = interp1(x,lambda1,xx);
% lambda05int = interp1(x,lambda0_5,xx);
% lambda025int = interp1(x,lambda0_25,xx);
% 
% plot(x,lambda1,'*')
% hold on
% plot(x,lambda0_5,'x')
% hold on
% plot(x,lambda0_25,'o')
% hold on
% plot(xx,lambda1int)
% hold on
% plot(xx,lambda05int)
% hold on
% plot(xx,lambda025int)
% hold on
% plot(AR,K_nacareport,'+')
% txt = ['\lambda = ', num2str(lambda) ];
% text(AR + 0.25,K_nacareport + 0.005,txt)
% title ('Fig. 2 NACA Report No. 823')
% xlabel('AR')
% ylabel('K','Rotation',0)
% ylim([.6 0.86])
% grid on
% legend('\lambda = 1', '\lambda = 0.5','\lambda = 0.25')
end