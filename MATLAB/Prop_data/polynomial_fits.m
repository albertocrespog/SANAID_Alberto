close all

%% Cp polynomial fit
% Modl: z = X*p_xy
% z : Ct
% x : J
% y : AoA
% z(x,y) = p00         + p10*x     + p01*y       + p20*x^2     + p11*x*y +  
%          p02*y^2     + p30*x^3   + p21*x^2*y   + p12*x*y^2   + p03*y^3 + 
%          p40*x^4     + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3   + p04*y^4 + 
%          p50*x^5     + p41*x^4*y + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 + 
%          p05*y^5
x = J_mat(:);
y = angle_mat(:);   

X = [x,         y,          x.^2,       x.*y,       y.^2,       x.^3, ...
     x.^2.*y,   x.*y.^2,    y.^3,       x.^4,       x.^3.*y,    x.^2.*y.^2, ...
     x.*y.^3,   y.^4,       x.^5,       x.^4.*y,    x.^3.*y.^2, x.^2.*y.^3, ...
     x.*y.^4, y.^5];
X = [ones(length(x),1), X];

p_Cp = (X'*X)\X'*Cp_mat(:);

% Evaluation
x_ev = linspace(min(x), max(x),100);
y_ev = linspace(min(y), max(y),100);

Cp_ev = X*p_Cp;
Cp_fit = figure;
figure(Cp_fit)
surf(J_mat,angle_mat*R2D,Cp_mat)
ylabel('\alpha (deg)')
xlabel('J (-)')
zlabel('Cp (-)')
hold on
plot3(J_mat(:),angle_mat(:)*R2D,Cp_ev,'ko')


%% Cq polynomial fit
% Modl: z = X*p_xy
% z : Ct
% x : J
% y : AoA
% z(x,y) = p00         + p10*x     + p01*y       + p20*x^2     + p11*x*y +  
%          p02*y^2     + p30*x^3   + p21*x^2*y   + p12*x*y^2   + p03*y^3 + 
%          p40*x^4     + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3   + p04*y^4 + 
%          p50*x^5     + p41*x^4*y + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 + 
%          p05*y^5
x = J_mat(:);
y = angle_mat(:);   

X = [x,         y,          x.^2,       x.*y,       y.^2,       x.^3, ...
     x.^2.*y,   x.*y.^2,    y.^3,       x.^4,       x.^3.*y,    x.^2.*y.^2, ...
     x.*y.^3,   y.^4,       x.^5,       x.^4.*y,    x.^3.*y.^2, x.^2.*y.^3, ...
     x.*y.^4, y.^5];
X = [ones(length(x),1), X];

p_Cq = (X'*X)\X'*Cq_mat(:);

% Evaluation
x_ev = linspace(min(x), max(x),100);
y_ev = linspace(min(y), max(y),100);

Cq_ev = X*p_Cq;
Cq_fit = figure;
figure(Cq_fit)
surf(J_mat,angle_mat*R2D,Cq_mat)
ylabel('\alpha (deg)')
xlabel('J (-)')
zlabel('Cq (-)')
hold on
plot3(J_mat(:),angle_mat(:)*R2D,Cq_ev,'ko')


%% Ct polynomial fit
% Modl: z = X*p_xy
% z : Ct
% x : J
% y : AoA
% z(x,y) = p00         + p10*x     + p01*y       + p20*x^2     + p11*x*y +  
%          p02*y^2     + p30*x^3   + p21*x^2*y   + p12*x*y^2   + p03*y^3 + 
%          p31*x^3*y   + p22*x^2*y^2 + p13*x*y^3   + p04*y^4 
x = J_mat(:);
y = angle_mat(:);   

X = [x,         y,          x.^2,       x.*y,       y.^2,       x.^3, ...
     x.^2.*y,   x.*y.^2,    y.^3,       x.^3.*y,    x.^2.*y.^2, ...
     x.*y.^3,   y.^4];
X = [ones(length(x),1), X];

p_Ct = (X'*X)\X'*Ct_mat(:);

% Evaluation
x_ev = linspace(min(x), max(x),100);
y_ev = linspace(min(y), max(y),100);

Ct_ev = X*p_Ct;
Ct_fit = figure;
figure(Ct_fit)
surf(J_mat,angle_mat*R2D,Ct_mat)
ylabel('\alpha (deg)')
xlabel('J (-)')
zlabel('Ct (-)')
hold on
plot3(J_mat(:),angle_mat(:)*R2D,Ct_ev,'ko')
