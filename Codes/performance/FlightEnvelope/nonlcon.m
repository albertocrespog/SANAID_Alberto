%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq] = nonlcon(x,S,C_D0,C_D1,C_D2,m,g,ROC_SC)
c = [];
ceq = propulsive_model(x(1),densityISA2ALT(x(2)),1) - (0.5*x(2)*x(1)^2*S*(C_D0+C_D1*(m*g/(0.5*x(2)*x(1)^2*S))+C_D2*(m*g/(0.5*x(2)*x(1)^2*S))^2)+m*g*ROC_SC/x(1));
