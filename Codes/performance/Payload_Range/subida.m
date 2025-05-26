%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = subida(h,y,V_EAS,S,C_D0,C_D1,C_D2,g,rho_SL)
m = y(1);

[rho,drho_dh] = densityISA(h); 

V     = V_EAS*sqrt(rho_SL/rho);
[T,c] = propulsive_model(V,h,1);

gamma = (T-.5*rho*V^2*S*(C_D0+C_D1*(m*g/(.5*rho*S*V^2))+C_D2*(m*g/(.5*rho*S*V^2))^2))/(m*g-.5*m*V^2/rho*drho_dh);

f = -c/(V*gamma);

