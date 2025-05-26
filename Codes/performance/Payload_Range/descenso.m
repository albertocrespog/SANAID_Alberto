
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = descenso(h,y,V_EAS,S,C_D0,C_D1,C_D2,g,rho_SL)
m = y(1);

[rho,drho_dh] = densityISA(h); 

V     = V_EAS*sqrt(rho_SL/rho);
T = 0;
c  = .01/60; %[kg/s] %Consumo con el motor en ralentí

gamma = (T-.5*rho*V^2*S*(C_D0+C_D1*(m*g/(.5*rho*S*V^2))+C_D2*(m*g/(.5*rho*S*V^2))^2))/(m*g-.5*m*V^2/rho*drho_dh);

f = -c/(V*gamma);

