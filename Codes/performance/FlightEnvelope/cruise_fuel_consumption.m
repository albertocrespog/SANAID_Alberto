%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = cruise_fuel_consumption(V,h,S,C_D0,C_D1,C_D2,m,g)
%Función que devuelve el consumo en crucero (es necesario determinar la 
%posición de palanca)
rho = densityISA(h);

[T,c,Vnmax,delta_p_max] = propulsive_model(V,h,1);

delta_p = fzero(@(delta_p) .5*rho*V^2*S*(C_D0+C_D1*(m*g/(.5*rho*S*V^2))+C_D2*(m*g/(.5*rho*S*V^2))^2) - propulsive_model(V,h,delta_p),[0 delta_p_max]);
[T,c] = propulsive_model(V,h,delta_p);
f = c;
