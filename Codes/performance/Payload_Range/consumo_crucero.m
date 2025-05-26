%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = consumo_crucero(V,h,S,C_D0,C_D1,C_D2,m,g)
%Función que devuelve el consumo en crucero (es necesario determinar la 
%posición de palanca)
rho = densityISA(h);

persistent delta_p
if isempty(delta_p)
    delta_p = 0.5;
end

delta_p = fzero(@(delta_p) .5*rho*V^2*S*(C_D0+C_D1*(m*g/(.5*rho*S*V^2))+C_D2*(m*g/(.5*rho*S*V^2))^2) - propulsive_model(V,h,delta_p),delta_p);
[T,c] = propulsive_model(V,h,delta_p);
f = c;