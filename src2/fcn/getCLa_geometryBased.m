function CLa = getCLa_geometryBased(AR,LAMc2,Minf,Cla)

beta = sqrt(1 - Minf^2);

if isempty(Cla)
    k       = 1;
else
    Cla_M   = Cla/beta; % Cla_0 Pendiente de la curva de sustentacion del perfil a angulo de ataque nulo
    k       = Cla_M/2/pi;
    if k>1
        k = 1;
    end
end

%CLa     = AR*CLa_AR_calc(AR,beta,LAMc2,k)  
CLa     = 2*pi*AR/(2 + sqrt((AR*beta/k)^2*(1+(tan(LAMc2)/beta)^2)+4));
end