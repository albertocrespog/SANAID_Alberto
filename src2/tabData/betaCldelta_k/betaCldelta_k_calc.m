%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   12 de Agosto de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo de betaCldelta_k: Aleiron Rolling Moment parameter. AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.4, Fig. 10.46 (pag 443, 2134 PDF))
function betaCldelta_k = betaCldelta_k_calc(y_b,lam,AR,LAMc4,Minf)

beta = sqrt(1-Minf^2);
betaAR = AR*beta;
betaAR_pos = [2, 4, 8];
[~,ii1] = min(abs(betaAR-betaAR_pos));

switch ii1
    case 1
        data = load('betaClb_k_AR2.dat');
    case 2
        data = load('betaClb_k_AR4.dat');
    case 3
        data = load('betaClb_k_AR8.dat');
end

lam_pos = [0, 0.5, 1];
[~,ii2] = min(abs(lam - lam_pos));
nlambdaPuntos = 33;
kk_i = (ii2-1)*nlambdaPuntos+1;
kk_f = kk_i + (nlambdaPuntos - 1);
data = data(kk_i:kk_f,:);

LAMbeta = atan(1/beta*tan(LAMc4));
LAMbeta_pos = [60, 0, 40]*pi/180;
[~,ii3] = min(abs(LAMbeta - LAMbeta_pos));
nLAMbetaPuntos = 11;
kk_i = (ii3-1)*nLAMbetaPuntos+1;
kk_f = kk_i + (nLAMbetaPuntos - 1);
data = data(kk_i:kk_f,:);

betaCldelta_k = interp1(data(:,1),data(:,2),y_b,'pchip');
end

