function Cl_da = getClda(modelo)
%(y0, yf, b, lam, AR, LAMc4, t_c, ca_c, Cla, Minf)
    %% DATOS REQUERIDOS DEL MODELO
    Minf    = modelo.general.Minf;
    
    yf      = modelo.ala.y1;
    y0      = modelo.ala.y0;
    b       = modelo.ala.b;
    TR      = modelo.ala.TR;
    LAMc4   = modelo.ala.LAMc4;
    ca_c    = modelo.ala.ca_c;
    t_c     = modelo.ala.t_c;
    Cla     = modelo.ala.Cla;

    %% CALCULO DE LA DERIVADA DE CONTROL
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 442, 2137 PDF)
    
    beta = sqrt(1-Minf^2);

    betaCldelta_kf  = betaCldelta_k_calc(yf,b,TR,AR,LAMc4,Minf);
    betaCldelta_k0  = betaCldelta_k_calc(y0,b,TR,AR,LAMc4,Minf);
    Cldelta_prima   = 1/beta*(betaCldelta_kf - betaCldelta_k0);

    Cldelta_int     = Cldelta_Cldeltatheory_calc(ca_c,Cla)*Cldeltatheory_calc(ca_c,t_c);

    alpha_delta     = Cldelta_int/Cla;

    Cl_da = Cldelta_prima*alpha_delta;
end


