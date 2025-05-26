
% Cálculo preliminar de la derivada de estabilidad C_l_da.

function Stab_Der = get_aileron_latdir_der(modelo,Stab_Der)
%(y0, yf, b, lam, AR, LAMc4, t_c, ca_c, Cla, Minf)
    %% DATOS REQUERIDOS DEL MODELO
    Minf    = modelo.general.Minf;
    
    yf_b2   = modelo.ala.y1_b2;
    y0_b2   = modelo.ala.y0_b2;
    TR      = modelo.ala.TR;
    LAMc4   = modelo.ala.LAMc4;
    ca_c    = modelo.ala.cm_c;
    t_c     = modelo.ala.t_c;
    Cla     = modelo.ala.Cla;
    AR      = modelo.ala.AR;

    CL      = W/qinf/Sref; 
    
    %% CALCULO DE LA DERIVADA DE CONTROL C_y_da
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 442, 2137 PDF)
    Cy_da = 0;
    Stab_Der.Cydeltaa = Cy_da;
    %% CALCULO DE LA DERIVADA DE CONTROL C_l_da
    %   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 442, 2137 PDF)
    
    beta = sqrt(1-Minf^2);

    betaCldelta_kf  = betaCldelta_k_calc(yf_b2,TR,AR,LAMc4,Minf);
    betaCldelta_k0  = betaCldelta_k_calc(y0_b2,TR,AR,LAMc4,Minf);
    
    airfoil_data = importdata('NACA_0012.dat');
    [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
    clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
    clalpha_w2 = clalpha_w2_M0/beta;
    Cla_theo_w2 = 2*pi + 5.0525*t_c;
    Cla_M0 = clalpha_w2_M0/Cla_theo_w2;   
    
    k = clalpha_w2_M0/(2*pi);
    
    Cldelta_prima   = k/beta*(betaCldelta_kf - betaCldelta_k0);
    
    Cldelta_int     = Cldelta_Cldeltatheory_calc(ca_c,Cla_M0)*Cldeltatheory_calc(ca_c,t_c);
    
    

    alpha_delta     = Cldelta_int/clalpha_w2;
    
    Cl_da = Cldelta_prima*alpha_delta;
    
    Stab_Der.Cldeltaa = Cl_da;
    
    %% CALCULO DE Cn_da
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 448, 2138 PDF)
Ka      = Ka_calc(yf_b2, AR, TR);

Cn_da   = 2*Ka*CL*Cl_da;

Stab_Der.Cndeltaa = Cn_da;
end


