%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   27 de Agosto de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Cálculo K_a: Correlation constant for yawing moment due to aleiron
% deflection.
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 Fig. 10.48 (pag 448, 2138 PDF)

function Ka = Ka_calc(y1_b2, AR, TR)
    TR_pos = [0.25 0.5 0.75 1];
    [~,ii1] = min(abs(TR-TR_pos));
    switch ii1
        case 1
            data = load('Ka_025.dat');
        case 2
            data = load('Ka_05.dat');
        case 3
            data = load('Ka_075.dat');
        case 4
            data = load('Ka_1.dat');
    end

    AR_pos      = [4 6 8];
    [~,ii2]     = min(abs(AR-AR_pos));
    nPuntos     = 6;
    
    kk_i = (ii2-1)*nPuntos+1;
    kk_f = kk_i + (nPuntos - 1);
    data = data(kk_i:kk_f,:);
    
    Ka = interp1(data(:,1),data(:,2), y1_b2, 'pchip');
end
