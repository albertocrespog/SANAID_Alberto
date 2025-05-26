function CNa_807 = CNa_807_calc(n_blades, beta_blade)
    % beta_blade in degrees!!
    n_pos   = [2 3 4 6];
    [~,jj]      = min(n_blades - n_pos);
    
    nPuntos = 10;
    jj_i = (jj-1)*nPuntos+1;
    jj_f = jj*nPuntos;
    data = load('CNa_807.dat');
    data = data(jj_i:jj_f,:);
    if beta_blade > data(end,1)
        beta_blade = data(end,1);
    elseif beta_blade < data(1,1)
        beta_blade = data(1,1);
    end
    CNa_807 = interp1(data(:,1), data(:,2), beta_blade, 'pchip');
end