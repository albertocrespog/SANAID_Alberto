function depsilon_dalpha = depsilon_dalpha_calc(x, crexp_w, XLE_w, l_fus, downwash,CLa_WB)
    [n,m] = size(x);
    depsilon_dalpha = zeros(n,m);
    
    XTE_w   = XLE_w + crexp_w;
    
for i = 1:max(n,m)    
    if x(i) <= 5*XLE_w/6
        data = load('depsilon_dalpha2.dat');
        depsilon_dalpha(i) = (1+interp1(data(:,1),data(:,2),x(i)/crexp_w,'pchip')*CLa_WB*pi/180/0.0785);
    elseif x(i) >= 5*XLE_w/6 && x(i) <= XLE_w 
        data = load('depsilon_dalpha.dat');
        depsilon_dalpha(i) = (1+interp1(data(:,1),data(:,2),x(i)/crexp_w,'pchip')*CLa_WB*pi/180/0.0785);
    elseif x(i) > XLE_w && x(i) < XTE_w
        depsilon_dalpha(i) = 0;
    elseif x(i) >= XTE_w
        depsilon_dalpha(i) = (x(i) - XTE_w)/(l_fus - XTE_w)*(1-downwash);
    end

end

