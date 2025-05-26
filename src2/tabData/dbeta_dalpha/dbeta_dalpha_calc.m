function dbeta_dalpha = dbeta_dalpha_calc(x1, d, tramo)
    
    data = load('dbeta_dalpha.dat');
    if tramo < 6
        data = data(1:8,:);
        dbeta_dalpha = interp1(data(:,1),data(:,2),x1/d,'pchip');
    elseif tramo == 6
        data = data(9:end,:);
        dbeta_dalpha = interp1(data(:,1),data(:,2),x1/d,'pchip');
%     elseif tramo > 6
%         dbeta_dalpha = x1/d*(1 - downw);
    end

end