function Ki = Ki_calc(z_w, D)

data = load('Ki.dat');
par = 2*z_w/D;
if par > 1
    par = 1;
elseif par < -1
    par  = -1;
end
Ki = interp1(data(:,1), data(:,2), par, 'linear');

end