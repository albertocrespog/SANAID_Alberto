function tau = controlEffec_calc(Smovil_Stotal)

data = load('deflec_effec.dat');
if Smovil_Stotal > max(data(:,1))
    Smovil_Stotal = max(data(:,1));
elseif Smovil_Stotal < min(data(:,1))
    Smovil_Stotal = min(data(:,1));
end

tau = interp1(data(:,1),data(:,2),Smovil_Stotal,'pchip');

end