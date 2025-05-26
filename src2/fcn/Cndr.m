function Cn_dr = Cndr(varargin)
switch length(varargin)
    case 7 % Calculo preliminar
        CLav    = varargin(1);
        Sv      = varargin(2);
        Srudd   = varargin(3);
        Sref    = varargin(4);
        b       = varargin(5);
        lv      = varargin(6);
        rend    = varargin(7);
        tau     = controlEffec_calc(Sv, Srudd);
        Cn_dr   = -CLav*rend*tau*Sv/Sref*lv/b;
    case 2

end
end