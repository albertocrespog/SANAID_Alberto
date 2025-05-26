%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   09 de Septiembre de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% sigma_beta_wb: DATCOM 7.4.4.4

function sigma_beta_wb = sigma_beta_wb_calc(Zca_b2, LAMc4, TR, W_b)
if W_b < 0.06
    W_b = 0.06;
elseif W_b > 0.24
    W_b = 0.24;
end

% LAMc4 en grados
    LAMc4_pos       = [0 35];
    W_b_pos         = [0.06 0.12 0.24];
    TR_pos          = [0 0.25 0.5 1];
    
    if LAMc4 > LAMc4_pos(2)
        LAMc4 = LAMc4_pos(2);
    elseif LAMc4 < LAMc4_pos(1)
        LAMc4 = LAMc4_pos(1);
    end
    
    data_006        = load('08_006.dat');
    val(1)          = interp1(data_006(:,1), data_006(:,2), Zca_b2, 'pchip');
    
    data_012        = load('08_012.dat');
    val(2)          = interp1(data_012(:,1), data_012(:,2), Zca_b2, 'pchip');
    
    
    data_024_0_0    = load('08_024_0_0.dat');
    val_024_0(1)    = interp1(data_024_0_0(:,1), data_024_0_0(:,2), Zca_b2, 'pchip');
    
    data_024_0_35   = load('08_024_0_35.dat');
    val_024_0(2)    = interp1(data_024_0_35(:,1), data_024_0_35(:,2), Zca_b2, 'pchip');
    
    val_024(1)      = interp1(LAMc4_pos,val_024_0,LAMc4,'pchip');
    
    
    data_024_025_0  = load('08_024_025_0.dat');
    val_024_025(1)  = interp1(data_024_025_0(:,1), data_024_025_0(:,2), Zca_b2, 'pchip');
    
    data_024_025_35 = load('08_024_025_0.dat');
    val_024_025(2)  = interp1(data_024_025_35(:,1), data_024_025_35(:,2), Zca_b2, 'pchip');
    
    val_024(2)      = interp1(LAMc4_pos,val_024_025,LAMc4,'pchip');
    
    
    data_024_05_0   = load('08_024_025_0.dat');
    val_024_05(1)   = interp1(data_024_05_0(:,1), data_024_05_0(:,2), Zca_b2, 'pchip');
    
    data_024_05_35  = load('08_024_025_0.dat');
    val_024_05(2)   = interp1(data_024_05_35(:,1), data_024_05_35(:,2), Zca_b2, 'pchip');
    
    val_024(3)      = interp1(LAMc4_pos,val_024_05,LAMc4,'pchip');
    
    
    data_024_1_0    = load('08_024_025_0.dat');
    val_024_1(1)    = interp1(data_024_1_0(:,1), data_024_1_0(:,2), Zca_b2, 'pchip');
    
    data_024_1_35   = load('08_024_025_0.dat');
    val_024_1(2)    = interp1(data_024_1_35(:,1), data_024_1_35(:,2), Zca_b2, 'pchip');
    
    val_024(4)      = interp1(LAMc4_pos,val_024_1,LAMc4,'linear');
    
    val(3)          = interp1(TR_pos, val_024, TR, 'linear');
    
    sigma_beta_wb   = interp1(W_b_pos, val, W_b, 'linear');
    
    
end