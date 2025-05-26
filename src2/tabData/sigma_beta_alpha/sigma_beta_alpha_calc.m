%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   09 de Septiembre de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% sigma_beta_alpha: DATCOM 7.4.4.4

function sigma_beta_alpha = sigma_beta_alpha_calc(Zca_b2, LAMc4, TR)
    % LAMc4 en grados
    LAMc4_pos       = [0 35];
    TR_pos          = [0 0.25 0.5 1];
    
    data_0_0        = load('08_0_0.dat');
    data_0_35       = load('08_0_35.dat');
    
    val_0(1)        = interp1(data_0_0(:,1), data_0_0(:,2), Zca_b2, 'pchip');
    val_0(2)        = interp1(data_0_35(:,1), data_0_35(:,2), Zca_b2, 'pchip');
    val(1)          = interp1(LAMc4_pos,val_0,LAMc4,'pchip');
    
    data_025_0      = load('08_025_0.dat');
    data_025_35     = load('08_025_35.dat');
    
    val_025(1)      = interp1(data_025_0(:,1), data_025_0(:,2), Zca_b2, 'pchip');
    val_025(2)      = interp1(data_025_35(:,1), data_025_35(:,2), Zca_b2, 'pchip');
    val(2)          = interp1(LAMc4_pos,val_025,LAMc4,'pchip');
    
    data_05_0       = load('08_05_0.dat');
    data_05_35      = load('08_05_35.dat');
    
    val_05(1)       = interp1(data_05_0(:,1), data_05_0(:,2), Zca_b2, 'pchip');
    val_05(2)       = interp1(data_05_35(:,1), data_05_35(:,2), Zca_b2, 'pchip');
    val(3)          = interp1(LAMc4_pos,val_05,LAMc4,'pchip');
    
    data_1_0        = load('08_05_0.dat');
    data_1_35       = load('08_05_35.dat');
    
    val_1(1)        = interp1(data_1_0(:,1), data_1_0(:,2), Zca_b2, 'pchip');
    val_1(2)        = interp1(data_1_35(:,1), data_1_35(:,2), Zca_b2, 'pchip');
    val(4)          = interp1(LAMc4_pos,val_1,LAMc4,'pchip');
       
    sigma_beta_alpha = interp1(TR_pos,val,TR,'linear'); 
      
end