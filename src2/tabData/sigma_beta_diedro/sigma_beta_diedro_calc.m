%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   09 de Septiembre de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% sigma_beta_diedro: DATCOM 7.4.4.4

function sigma_beta_diedro = sigma_beta_diedro_calc(Zca_b2, LAMc4, Minf)
    % LAMc4 en grados
    LAMc4_pos       = [0 35];
    Mpos            = [0.2 0.8];
    data_0_02       = load('02_0.dat');
    data_35_02      = load('02_35.dat');
       
    val_02(1)       = interp1(data_0_02(:,1), data_0_02(:,2), Zca_b2, 'pchip');
    val_02(2)       = interp1(data_35_02(:,1), data_35_02(:,2), Zca_b2, 'pchip');
    
    val(1)          = interp1(LAMc4_pos,val_02,LAMc4,'pchip');
    
    data_0_08       = load('02_0.dat');
    data_35_08      = load('02_35.dat');
       
    val_08(1)       = interp1(data_0_08(:,1), data_0_08(:,2), Zca_b2, 'pchip');
    val_08(2)       = interp1(data_35_08(:,1), data_35_08(:,2), Zca_b2, 'pchip');
        
    val(2)          = interp1(LAMc4_pos,val_08,LAMc4,'pchip');
    
    sigma_beta_diedro = interp1(Mpos,val,Minf,'pchip');   
end