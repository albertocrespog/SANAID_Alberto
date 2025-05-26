%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   08 de Septiembre de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% DATCOM 7.1.2.1 (PDF 2529)
function Cyp_CL_M0 = Cyp_CL_M0_calc(LAMc4, TR)
    
    data1   = load('Cyp_CL_M0_1.dat');
    val     = interp1(data1(:,1), data1(:,2), LAMc4, 'pchip');
    
    data2   = load('Cyp_CL_M0_2.dat');
    data2_1 = data2(1:2,:);
    data2_2 = data2(3:4,:);
    data2_3 = data2(5:6,:);
    
    TR_pos              = [1; 0.5; 0]; 
    Cyp_CL_M0_pos(1)    = interp1(data2_1(:,1), data2_1(:,2), val, 'linear');
    Cyp_CL_M0_pos(2)    = interp1(data2_2(:,1), data2_2(:,2), val, 'linear');
    Cyp_CL_M0_pos(3)    = interp1(data2_3(:,1), data2_3(:,2), val, 'linear');
    
    Cyp_CL_M0           = interp1(TR_pos, Cyp_CL_M0_pos, TR, 'linear');
    
end