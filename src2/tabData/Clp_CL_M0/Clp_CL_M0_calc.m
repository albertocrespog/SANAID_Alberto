%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   09 de Septiembre de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Clp_CL_CL0_M0: Performance, Stability, Dynamics and Control of Airplanes;
% Pamadi, Bandu N.; Fig 4.28, pagina 417

function Clp_CL_M0 = Clp_CL_M0_calc(AR, TR, LAMc4)
    
    data1       = load('Clp_CL_M0_1.dat');
    TR_pos      = [0 0.25 0.5 1];
    
    data1_1     = data1(1:10,:);
    data1_2     = data1(11:20,:);
    data1_3     = data1(21:30,:);
    data1_4     = data1(31:40,:);
    
    val1(1)     = interp1(data1_1(:,1), data1_1(:,2), AR, 'pchip');
    val1(2)     = interp1(data1_2(:,1), data1_2(:,2), AR, 'pchip');
    val1(3)     = interp1(data1_3(:,1), data1_3(:,2), AR, 'pchip');
    val1(4)     = interp1(data1_4(:,1), data1_4(:,2), AR, 'pchip');
        
    par1        = interp1(TR_pos, val1, TR, 'linear');
    
    LAMc4_pos   = [0 15 30 45 60];
    data2       = load('Cyp_CL_M0_2.dat');
    data2_1     = data2(1:2,:);
    data2_2     = data2(3:4,:);
    data2_3     = data2(5:6,:);
    data2_4     = data2(7:8,:);
    data2_5     = data2(9:10,:);
    
    val2(1)     = interp1(data2_1(:,1), data2_1(:,2), par1, 'linear');
    val2(2)     = interp1(data2_2(:,1), data2_2(:,2), par1, 'linear');
    val2(3)     = interp1(data2_3(:,1), data2_3(:,2), par1, 'linear');
    val2(4)     = interp1(data2_4(:,1), data2_4(:,2), par1, 'linear');
    val2(5)     = interp1(data2_5(:,1), data2_5(:,2), par1, 'linear');
    
    Clp_CL_M0   = interp1(LAMc4_pos, val2, LAMc4, 'linear');
    
end