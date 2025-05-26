function [Fig] = Generates_Plots_StabilityAnalysis_lat_VAR(Geo_tier,Plot_Options,Stab_Dyn_LatDir_var_V_XCG,OUTPUT_read_XLSX,Fig,filenameS)

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;

% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filenamec = 'Figs';
% filename_DATA = strcat(filename,filenamed);
% filename_Plots = strcat(filename,filenamec);
fname = filenameS.plots;

% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;
MATLAB_in = Plot_Options.MATLAB_in;

V_VAR = Plot_Options.V_VAR;
m_VAR = Plot_Options.W_VAR;

% Generates the stability analysis elements to plot
for i=1:length(V_VAR)
    for j=1:length(m_VAR)

        Pole1(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.pole1;
        re_p1(i,j)= real(Pole1(i,j));
        im_p1(i,j)= imag(Pole1(i,j));
        Pole2(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.pole2;
        re_p2(i,j)= real(Pole2(i,j));
        im_p2(i,j)= imag(Pole2(i,j));
        Pole3(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.pole3;
        re_p3(i,j)= real(Pole3(i,j));
        im_p3(i,j)= imag(Pole3(i,j));
        Pole4(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.pole4;
        re_p4(i,j)= real(Pole4(i,j));
        im_p4(i,j)= imag(Pole4(i,j));        
        Pole5(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.pole5;
        re_p5(i,j)= real(Pole5(i,j));
        im_p5(i,j)= imag(Pole5(i,j));        

        DR_poles1(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.DR_poles1;
        re_DR_poles1(i,j)= real(DR_poles1(i,j));
        im_DR_poles1(i,j)= imag(DR_poles1(i,j));
        DR_poles2(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.DR_poles2;
        re_DR_poles2(i,j)= real(DR_poles2(i,j));
        im_DR_poles2(i,j)= imag(DR_poles2(i,j));
        Spiral_pole(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.Spiral_pole;
        Roll_pole(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.Roll_pole;
        Yaw_pole(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.Yaw_pole;
        
        DAMP_DR(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.DR_damp;
        WN_DR(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.DR_wn;
        
        DR_T(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.DR_T;
        DR_T2(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.DR_T2;
        cycles_DR(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.cycles_DR;
        decrement_DR(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.decrement_DR;
        cycles_DR(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.cycles_DR;
        decrement_DR(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.decrement_DR;
        ROL_T2(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.ROL_T2;
        ESP_T2(i,j) = Stab_Dyn_LatDir_var_V_XCG{i,j}.lat.ESP_T2;
    end
end

Fig = Fig + 1;
figure(Fig)
subplot(4,3,1)
surf(m_VAR,V_VAR,re_p1)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Real Part Pole 1 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('re_{\lambda_1}')

subplot(4,3,4)
surf(m_VAR,V_VAR,im_p1)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Imaginary Part Pole 1 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('im_{\lambda_1}')

subplot(4,3,2)
surf(m_VAR,V_VAR,re_p2)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Real Part Pole 2 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('re_{\lambda_2}')

subplot(4,3,5)
surf(m_VAR,V_VAR,im_p2)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Imaginary Part Pole 2 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('im_{\lambda_2}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'poles_long_12_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end
% 
% Fig = Fig + 1;
% figure(Fig)
subplot(4,3,3)
surf(m_VAR,V_VAR,re_p3)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Real Part Pole 3 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('re_{\lambda_3}')

subplot(4,3,6)
surf(m_VAR,V_VAR,im_p3)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Imaginary Part Pole 3 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('im_{\lambda_3}')

subplot(4,3,7)
surf(m_VAR,V_VAR,re_p4)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Real Part Pole 4 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('re_{\lambda_4}')

subplot(4,3,10)
surf(m_VAR,V_VAR,im_p4)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Imaginary Part Pole 4 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('im_{\lambda_4}')

subplot(4,3,8)
surf(m_VAR,V_VAR,re_p5)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Real Part Pole 5 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('re_{\lambda_5}')

subplot(4,3,11)
surf(m_VAR,V_VAR,im_p5)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Imaginary Part Pole 5 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('im_{\lambda_5}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('properties_LatDir');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    saveas(gca, fullfile(fname, name), 'jpeg');
    saveas(gcf,fullfile(fname, name),'fig');
    saveas(gcf,fullfile(fname, name),'pdf');
    saveas(gcf,fullfile(fname, name),'bmp');
    saveas(gcf,fullfile(fname, name),'png');
end

Fig = Fig + 1;
figure(Fig)
subplot(2,3,1)
surf(m_VAR,V_VAR,re_DR_poles1)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Real Part Dutch Roll pole1 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('re_{\lambda_{1_{DR}}}')

subplot(2,3,4)
surf(m_VAR,V_VAR,im_DR_poles1)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Imaginary Part Dutch Roll Pole 1 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('im_{\lambda_{1_{DR}}}')

subplot(2,3,2)
surf(m_VAR,V_VAR,re_DR_poles2)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Real Part Dutch Roll Pole 2 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('re_{\lambda_{2_{DR}}}')

subplot(2,3,5)
surf(m_VAR,V_VAR,im_DR_poles2)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Imaginary Part Dutch Roll Pole 2 vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('im_{\lambda_{1_{DR}}}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'poles_DR_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end

subplot(2,3,3)
surf(m_VAR,V_VAR,WN_DR)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Dutch Roll natural frequency vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\omega_{n_{DR}}')

subplot(2,3,6)
surf(m_VAR,V_VAR,DAMP_DR)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Dutch Roll Damping Ratio vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\zeta_{DR}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('omega_chi_DR_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    saveas(gca, fullfile(fname, name), 'jpeg');
    saveas(gcf,fullfile(fname, name),'fig');
    saveas(gcf,fullfile(fname, name),'pdf');
    saveas(gcf,fullfile(fname, name),'bmp');
    saveas(gcf,fullfile(fname, name),'png');
end

Fig = Fig + 1;
figure(Fig)
subplot(3,1,1)
surf(m_VAR,V_VAR,Spiral_pole)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Spiral Pole vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\lambda_{spiral}')

subplot(3,1,2)
surf(m_VAR,V_VAR,Roll_pole)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Roll Pole vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\lambda_{roll}')
        
subplot(3,1,3)
surf(m_VAR,V_VAR,Yaw_pole)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Yaw Pole vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\lambda_{yaw}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('poles_lat2_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    saveas(gca, fullfile(fname, name), 'jpeg');
    saveas(gcf,fullfile(fname, name),'fig');
    saveas(gcf,fullfile(fname, name),'pdf');
    saveas(gcf,fullfile(fname, name),'bmp');
    saveas(gcf,fullfile(fname, name),'png');
end

% Fig = Fig + 1;
% figure(Fig)
% subplot(2,1,1)
% surf(m_VAR,V_VAR,WN_DR)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
% %    rotate3d on % allows you to rotate the graph with the mouse
% colorbar % adds the color bar to the right of the graph
% title('Dutch Roll natural frequency vs V and mass')
% ylabel('Velocity  (m/s)')
% xlabel('mass (kg)')
% zlabel('\omega_{n_{DR}}')
% 
% subplot(2,1,2)
% surf(m_VAR,V_VAR,DAMP_DR)
% shading interp
% colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
% %like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
% %    rotate3d on % allows you to rotate the graph with the mouse
% colorbar % adds the color bar to the right of the graph
% title('Dutch Roll Damping Ratio vs V and mass')
% ylabel('Velocity  (m/s)')
% xlabel('mass (kg)')
% zlabel('\zeta_{DR}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'omega_chi_DR_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end
%
Fig = Fig + 1;
figure(Fig)
subplot(3,2,1)
surf(m_VAR,V_VAR,DR_T)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Dutch Roll period vs V and mass')
ylabel('Velocity (m/s)')
xlabel('mass (kg)')
zlabel('T_{DR}')

subplot(3,2,4)
surf(m_VAR,V_VAR,DR_T2)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Dutch Roll time to half/double vs V and mass')
ylabel('Velocity (m/s)')
xlabel('mass (kg)')
zlabel('t_{1/2_{DR}}')

subplot(3,2,2)
surf(m_VAR,V_VAR,cycles_DR)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Dutch Roll cycles to half/double vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('cycles_{1/2_{DR}}')

subplot(3,2,5)
surf(m_VAR,V_VAR,decrement_DR)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Dutch Roll Logarithmic decrement vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('log_{\delta_{DR}}')

subplot(3,2,3)
surf(m_VAR,V_VAR,ROL_T2)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Roll Perido vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('cycles_{1/2_{DR}}')

subplot(3,2,6)
surf(m_VAR,V_VAR,ESP_T2)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Spiral Period vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('log_{\delta_{DR}}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('properties_DR');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
    saveas(gca, fullfile(fname, name), 'jpeg');
    saveas(gcf,fullfile(fname, name),'fig');
    saveas(gcf,fullfile(fname, name),'pdf');
    saveas(gcf,fullfile(fname, name),'bmp');
    saveas(gcf,fullfile(fname, name),'png');
end