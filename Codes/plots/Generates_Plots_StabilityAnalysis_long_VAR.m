function [Fig] = Generates_Plots_StabilityAnalysis_long_VAR(Geo_tier,Plot_Options,Stab_Dyn_Long_var_V_XCG,OUTPUT_read_XLSX,Fig,filenameS)

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
        Pole1(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.pole1;
        re_p1(i,j)= real(Pole1(i,j));
        im_p1(i,j)= imag(Pole1(i,j));
        Pole2(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.pole2;
        re_p2(i,j)= real(Pole2(i,j));
        im_p2(i,j)= imag(Pole2(i,j));
        Pole3(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.pole3;
        re_p3(i,j)= real(Pole3(i,j));
        im_p3(i,j)= imag(Pole3(i,j));
        Pole4(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.pole4;
        re_p4(i,j)= real(Pole4(i,j));
        im_p4(i,j)= imag(Pole4(i,j));        
                
        % Poles
        SP_poles(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.SP_poles;
        PH_poles(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.PH_poles;
        % Short Period Poles
        re_sp(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.re_sp;
        im_sp(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.im_sp;
        % Phugoid Poles
        re_ph(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.re_ph;
        im_ph(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.im_ph;

        WN_SP(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.SP_wn;
        DAMP_SP(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.SP_damp;
        WN_PH(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.PH_wn;
        DAMP_PH(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.PH_damp;        

        PH_T(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.PH_T;
        PH_T2(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.SP_T2;
        cycles_PH(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.cycles_PH;
        decrement_PH(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.decrement_PH;

        SP_T(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.SP_T;
        SP_T2(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.SP_T2;
        cycles_SP(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.cycles_SP;
        decrement_SP(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.decrement_SP;

        % Approximations
        WN_SP_approx(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.WN_SP_approx;
        DAMP_SP_approx(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.DAMP_SP_approx;
        re_sp_approx(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.re_sp_approx;
        im_sp_approx(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.im_sp_approx;
        % Approximations
        WN_PH_approx(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.WN_PH_approx;
        DAMP_PH_approx(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.DAMP_PH_approx;
        re_ph_approx(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.re_ph_approx;
        im_ph_approx(i,j) = Stab_Dyn_Long_var_V_XCG{i,j}.long.im_ph_approx;
        
    end
end

Fig = Fig + 1;
figure(Fig)
subplot(2,2,1)
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

subplot(2,2,2)
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

subplot(2,2,3)
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

subplot(2,2,4)
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
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('poles_long_12_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
subplot(2,2,1)
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

subplot(2,2,2)
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

subplot(2,2,3)
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

subplot(2,2,4)
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
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('poles_long_34_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
subplot(2,2,1)
surf(m_VAR,V_VAR,re_ph)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Real Part Phugoid Mode vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('re_{\lambda_{ph}}')

subplot(2,2,2)
surf(m_VAR,V_VAR,im_ph)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Imaginary Part Phugoid vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('im_{\lambda_{ph}}')

subplot(2,2,3)
surf(m_VAR,V_VAR,re_sp)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Real Part Short Period Mode vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('re_{\lambda_{sp}}')

subplot(2,2,4)
surf(m_VAR,V_VAR,im_sp)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Imaginary Part Short Period vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('im_{\lambda_{sp}}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('poles_sp_ph_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
subplot(2,2,1)
surf(m_VAR,V_VAR,WN_SP)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Short Period natural frequency vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\omega_{n_{SP}}')

subplot(2,2,2)
surf(m_VAR,V_VAR,DAMP_SP)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Short Period Damping Ratio vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\zeta_{SP}')

subplot(2,2,3)
surf(m_VAR,V_VAR,WN_PH)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Phugoid natural frequency vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\omega_{n_{PH}}')

subplot(2,2,4)
surf(m_VAR,V_VAR,DAMP_PH)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Phugoid vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\zeta_{PH}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('omega_chi_SPPH_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
subplot(2,2,1)
surf(m_VAR,V_VAR,WN_SP_approx)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Short Period approximation natural frequency vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\omega_{n_{SP_{approx}}}')

subplot(2,2,2)
surf(m_VAR,V_VAR,DAMP_SP_approx)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Short Period Damping Ratio approximation vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\zeta_{SP_{approx}}')

subplot(2,2,3)
surf(m_VAR,V_VAR,WN_PH_approx)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Phugoid natural frequency approximation vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\omega_{n_{PH_{approx}}}')

subplot(2,2,4)
surf(m_VAR,V_VAR,DAMP_PH_approx)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Phugois Damping Ratio approximation vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('\zeta_{PH_{approx}}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('omega_chi_SPPH_approx');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
subplot(2,2,1)
surf(m_VAR,V_VAR,SP_T)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Short Period period vs V and mass')
ylabel('Velocity (m/s)')
xlabel('mass (kg)')
zlabel('T_{SP}')

subplot(2,2,2)
surf(m_VAR,V_VAR,SP_T2)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Short Period time to half/double vs V and mass')
ylabel('Velocity (m/s)')
xlabel('mass (kg)')
zlabel('t_{1/2_{SP}}')

subplot(2,2,3)
surf(m_VAR,V_VAR,cycles_SP)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Short Period cycles to half/double vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('cycles_{1/2_{SP}}')

subplot(2,2,4)
surf(m_VAR,V_VAR,decrement_SP)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Short Period Logarithmic decrement vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('log_{\delta_{SP}}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('properties_SP');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

Fig = Fig + 1;
figure(Fig)
subplot(2,2,1)
surf(m_VAR,V_VAR,PH_T)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Phugoid period vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('T_{PH}')

subplot(2,2,2)
surf(m_VAR,V_VAR,PH_T2)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Phugoid time to half/double vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('t_{1/2_{PH}}')

subplot(2,2,3)
surf(m_VAR,V_VAR,cycles_PH)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Phugoid cycles to half/double vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('cycles_{1/2_{PH}}')

subplot(2,2,4)
surf(m_VAR,V_VAR,decrement_PH)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('Phugoid Logarithmic decrement vs V and mass')
ylabel('Velocity  (m/s)')
xlabel('mass (kg)')
zlabel('log_{\delta_{PH}}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('properties_PH');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end
