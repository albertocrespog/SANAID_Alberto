function [Fig] = Generates_Plots_Derivatives_VAR(Stab_Der_var_V_XCG,Geo_tier,Plot_Options,OUTPUT_read_XLSX,Fig,filenameS)

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
         
for i=1:length(V_VAR)
    for j=1:length(m_VAR)
        CL(i,j) = Stab_Der_var_V_XCG{i,j}.CL;
        CD(i,j) = Stab_Der_var_V_XCG{i,j}.CD;
        CM(i,j) = Stab_Der_var_V_XCG{i,j}.CM;

        % Alpha dot Derivatives
        CL_alpha_ac(i,j) = Stab_Der_var_V_XCG{i,j}.CL_alpha_ac;
        CD_alpha(i,j) = Stab_Der_var_V_XCG{i,j}.CD_alpha;
        CM_alpha_ac(i,j) = Stab_Der_var_V_XCG{i,j}.CM_alpha_ac;

        % Alpha dot Derivatives
        CLalphapunto(i,j) = Stab_Der_var_V_XCG{i,j}.CL_alphapunto;
        CXalfapunto(i,j) = Stab_Der_var_V_XCG{i,j}.CX_alphapunto;
        CMalphapunto(i,j) = Stab_Der_var_V_XCG{i,j}.CM_alphapunto;

        % Q Derivatives
        CLq(i,j) = Stab_Der_var_V_XCG{i,j}.CLq;
        CXq(i,j) = Stab_Der_var_V_XCG{i,j}.CXq;
        CMq(i,j) = Stab_Der_var_V_XCG{i,j}.CMq;

        % Theta Derivatives
        CZteta(i,j) = Stab_Der_var_V_XCG{i,j}.CZ_teta;
        CXteta(i,j) = Stab_Der_var_V_XCG{i,j}.CX_teta;
        CMteta(i,j) = Stab_Der_var_V_XCG{i,j}.CM_teta;

        % u Derivatives
        CLu(i,j) = Stab_Der_var_V_XCG{i,j}.CLu;
        CDu(i,j) = Stab_Der_var_V_XCG{i,j}.CDu;
        CMu(i,j) = Stab_Der_var_V_XCG{i,j}.CMu;
        
        % Y Derivatives
        Cyb(i,j) = Stab_Der_var_V_XCG{i,j}.Cyb;
        Clb(i,j) = Stab_Der_var_V_XCG{i,j}.Clb;
        Cnb(i,j) = Stab_Der_var_V_XCG{i,j}.Cnb;

        % P Derivatives
        Cyp(i,j) = Stab_Der_var_V_XCG{i,j}.Cyp;
        Clp(i,j) = Stab_Der_var_V_XCG{i,j}.Clp;
        Cnp(i,j) = Stab_Der_var_V_XCG{i,j}.Cnp;
        
        % R Derivatives
        Cyr(i,j) = Stab_Der_var_V_XCG{i,j}.Cyr;
        Clr(i,j) = Stab_Der_var_V_XCG{i,j}.Clr;
        Cnr(i,j) = Stab_Der_var_V_XCG{i,j}.Cnr;

        % Beta dot Derivatives
        Cybpunto(i,j) = Stab_Der_var_V_XCG{i,j}.Cybpunto;
        Clbpunto(i,j) = Stab_Der_var_V_XCG{i,j}.Clbpunto;
        Cnbpunto(i,j) = Stab_Der_var_V_XCG{i,j}.Cnbpunto;

        % delta_e Derivatives
        CL_delta_e(i,j) = Stab_Der_var_V_XCG{i,j}.CL_delta_e;
        CD_delta_e(i,j) = Stab_Der_var_V_XCG{i,j}.CD_delta_e;
        CM_delta_e(i,j) = Stab_Der_var_V_XCG{i,j}.CM_delta_e;

        % delta_a Derivatives
        Cydeltaa(i,j) = Stab_Der_var_V_XCG{i,j}.Cydeltaa;
        Cldeltaa(i,j) = Stab_Der_var_V_XCG{i,j}.Cldeltaa;
        Cndeltaa(i,j) = Stab_Der_var_V_XCG{i,j}.Cndeltaa;
        
        % delta_r Derivatives
        Cydeltar(i,j) = Stab_Der_var_V_XCG{i,j}.Cydeltar;
        Cldeltar(i,j) = Stab_Der_var_V_XCG{i,j}.Cldeltar;
        Cndeltar(i,j) = Stab_Der_var_V_XCG{i,j}.Cndeltar;
           
%         % delta_rv Derivatives
%         Cydeltarv(i,j) = Stab_Der_var_V_XCG{i,j}.Cydeltarv;
%         Cldeltarv(i,j) = Stab_Der_var_V_XCG{i,j}.Cldeltarv;
%         Cndeltarv(i,j) = Stab_Der_var_V_XCG{i,j}.Cndeltarv;

        %  Propulsive Derivatives
        CTx1(i,j) = Stab_Der_var_V_XCG{i,j}.CTx1;
        CMt1(i,j) = Stab_Der_var_V_XCG{i,j}.CMt1;

        CTxu(i,j) = Stab_Der_var_V_XCG{i,j}.CTxu;
        CMtu(i,j) = Stab_Der_var_V_XCG{i,j}.CMtu;
        
        CTxalpha(i,j) = Stab_Der_var_V_XCG{i,j}.CTxalpha;
        CMtalpha(i,j) = Stab_Der_var_V_XCG{i,j}.CMtalpha;
        
        CNTb(i,j) = Stab_Der_var_V_XCG{i,j}.CNTb;
        CyTb(i,j) = Stab_Der_var_V_XCG{i,j}.CyTb;
        
    end
end

Fig = Fig + 1;
figure(Fig)
subplot(2,3,1)
surf(m_VAR,V_VAR,CL)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{L} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{L}')

subplot(2,3,2)
surf(m_VAR,V_VAR,CD)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{D} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{D}')

subplot(2,3,3)
surf(m_VAR,V_VAR,CM)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M}')

% % u Derivatives
% CLu(i,j) = Stab_Der_var_V_XCG{i,j}.CLu;
% CDu(i,j) = Stab_Der_var_V_XCG{i,j}.CDu;
% CMu(i,j) = Stab_Der_var_V_XCG{i,j}.CMu;
% Fig = Fig + 1;
% figure(Fig)
subplot(2,3,4)
surf(m_VAR,V_VAR,CLu)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{L_{u}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{L_{u}}')

subplot(2,3,5)
surf(m_VAR,V_VAR,CDu)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{D_u} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{D_u}')

subplot(2,3,6)
surf(m_VAR,V_VAR,CMu)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M_u} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M_u}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('u_der');
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
subplot(2,3,1)
surf(m_VAR,V_VAR,CL_alpha_ac)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{L_{\alpha}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{L_{\alpha}}')

subplot(2,3,2)
surf(m_VAR,V_VAR,CD_alpha)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{D_\alpha} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{D_\alpha}')

subplot(2,3,3)
surf(m_VAR,V_VAR,CM_alpha_ac)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M_\alpha} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M_\alpha}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'alpha_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end
%
% % Alpha dot Derivatives
% CLalphapunto(i,j) = Stab_Der_var_V_XCG{i,j}.CLalphapunto;
% CXalfapunto(i,j) = Stab_Der_var_V_XCG{i,j}.CXalfapunto;
% CMalphapunto(i,j) = Stab_Der_var_V_XCG{i,j}.CMalphapunto;
% Fig = Fig + 1;
% figure(Fig)
subplot(2,3,4)
surf(m_VAR,V_VAR,CLalphapunto)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{L_{\alpha''}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{L_{\alpha''}}')

subplot(2,3,5)
surf(m_VAR,V_VAR,CXalfapunto)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{D_\alpha''} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{D_\alpha''}')

subplot(2,3,6)
surf(m_VAR,V_VAR,CMalphapunto)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M_\alpha''} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M_\alpha''}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('alpha_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end


% % Q Derivatives
% CLq(i,j) = Stab_Der_var_V_XCG{i,j}.CLq;
% CXq(i,j) = Stab_Der_var_V_XCG{i,j}.CXq;
% CMq(i,j) = Stab_Der_var_V_XCG{i,j}.CMq;
Fig = Fig + 1;
figure(Fig)
subplot(2,3,1)
surf(m_VAR,V_VAR,CLq)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{L_{q}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{L_{q}}')

subplot(2,3,2)
surf(m_VAR,V_VAR,CXq)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{D_q} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{D_q}')

subplot(2,3,3)
surf(m_VAR,V_VAR,CMq)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M_q} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M_q}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'q_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end
% % Theta Derivatives
% CZteta(i,j) = Stab_Der_var_V_XCG{i,j}.CZteta;
% CXteta(i,j) = Stab_Der_var_V_XCG{i,j}.CXteta;
% CMteta(i,j) = Stab_Der_var_V_XCG{i,j}.CMteta;
% Fig = Fig + 1;
% figure(Fig)
subplot(2,3,4)
surf(m_VAR,V_VAR,CZteta)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{L_{\theta}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{L_{\theta}}')

subplot(2,3,5)
surf(m_VAR,V_VAR,CXteta)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{D_\theta} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{D_\theta}')

subplot(2,3,6)
surf(m_VAR,V_VAR,CMteta)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M_\theta} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M_\theta}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('qtheta_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

% Y Derivatives
% Cyb(i,j) = Stab_Der_var_V_XCG{i,j}.Cyb;
% Clb(i,j) = Stab_Der_var_V_XCG{i,j}.Clb;
% Cnb(i,j) = Stab_Der_var_V_XCG{i,j}.Cnb;
Fig = Fig + 1;
figure(Fig)
subplot(2,3,1)
surf(m_VAR,V_VAR,Cyb)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{Y_{\beta}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{Y_{\beta}}')

subplot(2,3,2)
surf(m_VAR,V_VAR,Clb)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{l_\beta} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{l_\beta}')

subplot(2,3,3)
surf(m_VAR,V_VAR,Cnb)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{N_\beta} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{N_\beta''}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'beta_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end

% Beta dot Derivatives
% Cybpunto(i,j) = Stab_Der_var_V_XCG{i,j}.Cybpunto;
% Clbpunto(i,j) = Stab_Der_var_V_XCG{i,j}.Clbpunto;
% Cnbpunto(i,j) = Stab_Der_var_V_XCG{i,j}.Cnbpunto;
% Fig = Fig + 1;
% figure(Fig)
subplot(2,3,4)
surf(m_VAR,V_VAR,Cybpunto)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{Y_{\beta''}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{Y_{\beta''}}')

subplot(2,3,5)
surf(m_VAR,V_VAR,Clbpunto)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{l_{\beta''}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{l_{\beta''}}')

subplot(2,3,6)
surf(m_VAR,V_VAR,Cnbpunto)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{N_{\beta''}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{N_{\beta''}}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('beta_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

% % P Derivatives
% Cyp(i,j) = Stab_Der_var_V_XCG{i,j}.Cyp;
% Clp(i,j) = Stab_Der_var_V_XCG{i,j}.Clp;
% Cnp(i,j) = Stab_Der_var_V_XCG{i,j}.Cnp;
Fig = Fig + 1;
figure(Fig)
subplot(2,3,1)
surf(m_VAR,V_VAR,Cyp)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{Y_{p}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{Y_{p}}')

subplot(2,3,2)
surf(m_VAR,V_VAR,Clp)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{l_p} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{l_p}')

subplot(2,3,3)
surf(m_VAR,V_VAR,Cnp)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{N_p} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{n_p}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'p_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end

% % R Derivatives
% Cyr(i,j) = Stab_Der_var_V_XCG{i,j}.Cyr;
% Clr(i,j) = Stab_Der_var_V_XCG{i,j}.Clr;
% Cnr(i,j) = Stab_Der_var_V_XCG{i,j}.Cnr;
% Fig = Fig + 1;
% figure(Fig)
subplot(2,3,4)
surf(m_VAR,V_VAR,Cyr)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{Y_{r}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{Y_{r}}')

subplot(2,3,5)
surf(m_VAR,V_VAR,Clr)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{l_r} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{l_r}')

subplot(2,3,6)
surf(m_VAR,V_VAR,Cnr)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{N_r} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{n_r}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('pr_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

% % delta_e Derivatives
% CL_delta_e(i,j) = Stab_Der_var_V_XCG{i,j}.CL_delta_e;
% CD_delta_e(i,j) = Stab_Der_var_V_XCG{i,j}.CD_delta_e;
% CM_delta_e(i,j) = Stab_Der_var_V_XCG{i,j}.CM_delta_e;
Fig = Fig + 1;
figure(Fig)
subplot(3,3,1)
surf(m_VAR,V_VAR,CL_delta_e)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{L_{\delta_e}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{L_{\delta_e}}')

subplot(3,3,2)
surf(m_VAR,V_VAR,CD_delta_e)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{D_{\delta_e}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{D_{\delta_e}}')

subplot(3,3,3)
surf(m_VAR,V_VAR,CM_delta_e)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M_{\delta_e}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M_{\delta_e}}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'delta_e_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end
% 
% % delta_a Derivatives
% Cydeltaa(i,j) = Stab_Der_var_V_XCG{i,j}.Cydeltaa;
% Cldeltaa(i,j) = Stab_Der_var_V_XCG{i,j}.Cldeltaa;
% Cndeltaa(i,j) = Stab_Der_var_V_XCG{i,j}.Cndeltaa;
% Fig = Fig + 1;
% figure(Fig)
subplot(3,3,4)
surf(m_VAR,V_VAR,Cydeltaa)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{Y_{\delta_a}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{Y_{\delta_a}}')

subplot(3,3,5)
surf(m_VAR,V_VAR,Cldeltaa)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{l_{\delta_a}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{l_{\delta_a}}')

subplot(3,3,6)
surf(m_VAR,V_VAR,Cndeltaa)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{N_{\delta_a}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{N_{\delta_a}}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'delta_e_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end
% 
% % delta_r Derivatives
% Cydeltar(i,j) = Stab_Der_var_V_XCG{i,j}.Cydeltar;
% Cldeltar(i,j) = Stab_Der_var_V_XCG{i,j}.Cldeltar;
% Cndeltar(i,j) = Stab_Der_var_V_XCG{i,j}.Cndeltar;
% Fig = Fig + 1;
% figure(Fig)
subplot(3,3,7)
surf(m_VAR,V_VAR,Cydeltar)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{Y_{\delta_r}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{Y_{\delta_r}}')

subplot(3,3,8)
surf(m_VAR,V_VAR,Cldeltar)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{l_{\delta_r}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{l_{\delta_r}}')

subplot(3,3,9)
surf(m_VAR,V_VAR,Cndeltar)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{N_{\delta_r}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{N_{\delta_r}}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('delta_ear_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

% 
% % delta_rv Derivatives
% Cydeltarv(i,j) = Stab_Der_var_V_XCG{i,j}.Cydeltarv;
% Cldeltarv(i,j) = Stab_Der_var_V_XCG{i,j}.Cldeltarv;
% Cndeltarv(i,j) = Stab_Der_var_V_XCG{i,j}.Cndeltarv;
% 
%         %  Propulsive Derivatives
%         CTx1(i,j) = Stab_Der_var_V_XCG{i,j}.CTx1;
%         CMt1(i,j) = Stab_Der_var_V_XCG{i,j}.CMt1;
% 
%         CTxu(i,j) = Stab_Der_var_V_XCG{i,j}.CTxu;
%         CMtu(i,j) = Stab_Der_var_V_XCG{i,j}.CMtu;
%         
%         CTxalpha(i,j) = Stab_Der_var_V_XCG{i,j}.CTxalpha;
%         CMtalpha(i,j) = Stab_Der_var_V_XCG{i,j}.CMtalpha;
%         
%         CNTb(i,j) = Stab_Der_var_V_XCG{i,j}.CNTb;
%         CyTb(i,j) = Stab_Der_var_V_XCG{i,j}.CyTb;
Fig = Fig + 1;
figure(Fig)
subplot(3,3,1)
surf(m_VAR,V_VAR,CTx1)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{T_{x_1}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{T_{x_1}}')

subplot(3,3,4)
surf(m_VAR,V_VAR,CMt1)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M_{T_1}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M_{T_1}}')

subplot(3,3,2)
surf(m_VAR,V_VAR,CTxu)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{T_{x_u}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{T_{x_u}}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'delta_e_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end

% CMtu(i,j) = Stab_Der_var_V_XCG{i,j}.CMtu;
%         
%         CTxalpha(i,j) = Stab_Der_var_V_XCG{i,j}.CTxalpha;
%         CMtalpha(i,j) = Stab_Der_var_V_XCG{i,j}.CMtalpha;
% Fig = Fig + 1;
% figure(Fig)
subplot(3,3,5)
surf(m_VAR,V_VAR,CMtu)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M_{T_u}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M_{T_u}}')

subplot(3,3,3)
surf(m_VAR,V_VAR,CTxalpha)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{T_{x_\alpha}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{T_{x_\alpha}}')

subplot(3,3,6)
surf(m_VAR,V_VAR,CMtalpha)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{M_{T_\alpha}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{M_{T_\alpha}}')
% if SAVE_FIGS==1
%     name   = strcat(prefix,'delta_e_der');
%     saveas(gcf,name,'fig');
%     saveas(gca,name,'epsc');
%     saveas(gcf,name,'pdf');
%     saveas(gcf,name,'bmp');
%     saveas(gcf,name,'png');
% end

%         CNTb(i,j) = Stab_Der_var_V_XCG{i,j}.CNTb;
%         CyTb(i,j) = Stab_Der_var_V_XCG{i,j}.CyTb;
subplot(3,3,7)
surf(m_VAR,V_VAR,CyTb)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{Y_{T_\beta}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{Y_{T_\beta}}')

subplot(3,3,8)
surf(m_VAR,V_VAR,CNTb)
shading interp
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types
%like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
%    rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('C_{N_{T_\beta}} vs V and mass')
ylabel('Velocity  (km/h)')
xlabel('mass (kg)')
zlabel('C_{N_{T_\beta}}')
if SAVE_FIGS==1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('prop_der');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

