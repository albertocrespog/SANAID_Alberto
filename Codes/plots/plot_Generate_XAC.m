function [XAC,Fig] = plot_Generate_XAC(DATA_Ae,Aero,DATA_PL,CASOS,degrees_XAC,Fig,...
    XFLR5_DATA,V_Performance,X_OC,Plot_Options,VECTOR_XFLR5)

LS = Plot_Options.LS;
FS = Plot_Options.FS;
MATLAB_in = Plot_Options.MATLAB_in;
    
Casos_XAC_w1 = CASOS.Casos_XAC_w1;
Casos_XAC_w2 = CASOS.Casos_XAC_w2;

compare = VECTOR_XFLR5.compare;

switch compare
    case 1 % compare = 1 - w1 wing alone
        % W1 cases to be plotted
        VEC{1} = VECTOR_XFLR5.v1;
    case 2 % compare = 1 - w2 wing alone
        VEC{1} = VECTOR_XFLR5.v1;
        VEC{2} = VECTOR_XFLR5.v2;
    case 3 % compare = 3 - vtp wing alone
        VEC{1} = VECTOR_XFLR5.v1;
        VEC{2} = VECTOR_XFLR5.v2;
        VEC{3} = VECTOR_XFLR5.v3;
end

% Solo para una comparativa
VECTOR_XFLR5 = VECTOR_XFLR5.v1;;

rho = V_Performance.rho;
V = V_Performance.V;

% S_w1_e = Geo_tier.S_w1_e;
% cmac_w1 = Geo_tier.cmac_w1;

MAC_XFLR5 = XFLR5_DATA.MAC_XFLR5;
S_w1_XFLR5 = XFLR5_DATA.S_w1_XFLR5;
XNP_XFLR5 = XFLR5_DATA.XNP_XFLR5;
Q_dyn_pressure = 0.5*V^2*rho;

vector_X = linspace(-0.05,0.05,100);
% for i=1:length(VECTOR_XFLR5)
%     for j=1:length(degrees_XAC)
%         index_w1 = VECTOR_XFLR5(i);
%         CL_w1(i,j) = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CL,degrees_XAC(j),'spline');
%         CM_w1(i,j) = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).Cm,degrees_XAC(j),'spline');
%         XNP(i,j) = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).XCP,degrees_XAC(j),'spline');
%         for k=1:length(vector_X)
%             My_vec{i}(j,k) = Q_dyn_pressure*S_w1_XFLR5*MAC_XFLR5*CM_w1(i,j) + ...
%                 Q_dyn_pressure*S_w1_XFLR5*CL_w1(j)*(X_OC(index_w1) - vector_X(k));
%             X_vector_plot(i,k) = (X_OC(index_w1) - vector_X(k));
%         end
%     end
% end

for i=1:length(VECTOR_XFLR5)
    for j=1:length(degrees_XAC)
        index_w1 = VECTOR_XFLR5(i);
        CL_w1(j) = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CL,degrees_XAC(j),'spline');
        CM_w1(j) = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).Cm,degrees_XAC(j),'spline');
        XNP(j) = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).XCP,degrees_XAC(j),'spline');
        for k=1:length(vector_X)
            My_vec{i}(j,k) = Q_dyn_pressure*S_w1_XFLR5*MAC_XFLR5*CM_w1(j) + ...
                Q_dyn_pressure*S_w1_XFLR5*CL_w1(j)*(X_OC(index_w1) - vector_X(k));
            X_vector_plot(i,k) = (X_OC(index_w1) - vector_X(k));
        end
    end
end


% Plot of figures
Fig = Fig + 1;
figure(Fig)
surf(vector_X,degrees_XAC,My_vec{1})
shading interp 
colormap(jet) % colors paletts for the 3-D graph if you type help graph3d will give you a list of different color types 
              % like hot, or graybone, copper, pink, jet (the one used here) spring, winter, summer, autumn etc...
rotate3d on % allows you to rotate the graph with the mouse
colorbar % adds the color bar to the right of the graph
title('M_y')
xlabel('X (m)')
ylabel('\alpha (deg)')
zlabel('M_y')

for i=1:length(VECTOR_XFLR5)
    Fig = Fig + 1;
    figure(Fig)
    for j=1:length(degrees_XAC)
        plot_y = My_vec{i}(j,:);
        plot(X_vector_plot(i,:),plot_y,'LineWidth', LS); hold on
%         leg1{i} = strcat(mark_legend{VECTOR(i)});
    end
    title('M_y vs \alpha','FontSize',FS)
    xlabel('X (deg)','FontSize',FS)
    ylabel('C_L','FontSize',FS)
%     legend(leg1,0,'FontSize',LFS)
    hold off
    grid on
end
temp=1;
XAC = temp;