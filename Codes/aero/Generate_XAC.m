function [XAC,Fig] = Generate_XAC(DATA_Ae,Aero,DATA_PL,CASOS,degrees_XAC,Fig,...
    XFLR5_DATA,V_Performance,X_OC,Plot_Options,VECTOR_XFLR5)

LS = Plot_Options.LS;
FS = Plot_Options.FS;


Casos_XAC_w1 = CASOS.Casos_XAC_w1;
Casos_XAC_w2 = CASOS.Casos_XAC_w2;

VECTOR_XFLR5 = VECTOR_XFLR5.v1;

rho = V_Performance.rho;
V_CR = V_Performance.V;

MAC_XFLR5 = XFLR5_DATA.MAC_XFLR5;
S_w1_XFLR5 = XFLR5_DATA.S_w1_XFLR5;
XNP_XFLR5 = XFLR5_DATA.XNP_XFLR5;
Q_dyn_pressure = 0.5*V_CR^2*rho;

vector_X = linspace(-0.20,0.20,100);
for i=1:length(VECTOR_XFLR5)
    for j=1:length(degrees_XAC)
        index_w1 = VECTOR_XFLR5(i);
        CL_w1(i,j) = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).CL,degrees_XAC(j),'spline');
        CM_w1(i,j) = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).Cm,degrees_XAC(j),'spline');
        XNP(i,j) = interp1(DATA_Ae(index_w1).alpha,DATA_Ae(index_w1).XCP,degrees_XAC(j),'spline');
        for k=1:length(vector_X)
            My_vec{i}(j,k) = Q_dyn_pressure*S_w1_XFLR5*MAC_XFLR5*CM_w1(i,j) + ...
                Q_dyn_pressure*S_w1_XFLR5*CL_w1(j)*(X_OC(index_w1) - vector_X(k));
            X_vector_plot(i,k) = (X_OC(index_w1) - vector_X(k));
        end
    end
end
XAC = 1;


for i=1:length(VECTOR_XFLR5)
    Fig = Fig + 1;
    figure(Fig)
    for j=1:length(degrees_XAC)
        plot_y = My_vec{i}(j,:);
        plot(X_vector_plot(j,:),plot_y,'LineWidth', LS); hold on
%         leg1{i} = strcat(mark_legend{VECTOR(i)});
    end
    title('M_y vs \alpha','FontSize',FS)
    xlabel('alpha (deg)','FontSize',FS)
    ylabel('C_L','FontSize',FS)
%     legend(leg1,0,'FontSize',LFS)
    hold off
    grid on
end

