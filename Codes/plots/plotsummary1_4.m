function plotsummary1_4(RESULTS,POST,N_dimension,Engine,SAVE_FIG,fig)

close all

Fig = 0;
lv = N_dimension.lv;

WTO_est = RESULTS.WTO_est;

N_W_pl = N_dimension.N_W_pl;
N_s = N_dimension.N_s;
N_Vh = N_dimension.N_Vh;

W_pl = N_dimension.W_pl;
s = N_dimension.s;
Vh = N_dimension.Vh;

PP = N_dimension.PP;
SS = N_dimension.SS;
VH = N_dimension.VH;

Color=rand(15,3);

% Uses the VTOL prop for dimensioning the margins of the plots
D_prop_vec = POST.D_prop_vec_VTOL;

x_HTP_LE = RESULTS.x_HTP_LE;
x_w_LE_w1 = RESULTS.x_w_LE_w1;
x_VTP_LE = RESULTS.x_VTP_LE;
L_fus = RESULTS.L_fus_total;

% k=1;
% for k=1:N_Vh
%     Y_AXIS(k) = max(max(RESULTS.b_w1(:,:,k)));
%     X_AXIS(k) = max(max(RESULTS.L_fus_total(:,:,k)));
% end

for k=1:N_Vh
    n=1;
    figure(k)
    for i=1:N_s
        for j=1:N_W_pl
            PLOTTING_UAV = plot_Models_VTOL_ITER(RESULTS,POST,Engine,i,j,k);
            x_loc_w1(i,j,k) = PLOTTING_UAV.x_loc(1);
            x_loc_HTP(i,j,k) = PLOTTING_UAV.x_loc(2);
            x_loc_VTP(i,j,k) = PLOTTING_UAV.x_loc(3);
            
            y_loc_w1(i,j,k) = PLOTTING_UAV.Y_vec(1);
            y_loc_HTP(i,j,k) = PLOTTING_UAV.Y_vec(2);
            y_loc_VTP(i,j,k) = PLOTTING_UAV.Y_vec(3);
            
            z_loc_w1(i,j,k) = PLOTTING_UAV.Z_vec(1);
            z_loc_HTP(i,j,k) = PLOTTING_UAV.Z_vec(2);
            z_loc_VTP(i,j,k) = PLOTTING_UAV.Z_vec(3);
            
            z_max_fus(i,j,k) = PLOTTING_UAV.z_max_fus;
            z_min_fus(i,j,k) = PLOTTING_UAV.z_min_fus;
        end
    end
end

Delta_Plot = 0.05; % Margin to all sizes of plots
MAX_X = L_fus + Delta_Plot;
MAX_Y = RESULTS.b_w1 + 2*D_prop_vec + Delta_Plot + y_loc_w1/2;
MAX_Z = z_max_fus + RESULTS.b_VTP + Delta_Plot;

for k=1:N_Vh
    X_AXIS(k) = max(max(MAX_X(:,:,k)));
    Y_AXIS(k) = max(max(MAX_Y(:,:,k)));
    Z_AXIS(k) = max(max(MAX_Z(:,:,k)));
end

PLOT_DISK = PLOTTING_UAV.PLOT_DISK; % DATA Disk 

FS = 8; % fontsize
% for k=1:N_Vh
%     n=1;
%     figure(k)
%     for i=1:N_s
%         for j=1:N_W_pl
%             subplot(N_s,N_W_pl,n);
%             PLOTTING_UAV = plot_Models_VTOL_ITER(RESULTS,POST,Engine,i,j,k);
%             meshData = PLOTTING_UAV.meshData; % FUSELAJE
%             x_mesh_w_New = PLOTTING_UAV.x_mesh_w_New; % WING
%             x_mesh_HTP_New = PLOTTING_UAV.x_mesh_HTP_New; % HTP
%             x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
%             x_mesh_VTP2_New = PLOTTING_UAV.x_mesh_VTP2_New; % VTP 2
%             
%             y_mesh_w_New = PLOTTING_UAV.y_mesh_w_New; % WING
%             y_mesh_HTP_New = PLOTTING_UAV.y_mesh_HTP_New; % HTP
%             y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New; % VTP 1
%             y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New; % VTP 2
%             
%             z_mesh_w_New = PLOTTING_UAV.z_mesh_w_New; % WING
%             z_mesh_HTP_New = PLOTTING_UAV.z_mesh_HTP_New; % HTP
%             z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New; % VTP 1
%             z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New; % VTP 2
%             
%             C(:,:,1) = 0.2*ones(size(meshData{1}));
%             C(:,:,2) = 0.3*ones(size(meshData{1}));
%             C(:,:,3) = 0.7*ones(size(meshData{1}));
%             mesh(meshData{1},meshData{2},meshData{3},C);
%             st = strcat('UAV - R ',num2str(round(s(i)/1000)),' km, ',num2str(round(W_pl(j))),' kg, ',num2str(round(Vh(k))),' m/s, W_{TO}=',num2str(round(WTO_est(i,j,k))),'kg');
%             title(st,'fontsize',FS)
%             xlabel('y (m)')
%             ylabel('z (m)')
%             zlabel('z (m)')
%             hold on
%             mesh(x_mesh_w_New,y_mesh_w_New,z_mesh_w_New)
%             mesh(x_mesh_HTP_New,y_mesh_HTP_New,z_mesh_HTP_New)
%             mesh(x_mesh_VTP1_New,y_mesh_VTP1_New,z_mesh_VTP1_New)
%             mesh(x_mesh_VTP2_New,y_mesh_VTP2_New,z_mesh_VTP2_New)
%             hold off
%             axis equal
%             grid on
%             hold on
%             axis([0 X_AXIS(k) -Y_AXIS(k)/2 Y_AXIS(k)/2]);
%             n=n+1;
%         end
%     end
%     if SAVE_FIG==1
%         prefix = strcat('SingleWing_T_D_');
%         st = strcat(num2str(round(Vh(k),3,'significant')),'_m_s');
%         name   = strcat(prefix,st);
%         saveas(gcf,name,'fig');
%         %saveas(gcf,name,'pdf');
%         %saveas(gcf,name,'bmp');
%         saveas(gcf,name,'png');
%     end
%     pause
% end

FS = 8; % fontsize
for k=1:N_Vh
    n=1;
    Fig = Fig +1;
    figure(Fig)
    for i=1:N_s
        for j=1:N_W_pl
            subplot(N_s,N_W_pl,n);
            PLOTTING_UAV = plot_Models_VTOL_ITER(RESULTS,POST,Engine,i,j,k);
            
            % Engine
            Eng_L_x_w1 = PLOTTING_UAV.Eng_L_x_w1;
            Eng_L_y_w1 = PLOTTING_UAV.Eng_L_y_w1;
            Eng_L_z_w1 = PLOTTING_UAV.Eng_L_z_w1;
            Eng_L_x_w1_s = PLOTTING_UAV.Eng_L_x_w1_s;
            Eng_L_y_w1_s = PLOTTING_UAV.Eng_L_y_w1_s;
            Eng_L_z_w1_s = PLOTTING_UAV.Eng_L_z_w1_s;
            
            Eng_R_x_w1 = PLOTTING_UAV.Eng_R_x_w1;
            Eng_R_y_w1 = PLOTTING_UAV.Eng_R_y_w1;
            Eng_R_z_w1 = PLOTTING_UAV.Eng_R_z_w1;
            Eng_R_x_w1_s = PLOTTING_UAV.Eng_R_x_w1_s;
            Eng_R_y_w1_s = PLOTTING_UAV.Eng_R_y_w1_s;
            Eng_R_z_w1_s = PLOTTING_UAV.Eng_R_z_w1_s;

            % Prop
            Prop_L_x_w1 = PLOTTING_UAV.Prop_L_x_w1;
            Prop_L_y_w1 = PLOTTING_UAV.Prop_L_y_w1;
            Prop_L_z_w1 = PLOTTING_UAV.Prop_L_z_w1;
            Prop_L_x_w1_s = PLOTTING_UAV.Prop_L_x_w1_s;
            Prop_L_y_w1_s = PLOTTING_UAV.Prop_L_y_w1_s;
            Prop_L_z_w1_s = PLOTTING_UAV.Prop_L_z_w1_s;
            
            Prop_R_x_w1 = PLOTTING_UAV.Prop_R_x_w1;
            Prop_R_y_w1 = PLOTTING_UAV.Prop_R_y_w1;
            Prop_R_z_w1 = PLOTTING_UAV.Prop_R_z_w1;
            Prop_R_x_w1_s = PLOTTING_UAV.Prop_R_x_w1_s;
            Prop_R_y_w1_s = PLOTTING_UAV.Prop_R_y_w1_s;
            Prop_R_z_w1_s = PLOTTING_UAV.Prop_R_z_w1_s;

            meshData = PLOTTING_UAV.meshData; % FUSELAJE
            x_mesh_w_New = PLOTTING_UAV.x_mesh_w_New; % WING
            x_mesh_HTP_New = PLOTTING_UAV.x_mesh_HTP_New; % HTP
            x_mesh_VTP1_New = PLOTTING_UAV.x_mesh_VTP1_New; % VTP 1
            x_mesh_VTP2_New = PLOTTING_UAV.x_mesh_VTP2_New; % VTP 2
            
            y_mesh_w_New = PLOTTING_UAV.y_mesh_w_New; % WING
            y_mesh_HTP_New = PLOTTING_UAV.y_mesh_HTP_New; % HTP
            y_mesh_VTP1_New = PLOTTING_UAV.y_mesh_VTP1_New; % VTP 1
            y_mesh_VTP2_New = PLOTTING_UAV.y_mesh_VTP2_New; % VTP 2
            
            z_mesh_w_New = PLOTTING_UAV.z_mesh_w_New; % WING
            z_mesh_HTP_New = PLOTTING_UAV.z_mesh_HTP_New; % HTP
            z_mesh_VTP1_New = PLOTTING_UAV.z_mesh_VTP1_New; % VTP 1
            z_mesh_VTP2_New = PLOTTING_UAV.z_mesh_VTP2_New; % VTP 2
            
            C(:,:,1) = 0.2*ones(size(meshData{1}));
            C(:,:,2) = 0.3*ones(size(meshData{1}));
            C(:,:,3) = 0.7*ones(size(meshData{1}));
            mesh(meshData{1},meshData{2},meshData{3},C);
            st = strcat('UAV - R ',num2str(round(s(i)/1000)),' km, ',num2str(round(W_pl(j))),' kg, ',num2str(round(Vh(k))),' m/s, W_{TO}=',num2str(round(WTO_est(i,j,k))),'kg');
            title(st,'fontsize',FS)
            xlabel('y (m)')
            ylabel('z (m)')
            zlabel('z (m)')
            hold on
            mesh(x_mesh_w_New,y_mesh_w_New,z_mesh_w_New)
            mesh(x_mesh_HTP_New,y_mesh_HTP_New,z_mesh_HTP_New)
            mesh(x_mesh_VTP1_New,y_mesh_VTP1_New,z_mesh_VTP1_New)
            mesh(x_mesh_VTP2_New,y_mesh_VTP2_New,z_mesh_VTP2_New)
            
            % Plots disck actuators wing 2
            pL_w1 = patch(Prop_L_x_w1,Prop_L_y_w1,Prop_L_z_w1,'y');
            set(pL_w1,'FaceLighting','phong','EdgeLighting','phong');
            set(pL_w1,'EraseMode','normal');
            pL_w1_s = patch(Prop_L_x_w1_s,Prop_L_y_w1_s,Prop_L_z_w1_s,'y');
            set(pL_w1_s,'FaceLighting','phong','EdgeLighting','phong');
            set(pL_w1_s,'EraseMode','normal');
            
            pR_w1 = patch(Prop_R_x_w1,Prop_R_y_w1,Prop_R_z_w1,'y');
            set(pR_w1,'FaceLighting','phong','EdgeLighting','phong');
            set(pR_w1,'EraseMode','normal');
            pR_w1_s = patch(Prop_R_x_w1_s,Prop_R_y_w1_s,Prop_R_z_w1_s,'y');
            set(pR_w1_s,'FaceLighting','phong','EdgeLighting','phong');
            set(pR_w1_s,'EraseMode','normal');

            % Plots disck actuators Engine 2
            pL_w1 = patch(Eng_L_x_w1,Eng_L_y_w1,Eng_L_z_w1,'y');
            set(pL_w1,'FaceColor',[.1 .8 0],'FaceLighting','phong','EdgeLighting','phong');
            set(pL_w1,'EraseMode','normal');
            pL_w1_s = patch(Eng_L_x_w1_s,Eng_L_y_w1_s,Eng_L_z_w1_s,'y');
            set(pL_w1_s,'FaceColor',[.1 .8 0],'FaceLighting','phong','EdgeLighting','phong');
            set(pL_w1_s,'EraseMode','normal');
            
            pR_w1 = patch(Eng_R_x_w1,Eng_R_y_w1,Eng_R_z_w1,'y');
            set(pR_w1,'FaceColor',[.1 .8 0],'FaceLighting','phong','EdgeLighting','phong');
            set(pR_w1,'EraseMode','normal');
            pR_w1_s = patch(Eng_R_x_w1_s,Eng_R_y_w1_s,Eng_R_z_w1_s,'y');
            set(pR_w1_s,'FaceColor',[.1 .8 0],'FaceLighting','phong','EdgeLighting','phong');
            set(pR_w1_s,'EraseMode','normal');

            hold off
            axis equal
            grid on
            
            z_max_fus = PLOTTING_UAV.z_max_fus;
            z_min_fus = PLOTTING_UAV.z_min_fus;
            % Defines limits of plots axis
            x_min_PLOT = -Delta_Plot;
            x_max_PLOT = X_AXIS(k);
            y_min_PLOT = -Y_AXIS(k)/2 ;
            y_max_PLOT = Y_AXIS(k)/2 ;
            z_min_PLOT = (z_min_fus - D_prop_vec(i,j,k)); 
            z_max_PLOT = Z_AXIS(k);
            axis([x_min_PLOT x_max_PLOT y_min_PLOT y_max_PLOT z_min_PLOT z_max_PLOT]);
            
            axis equal
            grid on
            PLOTTING_DATA{i,j,k} = PLOTTING_UAV;
            n=n+1;
        end
    end
    if SAVE_FIG==1
        prefix = strcat('SingleWing_T_D_');
        st = strcat(num2str(round(Vh(k),3,'significant')),'_m_s');
        name   = strcat(prefix,st);
        saveas(gcf,name,'fig');
        %saveas(gcf,name,'pdf');
        %saveas(gcf,name,'bmp');
        saveas(gcf,name,'png');
    end

end
