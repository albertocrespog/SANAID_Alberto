function  [DATA_mass,Fig] = variation_mass_function(Prop_data,D_prop_vect,AC_CONFIGURATION,SF_prop,...
    SE_battery,RPMMAX_APC,Weight_tier,Fig,Prop_selection,N_contour_lines)
        
% Mass variation
m_min = 15;
m_max = 25;
N_mass = 100;
m_vect=linspace(m_min,m_max,N_mass);
rho = 1.225;
solve_w_fero = 0;
SAVE_FIGS = 1;

% actualizes the Prop geometry
Prop_data = Generation_Propulsion_Data(SF_prop,D_prop_vect(1),AC_CONFIGURATION,Prop_selection); % Defines Propulsion DATA
SP_motor = Weight_tier.M_ENERGY.SP_motor; % Motor Specific Power kW/kg - TIER 1
SP_ESC = Weight_tier.M_ENERGY.SP_ESC;
m_prop_din = Weight_tier.M_ENERGY.m_prop_din;
% Porcentaje de masa de baterias

n_eng_vec = [4 6 8 12];
m_payload = 0;
SF_estructure = 0.85;

% Identifies the plots that are conducted
plots_3D_mass_variation = 0;
plots_contour_mass_variation = 1;

for n =1:length(n_eng_vec)
    for j = 1:length(D_prop_vect)
        D_prop = D_prop_vect(j);
        D_inch = D_prop/2.54*100;
        rpm_max = SF_prop*RPMMAX_APC/D_inch;
        nps_max = rpm_max/60;
        for i = 1:length(m_vect)
            n_eng = n_eng_vec(n);
            V_Hov = 0;
            % Estimation of Desired Thrust
            m_TOW(i) = m_vect(i);
            WeW0(i) = SF_estructure*0.4666*(m_TOW(i)*9.81)^(-0.02);
            m_structure(i,j) = WeW0(i)*(m_TOW(i));
%             m_batt(i) = m_vect(i)*pp_bat_vect(n);
            Fdes_Hov = (m_TOW(i)*9.81); % For hover condition
            % Function that determines the rev per second for a given prop
            
            % Efficiencies
            eta_gear = 0.96;                   % Gear box efficiency.
            eta_m    = 0.88;                      % Engine efficiency (output/input).
            eta_esc  = 0.98;                    % Speed controller efficiency.
            eta_dist = 0.96;                   % Shaft efficiency.
            eta_ES   = eta_m*eta_esc*eta_dist*eta_gear; % Electric system efficiency

            M = 0.70; % Figure of merit
            S_disk = (pi*((D_prop/2)^2));
            
            % Induced velocity at Hover
            T_vpf = Fdes_Hov/n_eng;
            vi0 = sqrt((T_vpf)/(2*1.225*S_disk));
            P_ideal = n_eng*vi0*(T_vpf);
            Pe_plot_design_Hov(i,j) = P_ideal*(1/(M*eta_ES));
%             Pe_plot_design_Hov(i,j) = (1/(M*eta_ES))*(Fdes_Hov^(3/2))/sqrt(2*1.225*S_disk)
            Pe_eng_plot_design_Hov(i,j) = Pe_plot_design_Hov(i,j)/n_eng;
            % diameter
            n_Hov(i,j) = Selection_J_CT_F_design_mass(V_Hov,Prop_data,Fdes_Hov,rho,D_prop_vect(j),solve_w_fero,n_eng);
            
            delta_T_Hov(i,j) = n_Hov(i,j)/nps_max;
            RPM_Hov(i,j) = n_Hov(i,j)*60;
            Propulsion_Hov{i,j} = get_EngineProperties_design_mass(V_Hov,rho,n_Hov(i,j),Prop_data,AC_CONFIGURATION,D_prop_vect(j),n_eng);
            RPM_plot_design_Hov(i,j) = Propulsion_Hov{i,j}.RPM;
            T_eng_plot_design_Hov(i,j) = Propulsion_Hov{i,j}.Ti_eng;
            Qi_eng_Hov(i,j) = Propulsion_Hov{i,j}.Qi_eng;
%             Pe_plot_design_Hov(i,j) = Propulsion_Hov{i,j}.Pe;
%             Pe_eng_plot_design_Hov(i,j) = Propulsion_Hov{i,j}.Pe_eng;
            m_motor(i,j) = (Pe_eng_plot_design_Hov(i,j)/1000)/SP_motor;
            m_propeller(i,j) = n_eng*(m_prop_din*(D_prop_vect(j)*100/2.54));
            
            m_ESC(i,j) = (Pe_eng_plot_design_Hov(i,j)/1000)/SP_ESC;
            m_prop_system(i,j) = n_eng*(m_motor(i,j) + m_ESC(i,j))  + m_propeller(i,j);
            m_batt(i,j) =  m_TOW(i) - m_structure(i,j) - m_prop_system(i,j) - m_payload;
            pp_bat_vect(i,j) = m_batt(i,j)/m_TOW(i);
            Pi_plot_design_Hov(i,j) = Propulsion_Hov{i,j}.Pi;
            Pi_eng_plot_design_Hov(i,j) = Propulsion_Hov{i,j}.Pi_eng;
            etha_mp_plot_design_Hov(i,j) = Propulsion_Hov{i,j}.ethamp;
            etha_mp_total_plot_design_Hov(i,j) = Propulsion_Hov{i,j}.etha_emp;
            
            % Reserva remanente
            tau_reserve = 0.20;
            
            Energy_available_Hov(i,j) = (m_batt(i,j)*SE_battery)*(1-tau_reserve); % en W-h % Wh/kg
            time_Hov(i,j) = Energy_available_Hov(i,j)*3600/Pe_plot_design_Hov(i,j);
        end
    end
    DATA_mass = time_Hov;
    
    X = D_prop_vect*100/2.54;
    Y = m_vect;
   
    Z_time = time_Hov/60;
    min_time = min(time_Hov(:))/60;
    max_time = max(time_Hov(:))/60;
    vect_cc_time = linspace(min_time,max_time,N_contour_lines);
    
    Z_Pe_eng = Pe_eng_plot_design_Hov/1000;
    min_Pe_eng = min(Pe_eng_plot_design_Hov(:))/1000;
    max_Pe_eng = max(Pe_eng_plot_design_Hov(:))/1000;
    vect_cc_Pe_eng = linspace(min_Pe_eng,max_Pe_eng,N_contour_lines);
    
    Z_Pe = Pe_plot_design_Hov/1000;
    min_Pe = min(Pe_plot_design_Hov(:))/1000;
    max_Pe = max(Pe_plot_design_Hov(:))/1000;
    vect_cc_Pe = linspace(min_Pe,max_Pe,N_contour_lines);
    
    Z_Energy = Energy_available_Hov/1000;
    min_Energy = min(Energy_available_Hov(:))/1000;
    max_Energy = max(Energy_available_Hov(:))/1000;
    vect_cc_Energy = linspace(min_Energy,max_Energy,N_contour_lines);
    
    Z_m_prop = m_prop_system;
    min_m_prop = min(m_prop_system(:));
    max_m_prop = max(m_prop_system(:));
    vect_cc_m_prop = linspace(min_m_prop,max_m_prop,N_contour_lines);

    Z_m_structure = m_structure;
    min_m_structure = min(m_structure(:));
    max_m_structure = max(m_structure(:));
    vect_cc_m_structure = linspace(min_m_structure,max_m_structure,N_contour_lines);

    Z_m_propeller = m_propeller;
    min_m_propeller = min(m_propeller(:));
    max_m_propeller = max(m_propeller(:));
    vect_cc_m_propeller = linspace(min_m_propeller,max_m_propeller,N_contour_lines);

    Z_m_batt = m_batt;
    min_m_batt = min(m_batt(:));
    max_m_batt = max(m_batt(:));
    vect_cc_m_batt = linspace(min_m_batt,max_m_batt,N_contour_lines);
    
    Z_delta_T_Hov = delta_T_Hov;
    min_delta_T_Hov = min(delta_T_Hov(:));
    max_delta_T_Hov = max(delta_T_Hov(:));
    vect_cc_delta_T_Hov = linspace(min_delta_T_Hov,max_delta_T_Hov,N_contour_lines);


    if plots_3D_mass_variation == 1
        
        prefix = strcat('CERVERA_3D');
        % Legends
        name1   = strcat(prefix,'Time_vs_mbat_n_eng_',num2str(n_eng),'');
        name2   = strcat(prefix,'Pe_eng_vs_mbat_n_eng_',num2str(n_eng),'');
        name3   = strcat(prefix,'Pe_vs_mbat_n_eng_',num2str(n_eng),'');
        name4   = strcat(prefix,'E_vs_mbat_n_eng_',num2str(n_eng),'');
        name5   = strcat(prefix,'m_prop_vs_mbat_n_eng_',num2str(n_eng),'');
        name6   = strcat(prefix,'m_estructure_vs_mbat_n_eng_',num2str(n_eng),'');
        name7   = strcat(prefix,'m_propeller_vs_mbat_n_eng_',num2str(n_eng),'');
        name8   = strcat(prefix,'m_energy_vs_mbat_n_eng_',num2str(n_eng),'');
        
        figure(Fig)
        mesh(X,Y,Z_time)
        grid on
        colorbar
        Title1 = strcat('Time vs. D_{prop} & Battery Mass with n_{eng}=',num2str(n_eng));
        title(Title1)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        zlabel('Time (min)')
        if SAVE_FIGS==1
            saveas(gcf,name1,'fig');
            %         saveas(gca,name1,'epsc');
            %         saveas(gcf,name1,'pdf');
            %         saveas(gcf,name1,'bmp');
            saveas(gcf,name1,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        mesh(X,Y,Z_Pe_eng)
        grid on
        colorbar
        Title2 = strcat('Electric Power per engine vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title2)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        zlabel('P_e per engine (kW)')
        if SAVE_FIGS==1
            saveas(gcf,name2,'fig');
            %         saveas(gca,name2,'epsc');
            %         saveas(gcf,name2,'pdf');
            %         saveas(gcf,name2,'bmp');
            saveas(gcf,name2,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        mesh(X,Y,Z_Pe)
        grid on
        colorbar
        Title3 = strcat('Electric Power vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title3)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        zlabel('P_e (kW)')
        if SAVE_FIGS==1
            saveas(gcf,name3,'fig');
            %         saveas(gca,name3,'epsc');
            %         saveas(gcf,name3,'pdf');
            %         saveas(gcf,name3,'bmp');
            saveas(gcf,name3,'png');
        end
        Fig = Fig + 1;
   
        figure(Fig)
        mesh(X,Y,Z_Energy)
        grid on
        colorbar
        Title4 = strcat('Energy vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title4)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        zlabel('Energy (kW-h)')
        if SAVE_FIGS==1
            saveas(gcf,name4,'fig');
            %         saveas(gca,name4,'epsc');
            %         saveas(gcf,name4,'pdf');
            %         saveas(gcf,name4,'bmp');
            saveas(gcf,name4,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        mesh(X,Y,Z_m_prop)
        grid on
        colorbar
        Title5 = strcat('Propulsive System (Engine, ESC and Props) (kg) vs. D_{prop} & Drone Mass n_{eng}=',num2str(n_eng));
        title(Title5)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        zlabel('m_{prop} - Engine, ESC and Props (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name5,'fig');
            %         saveas(gca,name5,'epsc');
            %         saveas(gcf,name5,'pdf');
            %         saveas(gcf,name5,'bmp');
            saveas(gcf,name5,'png');
        end
        Fig = Fig + 1;
        

        figure(Fig)
        mesh(X,Y,Z_m_structure)
        grid on
        colorbar
        Title6 = strcat('Estructure (kg) vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title6)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        zlabel('m_{estructure} (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name6,'fig');
            %         saveas(gca,name6,'epsc');
            %         saveas(gcf,name6,'pdf');
            %         saveas(gcf,name6,'bmp');
            saveas(gcf,name6,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        mesh(X,Y,Z_m_propeller)
        grid on
        colorbar
        Title7 = strcat('Propellers mass (kg) vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title7)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        zlabel('m_{propeller} (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name7,'fig');
            %         saveas(gca,name6,'epsc');
            %         saveas(gcf,name6,'pdf');
            %         saveas(gcf,name6,'bmp');
            saveas(gcf,name7,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        mesh(X,Y,Z_m_batt)
        grid on
        colorbar
        Title8 = strcat('Energy System Mass (kg) vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title8)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        zlabel('m_{energy} (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name7,'fig');
            %         saveas(gca,name7,'epsc');
            %         saveas(gcf,name7,'pdf');
            %         saveas(gcf,name7,'bmp');
            saveas(gcf,name7,'png');
        end
        Fig = Fig + 1;
    end
    
    if plots_contour_mass_variation == 1
        prefix = strcat('CERVERA_contour');
        % Legends
        name1   = strcat(prefix,'Time_vs_mbat_n_eng_',num2str(n_eng),'');
        name2   = strcat(prefix,'Pe_eng_vs_mbat_n_eng_',num2str(n_eng),'');
        name3   = strcat(prefix,'Pe_vs_mbat_n_eng_',num2str(n_eng),'');
        name4   = strcat(prefix,'E_vs_mbat_n_eng_',num2str(n_eng),'');
        name5   = strcat(prefix,'m_prop_vs_mbat_n_eng_',num2str(n_eng),'');
        name6   = strcat(prefix,'m_estructure_vs_mbat_n_eng_',num2str(n_eng),'');
        name7   = strcat(prefix,'m_propeller_vs_mbat_n_eng_',num2str(n_eng),'');
        name8   = strcat(prefix,'m_energy_vs_mbat_n_eng_',num2str(n_eng),'');
        name9   = strcat(prefix,'delta_T_Hov_vs_mbat_n_eng_',num2str(n_eng),'');
        
        
        figure(Fig)
        [C_time,h_time] = contour(X,Y,Z_time,vect_cc_time');
        clabel(C_time,h_time)
        grid on
        Title1 = strcat('Time (min) vs. D_{prop} & Battery Mass with n_{eng}=',num2str(n_eng));
        title(Title1)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name1,'fig');
            %         saveas(gca,name1,'epsc');
            %         saveas(gcf,name1,'pdf');
            %         saveas(gcf,name1,'bmp');
            saveas(gcf,name1,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_Pe_eng,h_Pe_eng] = contour(X,Y,Z_Pe_eng,vect_cc_Pe_eng');
        clabel(C_Pe_eng,h_Pe_eng)
        grid on
        Title2b = strcat('Electric Power (kW) per engine vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title2b)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name2,'fig');
            %         saveas(gca,name2,'epsc');
            %         saveas(gcf,name2,'pdf');
            %         saveas(gcf,name2,'bmp');
            saveas(gcf,name2,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_Pe,h_Pe] = contour(X,Y,Z_Pe,vect_cc_Pe');
        clabel(C_Pe,h_Pe)
        grid on
        Title3b = strcat('Total Electric Power (kW) vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title3b)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name3,'fig');
            %         saveas(gca,name3,'epsc');
            %         saveas(gcf,name3,'pdf');
            %         saveas(gcf,name3,'bmp');
            saveas(gcf,name3,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        [C_Energy,h_Energy] = contour(X,Y,Z_Energy,vect_cc_Energy');
        clabel(C_Energy,h_Energy)
        grid on
        Title4b = strcat('Energy (kW-h) vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title4b)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name4,'fig');
            %         saveas(gca,name4,'epsc');
            %         saveas(gcf,name4,'pdf');
            %         saveas(gcf,name4,'bmp');
            saveas(gcf,name4,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_m_prop,h_m_prop] = contour(X,Y,Z_m_prop,vect_cc_m_prop');
        clabel(C_m_prop,h_m_prop)
        grid on
        Title5b = strcat('Propulsive System (Engine, ESC and Props) (kg) vs. D_{prop} & Drone Mass n_{eng}=',num2str(n_eng));
        title(Title5b)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name5,'fig');
            %         saveas(gca,name5,'epsc');
            %         saveas(gcf,name5,'pdf');
            %         saveas(gcf,name5,'bmp');
            saveas(gcf,name5,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_m_structure,h_m_structure] = contour(X,Y,Z_m_structure,vect_cc_m_structure');
        clabel(C_m_structure,h_m_structure)
        grid on
        Title6b = strcat('Estructure (kg) vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title6b)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name6,'fig');
            %         saveas(gca,name6,'epsc');
            %         saveas(gcf,name6,'pdf');
            %         saveas(gcf,name6,'bmp');
            saveas(gcf,name6,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_m_propeller,h_m_propeller] = contour(X,Y,Z_m_propeller,vect_cc_m_propeller');
        clabel(C_m_propeller,h_m_propeller)
        grid on
        Title7b = strcat('Propellers mass (kg) vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title7b)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name7,'fig');
            %         saveas(gca,name6,'epsc');
            %         saveas(gcf,name6,'pdf');
            %         saveas(gcf,name6,'bmp');
            saveas(gcf,name7,'png');
        end
        Fig = Fig + 1;
                
        figure(Fig)
        [C_m_batt,h_m_batt] = contour(X,Y,Z_m_batt,vect_cc_m_batt');
        clabel(C_m_batt,h_m_batt)
        grid on
        Title8b = strcat('Energy System Mass (kg) vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title8b)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name8,'fig');
            %         saveas(gca,name7,'epsc');
            %         saveas(gcf,name7,'pdf');
            %         saveas(gcf,name7,'bmp');
            saveas(gcf,name8,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_delta_T_Hov,h_delta_T_Hov] = contour(X,Y,Z_delta_T_Hov,vect_cc_delta_T_Hov');
        clabel(C_delta_T_Hov,h_delta_T_Hov)
        grid on
        Title9b = strcat('Throttle vs. D_{prop} & Drone Mass with n_{eng}=',num2str(n_eng));
        title(Title9b)
        xlabel('D_p_r_o_p (in)')
        ylabel('Drone mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name9,'fig');
            %         saveas(gca,name7,'epsc');
            %         saveas(gcf,name7,'pdf');
            %         saveas(gcf,name7,'bmp');
            saveas(gcf,name9,'png');
        end
        Fig = Fig + 1;
        
        
    
    end
% pause
end
