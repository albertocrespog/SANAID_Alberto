%%%%%% PRUEBA PLOTS %%%%%%
function Generates_Plots_Performance_v3(datos,OUTPUT_read_XLSX,Fig,Plot_Options,filenameS)
%close all

% Retrieves information for the plots
LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS;
% Fig = Plot_Options.Fig;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
mark_Type = Plot_Options.mark_Type;
mark_legend = Plot_Options.mark_legend;
MATLAB_in = Plot_Options.MATLAB_in;

% V_VAR = Plot_Options.V_VAR;
% m_VAR = Plot_Options.W_VAR;

% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filenamec = 'Figs';
% filename_DATA = strcat(filename,filenamed);
% filename_Plots = strcat(filename,filenamec);
fname = filenameS.plots;

time_zero = 0;
distance_zero = 0;
fuel_zero = 0;
if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
    mbatery_zero = 0;
    energy_zero = 0;
end
for i=1:length(datos)-1
    if strcmp(datos(i).nombre,'Climb') == 1
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+1)
            plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.mbaterias+mbatery_zero,'LineWidth',1.4)
            title('Batteries Mass Vs Time')
            xlabel('Time [min]')
            ylabel('Batteries Mass [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('BatteriesMass_vs_time');
                name   = strcat(prefix,st);
                % % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
                % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%                 saveas(gca, fullfile(fname, name), 'jpeg');
%                 saveas(gcf,fullfile(fname, name),'fig');
%                 saveas(gcf,fullfile(fname, name),'pdf');
%                 saveas(gcf,fullfile(fname, name),'bmp');
%                 saveas(gcf,fullfile(fname, name),'png');
            end
        else
            figure(Fig+1)
            plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.fuel+fuel_zero,'LineWidth',1.4)
            title('Fuel Vs Time')
            xlabel('Time [min]')
            ylabel('Fuel [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Fuel_vs_time');
                name   = strcat(prefix,st);
                % % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
                % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%                 saveas(gca, fullfile(fname, name), 'jpeg');
%                 saveas(gcf,fullfile(fname, name),'fig');
%                 saveas(gcf,fullfile(fname, name),'pdf');
%                 saveas(gcf,fullfile(fname, name),'bmp');
%                 saveas(gcf,fullfile(fname, name),'png');
            end
            
        end
        figure(Fig+2)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,(datos(i).segmento.distancia+distance_zero)/1000,'LineWidth',1.4)
        title('Distance Vs Time')
        xlabel('Time [min]')
        ylabel('Distance [km]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Distance_Vs_Time');
            name   = strcat(prefix,st);
            % % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%                 saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+3)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.palanca,'LineWidth',1.4)
        title('Throttle Vs Time')
        xlabel('Time [min]')
        ylabel('Throttle [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Throttle_Vs_Time');
            name   = strcat(prefix,st);
            % % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+4)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.empuje,'LineWidth',1.4)
        title('Thrust Vs Time')
        xlabel('Time [min]')
        ylabel('Thrust [N]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Thrust_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
                saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+5)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.veloc_vertical,'LineWidth',1.4)
        title('Vertical Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('VerticalVelocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+6)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.L_D,'LineWidth',1.4)
        title('L/D Vs Time')
        xlabel('Time [min]')
        ylabel('L/D [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('LD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+7)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.peso/9.80665,'LineWidth',1.4)
        title('Mass Vs Time')
        xlabel('Time [min]')
        ylabel('Mass [kg]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Mass_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+8)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.gamma*180/pi*ones(length(datos(i).segmento.tiempo(:,1)),1),'LineWidth',1.4)
        title('gamma Vs Time')
        xlabel('Time [min]')
        ylabel('gamma [º]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('gamma_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+9)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.velocidad,'LineWidth',1.4)
        title('Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Velocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+10)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.altura,'LineWidth',1.4)
        title('Height Vs Time')
        xlabel('Time [min]')
        ylabel('Height [m]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Height_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+11)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.CL,'LineWidth',1.4)
        title('CL Vs Time')
        xlabel('Time [min]')
        ylabel('CL [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CL_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+12)
        plot((datos(i).segmento.tiempo(:,1)+time_zero)/60,datos(i).segmento.CD,'LineWidth',1.4)
        title('CD Vs Time')
        xlabel('Time [min]')
        ylabel('CD [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+13)
        plot((datos(i).segmento.distancia+distance_zero)/1000,datos(i).segmento.altura,'LineWidth',1.4)
        title('Mission profile')
        xlabel('Distance [km]')
        ylabel('Height [m]')
        hold on
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+14)
            plot((datos(i).segmento.tiempo+time_zero)/60,(datos(i).segmento.energiav+energy_zero)/1000,'LineWidth',1.4)
            title('Energy Vs Time')
            xlabel('Time [min]')
            ylabel('Energy [KJ]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Energy_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
                % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%                 saveas(gca, fullfile(fname, name), 'jpeg');
%                 saveas(gcf,fullfile(fname, name),'fig');
%                 saveas(gcf,fullfile(fname, name),'pdf');
%                 saveas(gcf,fullfile(fname, name),'bmp');
%                 saveas(gcf,fullfile(fname, name),'png');
            end
            
        end
        
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('MissionProfile_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
    end
    if strcmp(datos(i).nombre,'VTOL Climb') == 1
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+1)
            plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.mbaterias+mbatery_zero,'LineWidth',1.4)
            title('Batteries Mass Vs Time')
            xlabel('Time [min]')
            ylabel('Batteries Mass [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('BatteriesMass_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
% Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%                 saveas(gca, fullfile(fname, name), 'jpeg');
%                 saveas(gcf,fullfile(fname, name),'fig');
%                 saveas(gcf,fullfile(fname, name),'pdf');
%                 saveas(gcf,fullfile(fname, name),'bmp');
%                 saveas(gcf,fullfile(fname, name),'png');
            end
            
        else
            figure(Fig+1)
            plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.fuel*ones(1,length(datos(i).segmento.tiempo))+fuel_zero,'LineWidth',1.4)
            title('Fuel Vs Time')
            xlabel('Time [min]')
            ylabel('Fuel [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Fuel_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
                % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%                 saveas(gca, fullfile(fname, name), 'jpeg');
%                 saveas(gcf,fullfile(fname, name),'fig');
%                 saveas(gcf,fullfile(fname, name),'pdf');
%                 saveas(gcf,fullfile(fname, name),'bmp');
%                 saveas(gcf,fullfile(fname, name),'png');
            end
            
        end
        figure(Fig+2)
        plot((datos(i).segmento.tiempo+time_zero)/60,zeros(length(datos(i).segmento.tiempo))+distance_zero,'LineWidth',1.4)
        title('Distance Vs Time')
        xlabel('Time [min]')
        ylabel('Distance [km]')
        hold on
        
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Distance_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+3)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.palanca,'LineWidth',1.4)
        title('Throttle Vs Time')
        xlabel('Time [min]')
        ylabel('Throttle [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Throttle_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
% Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%                 saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+4)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.empuje,'LineWidth',1.4)
        title('Thrust Vs Time')
        xlabel('Time [min]')
        ylabel('Thrust [N]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Thrust_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
% Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%                 saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+5)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.veloc_vertical,'LineWidth',1.4)
        title('Vertical Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('VerticalVelocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%                 saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+6)
        plot((datos(i).segmento.tiempo+time_zero)/60,zeros(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('L/D Vs Time')
        xlabel('Time [min]')
        ylabel('L/D [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('LD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
                SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+7)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.peso/9.80665*ones(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('Mass Vs Time')
        xlabel('Time [min]')
        ylabel('Mass [kg]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Mass_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+8)
        plot((datos(i).segmento.tiempo+time_zero)/60,90*ones(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('gamma Vs Time')
        xlabel('Time [min]')
        ylabel('gamma [º]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('gamma_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
                        % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+9)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.veloc_vertical,'LineWidth',1.4)
        title('Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Velocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+10)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.altura,'LineWidth',1.4)
        title('Height Vs Time')
        xlabel('Time [min]')
        ylabel('Height [m]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Height_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+11)
        plot((datos(i).segmento.tiempo+time_zero)/60,zeros(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('CL Vs Time')
        xlabel('Time [min]')
        ylabel('CL [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CL_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+12)
        plot((datos(i).segmento.tiempo+time_zero)/60,ones(1,length(datos(i).segmento.tiempo))*datos(i).segmento.CD,'LineWidth',1.4)
        title('CD Vs Time')
        xlabel('Time [min]')
        ylabel('CD [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+13)
        plot(zeros(length(datos(i).segmento.tiempo))+distance_zero,datos(i).segmento.altura,'LineWidth',1.4)
        title('Mission profile')
        xlabel('Distance [km]')
        ylabel('Height [m]')
        hold on
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+14)
            plot((datos(i).segmento.tiempo+time_zero)/60,(datos(i).segmento.energiav+energy_zero)/1000,'LineWidth',1.4)
            title('Energy Vs Time')
            xlabel('Time [min]')
            ylabel('Energy [KJ]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Energy_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
            end
            
        end
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('MissionProfile_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
    end
    if strcmp(datos(i).nombre,'Cruise') == 1
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+1)
            plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.mbaterias+mbatery_zero,'LineWidth',1.4)
            title('Batteries Mass Vs Time')
            xlabel('Time [min]')
            ylabel('Batteries Mass [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('BatteriesMass_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

            end
            
        else
            figure(Fig+1)
            plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.fuel+fuel_zero,'LineWidth',1.4)
            title('Fuel Vs Time')
            xlabel('Time [min]')
            ylabel('Fuel [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Fuel_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

            end
            
        end
        figure(Fig+2)
        plot((datos(i).segmento.tiempo+time_zero)/60,(datos(i).segmento.distancia+distance_zero)/1000,'LineWidth',1.4)
        title('Distance Vs Time')
        xlabel('Time [min]')
        ylabel('Distance [km]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Distance_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
                % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+3)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.palanca,'LineWidth',1.4)
        title('Throttle Vs Time')
        xlabel('Time [min]')
        ylabel('Throttle [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Throttle_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
              % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+4)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.empuje,'LineWidth',1.4)
        title('Thrust Vs Time')
        xlabel('Time [min]')
        ylabel('Thrust [N]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Thrust_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
                     % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+5)
        plot((datos(i).segmento.tiempo+time_zero)/60,zeros(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('Vertical Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('VerticalVelocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+6)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.L_D,'LineWidth',1.4)
        title('L/D Vs Time')
        xlabel('Time [min]')
        ylabel('L/D [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('LD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+7)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.peso/9.80665,'LineWidth',1.4)
        title('Mass Vs Time')
        xlabel('Time [min]')
        ylabel('Mass [kg]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Mass_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options
            SAVE_types(fname,name,gca,gcf);
%             saveas(gca, fullfile(fname, name), 'jpeg');
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+8)
        plot((datos(i).segmento.tiempo+time_zero)/60,zeros(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('gamma Vs Time')
        xlabel('Time [min]')
        ylabel('gamma [º]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('gamma_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+9)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.velocidad*ones(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Velocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+10)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.altura*ones(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('Height Vs Time')
        xlabel('Time [min]')
        ylabel('Height [m]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Height_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+11)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.CL,'LineWidth',1.4)
        title('CL Vs Time')
        xlabel('Time [min]')
        ylabel('CL [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CL_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+12)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.CD,'LineWidth',1.4)
        title('CD Vs Time')
        xlabel('Time [min]')
        ylabel('CD [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+13)
        plot((datos(i).segmento.distancia+distance_zero)/1000,datos(i).segmento.altura*ones(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('Mission profile')
        xlabel('Distance [km]')
        ylabel('Height [m]')
        hold on
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+14)
            plot((datos(i).segmento.tiempo+time_zero)/60,(datos(i).segmento.energiav+energy_zero)/1000,'LineWidth',1.4)
            title('Energy Vs Time')
            xlabel('Time [min]')
            ylabel('Energy [KJ]')
            hold on
        end
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Energy_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
    end
    if strcmp(datos(i).nombre,'Load Deployment') == 1
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+1)
            plot(0,0)
            title('Batteries Mass Vs Time')
            xlabel('Time [min]')
            ylabel('Batteries Mass [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('BatteriesMass_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
            end
            
        else
            figure(Fig+1)
            plot(0,0)
            title('Fuel Vs Time')
            xlabel('Time [min]')
            ylabel('Fuel [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Fuel_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
            end
            
        end
        figure(Fig+2)
        plot(0,0)
        title('Distance Vs Time')
        xlabel('Time [min]')
        ylabel('Distance [km]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Distance_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+3)
        plot(0,0)
        title('Throttle Vs Time')
        xlabel('Time [min]')
        ylabel('Throttle [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Throttle_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+4)
        plot(0,0)
        title('Thrust Vs Time')
        xlabel('Time [min]')
        ylabel('Thrust [N]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Thrust_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+5)
        plot(0,0)
        title('Vertical Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('VerticalVelocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+6)
        plot(0,0)
        title('L/D Vs Time')
        xlabel('Time [min]')
        ylabel('L/D [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('LD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+7)
        plot(0,0)
        title('Mass Vs Time')
        xlabel('Time [min]')
        ylabel('Mass [kg]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Mass_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+8)
        plot(0,0)
        title('gamma Vs Time')
        xlabel('Time [min]')
        ylabel('gamma [º]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('gamma_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+9)
        plot(0,0)
        title('Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Velocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            saveas(gca, fullfile(fname, name), 'jpeg');
            saveas(gcf,fullfile(fname, name),'fig');
            saveas(gcf,fullfile(fname, name),'pdf');
            saveas(gcf,fullfile(fname, name),'bmp');
            saveas(gcf,fullfile(fname, name),'png');
        end
        
        figure(Fig+10)
        plot(0,0)
        title('Height Vs Time')
        xlabel('Time [min]')
        ylabel('Height [m]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Height_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+11)
        plot(0,0)
        title('CL Vs Time')
        xlabel('Time [min]')
        ylabel('CL [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CL_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
               % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+12)
        plot(0,0)
        title('CD Vs Time')
        xlabel('Time [min]')
        ylabel('CD [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
               % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+13)
        plot(0,0)
        title('Mission profile')
        xlabel('Distance [km]')
        ylabel('Height [m]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('MissionProfile_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
               % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+14)
            plot(0,0)
            title('Energy Vs Time')
            xlabel('Time [min]')
            ylabel('Energy [KJ]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Energy_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
             % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

            end
            
        end
    end
    if strcmp(datos(i).nombre,'Descent') == 1
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+1)
            plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.mbaterias+mbatery_zero,'LineWidth',1.4)
            title('Batteries Mass Vs Time')
            xlabel('Time [min]')
            ylabel('Batteries Mass [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('BatteriesMass_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
               % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

            end
            
        else
            figure(Fig+1)
            plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.fuel+fuel_zero,'LineWidth',1.4)
            title('Fuel Vs Time')
            xlabel('Time [min]')
            ylabel('Fuel [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Fuel_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

            end
            
        end
        figure(Fig+2)
        plot((datos(i).segmento.tiempo+time_zero)/60,(datos(i).segmento.distancia+distance_zero)/1000,'LineWidth',1.4)
        title('Distance Vs Time')
        xlabel('Time [min]')
        ylabel('Distance [km]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Distance_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
                 % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+3)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.palanca,'LineWidth',1.4)
        title('Throttle Vs Time')
        xlabel('Time [min]')
        ylabel('Throttle [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Throttle_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
               % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+4)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.empuje,'LineWidth',1.4)
        title('Thrust Vs Time')
        xlabel('Time [min]')
        ylabel('Thrust [N]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Thrust_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+5)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.veloc_vertical,'LineWidth',1.4)
        title('Vertical Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('VerticalVelocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+6)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.L_D,'LineWidth',1.4)
        title('L/D Vs Time')
        xlabel('Time [min]')
        ylabel('L/D [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('LD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
             % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+7)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.peso/9.80665,'LineWidth',1.4)
        title('Mass Vs Time')
        xlabel('Time [min]')
        ylabel('Mass [kg]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Mass_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
               % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+8)
        % plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.gamma*180/pi*ones(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.gamma*180/pi,'LineWidth',1.4)
        title('gamma Vs Time')
        xlabel('Time [min]')
        ylabel('gamma [º]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('gamma_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+9)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.velocidad,'LineWidth',1.4)
        title('Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Velocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+10)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.altura,'LineWidth',1.4)
        title('Height Vs Time')
        xlabel('Time [min]')
        ylabel('Height [m]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Height_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+11)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.CL,'LineWidth',1.4)
        title('CL Vs Time')
        xlabel('Time [min]')
        ylabel('CL [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CL_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+12)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.CD,'LineWidth',1.4)
        title('CD Vs Time')
        xlabel('Time [min]')
        ylabel('CD [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+13)
        plot((datos(i).segmento.distancia+distance_zero)/1000,datos(i).segmento.altura,'LineWidth',1.4)
        title('Mission profile')
        xlabel('Distance [km]')
        ylabel('Height [m]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('MissionProfile_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+14)
            plot((datos(i).segmento.tiempo+time_zero)/60,(datos(i).segmento.energiav+energy_zero)/1000,'LineWidth',1.4)
            title('Energy Vs Time')
            xlabel('Time [min]')
            ylabel('Energy [KJ]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Energy_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
              % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

            end
            
        end
    end
    if strcmp(datos(i).nombre,'VTOL Descent') == 1
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+1)
            plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.mbaterias+mbatery_zero,'LineWidth',1.4)
            title('Batteries Mass Vs Time')
            xlabel('Time [min]')
            ylabel('Batteries Mass [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('BatteriesMass_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

            end
            
        else
            figure(Fig+1)
            plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.fuel*ones(1,length(datos(i).segmento.tiempo))+fuel_zero,'LineWidth',1.4)
            title('Fuel Vs Time')
            xlabel('Time [min]')
            ylabel('Fuel [kg]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Fuel_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

            end
            
        end
        figure(Fig+2)
        plot((datos(i).segmento.tiempo+time_zero)/60,zeros(length(datos(i).segmento.tiempo))+distance_zero,'LineWidth',1.4)
        title('Distance Vs Time')
        xlabel('Time [min]')
        ylabel('Distance [km]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Distance_Vs_Time');
            name   = strcat(prefix,st);
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+3)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.palanca,'LineWidth',1.4)
        title('Throttle Vs Time')
        xlabel('Time [min]')
        ylabel('Throttle [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Throttle_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+4)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.empuje,'LineWidth',1.4)
        title('Thrust Vs Time')
        xlabel('Time [min]')
        ylabel('Thrust [N]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Thrust_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+5)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.veloc_vertical,'LineWidth',1.4)
        title('Vertical Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('VerticalVelocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+6)
        plot((datos(i).segmento.tiempo+time_zero)/60,zeros(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('L/D Vs Time')
        xlabel('Time [min]')
        ylabel('L/D [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('LD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+7)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.peso/9.80665*ones(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('Mass Vs Time')
        xlabel('Time [min]')
        ylabel('Mass [kg]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Mass_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+8)
        plot((datos(i).segmento.tiempo+time_zero)/60,-90*ones(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('gamma Vs Time')
        xlabel('Time [min]')
        ylabel('gamma [º]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('gamma_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
                % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+9)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.veloc_vertical,'LineWidth',1.4)
        title('Velocity Vs Time')
        xlabel('Time [min]')
        ylabel('Velocity [m/s]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Velocity_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+10)
        plot((datos(i).segmento.tiempo+time_zero)/60,datos(i).segmento.altura,'LineWidth',1.4)
        title('Height Vs Time')
        xlabel('Time [min]')
        ylabel('Height [m]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('Height_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+11)
        plot((datos(i).segmento.tiempo+time_zero)/60,zeros(1,length(datos(i).segmento.tiempo)),'LineWidth',1.4)
        title('CL Vs Time')
        xlabel('Time [min]')
        ylabel('CL [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CL_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+12)
        plot((datos(i).segmento.tiempo+time_zero)/60,ones(1,length(datos(i).segmento.tiempo))*datos(i).segmento.CD,'LineWidth',1.4)
        title('CD Vs Time')
        xlabel('Time [min]')
        ylabel('CD [-]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('CD_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        figure(Fig+13)
        plot(zeros(length(datos(i).segmento.tiempo))+distance_zero,datos(i).segmento.altura,'LineWidth',1.4)
        title('Mission profile')
        xlabel('Distance [km]')
        ylabel('Height [m]')
        hold on
        if SAVE_FIGS==1
            prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
            st = strcat('MissionProfile_Vs_Time');
            name   = strcat(prefix,st);
            % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

        end
        
        if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
            figure(Fig+14)
            plot((datos(i).segmento.tiempo+time_zero)/60,(datos(i).segmento.energiav+energy_zero)/1000,'LineWidth',1.4)
            title('Energy Vs Time')
            xlabel('Time [min]')
            ylabel('Energy [KJ]')
            hold on
            if SAVE_FIGS==1
                prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
                st = strcat('Energy_Vs_Time');
                name   = strcat(prefix,st);
                % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

            end
            
        end
    end
    time_zero = time_zero+datos(i).segmento.tiempo(end);
    distance_zero = (distance_zero+datos(i).segmento.distancia(end))/1000;
    if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
        mbatery_zero = mbatery_zero+datos(i).segmento.mbaterias(end);
        energy_zero = (energy_zero+datos(i).segmento.energiav(end))/1000;
    else
        fuel_zero = fuel_zero+datos(i).segmento.fuel(end);
    end
end
if OUTPUT_read_XLSX.Propulsive_flags.propul(1) == 4
    for i=(Fig+1):(Fig+14)
        figure(i)
        grid on
        legend(datos(1:(end-1)).nombre)
    end
else
    for i=(Fig+1):(Fig+13)
        figure(i)
        grid on
        legend(datos(1:(end-1)).nombre)
    end
end