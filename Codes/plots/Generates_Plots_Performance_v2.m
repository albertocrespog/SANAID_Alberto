%%%%%% PRUEBA PLOTS %%%%%%
function [Fig] =  Generates_Plots_Performance_v2(datos,Fig)
close all
time_zero = 0;
distance_zero = 0;


for i=1:length(datos)-1
    if strcmp(datos(i).nombre,'Climb') == 1
        Fig = Fig + 1;
        Fig0 = Fig;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.mbaterias*ones(1,length(datos(i).segmento.tiempo)))
        title('Batteries Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Batteries Mass [kg]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.distancia+distance_zero+distance_zero)
        title('Distance Vs Time')
        xlabel('Time [s]')
        ylabel('Distance [m]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.palanca)
        title('Throttle Vs Time')
        xlabel('Time [s]')
        ylabel('Throttle [-]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.empuje)
        title('Thrust Vs Time')
        xlabel('Time [s]')
        ylabel('Thrust [N]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.veloc_vertical)
        title('Vertical Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.L_D)
        title('L/D Vs Time')
        xlabel('Time [s]')
        ylabel('L/D [-]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.peso)
        title('Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Mass [kg]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.gamma*180/pi*ones(length(datos(i).segmento.tiempo)))
        title('gamma Vs Time')
        xlabel('Time [s]')
        ylabel('gamma [º]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.velocidad)
        title('Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.altura)
        title('Height Vs Time')
        xlabel('Time [s]')
        ylabel('Height [m]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.CL)
        title('CL Vs Time')
        xlabel('Time [s]')
        ylabel('CL [-]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.CD)
        title('CD Vs Time')
        xlabel('Time [s]')
        ylabel('CD [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.distancia+distance_zero,datos(i).segmento.altura)
        title('Mission profile')
        xlabel('Distance [m]')
        ylabel('Height [m]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.energiav)
        title('Energy Vs Time')
        xlabel('Time [s]')
        ylabel('Energy [J]')
        hold on
    end
    if strcmp(datos(i).nombre,'VTOL Climb') == 1
        Fig = Fig + 1;
        Fig0 = Fig;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.mbaterias*ones(1,length(datos(i).segmento.tiempo)))
        title('Batteries Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Batteries Mass [kg]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,zeros(length(datos(i).segmento.tiempo))+distance_zero)
        title('Distance Vs Time')
        xlabel('Time [s]')
        ylabel('Distance [m]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.palanca*ones(1,length(datos(i).segmento.tiempo)))
        title('Throttle Vs Time')
        xlabel('Time [s]')
        ylabel('Throttle [-]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.empuje)
        title('Thrust Vs Time')
        xlabel('Time [s]')
        ylabel('Thrust [N]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.veloc_vertical)
        title('Vertical Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,zeros(1,length(datos(i).segmento.tiempo)))
        title('L/D Vs Time')
        xlabel('Time [s]')
        ylabel('L/D [-]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.peso*ones(1,length(datos(i).segmento.tiempo)))
        title('Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Mass [kg]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,90*ones(1,length(datos(i).segmento.tiempo)))
        title('gamma Vs Time')
        xlabel('Time [s]')
        ylabel('gamma [º]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.veloc_vertical)
        title('Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.altura)
        title('Height Vs Time')
        xlabel('Time [s]')
        ylabel('Height [m]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,zeros(1,length(datos(i).segmento.tiempo)))
        title('CL Vs Time')
        xlabel('Time [s]')
        ylabel('CL [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,ones(1,length(datos(i).segmento.tiempo))*datos(i).segmento.CD)
        title('CD Vs Time')
        xlabel('Time [s]')
        ylabel('CD [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(zeros(length(datos(i).segmento.tiempo))+distance_zero,datos(i).segmento.altura)
        title('Mission profile')
        xlabel('Distance [m]')
        ylabel('Height [m]')
        hold on

        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.energiav)
        title('Energy Vs Time')
        xlabel('Time [s]')
        ylabel('Energy [J]')
        hold on
    end
    if strcmp(datos(i).nombre,'Cruise') == 1
        Fig = Fig + 1;
        Fig0 = Fig;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.mbateria*ones(1,length(datos(i).segmento.tiempo)))
        title('Batteries Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Batteries Mass [kg]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.distancia+distance_zero)
        title('Distance Vs Time')
        xlabel('Time [s]')
        ylabel('Distance [m]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.palanca)
        title('Throttle Vs Time')
        xlabel('Time [s]')
        ylabel('Throttle [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.empuje)
        title('Thrust Vs Time')
        xlabel('Time [s]')
        ylabel('Thrust [N]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,zeros(1,length(datos(i).segmento.tiempo)))
        title('Vertical Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.L_D)
        title('L/D Vs Time')
        xlabel('Time [s]')
        ylabel('L/D [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.peso)
        title('Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Mass [kg]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,zeros(1,length(datos(i).segmento.tiempo)))
        title('gamma Vs Time')
        xlabel('Time [s]')
        ylabel('gamma [º]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.velocidad*ones(1,length(datos(i).segmento.tiempo)))
        title('Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.altura*ones(1,length(datos(i).segmento.tiempo)))
        title('Height Vs Time')
        xlabel('Time [s]')
        ylabel('Height [m]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.CL)
        title('CL Vs Time')
        xlabel('Time [s]')
        ylabel('CL [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.CD)
        title('CD Vs Time')
        xlabel('Time [s]')
        ylabel('CD [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.distancia+distance_zero,datos(i).segmento.altura*ones(length(datos(i).segmento.tiempo)))
        title('Mission profile')
        xlabel('Distance [m]')
        ylabel('Height [m]')
        hold on
        
%         Fig = Fig + 1;
%         figure(Fig)
%         plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.energiav)
%         title('Energy Vs Time')
%         xlabel('Time [s]')
%         ylabel('Energy [J]')
%         hold on
    end
    if strcmp(datos(i).nombre,'Descent') == 1
        Fig = Fig + 1;
        Fig0 = Fig;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.mbaterias*ones(1,length(datos(i).segmento.tiempo)))
        title('Batteries Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Batteries Mass [kg]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.distancia+distance_zero)
        title('Distance Vs Time')
        xlabel('Time [s]')
        ylabel('Distance [m]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.palanca)
        title('Throttle Vs Time')
        xlabel('Time [s]')
        ylabel('Throttle [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.empuje)
        title('Thrust Vs Time')
        xlabel('Time [s]')
        ylabel('Thrust [N]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.veloc_vertical)
        title('Vertical Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.L_D)
        title('L/D Vs Time')
        xlabel('Time [s]')
        ylabel('L/D [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.peso)
        title('Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Mass [kg]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.gamma*180/pi*ones(length(datos(i).segmento.tiempo)))
        title('gamma Vs Time')
        xlabel('Time [s]')
        ylabel('gamma [º]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.velocidad)
        title('Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.altura)
        title('Height Vs Time')
        xlabel('Time [s]')
        ylabel('Height [m]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.CL)
        title('CL Vs Time')
        xlabel('Time [s]')
        ylabel('CL [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.CD)
        title('CD Vs Time')
        xlabel('Time [s]')
        ylabel('CD [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.distancia+distance_zero,datos(i).segmento.altura)
        title('Mission profile')
        xlabel('Distance [m]')
        ylabel('Height [m]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.energiav)
        title('Energy Vs Time')
        xlabel('Time [s]')
        ylabel('Energy [J]')
        hold on
    end
    if strcmp(datos(i).nombre,'VTOL Descent') == 1
        Fig = Fig + 1;
        Fig0 = Fig;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.mbaterias*ones(length(datos(i).segmento.tiempo),1))
        title('Batteries Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Batteries Mass [kg]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,zeros(length(datos(i).segmento.tiempo))+distance_zero)
        title('Distance Vs Time')
        xlabel('Time [s]')
        ylabel('Distance [m]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.palanca)
        title('Throttle Vs Time')
        xlabel('Time [s]')
        ylabel('Throttle [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.empuje)
        title('Thrust Vs Time')
        xlabel('Time [s]')
        ylabel('Thrust [N]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.veloc_vertical)
        title('Vertical Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Vertical Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,zeros(1,length(datos(i).segmento.tiempo)))
        title('L/D Vs Time')
        xlabel('Time [s]')
        ylabel('L/D [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.peso*ones(length(datos(i).segmento.tiempo)))
        title('Mass Vs Time')
        xlabel('Time [s]')
        ylabel('Mass [kg]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,-90*ones(1,length(datos(i).segmento.tiempo)))
        title('gamma Vs Time')
        xlabel('Time [s]')
        ylabel('gamma [º]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.veloc_vertical)
        title('Velocity Vs Time')
        xlabel('Time [s]')
        ylabel('Velocity [m/s]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.altura)
        title('Height Vs Time')
        xlabel('Time [s]')
        ylabel('Height [m]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,zeros(1,length(datos(i).segmento.tiempo)))
        title('CL Vs Time')
        xlabel('Time [s]')
        ylabel('CL [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,ones(1,length(datos(i).segmento.tiempo))*datos(i).segmento.CD)
        title('CD Vs Time')
        xlabel('Time [s]')
        ylabel('CD [-]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(zeros(length(datos(i).segmento.tiempo))+distance_zero,datos(i).segmento.altura)
        title('Mission profile')
        xlabel('Distance [m]')
        ylabel('Height [m]')
        hold on
        
        Fig = Fig + 1;
        figure(Fig)
        plot(datos(i).segmento.tiempo+time_zero,datos(i).segmento.energiav)
        title('Energy Vs Time')
        xlabel('Time [s]')
        ylabel('Energy [J]')
        hold on
    end
    time_zero = time_zero+datos(i).segmento.tiempo(end);
    distance_zero = distance_zero+datos(i).segmento.distancia(end);
end

% Reset of Plot for the legend
Fig0 = Fig-1;
for i=1:14 
    Fig = Fig0 + 1;
    figure(Fig)
%     figure(i)
    grid on
    legend(datos(1:(end-1)).nombre)
end