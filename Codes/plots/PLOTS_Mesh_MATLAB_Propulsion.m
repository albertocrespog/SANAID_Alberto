function PLOTS_Mesh_MATLAB_Propulsion(OUTPUT_read_XLSX,PLOTTING_UAV,COLOR_scheme,AC_type,Engine_loc,n_eng)

%% Defines colors RGB
color_fus = COLOR_scheme.color_fus;
color_w1 = COLOR_scheme.color_w1;
color_HTP = COLOR_scheme.color_HTP;
color_vee = COLOR_scheme.color_vee;
color_vee2 = COLOR_scheme.color_vee2;
color_can = COLOR_scheme.color_can;
color_prop = COLOR_scheme.color_prop;
color_vtp = COLOR_scheme.color_vtp;
color_eng = COLOR_scheme.color_eng;
color_nac = COLOR_scheme.color_nac;
color_ac = COLOR_scheme.color_ac;

% Plots Fuselage
% Engine Configuration
% Engine location
switch Engine_loc    
    case 1 % Engine_loc = 1 - under wings
        if n_eng < 2
            disp('For engines at the wing tips 2 engines need to be defined. Please do correct Input Data');
            pause
        end
       if n_eng == 2
           % Engine
           plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
        elseif n_eng == 4
            % Engine and Nacelle first pair
           % Engine
           plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
           % Engine and Nacelle second pair
           % Engine
           plot_2engine2(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_2nacelle2(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_2prop2(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop);
        end
    case 2 % Engine_loc = 2 - fuselage front

           % Engine
           plot_1engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_1nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_1prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)

    case 3 % Engine_loc = 3 - fuselage rear

           % Engine
           plot_1engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_1nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_1prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)

    case 4 % Engine_loc = 4 - wingtips
        if n_eng < 2
            disp('For engines at the wing tips 2 engines need to ve defined. Please do correct Input Data');
            pause
        end
        % Engine
        plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
        % Nacelle
        plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
        % Propeller disk
        plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)

    case 5 % Engine_loc = 5 - wingtips 2 pairs
        if n_eng < 4
            disp('For engines at the wing tips 2 engines need to be defined. Please do correct Input Data');
            pause
        end
        if n_eng == 2
            % Engine
            plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
            % Nacelle
            plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
            % Propeller disk
            plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
        elseif n_eng == 4
            % Engine and Nacelle first pair
            % Engine
            plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
            % Nacelle
            plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
            % Propeller disk
            plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
            % Engine and Nacelle second pair
            % Engine
            plot_2engine2(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
            % Nacelle
            plot_2nacelle2(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
            % Propeller disk
            plot_2prop2(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop);
        end
end