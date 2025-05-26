function Conduct_Prop_Optimzation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Determine Minim Prop Diameter %%%%%%%%%%%%%%%
    %%
    %% Lets user select if wants to conduct study for the propellers
    answer = questdlg('Would you like to conduct PROPELLER''s STUDY?', ...
        'Prop Limits', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            %         disp([answer ' coming right up.'])
            Study_Prop_Limits = 1;
        case 'No'
            %         disp([answer ' coming right up.'])
            Study_Prop_Limits = 0;
    end
    
    if Study_Prop_Limits == 1
        % Determination
        switch case_AC
            case 1 % case_AC = 1 - EMERGENTIA 1:1
                pp = 0.90; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 1; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1;
            case 2 % case_AC = 1 - EMERGENTIA 1:2
                pp = 0.80; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 0.01; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1;
            case 3 % case_AC = 3 - PEPIÑO XXL
                pp = 0.95; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 0;
                Study_horizontal = 0;
            case 6 % case_AC = 6 - CERVERA
                pp = 1; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 1; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1    ;
            case 12 % case_AC = 12 - ALO
                pp = 1; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 1; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1    ;
            case 13 % case_AC = 13 - ALO Fuel Cell
                pp = 1; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                %% Vertical Analysis
                % Range of vertical speeds
                Vv_min = 1; % Vertical climb speed
                Vv_max = 10; % Vertical climb speed
                N_Vv = 100; % Number of points
                Vv_vect = linspace(Vv_min, Vv_max,N_Vv);
                
                % Range of horizontal speeds
                Vh_min = Performance.V_min; % Vertical climb speed
                Vh_max = Vh_min*1.2; % Vertical climb speed
                N_Vh = 100; % Number of points
                Vh_vect = linspace(Vh_min, Vh_max,N_Vh);
                
                % Flag to determine if vertical flight & horizontal flight
                % prop limits are calculated
                Study_vertical = 1;
                Study_horizontal = 1    ;
        end
        % uses fsolve to find associated rpms
        solve_w_fero = 0;
        
        %% Conducts study of sensitivity analysis of the Performance
        %% Lets user select if wants to conduct study for the propellers
        answer = questdlg('Would you like to conduct Propeller Performance?', ...
            'Prop Limits', ...
            'Yes','No','No');
        % Handle response
        switch answer
            case 'Yes'
                % Selects if 3D Plots are shown
                answer1 = questdlg('Would you like to conduct horizontal flight limits?', ...
                    'Horizontal limits', ...
                    'Yes','No','No');
                % Handle response
                switch answer1
                    case 'Yes'
                        Study_vertical = 1;
                    case 'No'
                        Study_vertical = 0;
                end
                
                % Selects if Contour Plots are shown
                answer2 = questdlg('Would you like to conduct vertical flight limits?', ...
                    '3D plots', ...
                    'Yes','No','No');
                % Handle response
                switch answer2
                    case 'Yes'
                        Study_horizontal = 1;
                    case 'No'
                        Study_horizontal = 0;
                end
                flight_minits_study = 1;
            case 'No'
                %         disp([answer ' coming right up.'])
                flight_minits_study = 0;
        end
        
        if flight_minits_study == 1
            % Flag to determine if vertical flight & horizontal flight
            % prop limits are calculated
            if Study_vertical == 1;
                % Selects to plot Get Vertical Flight Limits that determines the
                % minimum prop diameter for the different maneuvers
                PLOT_Get_Vertical_Flight_Limits = 1;
                [Fig] = Get_Vertical_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Performance,Geo_tier,...
                    AC_CONFIGURATION,N_Vv,Vv_vect,Fig,PLOT_Get_Vertical_Flight_Limits,solve_w_fero,Plot_Options);
                Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
                disp(Warning)
                pause
            end
            
            % Flag to determine if vertical flight & horizontal flight
            % prop limits are calculated
            if Study_horizontal == 1
                % Selects to plot Get Vertical Flight Limits that determines the
                % minimum prop diameter for the different maneuvers
                PLOT_Get_Horizontal_Flight_Limits = 1;
                [Fig] = Get_Horizontal_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Performance,Geo_tier,...
                    AC_CONFIGURATION,N_Vh,Vh_vect,Fig,PLOT_Get_Horizontal_Flight_Limits,solve_w_fero,Plot_Options);
                Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
                disp(Warning)
                pause
            end
        end
    end
    
    %% Conducts study of sensitivity analysis of the Performance
    %% Lets user select if wants to conduct study for the propellers
    answer = questdlg('Would you like to conduct Sensitivity Study To Determine the Propeller Diameter?', ...
        'Prop Limits', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            
            % Selects if 3D Plots are shown
            answer1 = questdlg('Would you like to show Hover Performence vs. D_{prop}?', ...
                'Hover plots', ...
                'Yes','No','No');
            % Handle response
            switch answer1
                case 'Yes'
                    Hover_PLOTS = 1;
                case 'No'
                    Hover_PLOTS = 0;
            end
            
            % Selects if 3D Plots are shown
            answer2 = questdlg('Would you like to show 3D plots?', ...
                '3D plots', ...
                'Yes','No','No');
            % Handle response
            switch answer2
                case 'Yes'
                    plots_3d_sensitivity = 1;
                case 'No'
                    plots_3d_sensitivity = 0;
            end
            
            % Selects if Contour Plots are shown
            answer3 = questdlg('Would you like to show contour plots?', ...
                '3D plots', ...
                'Yes','No','No');
            % Handle response
            switch answer3
                case 'Yes'
                    plots_contour_sensitivity = 1;
                case 'No'
                    plots_contour_sensitivity = 0;
            end
            
            % Selects the minimum Energy plots
            answer4 = questdlg('Would you like to show plots searching for minimum Energy?', ...
                'Minimum Energy', ...
                'Yes','No','No');
            % Handle response
            switch answer4
                case 'Yes'
                    special_PLOTS = 1; % represents the plots searching for minimum Energy
                case 'No'
                    special_PLOTS = 0; % represents the plots searching for minimum Energy
            end
            
            % Selects the minimum Energy plots
            answer5 = questdlg('Would you like to show plots for Performance Study vs battery mass variation for hover flight?', ...
                'Mass Variation', ...
                'Yes','No','No');
            % Handle response
            switch answer5
                case 'Yes'
                    variation_mass = 1; % represents the plots searching for minimum Energy
                case 'No'
                    variation_mass = 0; % represents the plots searching for minimum Energy
            end
            
            % Selects if Intersection curves between Vertical and
            % horizontal Flight are shown
            
            answer3 =questdlg ('Would you like to perform the Diameter optimization?',...
                '3D plots',...
                'Yes','No','No');
            %Handle response
            switch answer3
                case 'Yes'
                    diameter_optim=1;
                case 'No'
                    diameter_optim=0;
            end
            
            PERFORMANCE_sensitivity = 1;
            % Saves selection of plotting options
            PLOTS_Prop_Performance.Hover_PLOTS = Hover_PLOTS;
            PLOTS_Prop_Performance.PERFORMANCE_sensitivity = PERFORMANCE_sensitivity;
            PLOTS_Prop_Performance.plots_contour_sensitivity = plots_contour_sensitivity;
            PLOTS_Prop_Performance.plots_3d_sensitivity = plots_3d_sensitivity;
            PLOTS_Prop_Performance.special_PLOTS = special_PLOTS;
            PLOTS_Prop_Performance.variation_mass = variation_mass;
            
        case 'No'
            %         disp([answer ' coming right up.'])
            PERFORMANCE_sensitivity = 0;
            % Saves selection of plotting options
            PLOTS_Prop_Performance.Hover_PLOTS = 0;
            PLOTS_Prop_Performance.PERFORMANCE_sensitivity = 0;
            PLOTS_Prop_Performance.plots_contour_sensitivity = 0;
            PLOTS_Prop_Performance.plots_3d_sensitivity = 0;
            PLOTS_Prop_Performance.special_PLOTS = 0;
            
    end
    
    PLOTS_SENSITIVITY_VERTICAL = 1;
    PLOTS_SENSITIVITY_HORIZONTAL = 1;
    % Uses Fzero which takes longer time
    solve_w_fero = 0;
    
    if PERFORMANCE_sensitivity == 1
        % Selects the Engine-Prop configuration according to the Limit Study
        N_prop = 500; % number of prop elements to analyze between the limits
        pp_D_prop_min = 0.75;
        pp_D_prop_max = 1.5;
        N_contour_lines = 25; % number of contour lines
        
        switch case_AC
            case 1 % case_AC = 1 - EMERGENTIA 1:1
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 2 % case_AC = 1 - EMERGENTIA 1:2
                D_prop = 22*2.54/100;
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 3 % case_AC = 3 - PEPIÑO XXL
                D_prop = 32*2.54/100;
                % Actualizes Prop Data
                SF_prop = 0.85;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Thin Electric (E) Propellers Maximum RPM=150,000/prop diameter (inches)
                % Flag to determine if vertical flight
                sensitivity_vertical = 0;
                sensitivity_horizontal = 1;
            case 6 % case_AC = 6 - CERVERA
                D_prop = 30*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = pp; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 9 % case_AC = 1 - EMERGENTIA Wind Tunnel w2
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel w1 (sweep)
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 11 % case_AC = 11 - EMERGENTIA Manufactured
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = 0.95;
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 12 % case_AC = 12 - ALO
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = pp; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
            case 13 % case_AC = 13 - ALO Fuel Cell
                D_prop = 28*2.54/100; % Prop diameter around the center for the sensitivity study ()
                % Actualizes Prop Data
                SF_prop = pp; % Percentaje of throttle
                RPMMAX_APC = 150000; % From APC https://www.apcprop.com/technical-information/rpm-limits/
                % Flag to determine if vertical flight
                sensitivity_vertical = 1;
                sensitivity_horizontal = 1;
        end
        
        % Loop that conducts sensitivity study for vertical flight
        if sensitivity_vertical == 1;
            % Sensitivity analysis for the Vertical Performance
            [Xv,Yv,Zv_study,Fig] = Sensitivity_Study_Vertical_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,...
                Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_VERTICAL,Performance,solve_w_fero,...
                N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines,Vv_min,PLOTS_Prop_Performance,Prop_selection,prefix,Plot_Options);
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            query_properties = 1;
            while query_properties == 1
                % Lets user select if wants to conduct study
                answer_prop = questdlg('Would you like to calculate properties for a given prop diameter and vertical speed?', ...
                    'Prop Limits', ...
                    'Yes','No','Yes');
                % Handle response
                switch answer_prop
                    case 'Yes'
                        prompt = {'Enter the prop diameter (inches):','Enter vertical VTOl speed (m/s):'};
                        dlgtitle = 'Performance properties';
                        dims = [1 35];
                        definput = {'30','5'};
                        answer_sv = inputdlg(prompt,dlgtitle,dims,definput);
                        
                        D_propin_q = str2num(answer_sv{1});
                        Vv_q = str2num(answer_sv{2});
                        
                        z_Tv = interp2(Xv,Yv,Zv_study.Z_T,D_propin_q,Vv_q,'cubic');
                        z_Pev = interp2(Xv,Yv,Zv_study.Z_Pe,D_propin_q,Vv_q,'cubic');
                        %                         z_Pe_CHv = interp2(Xv,Yv,Zv_study.Z_Pe_CH,D_propin_q,Vv_q,'cubic');
                        z_Energyv = interp2(Xv,Yv,Zv_study.Z_Energy,D_propin_q,Vv_q,'cubic');
                        z_Energy_CHv = interp2(Xv,Yv,Zv_study.Z_Energy_CH,D_propin_q,Vv_q,'cubic');
                        z_Qv = interp2(Xv,Yv,Zv_study.Z_Q,D_propin_q,Vv_q,'cubic');
                        z_deltav = interp2(Xv,Yv,Zv_study.Z_delta,D_propin_q,Vv_q,'cubic');
                        z_RPMv = interp2(Xv,Yv,Zv_study.Z_RPM,D_propin_q,Vv_q,'cubic');
                        z_ethav = interp2(Xv,Yv,Zv_study.Z_etha,D_propin_q,Vv_q,'cubic');
                        z_Battery_massv = interp2(Xv,Yv,Zv_study.Z_Battery_mass,D_propin_q,Vv_q,'cubic');
                        z_Battery_mass_CHv = interp2(Xv,Yv,Zv_study.Z_Battery_mass_CH,D_propin_q,Vv_q,'cubic');
                        z_Pe_CHv = 0;
                        
                        uiwait(msgbox({['Total Thrust (climb) = ',num2str(z_Tv),' (N)']...
                            ['Total Electric Power (climb) = ',num2str(z_Pev),' (W)']...
                            ['Total Electric Power (climb & hover) = ',num2str(z_Pe_CHv),' (W)']...
                            ['Total Energy (climb) = ',num2str(z_Energyv),' (kW-h)']...
                            ['Total Energy (climb & hover) = ',num2str(z_Energy_CHv),' (kW-h)']...
                            ['Torque = ',num2str(z_Qv),' (Nm)']...
                            ['Throttle = ',num2str(z_deltav),' (%)']...
                            ['RPM = ',num2str(z_RPMv),' (RPM)']...
                            ['\eta = ',num2str(z_ethav)]...
                            ['m_{bat} (climb) = ',num2str(z_Battery_massv),' (kg)']...
                            ['m_{bat} (climb & ) hover= ',num2str(z_Battery_mass_CHv),' (kg)']}));
                        query_properties = 1;
                    case 'No'
                        query_properties = 0;
                end
            end
        end
        
        % Loop that conducts sensitivity study for horizontal flight
        if sensitivity_horizontal == 1
            [Xh,Yh,Zh_study,Fig] = Sensitivity_Study_Horizontal_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,...
                Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_VERTICAL,Performance,solve_w_fero,...
                N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines,PLOTS_Prop_Performance,Prop_selection,prefix,Plot_Options);
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            query_properties = 1;
            Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
            disp(Warning)
            pause

            while query_properties == 1
                % Lets user select if wants to conduct study
                answer_prop = questdlg('Would you like to calculate properties for a given prop diameter and horizontal speed?', ...
                    'Prop Limits', ...
                    'Yes','No','Yes');
                % Handle response
                switch answer_prop
                    case 'Yes'
                        prompt = {'Enter the prop diameter (inches):','Enter horizontal speed (m/s):'};
                        dlgtitle = 'Performance properties';
                        dims = [1 35];
                        definput = {'30','30'};
                        answer_sv = inputdlg(prompt,dlgtitle,dims,definput);
                        
                        D_propin_q = str2num(answer_sv{1});
                        Vv_q = str2num(answer_sv{2});
                        
                        z_Th = interp2(Xh,Yh,Zh_study.Z_T,D_propin_q,Vv_q,'cubic');
                        z_Peh = interp2(Xh,Yh,Zh_study.Z_Pe,D_propin_q,Vv_q,'cubic');
                        z_Energyh = interp2(Xh,Yh,Zh_study.Z_Energy,D_propin_q,Vv_q,'cubic');
                        z_Qh = interp2(Xh,Yh,Zh_study.Z_Q,D_propin_q,Vv_q,'cubic');
                        z_deltah = interp2(Xh,Yh,Zh_study.Z_delta,D_propin_q,Vv_q,'cubic');
                        z_RPMh = interp2(Xh,Yh,Zh_study.Z_RPM,D_propin_q,Vv_q,'cubic');
                        z_ethah = interp2(Xh,Yh,Zh_study.Z_etha,D_propin_q,Vv_q,'cubic');
                        z_Battery_massh = interp2(Xh,Yh,Zh_study.Z_Battery_mass,D_propin_q,Vv_q,'cubic');
                        
                        uiwait(msgbox({['Total Thrust = ',num2str(z_Th),' (N)']...
                            ['Total Electric Power = ',num2str(z_Peh),' (W)']...
                            ['Total Energy = ',num2str(z_Energyh),' (kW-h)']...
                            ['Torque = ',num2str(z_Qh),' (Nm)']...
                            ['Throttle = ',num2str(z_deltah),' (%)']...
                            ['RPM = ',num2str(z_RPMh),' (RPM)']...
                            ['\eta = ',num2str(z_ethah)]...
                            ['m_{bat} (climb) = ',num2str(z_Battery_massh),' (kg)']}));
                        query_properties = 1;
                    case 'No'
                        query_properties = 0;
                end
            end
        end
        
        if sensitivity_horizontal == 1 && sensitivity_vertical == 1 && diameter_optim== 1
            
            [Vvopt_mc,Vhopt_mc,Dopt_mc] = diameteroptimization(Xh,Yv,Yh,Zh_study.Z_T,Zh_study.Z_Pe,Zh_study.Z_Energy,Zh_study.Z_etha,...
                Zv_study.Z_T,Zv_study.Z_Pe,Zv_study.Z_Energy,Zv_study.Z_etha);
            
        end
        
    end