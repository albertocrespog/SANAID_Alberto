function [Body_Geo,meshData] = Read_Fuselage_Data(Geo_tier,XFLR5_DATA,ESCALADO,XFLR5_file,STL_files,OUTPUT_read_XLSX)
% (Geo_tier,OUTPUT_read_XLSX,XFLR5_DATA,XFLR5_file,STL_file)  % Defines Fuselage DATA

% Retieving data
SF = OUTPUT_read_XLSX.AC_Data_flags.SF;
CASE_fuse = OUTPUT_read_XLSX.AC_Data_flags.CASE_fuse;
ESCALADO = OUTPUT_read_XLSX.AC_Data_flags.ESCALADO;
% XFLR5_DATA = OUTPUT_read_XLSX.Fuselage_flags.XFLR5_DATA;
% XFLR5_file = OUTPUT_read_XLSX.Fuselage_flags.XFLR5_file;
% STL_file = OUTPUT_read_XLSX.Fuselage_flags.STL_file;

% XFLR5_file = strcat(OUTPUT_read_XLSX.Fuselage_flags.XFLR5_file);
% STL_files.STL_ac =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_ac);
% STL_files.STL_fus =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_fus);
% STL_files.STL_wing =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_wing);
% STL_files.STL_canard =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_canard);
% STL_files.STL_HTP =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_HTP);
% STL_files.STL_VTP =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_VTP);
% STL_files.STL_Vee =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_Vee);
% STL_files.STL_engine =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_engine);
% STL_files.STL_nacelle =  strcat(OUTPUT_read_XLSX.Fuselage_flags.STL_nacelle);
STL_ac = STL_files.STL_ac;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XFLR5_file = strcat('fuse_pepi_desplazado.txt');
% XFLR5_file = strcat('Fuselage_DATA_scaled.txt');
% STL_file = strcat('Fus_ProVant_v4.0.stl');

%--------------------------- CALCULATES GEOMETRY DATA ---------------------------------
% Data to read XFLR5 files
nSections = XFLR5_DATA.nSections;
nPoints = XFLR5_DATA.nPoints;
lecture_Geo = XFLR5_DATA.lecture_Geo;
STL_PLOT = XFLR5_DATA.STL_PLOT;

% % Scale from fuselage in XFLR5 to CAD
% l_fus1 = 0.8799
% l_fus2 = Geo_tier.l_fus
% Ratio_SCALE = l_fus2/l_fus1
% pause
% Reads to extract length of fuselage
% assumes Ratio_SCALE = 1
Ratio_SCALE_l_fus_tmp = 1;
Ratio_SCALE_w_fus_tmp = 1;
Ratio_SCALE_h_fus_tmp = 1;
Ratio_SCALE_tmp.Ratio_SCALE_l_fus = Ratio_SCALE_l_fus_tmp;
Ratio_SCALE_tmp.Ratio_SCALE_w_fus = Ratio_SCALE_w_fus_tmp;
Ratio_SCALE_tmp.Ratio_SCALE_h_fus = Ratio_SCALE_h_fus_tmp;

CASE_fuse = 3;
switch CASE_fuse
    case 1 % STL of Fuselaje
        
    case 2
        % sTL of complete aircraft
        
        [tri, fileform, A, S] = stlread(STL_ac)
        
    case 3 % Read File from XFLR5
        
        Ratio_SCALE_l_fus = Ratio_SCALE_tmp.Ratio_SCALE_l_fus;
        Ratio_SCALE_w_fus = Ratio_SCALE_tmp.Ratio_SCALE_w_fus;
        Ratio_SCALE_h_fus = Ratio_SCALE_tmp.Ratio_SCALE_h_fus;

        n_ele_geo = nPoints; % elementos que se han de leer
        n_sections_body = nSections; % Elemtos del fuselaje

        for i=1:n_sections_body
            lecture_Geo = lecture_Geo;
            % Extracción de datos VLM
            file_id=fopen(XFLR5_file);
            aux = textscan(file_id,'%n%n%n','Headerlines',lecture_Geo);
            %      size(aux{i})
            %      pause
            Data_P(i).raw_data= aux;
            Data_P(i).label = XFLR5_file;
            Data_P(i).x = aux{1}*Ratio_SCALE_l_fus;
            Data_P(i).y = aux{2}*Ratio_SCALE_w_fus;
            Data_P(i).z = aux{3}*Ratio_SCALE_h_fus;
            x(:,i)      = aux{1}*Ratio_SCALE_l_fus;
            y(:,i)      = aux{2}*Ratio_SCALE_w_fus;
            z(:,i)      = aux{3}*Ratio_SCALE_h_fus;
            x_max(i) = max(aux{1}*Ratio_SCALE_l_fus);
            y_max(i) = max(aux{2}*Ratio_SCALE_w_fus);
            z_max(i) = max(aux{3}*Ratio_SCALE_h_fus);
            fclose(file_id);
            lecture_Geo = lecture_Geo + n_ele_geo  + 2;
        end      
        
%         model = createpde(3);
%         gm = importGeometry(model,STL_ac);
%         pdegplot(model,'FaceLabels','on');
%         model.Geometry;

%         [tri, fileform, A, S] = stlread(STL_ac);
        
%         figure(1)
%         subplot(2,3,1)
%         plot(tri(:,1),tri(:,2))
%         grid on
%         subplot(2,3,2)
%         plot(tri(:,1),tri(:,3))
%         grid on
%         subplot(2,3,3)
%         plot(tri(:,2),tri(:,3))
%         grid on
%         
%         figure(2)
%         subplot(1,3,1)
%         plot(x,y)
%         grid on
%         subplot(1,3,2)
%         plot(x,z)
%         grid on
%         subplot(1,3,3)
%         plot(y,z)
%         grid on
        
%         modela = stlread(STL_ac);
%         modelb = stl2matlab(STL_ac);

        plotData{1}  = [x(:,end:(-1):1),x];
        plotData{2}  = [-y(:,end:(-1):1), y];
        plotData{3} = [z(:,end:(-1):1), z];

        size(plotData{1})
        size(plotData{2})
        size(plotData{3})
        pause
        
        figure(2)
        subplot(1,3,1)
        plot(plotData{1},plotData{2})
        grid on
        subplot(1,3,2)
        plot(plotData{1},plotData{3})
        grid on
        subplot(1,3,3)
        plot(plotData{2},plotData{3})
        grid on

        figure(3)
        mesh(plotData{1},plotData{2},plotData{3});
        
        % Loads STL file into a mesh
        model = stl2matlab(STL_ac);
        ScaleFactor = 1/1000; % CAD designs are in mm, and changing to m
        plotData{1} = model{1}*ScaleFactor;
        plotData{2} = model{2}*ScaleFactor;
        plotData{3} = model{3}*ScaleFactor;
        
        size(plotData{1})
        size(plotData{2})
        size(plotData{3})
        pause
        
        figure(4)
        subplot(1,3,1)
        plot(plotData{1},plotData{2})
        grid on
        subplot(1,3,2)
        plot(plotData{1},plotData{3})
        grid on
        subplot(1,3,3)
        plot(plotData{2},plotData{3})
        grid on

        figure(5)
        mesh(plotData{1},plotData{2},plotData{3});
        pause
        
        % Rewrites data in the format of the XFLR5 lecture files in order to be used the
        x_rec(:,1) = plotData{1}(:,3);
        x_rec(:,2) = plotData{1}(:,4);
        y_rec(:,1) = plotData{2}(:,3);
        y_rec(:,2) = plotData{2}(:,4);
        z_rec(:,1) = plotData{3}(:,3);
        z_rec(:,2) = plotData{3}(:,4);
        pause
        
        figure(6)
        subplot(1,3,1)
        plot(plotData{1},plotData{2})
        grid on
        subplot(1,3,2)
        plot(plotData{1},plotData{3})
        grid on
        subplot(1,3,3)
        plot(plotData{2},plotData{3})
        grid on
                              
end

size(plotData)
size(x)
pause





[Body_Geo_XFLR5,meshData_XFLR5] = Calc_Body_Geometry_Feb2021(XFLR5_file,STL_file,Geo_tier, nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE_tmp,SF,plotData);
% Scale from fuselage in XFLR5 to CAD
l_fus1 = Body_Geo_XFLR5.l_fus; % Length of original XFLR5
w_fus1 = Body_Geo_XFLR5.w_fus; % width
h_fus1 = Body_Geo_XFLR5.h_fus; % height

% Desired Lenth from input data
l_fus2 = Geo_tier.l_fus; % Desired Length
w_fus2 = Geo_tier.w_fus; % width
h_fus2 = Geo_tier.h_fus; % height
% Determine the type of scaling
% ESCALADO = 0; % no scaling Ratio =1
% ESCALADO = 1; % scaling acording to desired for each 3 axix
% ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
if ESCALADO == 0
    % Maintaining the ratio
    Ratio_SCALE_l_fus = 1;
    Ratio_SCALE_w_fus = 1;
    Ratio_SCALE_h_fus = 1;
elseif ESCALADO == 1
    % Applies scaling factor according to desired geometry
    Ratio_SCALE_l_fus = l_fus2/l_fus1;
    Ratio_SCALE_w_fus = w_fus2/w_fus1;
    Ratio_SCALE_h_fus = h_fus2/h_fus1;
elseif ESCALADO == 2
    % Maintaining the ratio
    Ratio_SCALE_l_fus = l_fus2/l_fus1;
    Ratio_SCALE_w_fus = l_fus2/l_fus1;
    Ratio_SCALE_h_fus = l_fus2/l_fus1;
end

% Variable that define the scaling factor
Ratio_SCALE.Ratio_SCALE_l_fus = Ratio_SCALE_l_fus;
Ratio_SCALE.Ratio_SCALE_w_fus = Ratio_SCALE_w_fus;
Ratio_SCALE.Ratio_SCALE_h_fus = Ratio_SCALE_h_fus;

% Being used only if STL model used to determine geometry of fuselage
% Nth = size(plotData{1});
% STLnSections = Nth(1);% number of sectins (x coordinate) (eliminated 1 for convergence)
% STLnPoints = Nth(2); % number of elements per section (y and z coordinates)
% Determines fuselage information
[Body_Geo,meshData] = Calc_Body_Geometry_Dec2018(XFLR5_file,STL_file,Geo_tier, nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE,SF);
% Identifies the fuselage mode, between data from XFLR5 and STL files (CAD)
if STL_PLOT == 1
    switch CASE_fuse
        case 1 % STL of Fuselaje
            model = stl2matlab('ProVant.stl');
            ScaleFactor = SF*1/1000; % CAD designs are in mm, and changing to m
            plotData{1} = model{1}*ScaleFactor + Geo_tier.x_offset_CAD;
            plotData{2} = model{2}*ScaleFactor;
            plotData{3} = model{3}*ScaleFactor + Geo_tier.z_offset_CAD;
            meshData = plotData;
        case 2
            % sTL of complete aircraft
            model = stl2matlab('VTOL_Single.stl');
            model = stl2matlab('VTP_WIG.stl');
            model = stl2matlab('VTP_WIG.stl');            
            ScaleFactor = SF; % CAD designs are in mm, and changing to m
            plotData{1} = model{1}*ScaleFactor + Geo_tier.x_offset_CAD;
            plotData{2} = model{2}*ScaleFactor;
            plotData{3} = model{3}*ScaleFactor + Geo_tier.z_offset_CAD;
            meshData = plotData;

            model = stl2matlab('HTP_WIG.stl');            
            ScaleFactor = SF; % CAD designs are in mm, and changing to m
            plotData{1} = model{1}*ScaleFactor + Geo_tier.x_offset_CAD;
            plotData{2} = model{2}*ScaleFactor;
            plotData{3} = model{3}*ScaleFactor + Geo_tier.z_offset_CAD;
            meshData = plotData;
        case 3
            meshData = meshData;    
    end    
end