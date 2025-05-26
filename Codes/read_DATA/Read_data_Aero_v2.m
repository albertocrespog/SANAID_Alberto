%% Reads the data from the txt files gnerated from XFLR5
function [Data_P,casos,prefix,mark_legend,X_OC] = Read_data_Aero(Performance_preliminar,case_AC,OUTPUT_read_XLSX)

%% Selects the file to red aero files
switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        [casos,prefix,mark_legend,X_OC] = read_aero_files_March2022_EMERGENTIA(Performance_preliminar);
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        [casos,prefix,mark_legend,X_OC] = read_aero_files_March2022_EMERGENTIA(Performance_preliminar);
    case 3 % case_AC = 3 - PEPIÑO XXL
        [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_PEPINOXXL(Performance_preliminar);
    case 4 % case_AC = 4 - A400
        [casos,prefix,mark_legend,X_OC] = read_aero_files_A400(Performance_preliminar);
    case 5 % case_AC = 5 - WIGL
        [casos,prefix,mark_legend,X_OC] = read_aero_files_Feb2020_v3_WIG(Performance_preliminar);
    case 6 % case_AC = 6 - MILVUS
        % [casos prefix mark_legend,X_OC] = read_aero_files_MILVUS(Performance_preliminar);
        % [casos prefix mark_legend,X_OC] = read_aero_files_MILVUS_v2(Performance_preliminar);
        [casos prefix mark_legend,X_OC] = read_aero_files_MILVUS_v3(Performance_preliminar);
    case 7 % Existing Aircraft = 7
        [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
    case 8 % case_AC = 8 - TAMIZ
        [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
        [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
        [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
    case 11 % case_AC = 11 - EMERGENTIA Manufacturing
%         [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        [casos,prefix,mark_legend,X_OC] = read_aero_files_March2022_EMERGENTIA(Performance_preliminar);
    case 12 % case_AC = 12 - ALO
        [casos,prefix,mark_legend,X_OC] = read_aero_files_ALO_v1(Performance_preliminar);
    case 13 % case_AC = 13 - ALO Fuel Cell
        [casos,prefix,mark_legend,X_OC] = read_aero_files_ALO_v1(Performance_preliminar);
    case 14 % case_AC = 14 - A320
        [casos,prefix,mark_legend,X_OC] = read_aero_files_A320(Performance_preliminar);
    case 15 % case_AC = 15 - EMERGENTIA Tandem Wing
%         [casos,prefix,mark_legend,X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance_preliminar);
        [casos,prefix,mark_legend,X_OC] = read_aero_files_March2022_EMERGENTIA(Performance_preliminar);       
    case 16 % case_AC = 16 - Tarsis T75
        [casos,prefix,mark_legend,X_OC] = read_aero_files_March2022_TARSIS(Performance_preliminar);       
    case 17 % case_AC = 17 - Tarsis T120
        [casos,prefix,mark_legend,X_OC] = read_aero_files_March2022_TARSIS(Performance_preliminar);       
    case 18 % case_AC = 18 - BAT
        [casos,prefix,mark_legend,X_OC] = read_aero_files_BAT(Performance_preliminar);       
    case 19 % case_AC = 19 - FALCON2000
        [casos,prefix,mark_legend,X_OC] = read_aero_files_FALCON2000(Performance_preliminar);       
    case 20 % case_AC = 20 - SOLARtii T120
        % [casos,prefix,mark_legend,X_OC] = read_aero_files_April2024_SOLARTII(Performance_preliminar);   
        [casos,prefix,mark_legend,X_OC] = read_aero_files_April2024_SOLARTII_v1(Performance_preliminar);   
    case 21 % case_AC = 21 - VANTUS
        [casos,prefix,mark_legend,X_OC] = read_aero_files_June2024_VANTUS_v1(Performance_preliminar);       
    case 22 % case_AC = 22 - Cessna 208
        [casos,prefix,mark_legend,X_OC] = read_aero_files_July2024_Cessna208_v1(Performance_preliminar);       
    case 23 % case_AC = 23 - King Air 350
        [casos,prefix,mark_legend,X_OC] = read_aero_files_July2024_KingAir350_v1(Performance_preliminar);       
    case 24 % case_AC = 24 - Future 1
        [casos,prefix,mark_legend,X_OC] = read_aero_files_Oct2024_MK84_v1(Performance_preliminar);       
    case 25 % case_AC = 25 - Future 2
        [casos,prefix,mark_legend,X_OC] = read_aero_files_Oct2024_MK84_v1(Performance_preliminar);       
    case 26 % case_AC = 26 - Future 3
        [casos,prefix,mark_legend,X_OC] = read_aero_files_Nov2024_HA_v1(Performance_preliminar);       
    case 27 % case_AC = 27 - Future 4
        [casos,prefix,mark_legend,X_OC] = read_aero_files_Nov2024_FR1_v1(Performance_preliminar);       
    case 28 % case_AC = 28 - Future 5
        [casos,prefix,mark_legend,X_OC] = read_aero_files_Dec2024_H2_v1(Performance_preliminar);       
    case 29 % case_AC = 29 - Future 6
        [casos prefix mark_legend,X_OC] = read_aero_files_MILVUS2_v3(Performance_preliminar);       
    case 30 % case_AC = 30 - Future 7
        [casos,prefix,mark_legend,X_OC] = read_aero_files_July2024_Future7_v1(Performance_preliminar);       
end

prefix = OUTPUT_read_XLSX.PLOT_flags.prefix;

read_XFLR5 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.read_XFLR5; %
read_FLOW5 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.read_FLOW5; %

%% Reads files
% Number of cases
% num_cases = length(casos);
if read_XFLR5
    for i =1:length(casos)
        conv_w = 1;
        % Extracción de datos XFLR5
        file_id = fopen(casos{i});
        aux = textscan(file_id,'%n%n%n%n%n%n%n%n%n%n%n%n%n','Headerlines',8);
        % New results 13 entries
        Data_P(i).raw_data= aux;
        Data_P(i).label = casos{i};
        Data_P(i).alpha = aux{1};
        Data_P(i).beta = aux{2};
        Data_P(i).CL    = aux{3}*conv_w;
        Data_P(i).CDi   = aux{4}*conv_w;
        Data_P(i).CDv   = aux{5}*conv_w;
        Data_P(i).CD    = aux{6}*conv_w;
        Data_P(i).CY    = aux{7}*conv_w;
        Data_P(i).Cl    = aux{8}*conv_w;
        Data_P(i).Cm    = aux{9}*conv_w;
        Data_P(i).Cn    = aux{10}*conv_w;
        Data_P(i).Cni   = aux{11}*conv_w;
        Data_P(i).Qinf  = aux{12};
        Data_P(i).XCP   = aux{13};
        fclose(file_id);
    end
end

% Number of cases
num_cases = length(casos);

% Define format specifier for textscan once
formatSpec = repmat('%n', 1, 56); % Adjust based on the expected number of columns
headerLines = 6;

% Preallocate the structure with the expected fields
Data_P = repmat(struct( ...
    'raw_data', [], ...
    'label', '', ...
    'alpha', [], ...
    'beta', [], ...
    'CL', [], 'CD', [], 'CDv', [], 'CDi', [], 'CY', [], ...
    'Cm', [], 'Cmv', [], 'Cmp', [], 'Cl', [], 'Cn', [], 'Cnv', [], 'Cnp', [], ...
    'CLCD', [], 'CL32CD', [], 'inv_sqrt_CL', [], ...
    'L', [], 'D', [], 'Fx_FF', [], 'Fy_FF', [], 'Fz_FF', [], ...
    'Fx_sum', [], 'Fy_sum', [], 'Fz_sum', [], ...
    'Extra_drag', [], 'Fuse_drag', [], 'Cf_Fuse', [], ...
    'Vx', [], 'Vz', [], 'V', [], 'Gamma', [], ...
    'M', [], 'N', [], ...
    'CPx', [], 'CPy', [], 'CPz', [], ...
    'BM', [], 'm', struct('g', struct('Vz', [])), ...
    'Drag_x_V', [], 'Efficiency', [], ...
    'XCp', struct('Cl', []), ...
    'XNP', [], ...
    'Phugoid_Freq', [], 'Phugoid_Damping', [], ...
    'Short_Period_Freq', [], 'Short_Period_Damping_Ratio', [], ...
    'Dutch_Roll_Freq', [], 'Dutch_Roll_Damping', [], ...
    'Roll_Damping', [], 'Spiral_Damping', [], ...
    'Mass', [], 'CoG_x', [], 'CoG_z', [] ...
), num_cases, 1);

%% Reads files
if read_FLOW5% Loop through each case and populate Data_P
    for i = 1:num_cases
        conv_w = 1; % Conversion factor

        % Open file and read data
        file_id = fopen(casos{i});
        aux = textscan(file_id, formatSpec, 'HeaderLines', headerLines);

        % Populate the structure fields
        Data_P(i).raw_data = aux;
        Data_P(i).label = casos{i};

        e=2; % starting line
        Data_P(i).alpha = aux{e}; e = e+1;
        Data_P(i).beta = aux{e}; e = e+1;
        % Dimensionless Forces
        Data_P(i).CL    = aux{e}*conv_w; e = e+1;
        Data_P(i).CD    = aux{e}*conv_w; e = e+1;
        Data_P(i).CDv   = aux{e}*conv_w; e = e+1;
        Data_P(i).CDi   = aux{e}*conv_w; e = e+1;
        Data_P(i).CY    = aux{e}*conv_w; e = e+1;
        % Dimensionless Moments
        Data_P(i).Cm    = aux{e}*conv_w; e = e+1;
        Data_P(i).Cmv    = aux{e}*conv_w; e = e+1;
        Data_P(i).Cmp    = aux{e}*conv_w; e = e+1;
        Data_P(i).Cl    = aux{e}*conv_w; e = e+1;
        Data_P(i).Cn    = aux{e}*conv_w; e = e+1;
        Data_P(i).Cnv   = aux{e}*conv_w; e = e+1;
        Data_P(i).Cnp   = aux{e}*conv_w; e = e+1;
        % Aero
        Data_P(i).CLCD    = aux{e}*conv_w; e = e+1;
        Data_P(i).CL32CD    = aux{e}*conv_w; e = e+1;
        Data_P(i).inv_sqrt_CL    = aux{e}*conv_w; e = e+1;
        % Forces
        Data_P(i).L    = aux{e}; e = e+1;
        Data_P(i).D    = aux{e}; e = e+1;
        Data_P(i).Fx_FF    = aux{e}; e = e+1;
        Data_P(i).Fy_FF    = aux{e}; e = e+1;
        Data_P(i).Fz_FF    = aux{e}; e = e+1;
        Data_P(i).Fx_sum    = aux{e}; e = e+1;
        Data_P(i).Fy_sum    = aux{e}; e = e+1;
        Data_P(i).Fz_sum    = aux{e}; e = e+1;
        Data_P(i).Extra_drag = aux{e}; e = e+1;
        Data_P(i).Fuse_drag = aux{e}; e = e+1;
        Data_P(i).Cf_Fuse = aux{e}; e = e+1;
        % Velocities
        Data_P(i).Vx = aux{e}; e = e+1;
        Data_P(i).Vz = aux{e}; e = e+1;
        Data_P(i).V = aux{e}; e = e+1;
        Data_P(i).Gamma = aux{e}; e = e+1;
        % Moments
        Data_P(i).L = aux{e}; e = e+1;
        Data_P(i).M = aux{e}; e = e+1;
        Data_P(i).N = aux{e}; e = e+1;
        % Pressude Coefficients
        Data_P(i).CPx = aux{e}*conv_w; e = e+1;
        Data_P(i).CPy = aux{e}*conv_w; e = e+1;
        Data_P(i).CPz = aux{e}*conv_w; e = e+1;
        % Bending Moment
        Data_P(i).BM = aux{e}; e = e+1;
        Data_P(i).m.g.Vz = aux{e}; e = e+1;
        Data_P(i).Drag_x_V = aux{e}; e = e+1;
        Data_P(i).Efficiency = aux{e};  e = e+1;
        Data_P(i).XCp.Cl = aux{e}; e = e+1;
        Data_P(i).XNP = aux{e}; e = e+1;
        % Stability
        Data_P(i).Phugoid_Freq = aux{e}; e = e+1;
        Data_P(i).Phugoid_Damping = aux{e}; e = e+1;
        Data_P(i).Short_Period_Freq = aux{e}; e = e+1;
        Data_P(i).Short_Period_Damping_Ratio = aux{e}; e = e+1;
        Data_P(i).Dutch_Roll_Freq = aux{e}; e = e+1;
        Data_P(i).Dutch_Roll_Damping = aux{e}; e = e+1;
        Data_P(i).Roll_Damping = aux{e}; e = e+1;
        Data_P(i).Spiral_Damping = aux{e}; e = e+1;
        Data_P(i).Mass = aux{e}; e = e+1;
        Data_P(i).CoG_x = aux{e}; e = e+1;
        Data_P(i).CoG_z = aux{e}; e = e+1;

        % Close the file
        fclose(file_id);
    end
end

% 
% %% Reads files
% if read_FLOW5
%     % Preallocating data 
%     for i =1:length(casos)
%         conv_w = 1;
%         % Extracción de datos XFLR5
%         file_id = fopen(casos{i});
%         aux = textscan(file_id,'%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n','Headerlines',6);
%         % New results 13 entries
%         Data_P(i).raw_data= aux;
%         Data_P(i).label = casos{i};
%         e=2; % starting line
%         Data_P(i).alpha = aux{e}; e = e+1;
%         Data_P(i).beta = aux{e}; e = e+1;
%         % Dimensionless Forces
%         Data_P(i).CL    = aux{e}*conv_w; e = e+1;
%         Data_P(i).CD    = aux{e}*conv_w; e = e+1;
%         Data_P(i).CDv   = aux{e}*conv_w; e = e+1;
%         Data_P(i).CDi   = aux{e}*conv_w; e = e+1;
%         Data_P(i).CY    = aux{e}*conv_w; e = e+1;
%         % Dimensionless Moments
%         Data_P(i).Cm    = aux{e}*conv_w; e = e+1;
%         Data_P(i).Cmv    = aux{e}*conv_w; e = e+1;
%         Data_P(i).Cmp    = aux{e}*conv_w; e = e+1;
%         Data_P(i).Cl    = aux{e}*conv_w; e = e+1;
%         Data_P(i).Cn    = aux{e}*conv_w; e = e+1;
%         Data_P(i).Cnv   = aux{e}*conv_w; e = e+1;
%         Data_P(i).Cnp   = aux{e}*conv_w; e = e+1;
%         % Aero
%         Data_P(i).CLCD    = aux{e}*conv_w; e = e+1;
%         Data_P(i).CL32CD    = aux{e}*conv_w; e = e+1;
%         Data_P(i).inv_sqrt_CL    = aux{e}*conv_w; e = e+1;
%         % Forces
%         Data_P(i).L    = aux{e}; e = e+1;
%         Data_P(i).D    = aux{e}; e = e+1;
%         Data_P(i).Fx_FF    = aux{e}; e = e+1;
%         Data_P(i).Fy_FF    = aux{e}; e = e+1;
%         Data_P(i).Fz_FF    = aux{e}; e = e+1;
%         Data_P(i).Fx_sum    = aux{e}; e = e+1;
%         Data_P(i).Fy_sum    = aux{e}; e = e+1;
%         Data_P(i).Fz_sum    = aux{e}; e = e+1;
%         Data_P(i).Extra_drag = aux{e}; e = e+1;
%         Data_P(i).Fuse_drag = aux{e}; e = e+1;  
%         Data_P(i).Cf_Fuse = aux{e}; e = e+1;
%         % Velocities
%         Data_P(i).Vx = aux{e}; e = e+1;
%         Data_P(i).Vz = aux{e}; e = e+1;     
%         Data_P(i).V = aux{e}; e = e+1;      
%         Data_P(i).Gamma = aux{e}; e = e+1;       
%         % Moments
%         Data_P(i).L = aux{e}; e = e+1;
%         Data_P(i).M = aux{e}; e = e+1; 
%         Data_P(i).N = aux{e}; e = e+1;
%         % Pressude Coefficients
%         Data_P(i).CPx = aux{e}*conv_w; e = e+1;     
%         Data_P(i).CPy = aux{e}*conv_w; e = e+1;      
%         Data_P(i).CPz = aux{e}*conv_w; e = e+1;     
%         % Bending Moment
%         Data_P(i).BM = aux{e}; e = e+1;    
%         Data_P(i).m.g.Vz = aux{e}; e = e+1;  
%         Data_P(i).Drag_x_V = aux{e}; e = e+1; 
%         Data_P(i).Efficiency = aux{e};  e = e+1; 
%         Data_P(i).XCp.Cl = aux{e}; e = e+1;      
%         Data_P(i).XNP = aux{e}; e = e+1;     
%         % Stability
%         Data_P(i).Phugoid_Freq = aux{e}; e = e+1; 
%         Data_P(i).Phugoid_Damping = aux{e}; e = e+1; 
%         Data_P(i).Short_Period_Freq = aux{e}; e = e+1; 
%         Data_P(i).Short_Period_Damping_Ratio = aux{e}; e = e+1; 
%         Data_P(i).Dutch_Roll_Freq = aux{e}; e = e+1; 
%         Data_P(i).Dutch_Roll_Damping = aux{e}; e = e+1; 
%         Data_P(i).Roll_Damping = aux{e}; e = e+1; 
%         Data_P(i).Spiral_Damping = aux{e}; e = e+1;
%         Data_P(i).Mass = aux{e}; e = e+1;
%         Data_P(i).CoG_x = aux{e}; e = e+1;
%         Data_P(i).CoG_z = aux{e}; e = e+1;
%         % Close file
%         fclose(file_id);
%     end
% end
% toc

