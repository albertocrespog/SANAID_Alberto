function [Body_Geo, plotData] = getAdimFus(XFLR5_file, nPoints, nSections)
    
    % nSections = Número de secciones del fuselaje
    % nPoints   = Número de puntos de cada sección
    lecture_Geo = 5;                % lineas a partir desde donde empieza a leer
    skip        = 2;
    
   
    for i=1:nSections
        % Extracción de datos VLM
        file_id=fopen(XFLR5_file);
        aux = textscan(file_id,'%n%n%n','Headerlines',lecture_Geo);
        Data_P(i).raw_data  = aux;
        Data_P(i).label     = XFLR5_file;
        Data_P(i).x         = aux{1};
        Data_P(i).y         = aux{2};
        Data_P(i).z         = aux{3};
        x(:,i)              = aux{1};
        y(:,i)              = aux{2};
        z(:,i)              = aux{3};
        x_max(i)            = max(aux{1});
        y_max(i)            = max(aux{2});
        z_max(i)            = max(aux{3});
        z_min(i)            = min(aux{3});
        fclose(file_id);
        lecture_Geo = lecture_Geo + nPoints  + skip;
    end
    
    Body_Geo.x      = x_max;
    Body_Geo.l      = max(Body_Geo.x);
    
    Data_P(i).x     = Data_P.x/Body_Geo.l;
    Data_P(i).y     = Data_P.y/Body_Geo.l;
    Data_P(i).z     = Data_P.z/Body_Geo.l;
    
    x               = x/Body_Geo.l;
    y               = y/Body_Geo.l;
    z               = z/Body_Geo.l;
    
    x_max           = x_max/Body_Geo.l;
    y_max           = y_max/Body_Geo.l;
    z_max           = z_max/Body_Geo.l;
    
    
    plotData{1}     = [x(:,end:(-1):1),x];
    plotData{2}     = [-y(:,end:(-1):1), y];
    plotData{3}     = [z(:,end:(-1):1), z];
    
    % Calculo de los perímetros de cada sección
    for i = 1:nSections
        Per(i).y = [Data_P(i).y; -Data_P(i).y; Data_P(i).y(1)];
        Per(i).z = [Data_P(i).z; fliplr(Data_P(i).z')';Data_P(i).z(1)];
    end
    
    
    % [ GEOM, INER, CPMO ] = polygeom( X, Y ) returns
    % area, centroid, perimeter and area moments of 
    % inertia for the polygon.
    
    % función que permite calcular el área, la posición x e y del centroide y
    % el perímetro de cada sección: obvia el primer y último punto del
    % fuselaje
    
    for i=2:(nSections)
        [geom,~,~]      = polygeom(Per(i).y,Per(i).z);
        Sfront_x(i)     = geom(1);
        X_cen(i)        = geom(2);
        Z_cen(i)        = geom(3);
        Perimeter(i)    = geom(4);
    end
    
    Body_Geo.S_x    = [0, Sfront_x, 0];
    
    
    Delta_h             = Body_Geo.x(2:end) - Body_Geo.x(1:(end-1));
    Delta_Zcent         = Z_cen(2:end)-Z_cen(1:(end-1));
    alpha               = -atan(Delta_h./Delta_Zcent);
    Body_Geo.alpha_x    = [alpha,alpha(end)];

%     % Metodo alternativo aproximando con troncos de cono
%     r_eq = Perimeter/2/pi;              % Radio equivalente de cada sección
%     g = sqrt(Delta_h.^2 + (r_eq(2:end)-r_eq(1:(end-1))).^2); % Generatriz 
%     Surface = pi*g.*(r_eq(2:end) + r_eq(1:(end-1)));   % Area de cada tronco
% 
% 
%     % Surface total geometría
%     Surf_TOT = sum(Surface);
%     Body_Geo.Slat = Surf_TOT;
   
    % Calculo del volumen (aproximación)    
    Vol = Delta_h.*(Body_Geo.S_x(2:end) + Body_Geo.S_x(1:(end-1)) + sqrt(Body_Geo.S_x(2:end).*Body_Geo.S_x(1:(end-1))))/3;


    % Volumen total geometría
    Vol_TOT = sum(Vol);
    Body_Geo.vol = Vol_TOT;


    dSdX                = (Body_Geo.S_x(2) - Body_Geo.S_x(1))/(Body_Geo.x(2) - Body_Geo.x(1));
    dSdX                = [dSdX,(Body_Geo.S_x(3:end) - Body_Geo.S_x(1:(end-2)))./(Body_Geo.x(3:end) - Body_Geo.x(1:(end-2)))];
    dSdX                = [dSdX, (Body_Geo.S_x(end) - Body_Geo.S_x(end-1))/(Body_Geo.x(end) - Body_Geo.x(end-1))];
    Body_Geo.dSdX       = dSdX; 
    
    i_SfrontMax = 1;
    while dSdX(i_SfrontMax) > 0 && i_SfrontMax <= nSections
    i_SfrontMax=i_SfrontMax+1;
    end
    
    
    Body_Geo.Sfront = Body_Geo.S_x(i_SfrontMax);
    Body_Geo.xSmax  = x_max(i_SfrontMax);
    w_Area_b_max    = 2*y_max(i_SfrontMax);
    h_Area_b_max    = 2*z_max(i_SfrontMax);
    
    % Calculo del Area proyectada en plano x-y
    Body_Geo.Stop  = 2*polyarea(x_max,y_max);
    % Calculo del Area proyectada en plano x-z (se asume seccion simetrica)
    Body_Geo.Sside = 2*polyarea(x_max,z_max);

    %length_x_position = x_max;
    Body_Geo.W_x = 2*y_max;
    Body_Geo.W   = max(Body_Geo.W_x);
    Body_Geo.D_x = z_max - z_min;
    Body_Geo.D   = max(Body_Geo.D_x);    
end
