function [Body_Geo] = getAdimFus(XFLR5_file, nPoints, nSections, skip)
    
    % nSections = Número de secciones del fuselaje
    % nPoints   = Número de puntos de cada sección
    % skip      = Salto de lineas entre sección y sección
    lecture_Geo = 5;                % lineas a partir desde donde empieza a leer
    
   
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
    
    l_fus           = max(x_max);
    % Check if there are repeated elements in x_max (several sections for the same x coordinate)
    rep = find(diff(x_max) == 0);
    
    for k = rep
        if k ~= 1
            x_max(k) = (x_max(k-1) + x_max(k+1))/2;
            y_max(k) = (y_max(k-1) + y_max(k+1))/2;
            z_max(k) = (z_max(k-1) + z_max(k+1))/2;
            z_min(k) = (z_min(k-1) + z_min(k+1))/2;
            
            Data_P(k).x = (Data_P(k-1).x + Data_P(k+1).x)/2;
            Data_P(k).y = (Data_P(k-1).y + Data_P(k+1).y)/2;
            Data_P(k).z = (Data_P(k-1).z + Data_P(k+1).z)/2;
            
            x(:,k) = (x(:,k-1) + x(:,k+1))/2;
            y(:,k) = (y(:,k-1) + y(:,k+1))/2;
            z(:,k) = (z(:,k-1) + z(:,k+1))/2;
        else
            x_max(k+1) = (x_max(k) + x_max(k+2))/2;
            y_max(k+1) = (y_max(k) + y_max(k+2))/2;
            z_max(k+1) = (z_max(k) + z_max(k+2))/2;
            z_min(k+1) = (z_min(k) + z_min(k+2))/2;
            
            Data_P(k+1).x = (Data_P(k).x + Data_P(k+2).x)/2;
            Data_P(k+1).y = (Data_P(k).y + Data_P(k+2).y)/2;
            Data_P(k+1).z = (Data_P(k).z + Data_P(k+2).z)/2;
            
            x(:,k+1) = (x(:,k) + x(:,k+2))/2;
            y(:,k+1) = (y(:,k) + y(:,k+2))/2;
            z(:,k+1) = (z(:,k) + z(:,k+2))/2;
        end
    end
    
    Body_Geo.x      = x_max/l_fus;
    Body_Geo.l      = 1;
    
    for i = 1:nSections 
        Data_P(i).x     = Data_P(i).x/l_fus;
        Data_P(i).y     = Data_P(i).y/l_fus;
        Data_P(i).z     = Data_P(i).z/l_fus;
    end
    
    x               = x/l_fus;
    y               = y/l_fus;
    z               = z/l_fus;
    
    x_max           = x_max/l_fus;
    y_max           = y_max/l_fus;
    z_max           = z_max/l_fus;
    z_min           = z_min/l_fus;
    
    
    
    Body_Geo.meshData{1}     = [x(:,end:(-1):1),x];
    Body_Geo.meshData{2}     = [-y(:,end:(-1):1), y];
    Body_Geo.meshData{3}     = [z(:,end:(-1):1), z];
    %plot(Body_Geo.meshData{2}(end,:),Body_Geo.meshData{3}(end,:))
    
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
    
    for i=2:(nSections-1)
        [geom,~,~]      = polygeom(Per(i).y,Per(i).z);
        Sfront_x(i)     = geom(1);
        X_cen(i)        = geom(2);
        Z_cen(i)        = geom(3);
        Perimeter(i)    = geom(4);
    end
    
    Sfront_x(1)         = 0;
    X_cen(1)            = 0;
    Z_cen(1)            = 0;
    Perimeter(1)        = 0;

    Sfront_x(nSections)         = 0;
    X_cen(nSections)            = 0;
    Z_cen(nSections)            = 0;
    Perimeter(nSections)        = 0;
    
    Body_Geo.S_x        = Sfront_x;
        
    Delta_h             = Body_Geo.x(2:end) - Body_Geo.x(1:(end-1));
    Delta_Zcent         = Z_cen(2:end)-Z_cen(1:(end-1));
    alpha               = -atan(Delta_h./Delta_Zcent);
    Body_Geo.alpha_x    = [alpha,alpha(end)];
   
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
