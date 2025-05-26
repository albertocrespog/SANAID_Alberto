    function [sustSurface,centers] = initGraph(model,dir)
    cla;
    c_w = model.ala.MAC;
    c_t = model.horizontal.MAC;
    c_c = model.canard.MAC;
    L   = model.general.L;
    hold on;
    axis equal;
    grid on;
    ax = gca;
    xlabel('x/L');
    
    
%     file_id=fopen('./NACA0015.txt');
%     aux = textscan(file_id,'%n%n','Headerlines',1);
%     x = aux{1};
%     y = aux{2};
    
    load([dir.images,'airfoil.mat']);
    w_foil = fill(1.4*c_w/L*(x-1/4),1.4*c_w/L*y,[0 0.5 1],'LineWidth',2,'EdgeColor',[0 0.5 1]*0.6);
    %w_foil = plot(x-1/4,y);
    wing = hgtransform('Parent',ax);
    set(w_foil,'Parent',wing);
    sustSurface(1) = wing;
    
    t_foil = fill(1.4*c_t/L*(x-1/4),1.4*c_t/L*y,[0 1 0.5],'LineWidth',2,'EdgeColor',[0 1 0.5]*0.6);
    tail = hgtransform('Parent', ax);
    %t_foil = plot(c_t/L*x-c_t/L/4, c_t/L*y);
    set(t_foil,'Parent',tail);
    sustSurface(2) = tail;
            
    c_foil = fill(1.4*c_c/L*(x-1/4),1.4*c_c/L*y,[1 0.2 0.2],'LineWidth',2,'EdgeColor',[1 0.2 0.2]*0.6);
    %c_foil = fill(c_c/L*(x-1/4),c_c/L*y,[0 1 0],'LineWidth',2,'EdgeColor',[0 1 0]*0.6);
    canard = hgtransform('Parent', ax);
    %c_foil = plot(c_c/L*x-c_c/L/4, c_c/L*y);
    set(c_foil,'Parent',canard);
    sustSurface(3) = canard;    
    
    switch model.conf
        case 'convencional'
            set(sustSurface(3), 'Visible', 'off');
            
        case 'ala_canard'
            set(sustSurface(2), 'Visible', 'off');

        case 'flyWing'
            set(sustSurface(2), 'Visible', 'off');
            set(sustSurface(3), 'Visible', 'off');
    end
    
    Xna = hgtransform('Parent',ax);
    xna_point = scatter(0,0,80,'filled','b');
    set(xna_point,'Parent',Xna);
    centers(1) = Xna;
    
    Xcg = hgtransform('Parent',ax);
    xcg_point = scatter(0,0,80,'filled','r');
    set(xcg_point,'Parent',Xcg);
    centers(2) = Xcg;
    
    drawnow;
        
end