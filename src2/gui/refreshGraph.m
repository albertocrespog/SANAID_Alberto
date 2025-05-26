function refreshGraph(surfObj,centers,model)
    
    L = model.general.L;
    X_cg = model.general.Xcg;
    
    X_na = model.general.Xna;
    
    
    Xw = model.ala.Xca;
    if isempty(model.ala.i)    
        iw = 0;
    else
        iw = model.ala.i;
    end
    c_w = model.ala.MAC;
    
    
    if isempty(model.horizontal.Xca)
        Xh = 0;
    else
        Xh = model.horizontal.Xca;
    end
    
    if isempty(model.horizontal.i)
        ih = 0;
    else
        ih = model.horizontal.i;
    end
    
    if isempty(model.canard.Xca)
        Xc = 0;
    else
        Xc = model.canard.Xca;
    end
    if isempty(model.canard.i)
        ic = 0;
    else
        ic = model.canard.i;
    end
        
    
    Mw = makehgtform('translate',Xw/L,0,0)*makehgtform('zrotate',-iw);
    set(surfObj(1),'Matrix',Mw);
    
    M1 = makehgtform('translate',Xh/L,0,0)*makehgtform('zrotate',-ih);
    set(surfObj(2), 'Matrix',M1, 'Visible', 'off');
    
    M2 = makehgtform('translate',Xc/L,0,0)*makehgtform('zrotate',-ic);
    set(surfObj(3),'Matrix',M2, 'Visible', 'off');
    
    switch model.conf
        case 'convencional'
            set(surfObj(2), 'Visible', 'on')
        case 'canard'
            set(surfObj(3), 'Visible', 'on')
        case 'convencional_canard'
            set(surfObj(2), 'Visible', 'on')
            set(surfObj(3), 'Visible', 'on')
        case 'flyWing'
            set(surfObj(2), 'Visible', 'off')
            set(surfObj(3), 'Visible', 'off')
    end
    
    
    Mna = makehgtform('translate',X_na/L,0,0);
    set(centers(1),'Matrix',Mna); 
    Mcg = makehgtform('translate',X_cg/L,0,0);
    set(centers(2),'Matrix',Mcg);
    
    c_w = model.ala.MAC;
    xlim([-1/4 5/4]);
    ylim([-0.7, 0.3]*c_w/L);
    
    drawnow;
end
