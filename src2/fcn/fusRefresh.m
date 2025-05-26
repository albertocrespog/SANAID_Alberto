function fus = fusRefresh(fus, lnew, Sref)
    rel = lnew/fus.l;
% Geometry Update
    fus.l       = fus.l*rel;
    fus.D       = fus.D*rel;
    fus.W       = fus.W*rel;
    fus.Sside   = fus.Sside*rel^2;
    fus.Sfront  = fus.Sfront*rel^2;
    fus.Stop    = fus.Stop*rel^2;
    fus.x       = fus.x*rel;
    fus.D_x     = fus.D_x*rel;
    fus.W_x     = fus.W_x*rel;
    fus.S_x     = fus.S_x*rel;
    fus.dSdX    = fus.dSdX*rel;
    fus.vol     = fus.vol*rel^3;

% Lift Slope Coefficient Update
    % Fuselage apparent mass coefficient. Pamadi, Figure 3.6
    f_k2_k1         = k2_k1_calc(fus.l, fus.D);
    % Nose Lift Slope Coefficient
    fus.CLa    = 2*f_k2_k1*(fus.Sfront/Sref);
    
% Mesh Update
    fus.meshData{1} = fus.meshData{1}*rel;
    fus.meshData{2} = fus.meshData{2}*rel;
    fus.meshData{3} = fus.meshData{3}*rel;
    
    
    
    



