function Get_Prop_dw_Surface = Get_Prop_DownWash_Surface(Geo_tier)

x_prop_cF = Geo_tier.x_prop_cF;
y_prop_cF = Geo_tier.y_prop_cF;
y_y2_prop_dw = Geo_tier.y_y2_prop_dw; % inner location of prop disk

y_y1_prop_dw = Geo_tier.y_y1_prop_dw; % inner portion of the prop-wash

y_cT_w1_LE = Geo_tier.y_cT_w1_LE; % y-location of tip chord
y_cT_w1_TE = Geo_tier.y_cT_w1_TE; % y-location of tip chord
y_cT_w2_LE = Geo_tier.y_cT_w2_LE; % y-location of tip chord
y_cT_w2_TE = Geo_tier.y_cT_w2_TE; % y-location of tip chord

x_1R_y1_w1 = Geo_tier.x_1R_y1_w1;
x_2R_y1_w1 = Geo_tier.x_2R_y1_w1;
x_1R_y1_w2 = Geo_tier.x_1R_y1_w2;
x_2R_y1_w2 = Geo_tier.x_2R_y1_w2;

x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
x_cT_w2_LE = Geo_tier.x_cT_w2_LE;
x_cT_w2_TE = Geo_tier.x_cT_w2_TE;

y_offset_w1 = Geo_tier.y_offset_w1;
Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
Lambda_TE_w1 = Geo_tier.Lambda_TE_w1;

y_offset_w2 = Geo_tier.y_offset_w2;
Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
Lambda_TE_w2 = Geo_tier.Lambda_TE_w2;

% Checks if there is prop wash affecting w1
if y_y1_prop_dw < y_cT_w1_LE
    
    % Estimation of the corners of the affected w1 surface affected by prop
    % wash
    x_1R_y1_w1_pw = x_1R_y1_w1 + (y_y1_prop_dw - y_offset_w1)*tan(Lambda_LE_w1); % defines inner position of wing affeted by proopwash 
    x_2R_y1_w1_pw = x_2R_y1_w1 + (y_y1_prop_dw - y_offset_w1)*tan(Lambda_TE_w1); % defines inner position of wing affeted by proopwash
    x_1R_y2_w1_pw = x_cT_w1_LE;
    x_2R_y2_w1_pw = x_cT_w1_TE;
    % Chord
    c_y1_w1_pw = x_2R_y1_w1_pw - x_1R_y1_w1_pw; % defines chord of wing at beginning of prop wash surface
    c_y2_w1_pw = x_2R_y2_w1_pw - x_1R_y2_w1_pw; % defines chord of wing at beginning of prop wash surface
    % Storing DATA
    Get_Prop_dw_Surface.x_1R_y1_w1_pw = x_1R_y1_w1_pw;
    Get_Prop_dw_Surface.x_2R_y1_w1_pw = x_2R_y1_w1_pw;
    Get_Prop_dw_Surface.x_1R_y2_w1_pw = x_1R_y2_w1_pw;
    Get_Prop_dw_Surface.x_2R_y2_w1_pw = x_2R_y2_w1_pw;
    Get_Prop_dw_Surface.c_y1_w1_pw = c_y1_w1_pw;
    Get_Prop_dw_Surface.c_y2_w1_pw = c_y2_w1_pw;
    % y-location
    y_1R_y1_w1_pw = y_y1_prop_dw;
    y_2R_y1_w1_pw = y_y1_prop_dw;
    y_1R_y2_w1_pw = y_cT_w1_LE;
    y_2R_y2_w1_pw = y_cT_w1_TE;
    % Storing DATA
    Get_Prop_dw_Surface.y_1R_y1_w1_pw = y_1R_y1_w1_pw;
    Get_Prop_dw_Surface.y_2R_y1_w1_pw = y_2R_y1_w1_pw;
    Get_Prop_dw_Surface.y_1R_y2_w1_pw = y_1R_y2_w1_pw;
    Get_Prop_dw_Surface.y_2R_y2_w1_pw = y_2R_y2_w1_pw;
    
    % Geometry
    lambda_w1_pw = c_y2_w1_pw/c_y1_w1_pw; % w1 prop wash taper ratio
    b_w1_pw = y_1R_y2_w1_pw - y_1R_y1_w1_pw; % same span
    S_w1_pw_eng = b_w1_pw*((c_y1_w1_pw + c_y2_w1_pw)/2); % total area effective surface
    S_w1_pw = 2*S_w1_pw_eng; % total area effective surface for both engines
    AR_w1_pw = b_w1_pw^2/S_w1_pw; % Aspect Ratio
    % Storing DATA
    Get_Prop_dw_Surface.lambda_w1_pw = lambda_w1_pw;
    Get_Prop_dw_Surface.b_w1_pw = b_w1_pw;
    Get_Prop_dw_Surface.S_w1_pw = S_w1_pw;
    Get_Prop_dw_Surface.AR_w1_pw = AR_w1_pw;    
else
    S_w1_pw = 0;
    Get_Prop_dw_Surface.S_w1_pw = S_w1_pw;
    
end

% % Checks if there is prop wash affecting w1
if y_y1_prop_dw < y_cT_w2_LE
    
    % Estimation of the corners of the affected w2 surface affected by prop
    % wash
    x_1R_y1_w2_pw = x_1R_y1_w2 + (y_y1_prop_dw - y_offset_w2)*tan(Lambda_LE_w2); % defines inner position of wing affeted by proopwash 
    x_2R_y1_w2_pw = x_2R_y1_w2 + (y_y1_prop_dw - y_offset_w2)*tan(Lambda_TE_w2); % defines inner position of wing affeted by proopwash
    x_1R_y2_w2_pw = x_cT_w2_LE;
    x_2R_y2_w2_pw = x_cT_w2_TE;
    % Chord
    c_y1_w2_pw = x_2R_y1_w2_pw - x_1R_y1_w2_pw; % defines chord of wing at beginning of prop wash surface
    c_y2_w2_pw = x_2R_y2_w2_pw - x_1R_y2_w2_pw; % defines chord of wing at beginning of prop wash surface
    % Storing DATA
    Get_Prop_dw_Surface.x_1R_y1_w2_pw = x_1R_y1_w2_pw;
    Get_Prop_dw_Surface.x_2R_y1_w2_pw = x_2R_y1_w2_pw;
    Get_Prop_dw_Surface.x_1R_y2_w2_pw = x_1R_y2_w2_pw;
    Get_Prop_dw_Surface.x_2R_y2_w2_pw = x_2R_y2_w2_pw;
    Get_Prop_dw_Surface.c_y1_w2_pw = c_y1_w2_pw;
    Get_Prop_dw_Surface.c_y2_w2_pw = c_y2_w2_pw;
    % y-location
    y_1R_y1_w2_pw = y_y1_prop_dw;
    y_2R_y1_w2_pw = y_y1_prop_dw;
    y_1R_y2_w2_pw = y_cT_w2_LE;
    y_2R_y2_w2_pw = y_cT_w2_TE;
    % Storing DATA
    Get_Prop_dw_Surface.y_1R_y1_w2_pw = y_1R_y1_w2_pw;
    Get_Prop_dw_Surface.y_2R_y1_w2_pw = y_2R_y1_w2_pw;
    Get_Prop_dw_Surface.y_1R_y2_w2_pw = y_1R_y2_w2_pw;
    Get_Prop_dw_Surface.y_2R_y2_w2_pw = y_2R_y2_w2_pw;
    
    % Geometry
    lambda_w2_pw = c_y2_w2_pw/c_y1_w2_pw; % w2 prop wash taper ratio
    b_w2_pw = y_1R_y2_w2_pw - y_1R_y1_w2_pw; % same span
    S_w2_pw_eng = b_w2_pw*((c_y1_w2_pw + c_y2_w2_pw)/2); % total area effective surface
    S_w2_pw = 2*S_w2_pw_eng; % total area effective surface for both engines
    AR_w2_pw = b_w2_pw^2/S_w2_pw; % Aspect Ratio
    % Storing DATA
    Get_Prop_dw_Surface.lambda_w2_pw = lambda_w2_pw;
    Get_Prop_dw_Surface.b_w2_pw = b_w2_pw;
    Get_Prop_dw_Surface.S_w2_pw = S_w2_pw;
    Get_Prop_dw_Surface.AR_w2_pw = AR_w2_pw;    
else
    S_w2_pw = 0;
    Get_Prop_dw_Surface.S_w2_pw = S_w2_pw;
    
end

