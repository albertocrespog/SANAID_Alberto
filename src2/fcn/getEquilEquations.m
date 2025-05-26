%%  -----------------------------------------------------------------------
%   Academic Stability Pro
%   Developed by Álvaro Fernández Cobo
%   On January 4th, 2016

%%   ----------------------------------------------------------------------
% This function returns the equilibrium system of equations, composed of
% forces and moments equations

function [eq_F,eq_M] = getEquilEquations()
    
syms Sref Sw CL0w Sw c_w Xw CLaw CMw St c_t Xt CL0t CLat CMt eta_t E0t Sc c_c Xc CL0c CLac CMc eta_c E0c W q Xcg iw it ic;
    eq_F    (Sw,CL0w,CLaw,St,CL0t,CLat,eta_t,E0t,Sc,CL0c,CLac,eta_c,...
            E0c,Sref,W,q,iw,it,ic)    =     Sw/Sref*CL0w + Sw/Sref*CLaw*iw + (CL0t*St*eta_t)/Sref + ...
                                            (CLat*St*eta_t*(it - E0t))/Sref + eta_c*Sc/Sref*CL0c + ...
                                            eta_c*Sc/Sref*CLac*(ic + E0c) == W/(Sref*q);
        
    eq_M    (Sw,c_w,Xw,CL0w,CLaw,CMw,St,c_t,Xt,CL0t,CLat,CMt,eta_t,...
            E0t,Sc,c_c,Xc,CL0c,CLac,CMc,eta_c,E0c,Sref,W,q,Xcg,iw,it,ic) =  eta_c*Sc/Sref*c_c/c_w*CMc + eta_c*Sc/Sref/c_w*(Xcg-Xc)*...
                                                                            (CL0c + CLac*(ic + E0c)) + Sw/Sref*CMw + Sw/Sref/c_w*...
                                                                            (Xcg-Xw)*(CL0w+CLaw*iw)+eta_t*St/Sref*c_t/c_w*CMt...
                                                                            + eta_t*St/Sref*(Xcg-Xt)*(CL0t+CLat*(it - E0t))==0;
end




                                                                    