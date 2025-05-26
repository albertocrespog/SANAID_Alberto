function [A, B, C] = getPropModel(Psl, delta_T, M, h)
    gamma   = 1.4;
    [~, ~, p_SL, ~] = atmos_inter(0);
    [~, ~, pinf, ainf] = atmos_inter(h);
    Vinf    = ainf*M;
    V       = (Vinf-150):0.1:(Vinf+150);
    P_eng   = delta_T*Psl*(1+(gamma-1)/2*(V/ainf).^2).^((gamma)/(gamma-1))*pinf/p_SL;
    p       = polyfit(V,P_eng,2);
    A       = p(1);
    B       = p(2);
    C       = p(3);
end
