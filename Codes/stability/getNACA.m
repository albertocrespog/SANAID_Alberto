function [x,z] = getNACA(t_c)
    x                   = (logspace(0,log10(11))-1)/10;
    n                   = length(x);
    z(1:n)              = t_c/0.2*(0.2969*sqrt(x) - 0.126*x - 0.3516*x.^2 + 0.2843*x.^3-0.1015*x.^4);
    z(2*n:(-1):(n+1))   = -t_c/0.2*(0.2969*sqrt(x) - 0.126*x - 0.3516*x.^2 + 0.2843*x.^3-0.1015*x.^4);
    x                   = [x,x(end:(-1):1)];
end