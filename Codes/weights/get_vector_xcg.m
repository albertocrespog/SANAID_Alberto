function [x] = get_vector_xcg(a,b,N,val)
n1 = round((val-a)/(b-a)*N);
n2 = N - n1;
salto1 = (val-a)/(n1-1);
salto2 = (b-val)/(n2);
x1 = a:salto1:val;
x2 = val:salto2:b;
x = [x1,x2(2:end)]
end