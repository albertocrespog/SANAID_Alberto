
function [F]=ajuste_cE_F(deltae,n_ref,ABC)

n_aux=[1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000]; % valores de 'x'
n=n_aux./60; % pasamos de rpm a rev/s


if deltae==100
    ce_aux=[889,769,708,697,618,554,517,571,643,641,630]; % consumo específico en g/kWh
elseif deltae==90
    ce_aux=[919,762,700,702,624,561,510,573,626,641,632];
elseif deltae==85
    ce_aux=[902,748,697,722,615,546,499,557,616,641,640];
elseif deltae==80
    ce_aux=[902,759,680,693,604,523,492,527,593,644,643];
elseif deltae==75
    ce_aux=[884,752,678,659,565,510,469,499,588,636,635];
elseif deltae==70
    ce_aux=[884,823,688,665,568,484,457,490,572,614,624];
elseif deltae==65
    ce_aux=[906,840,672,672,555,461,449,466,551,599,610];
elseif deltae==60
    ce_aux=[879,841,646,606,521,435,427,468,556,601,580];
elseif deltae==55
    ce_aux=[862,823,586,580,513,400,388,446,556,609,577];
elseif deltae==50
    ce_aux=[988,804,563,561,484,396,377,434,544,554,540];
elseif deltae==40
    ce_aux=[1182,811,589,541,500,397,383,431,542,559,553];
elseif deltae==30
    ce_aux=[1461,972,811,648,517,457,431,412,500,522,512];
elseif deltae==20
    ce_aux=[2043,1674,1360,1143,769,669,591,422,448,523,534];
elseif deltae==15
    ce_aux=[4283,2798,2196,1475,960,736,588,400,401,459,446];
elseif deltae==10
    ce_aux=[7076,3457,2845,1810,1214,988,611,367,347,369,336];
end

ce=ce_aux/(1000); % pasamos de g/kWh a kg/kWh


% Defino matrices y funciones para resolucion por minimos cuadrados -->
% cálculo de F (que es función de deltae) habiendo determinado ya C1, C2, y C3

d= @(f1) n_ref.^2;
f= @(f1) f1*(ABC(1)*n.^2+ABC(2)*n*n_ref+ABC(3)*n_ref^2)./d(f1);
M= @(f1) sum((ce-f(f1)).*((ABC(1)*n.^2+ABC(2)*n*n_ref+ABC(3)*n_ref^2)./(n_ref.^2)));
M(1)

F=fsolve(M,0.5);

f=@(f1,n) f1*(ABC(1)*n.^2+ABC(2)*n*n_ref+ABC(3)*n_ref^2)./d(f1);

figure

plot(n,ce,'or')
xlim([10 100])
ylim([0 1])
hold on
xx=16:100;
yy=f(F,xx);
plot(xx,yy,'blue')
hold off

end