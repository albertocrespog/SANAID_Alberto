function [cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c_d)

if isstruct(airfoil_data) == 1
    airfoil_data = airfoil_data.data;
else
end

airfoil_data(:,1) = airfoil_data(:,1) - min(airfoil_data(:,1));

xx = airfoil_data(:,1)';
yy = airfoil_data(:,2)';
x_ie = find(xx==0);
xx_e = xx(1:x_ie);
yy_e = yy(1:x_ie);
yy_i = yy(x_ie+1:end);
xx_i = xx(x_ie+1:end);
y90_i = interp1(xx_i,yy_i,0.9);
y90_e = interp1(xx_e,yy_e,0.9);
y99_e = interp1(xx_e,yy_e,0.99);
y99_i = interp1(xx_i,yy_i,0.99);
y90 = (y90_e-y90_i)*100;
y99 = (y99_e-y99_i)*100;
tan_1_2_LE = (0.5*y90-0.5*y99)/9;

%Fig 4.1.1.2-8a 
  tan_12_LE = [0., 0.02,  0.04,  0.06,  0.08,  0.10,  0.12,  0.14,  0.16,...
          0.18,  0.20 ];
   
  claclatheo =  [ 0.900,  0.878,  0.858,  0.836,  0.815,  0.794,  0.772,  0.750,...
      0.728,  0.708,  0.685];
    
    if tan_1_2_LE > 0.2
        tan_1_2_LE = 0.2;
    end
  cla_clatheo = interp1(tan_12_LE,claclatheo,tan_1_2_LE);

 %Fig 4.1.1.2-8b 
 
 tdc = [0,0.02,0.06,0.1,0.14,0.18,0.2];
 clatheory = [6.3, 6.4, 6.6, 6.8, 7, 7.2,7.3];
if t_c_d >0.2
    t_c_d = 0.2;
end
 clatheo = interp1(tdc, clatheory, t_c_d);
end