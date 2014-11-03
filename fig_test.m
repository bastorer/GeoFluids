[x,y,z] = meshgrid(-1:0.01:1, -1:0.01:1, -1:0.01:1);
% C = sin(pi*x).*sin(pi*y).*sin(pi*z);
% C = exp(-20*(x.^2+y.^2+z.^2));

C = exp(z.^2).*x.*cos(pi*y/2);

figure;
VolumePlot(x,y,z,C)
