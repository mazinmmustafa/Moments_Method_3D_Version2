close all; clear; clc;

data  =  load("../data/near_field_2d_data.txt");
x  =   load("../data/near_field_2d_x.txt");
y  =   load("../data/near_field_2d_y.txt");

figure()
pcolor(x, y, 20*log10(data))
hold on
plot(cos(linspace(0, 2*pi))*0.5, sin(linspace(0, 2*pi))*0.5, '-k')
hold off
colormap(jet)
shading flat
axis equal
colorbar
##caxis([-20 30])
caxis([0 2])
axis off