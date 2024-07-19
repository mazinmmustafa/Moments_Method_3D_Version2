close all; clear; clc;
##
data  =   load("../data/test.dat");

hold on
plot(data(:,1), data(:,2))
plot(data(:,1), data(:,3))
hold off
set(gca, "FontSize", 22)
grid on 
grid minor
##