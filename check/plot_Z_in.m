close all; clear; clc;

data  =   load("../data/Z_in.txt");

figure()
hold on
plot(data(:, 1), data(:, 2)/1E3, '-')
plot(data(:, 1), data(:, 3)/1E3, '-')
hold off

ylim([-5 +5])