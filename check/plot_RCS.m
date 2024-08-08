close all; clear; clc;

data  =   load("../data/RCS.txt");

figure()
hold on
plot(data(:, 1), 10*log10(data(:, 2)))
plot(data(:, 1), 10*log10(data(:, 3)))
hold off
xlim([0 4])
ylim([-30 +10])