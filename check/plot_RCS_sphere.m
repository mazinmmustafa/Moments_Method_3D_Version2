close all; clear; clc;

data1  =   load("../data/RCS_1.txt");
data2  =   load("../data/RCS_2.txt");

figure()
hold on
plot(data1(:, 1), 10*log10(data1(:, 2)))
plot(data1(:, 1), 10*log10(data1(:, 3)))
hold off
xlim([0 180])
ylim([-20 +30])

figure()
hold on
plot(data2(:, 1), 10*log10(data2(:, 2)))
plot(data2(:, 1), 10*log10(data2(:, 3)))
hold off
xlim([0 180])
ylim([-20 +30])