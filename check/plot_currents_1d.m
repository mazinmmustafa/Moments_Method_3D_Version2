close all; clear; clc;

data = load("../data/currents.txt");
data = unique(data, 'rows');

I1 = data(:, 4)(data(:, 2)==-0.25);
x = data(:,1)(data(:, 2)==-0.25)

I2 = data(:, 4)(data(:, 2)==+0.25);
y = data(:,1)(data(:, 2)==+0.25)

figure()
hold on
plot(x, I1*1E3)
plot(y, I2*1E3, '--')
hold off
xlim([-1 +1]*0.75)
ylim([0 3])