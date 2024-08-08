close all; clear; clc;

data  =   load("../data/near_field.txt");
Ex = data(:, 2)+1j*data(:, 3);
Ey = data(:, 4)+1j*data(:, 5);
Ez = data(:, 6)+1j*data(:, 7);

figure()
hold on
plot(data(:, 1), abs(Ex))
plot(data(:, 1), abs(Ey), '--')
plot(data(:, 1), abs(Ez))
hold off
