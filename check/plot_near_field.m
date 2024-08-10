close all; clear; clc;

data  =   load("../data/near_field.txt");
Ex = data(:, 2)+1j*data(:, 3);
Ey = data(:, 4)+1j*data(:, 5);
Ez = data(:, 6)+1j*data(:, 7);

figure()
hold on
plot(data(:, 1), abs(Ex))
plot(data(:, 1), abs(Ey), '--')
plot(data(:, 1), abs(Ez), '-k')
hold off

Hx = data(:, 8)+1j*data(:, 9);
Hy = data(:, 10)+1j*data(:, 11);
Hz = data(:, 12)+1j*data(:, 13);


figure()
hold on
plot(data(:, 1), abs(Hx))
plot(data(:, 1), abs(Hy), '--')
plot(data(:, 1), abs(Hz), '-k')
hold off
