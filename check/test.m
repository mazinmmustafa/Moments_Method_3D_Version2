close all; clear; clc;
%%
x = linspace(0, 0.1, 101);
a = 1.0E-3;

figure()
plot(x, log((sqrt(a^2+x.^2)+x)./(sqrt(a^2+x.^2)-x)))

figure()
plot(x, log((sqrt(a^2+x.^2)-x)./(sqrt(a^2+x.^2)+x)))

figure()
plot(x, log((sqrt(a^2+x.^2)-x)./(sqrt(a^2+x.^2)-x)))