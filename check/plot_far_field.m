close all; clear; clc;

data  =   load("../data/far_field.txt");
E_far   =   data(:, 3).^2;
E_far   =   E_far/max(E_far);
E_far   =   10*log10(E_far);

figure()
range = [-30 0];
E_far(isnan(E_far))       =   min(range);
E_far(E_far < min(range))	=   min(range);
E_far = E_far-min(range);
polar(data(:, 1)*pi/180,E_far)
axis equal
ylim([0 30])
