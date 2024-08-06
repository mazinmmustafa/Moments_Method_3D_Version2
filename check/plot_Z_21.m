close all; clear; clc;

N = 2;
Z0 = 50;
data  =   load("../data/S_matrix.txt");

[Ns, ~] = size(data);

S = zeros(Ns, N, N);
for i=1:Ns
  counter = 1;
  for m=1:N
    for n=1:N
      counter = counter+1;
      S(i, m, n) = data(i, counter);
      counter = counter+1;
      S(i, m, n) = S(i, m, n)+1j*data(i, counter);
    end
  end
  
end

Z = zeros(Ns, N, N);
I = eye(N);
for i=1:Ns
  Z(i, :, :) = (I+reshape(S(i, :, :), [2 2]))*inv(I-reshape(S(i, :, :), [2 2]))*I*Z0;
end
 

figure()
hold on
plot(data(:, 1), real(Z(:, 2, 1)), '-')
plot(data(:, 1), imag(Z(:, 2, 1)), '-')
plot([0 max(data(:, 1))], [0 0], '-k')
hold off

##figure()
##hold on
##plot(data(:, 1), 20*log10(abs(S(:, 1, 1))))
##plot(data(:, 1), 20*log10(abs(S(:, 2, 1))))
##hold off