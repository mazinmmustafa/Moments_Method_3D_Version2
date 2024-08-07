close all; clear; clc;

N = 2;
Z0 = 50;
data  =   load("../data/S_matrix.txt");

[Ns, ~] = size(data);

for i=1:Ns
  fprintf("%f\n", abs(data(i, 4)+1j*data(i, 5)));
  
  end
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

figure()
hold on
plot(data(:, 1), 20*log10(abs(S(:, 1, 1))))
plot(data(:, 1), 20*log10(abs(S(:, 2, 1))))
hold off