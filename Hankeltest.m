% Parameters
r = linspace(0, 10, 100); % Radial positions
K = 10; % Maximum frequency
n = 0; % Transform order

% Define a test function h(r)
h_original = exp(-r.^2);

% Perform Hankel Transform
[H, k] = fht(h_original, K, r, n);

% Perform Inverse Hankel Transform
h_recovered = ifht(H, k, r, n);

% Plot results
figure;
subplot(1,2,1);
plot(r, h_original);
title('Original Function h(r)');
xlabel('r');
ylabel('h(r)');

subplot(1,2,2);
plot(r, h_recovered);
title('Recovered Function h(r)');
xlabel('r');
ylabel('h(r)');

legend('Original', 'Recovered');
