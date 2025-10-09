% ------------------------------
% File name: M370Exam1.m
% Name: Troy Chin
% Date: 10/02/2025
% ------------------------------

% Define data arrays
x = [0.0010 0.3832 0.7655 0.9566 1.1477 1.5299 1.7211 1.9122 2.2944 2.4855 2.8678 3.2500]';
y = [0.0017 0.3057 0.2829 0.2407 0.1966 0.1214 0.0930 0.0703 0.0391 0.0288 0.0154 0.0081]';

% Define coefficients
X = x; Y = y;
[A,B] = mylsfit(X, Y);
C = A; D = B; LY = A*X + B; % Line of best fit/Least Squares fit
R_squared = myRsqfit(y, LY); % Calculate R^2
FF = D*x*exp(C*x.^2);
MSE = mean((F - FF).^2); % Mean squared error

% Display results
fprintf('Coefficients: C = %.4f, D = %.4f\n', C, D);
fprintf('R-squared: %.4f\n', R_squared);
fprintf('Mean Squared Error: %.4f\n', MSE);

% Plot results
figure(2);clf;
hold on;
plot(n, F, '*', 'Color', 'r', 'Linewidth', 2.5);
plot(n, FF, 'Color', 'b', 'Linestyle', '-', 'Linewidth', 2.5); grid on;
legend('Data', ...
    ['Fitted Model with MSE = ', num2str(MSE), ' and R^2 = ', num2str(R_squared)], ...
    'Location', 'NW');
title(['Troy Chin MSE = ', num2str(MSE), ' and R^2 = ', num2str(R_squared)]); grid on;




