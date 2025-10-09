% ------------------------------------------- %

% Homework 4
% NAME: Troy Chin
% DATE: 09/30/2025

% ------------------------------------------- %

% Question 2
tol = 1e-6; a = 1; b = 5;
pi_val = pi;

while (b-a)/2 > tol
    c = (b+a)/2;
    A = hilb(5);
    A(1,1) = c;
    lam_max = max(eig(A));
    if lam_max - pi_val > 0
        b = c;
    else
        a = c;
    end
end

% Output first entry of A using Bisection method
a11 = (a+b)/2;
fprintf('First entry of A to six decimal places: %.6f\n', a11);

% ------------------------------------------- %

% Question 4.10
J0 = @(x)besselj(0,x);
Y0 = @(x)bessely(0,x);
xx = 0:.01:10; close all;
xlabel('x');
ylabel('y');
title('Newton Method Iterations');
plot(xx,J0(xx), xx, Y0(xx), 'linewidth', 2);
grid on;

