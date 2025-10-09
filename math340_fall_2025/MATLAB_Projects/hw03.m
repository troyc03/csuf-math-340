% ------------------------------------------- %

% Homework 3
% NAME: Troy Chin
% DATE: 09/27/2025

% ------------------------------------------- %

% Question 1

x = 0:6; 
y = [15 18 21 20 17 14 11];
u = linspace(min(x), max(x), 200);
v_linear = interp1(x,y,u,'linear');
v_poly = polyval(polyfit(x,y,length(x)-1), u);
v_spline = spline(x,y,u);
v_pchip = pchip(x,y,u);

figure;
subplot(2,2,1);
plot(x,y,'o', u,v_linear,'g','LineWidth', 1.5);
title('Linear Interpolant');
legend('Data','Linear','Location','Best');
grid on;

subplot(2,2,2);
plot(x,y,'o',u,v_poly,'r','LineWidth', 1.5);
title('Polynomial Interpolant');
legend('Data','Poly', 'Location', 'Best');
grid on;

subplot(2,2,3);
plot(x,y,'o',u,v_spline, 'c', 'LineWidth', 1.5);
title('Cubic Spline Interpolant');
legend('Data', 'Spline','Location','Best');
grid on;

subplot(2,2,4);
plot(x,y,'o',u,v_pchip,'m','Linewidth', 1.5);
title('PCHIP');
legend('Data','PCHIP','Location','Best');
grid on;

sgtitle('Comparison of Interpolants');

% ------------------------------------------- %

% Question 3

W = [10 27 2001 5 10 4 8
    11 19 2001 7 4 5 11
    12 03 2001 8 12 6 4
    12 20 2001 10 14 8 7
    01 09 2002 12 13 10 3
    01 23 2002 14 8 12 0
    03 06 2002 16 10 13 10];

t = datetime(W(:, [3,1,2]));

Tom = W(:, 4) + W(:, 5)/16;
Ben = W(:, 6) + W(:, 7)/16;

figure(2); clf;
plot(t, Tom, 'o-', 'LineWidth', 1.5);
hold on;
plot(t, Ben, 's--', 'LineWidth', 1.5);
title('Tom and Ben Data Over Time');
xlabel('Date');
ylabel('Value');
legend('Tom', 'Ben', 'Location', 'Best');
grid on;
hold off;

% Define interpolant data
x = [-1 -0.5 0 0.5 1];
y = [0 -0.5 1 0.5 0];
xx = linspace(-1,1,200);

% Interpolants
yy_piecelin = interp1(x,y,xx,'linear');     % Piecewise linear
p = polyfit(x,y,length(x)-1);               % Polynomial coefficients
yy_poly = polyval(p,xx);                    % Polynomial values
yy_spline = spline(x,y,xx);                 % Cubic spline
yy_pchip  = pchip(x,y,xx);                  % PCHIP

% Plot results
figure(3); clf;

subplot(2,2,1)
plot(x,y,'o',xx,yy_piecelin,'g','LineWidth',1.5)
title('Piecewise Linear')
legend('Data','piecelin','Location','Best')
grid on

subplot(2,2,2)
plot(x,y,'o',xx,yy_poly,'r','LineWidth',1.5)
title('Polynomial')
legend('Data','polyinterp','Location','Best')
grid on

subplot(2,2,3)
plot(x,y,'o',xx,yy_spline,'c','LineWidth',1.5)
title('Cubic Spline')
legend('Data','splinetx','Location','Best')
grid on

subplot(2,2,4)
plot(x,y,'o',xx,yy_pchip,'m','LineWidth',1.5)
title('PCHIP')
legend('Data','pchiptx','Location','Best')
grid on

sgtitle('Comparison of Interpolants for [-1,1]')

% ------------------------------------------- %

% Question 4
figure('position',get(0,'screensize'))
axes('position',[0 0 1 1])
[x,y] = ginput;
n = length(x);
s = (1:n)';
t = (1:.05:n)';
u = spline(s,x,t);
clf reset;
plot(x,y,'.',u,v,'-');

% ------------------------------------------- %

% Question 7
% This problem was solved on paper. 

% PROOF:
% simplified proof for 3.7: 

% Suppose P(x) and Q(x) are degree (n-1) polynomials and agree at n distinct values of x.
% This implies that the polynomial (P(x) - Q(x)) is also degree (n-1), 
% and has zeros at all of those values of x. However, a degree (n-1) polynomial cannot have 
% (n) zeros, unless it is identically 0. Therefore, P(x) - Q(x) = 0, and P(x) = Q(x).

function rungeinterp(arg)
%RUNGEINTERP  Runge's polynomial interpolation example.
%   F(x) = 1/(1+25*x^2)
%   Polynomial interpolation at equally spaced points, -1 <= x <= 1.
%   Does interpolant converge as number of points is increased?

%   Copyright 2014 Cleve Moler
%   Copyright 2014 The MathWorks, Inc.


if nargin == 0

   % Initialize plot and uicontrols

   shg
   clf reset
   set(gcf,'numbertitle','off','menu','none', ...
       'name','Runge''s interpolation example')
   n = 1;
   u = -1.1:.01:1.1;
   z = rungerat(u);
   h.plot = plot(u,z,'-', 0,1,'o', u,z,'-');
   set(h.plot(1),'color',[.6 .6 .6]);
   set(h.plot(2),'color','blue');
   set(h.plot(3),'color',[0 2/3 0]);
   axis([-1.1 1.1 -0.1 1.1])
   title('1/(1+25*x^2)','interpreter','none')

   h.minus = uicontrol('units','norm','pos',[.38 .01 .06 .05], ...
          'fontsize',12,'string','<','callback','rungeinterp(''n--'')');
   h.n = uicontrol('units','norm','pos',[.46 .01 .12 .05], ...
          'fontsize',12,'userdata',n,'callback','rungeinterp(''n=1'')');
   h.plus = uicontrol('units','norm','pos',[.60 .01 .06 .05], ...
          'fontsize',12,'string','>','callback','rungeinterp(''n++'')');
   h.close = uicontrol('units','norm','pos',[.80 .01 .10 .05], ...
          'fontsize',12,'string','close','callback','close');

   set(gcf,'userdata',h)
   arg = 'n=1';
end

% Update plot.

h = get(gcf,'userdata');

% Number of interpolation points.

n = get(h.n,'userdata');
switch arg
   case 'n--', n = n-2;
   case 'n++', n = n+2;
   case 'n=1', n = 1;
end
set(h.n,'string',['n = ' num2str(n)],'userdata',n);
if n==1
   set(h.minus,'enable','off');
else
   set(h.minus,'enable','on');
end

if n == 1;
   x = 0;
else
   k = 1:n;
   x = cos((2*k - 1) * pi / (2*n));
end
y = rungerat(x);
u = get(h.plot(1),'xdata');
v = polyinterp(x,y,u);
set(h.plot(2),'xdata',x,'ydata',y);
set(h.plot(3),'xdata',u,'ydata',v);

end

% ------------------------

function y = rungerat(x);
y = 1./(1+25*x.^2);
end

% ------------------------------------------- %

% Question 9
% The function f(x) = 1/(25+x^2) will not converge at
% equally spaced points on the interval [-1,1]. However,
% if we used Chebyshev nodes then the interpolant function
%  P_n(x) will approach f(x) uniformly in the interval [-1,1].

% ------------------------------------------- %

% Question 11

function [v, v_prime] = splinetx(x,y,u)
%SPLINETX  Textbook spline function.
%  v = splinetx(x,y,u) finds the piecewise cubic interpolatory
%  spline S(x), with S(x(j)) = y(j), and returns v(k) = S(u(k)).
%
%  See SPLINE, PCHIPTX.

%   Copyright 2014 Cleve Moler
%   Copyright 2014 The MathWorks, Inc.

%  First derivatives

   h = diff(x);
   delta = diff(y)./h;
   d = splineslopes(h,delta);

%  Piecewise polynomial coefficients

   n = length(x);
   c = (3*delta - 2*d(1:n-1) - d(2:n))./h;
   b = (d(1:n-1) - 2*delta + d(2:n))./h.^2;

%  Find subinterval indices k so that x(k) <= u < x(k+1)

   k = ones(size(u));
   for j = 2:n-1
      k(x(j) <= u) = j;
   end

%  Evaluate spline

   s = u - x(k);
   v = y(k) + s.*(d(k) + s.*(c(k) + s.*b(k)));
   v_prime = d(k) + 2*c(k).*s + 3*b(k).*s.^2;


% -------------------------------------------------------

function d = splineslopes(h,delta)
%  SPLINESLOPES  Slopes for cubic spline interpolation.
%  splineslopes(h,delta) computes d(k) = S'(x(k)).
%  Uses not-a-knot end conditions.

%  Diagonals of tridiagonal system

   n = length(h)+1;
   a = zeros(size(h)); b = a; c = a; r = a;
   a(1:n-2) = h(2:n-1);
   a(n-1) = h(n-2)+h(n-1);
   b(1) = h(2);
   b(2:n-1) = 2*(h(2:n-1)+h(1:n-2));
   b(n) = h(n-2);
   c(1) = h(1)+h(2);
   c(2:n-1) = h(1:n-2);

%  Right-hand side

   r(1) = ((h(1)+2*c(1))*h(2)*delta(1)+ ...
          h(1)^2*delta(2))/c(1);
   r(2:n-1) = 3*(h(2:n-1).*delta(1:n-2)+ ...
              h(1:n-2).*delta(2:n-1));
   r(n) = (h(n-1)^2*delta(n-2)+ ...
          (2*a(n-1)+h(n-1))*h(n-2)*delta(n-1))/a(n-1);

%  Solve tridiagonal linear system

   d = tridisolve(a,b,c,r);

% -------------------------------------------------------

function [v, v_prime] = pchiptx(x,y,u)
%PCHIPTX  Textbook piecewise cubic Hermite interpolation.
%  v = pchiptx(x,y,u) finds the shape-preserving piecewise cubic
%  interpolant P(x), with P(x(j)) = y(j), and returns v(k) = P(u(k)).
%
%  See PCHIP, SPLINETX.
 
%   Copyright 2014 Cleve Moler
%   Copyright 2014 The MathWorks, Inc.

%  First derivatives
 
   h = diff(x);
   delta = diff(y)./h;
   d = pchipslopes(h,delta);

%  Piecewise polynomial coefficients

   n = length(x);
   c = (3*delta - 2*d(1:n-1) - d(2:n))./h;
   b = (d(1:n-1) - 2*delta + d(2:n))./h.^2;

%  Find subinterval indices k so that x(k) <= u < x(k+1)

   k = ones(size(u));
   for j = 2:n-1
      k(x(j) <= u) = j;
   end

%  Evaluate interpolant

   s = u - x(k);
   v = y(k) + s.*(d(k) + s.*(c(k) + s.*b(k)));
   v_prime = d(k) + 2*c(k).*s + 3*b(k).*s.^2;


% -------------------------------------------------------

function d = pchipslopes(h,delta)
%  PCHIPSLOPES  Slopes for shape-preserving Hermite cubic
%  pchipslopes(h,delta) computes d(k) = P'(x(k)).

%  Slopes at interior points
%  delta = diff(y)./diff(x).
%  d(k) = 0 if delta(k-1) and delta(k) have opposites
%         signs or either is zero.
%  d(k) = weighted harmonic mean of delta(k-1) and
%         delta(k) if they have the same sign.

   n = length(h)+1;
   d = zeros(size(h));
   k = find(sign(delta(1:n-2)).*sign(delta(2:n-1))>0)+1;
   w1 = 2*h(k)+h(k-1);
   w2 = h(k)+2*h(k-1);
   d(k) = (w1+w2)./(w1./delta(k-1) + w2./delta(k));

%  Slopes at endpoints

   d(1) = pchipend(h(1),h(2),delta(1),delta(2));
   d(n) = pchipend(h(n-1),h(n-2),delta(n-1),delta(n-2));

% -------------------------------------------------------

function d = pchipend(h1,h2,del1,del2)
%  Noncentered, shape-preserving, three-point formula.
   d = ((2*h1+h2)*del1 - h1*del2)/(h1+h2);
   if sign(d) ~= sign(del1)
      d = 0;
   elseif (sign(del1)~=sign(del2))&(abs(d)>abs(3*del1))
      d = 3*del1;
   end
% -------------------------------------------------------

% Question 18
% Using a Vandermonde matrix is generally a bad idea because
% the higher the degree of your interpolant polynomial the more
% unstable your system will be. The Basic Fitting tool will
% plot a line of best fit of the census data. We must choose
% mu = mean(t) and sigma = std(t) so that we can center and
% scale the data. 

% -------------------------------------------------------

