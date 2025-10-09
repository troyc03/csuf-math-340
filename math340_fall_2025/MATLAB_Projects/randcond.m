% RANDNCOND  Condition of random matrices with output estimates

nmax = 100;
nsamp = 200;   % number of random matrices per size
sizes = 5:5:nmax;

p = 2;         % initial guess for bounding lines
c1 = 0.1;
c2 = 10;

% Parallel bounding lines
n = 2:nmax;
kappalo = c1*n.^p;
kappahi = c2*n.^p;

shg
clf reset
h = loglog(n,[kappalo; kappahi],'-',nmax,NaN,'.');
set(h(1:2),'color',[0 2/3 0]);   % green bounds
set(h(3),'color','blue')
set(gca,'xtick',[2:2:10 20:20:nmax])
kappamax = 1.e6;
axis([2 nmax 2 kappamax])
hold on

% ---------------------------
% Collect condition numbers
% ---------------------------
all_n = [];
all_kappa = [];

for n = sizes
    for t = 1:nsamp
        A = randn(n,n);
        kappa = condest(A,1);
        loglog(n,kappa,'.','color','blue')
        all_n(end+1,1) = n;
        all_kappa(end+1,1) = kappa;
    end
    drawnow
end

hold off

% ---------------------------
% Fit log-log regression
% ---------------------------
logn = log(all_n);
logk = log(all_kappa);

coeffs = polyfit(logn,logk,1);
p_est = coeffs(1);
c_est = exp(coeffs(2));

% Estimate c1 and c2 from percentiles
kappa_p5 = zeros(size(sizes));
kappa_p95 = zeros(size(sizes));
for j = 1:length(sizes)
    k = all_kappa(all_n==sizes(j));
    kappa_p5(j) = prctile(k,5);
    kappa_p95(j) = prctile(k,95);
end
c1_est = median(kappa_p5 ./ (sizes'.^p_est));
c2_est = median(kappa_p95 ./ (sizes'.^p_est));

% ---------------------------
% Print results
% ---------------------------
fprintf('Estimated exponent p ≈ %.2f\n',p_est);
fprintf('Estimated c (for median fit) ≈ %.2f\n',c_est);
fprintf('Estimated c1 ≈ %.2f\n',c1_est);
fprintf('Estimated c2 ≈ %.2f\n',c2_est);
