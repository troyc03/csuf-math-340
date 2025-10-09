A = randn(8);
B = randn(8);
A = A + A';
B = B + B';

t=0:0.01:1;
e = [];
for tt = t
    e = [e eig(A + tt*B)];
end

figure(11)
plot(t, e)