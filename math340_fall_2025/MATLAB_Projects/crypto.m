function y = crypto(x)
% Performs some cryptographic operation y = crypto(x) on an input x
% Choose two characters c1, c2 and a modulo operator
p = 97; % Modulo operator
c1 = char(169); % character 1
c2 = char(174); % character 2
x(x==c1) = 127;
x(x==c2) = 128;
x = mod(real(x-32), p); % ASCII to numerical conversion
n = 2*floor(length(x)/2); % Matrix-vector product
X = reshape(x(1:n), 2, n/2);
% Encode with matrix multiplication modulo p.
A = [71 2; 2 26];
Y = mod(A*X,p);
y = reshape(Y,1,n);
% If the length of the input is odd, encode the last character
if length(x) > n
    y(n+1) = mod((p-1)*x(n+1),p);
end
y = char(y+32);
y(y==127) = c1;
y(y==128) = c2;

