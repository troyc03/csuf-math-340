function e = qr1(A) 

% Computes the eigenvalues of A, if the eigenvaleus are real numbers. 

n = size(A,1); 
I = eye(n); 
e = []; 

for k = 1:n
    while max(abs(A(n, 1:(n-1)))) > 1e-8 % 10 to the minus 8
        % repeats operations until the last row is all zero except for A(n,n) 
        % checks the last row to be all zero, excluding the last element
        s = A(n,n); 
        [Q, R] = qr(A - s*I);
        A = R*Q + s*I; 
    end
        e = [e A(n,n)];
        
        A = A(1:(n-1), 1:(n-1));
        n = n-1;
        I = eye(n); 
end


% A_k = Q_k*R_k + s*I
% A_(k+1) = R_k*Q_k + s*I 
% Q_(k+1)*R_(k+1) = A_(k+1)-s*I