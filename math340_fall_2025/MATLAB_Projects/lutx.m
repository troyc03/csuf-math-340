function [L,U,p,sig] = lutx(A)
%LUTX  Triangular factorization with permutation sign
%   [L,U,p,sig] = lutx(A) produces a unit lower triangular matrix L,
%   an upper triangular matrix U, a permutation vector p, and a scalar sig,
%   so that L*U = A(p,:) and sig = +1 or -1 depending on the parity of p.

[n,n] = size(A);
p = (1:n)';   % permutation vector
sig = 1;      % start with +1

for k = 1:n-1
   % Find index of largest element below diagonal in k-th column
   [r,m] = max(abs(A(k:n,k)));
   m = m+k-1;

   % Skip elimination if column is zero
   if (A(m,k) ~= 0)
      % Swap pivot row
      if (m ~= k)
         A([k m],:) = A([m k],:);
         p([k m]) = p([m k]);
         sig = -sig;   % flip sign for each row swap
      end

      % Compute multipliers
      i = k+1:n;
      A(i,k) = A(i,k)/A(k,k);

      % Update remainder of matrix
      j = k+1:n;
      A(i,j) = A(i,j) - A(i,k)*A(k,j); 
   end
end

% Separate result
L = tril(A,-1) + eye(n,n);
U = triu(A);
end

function d = mydet(A)
%MYDET  Determinant using LU factorization with row pivoting
%   d = mydet(A) uses lutx to compute the determinant.

A = [2 1 1; 4 -6 0; -2 7 2];
[L,U,p,sig] = lutx(A);

% Product of diagonal entries of U
A = [2 1 1; 4 -6 0; -2 7 2];
d = sig * prod(diag(U));
end

