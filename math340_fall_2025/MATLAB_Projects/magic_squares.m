function magic

% Initialize a 3x3 magic square
A = magic(3);
for k = 0:3
    rot90(A,k)
    rot90(A', k)
end