function L = compute_collatz(n)
    % Initialize length of the sequence
    L = 1;
    while L ~= 1
        if mod(n, 2) == 0
            n = n / 2;
        else
            n = 3 * n + 1;
        end
        L = L + 1; % Increment length for each term
    end
end
% Compute L(n) for [1, 1000]
n_values = 1:1000;
L_values = arrayfun(@compute_collatz, n_values);
% Plotting sequence
figure;
xlabel('n');
ylabel('L(n)');
title("Collatz Sequence for [1,1000]");
% Finding the maximum value of L(n) and the corresponding n
[max_L, max_n] = max(L_values);
fprintf("The maximum value of L(n) is %d and it occurs at n = %d\n", max_L, max_n);

function collatz_sequence(n)
% Initialize the sequence array
sequence = n;
    while n ~= 1
        if mod(n, 2) == 0
            n = n/2;
        else
            n = 3*n + 1;
        end
        sequence = [sequence, n];
    end
    
    % Plot sequence 
    plot(sequence, '-o');
    title(['Collatz Sequence Starting at ', num2str(sequence(1))]);
    xlabel('Step');
    ylabel('Value');
    grid on;
end

% Example usage
collatz_sequence(5);
collatz_sequence(10);
collatz_sequence(20);
collatz_sequence(40);



