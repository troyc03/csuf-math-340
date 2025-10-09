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



