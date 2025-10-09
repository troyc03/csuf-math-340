function fractal_fern
    % Clear the command window and figure
    shg;
    clf reset;

    % Set up the figure window properties
    set(gcf, 'color', 'white', ...
        'name', 'Fractal Fern', ...
        'numbertitle', 'off', ...
        'menubar', 'none');

    % Define the four affine transformations for the fern
    % Each matrix represents a transformation: [ a b; c d ] * x + [ e; f ]
    % The probabilities determine how often each transformation is chosen.
    A1 = [0.85 0.04; -0.04 0.85];
    b1 = [0; 1.6];
    p1 = 0.85;

    A2 = [0.20 -0.26; 0.23 0.22];
    b2 = [0; 1.6];
    p2 = 0.07;

    A3 = [-0.15 0.28; 0.26 0.24];
    b3 = [0; 0.44];
    p3 = 0.07;

    A4 = [0 0; 0 0.16];
    b4 = [0; 0];
    p4 = 0.01;

    % Initialize the starting point
    x = [0; 0];
    
    % Use a pre-allocated array for better performance
    N = 100000; % Number of points to generate
    points = zeros(2, N);

    % Create the plot handle
    darkgreen = [0 2/3 0];
    h = plot(NaN, NaN, '.', 'markersize', 1, 'color', darkgreen);
    
    % Set the axis limits and turn off the axis labels
    axis([-3 3 0 10]);
    axis off;
    
    % Create a stop button
    stop = uicontrol('style', 'toggle', 'string', 'Stop', 'background', 'white');
    
    % Main loop to generate the points
    for i = 1:N
        if get(stop, 'Value')
            break; % Break the loop if the stop button is pressed
        end
        
        % Choose a transformation based on the probabilities
        r = rand;
        if r < p1
            x = A1 * x + b1;
        elseif r < p1 + p2
            x = A2 * x + b2;
        elseif r < p1 + p2 + p3
            x = A3 * x + b3;
        else
            x = A4 * x + b4;
        end
        
        % Store the new point
        points(:, i) = x;
        
        % Update the plot every 100 iterations for efficiency
        if mod(i, 100) == 0
            set(h, 'XData', points(1, 1:i), 'YData', points(2, 1:i));
            drawnow limitrate; % Update the figure window
        end
    end

    % Plot all points at the end for the final render
    set(h, 'XData', points(1, 1:i), 'YData', points(2, 1:i));
    drawnow;
    
    % Reset the stop button
    set(stop, 'Value', 0);
end
