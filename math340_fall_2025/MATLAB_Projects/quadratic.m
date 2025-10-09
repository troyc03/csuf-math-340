for k = 1:3
    disp('For the equation ax^2 + bx + c')
    a = input('Enter the value of a: ');
    b = input('Enter the value of b: ');
    c = input('Enter the value of c: ');
    D = b^2-4*a*c; % Discriminant
    if D < 0
        disp('The equation has no real roots.')
    elseif D==0
        disp('The equation has one real root.')
        root = -b / (2*a);
        disp(['Root: ', num2str(root)]);
    else
        root1 = (-b + sqrt(D)) / (2*a);
        root2 = (-b - sqrt(D)) / (2*a);
        disp(['Roots: ', num2str(root1), ' and ', num2str(root2)]);
    end
end

