    % Coefficient Matrix
function F = problem2(t,x)
%PROB1
    F = [0; 0];
    F(1) = 4*x(1) - 3*x(2);
    F(2) = 6*x(1) - 5*x(2);
end