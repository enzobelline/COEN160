function F = problem1(t,x)
%PROB1
    F = [0; 0];
    F(1) = -2*x(1) + 2*x(2);
    F(2) = -2*x(1) - 7*x(2);
end