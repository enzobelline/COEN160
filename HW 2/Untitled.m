%% Part 1
%A)
X0 = [-1 0 5];
xs = [];
for x0=X0
    x2 = prob1a(x0);
    x1 = x2 - 3;
    fprintf("Equlibrium:\n   x1=%f\n   x2=%f\n",x1, x2);
    xs = [xs, [x1;x2]];
    snapnow();
end

%B)
for x = xs
    J = [-1 1; 2*x(1)+x(2)^2 2*x(1)*x(2)];
    v = eig(J);
    fprintf("Equlibrium:\n   x1=%f\n   x2=%f\n",x(1), x(2));
    fprintf("with eigen values at %f and %f so it is", v(1), v(2));
    if v(1) < 0 && v(2) < 0
        fprintf("stable\n");
    else
        fprintf("unstable\n");
    end
    fprintf('\n');
end

%C)
X0 = [[0;1] [-3;0] [-1;3] [-4;-1] [-2;3] [-2;-1]];
tend = 5;
figure();
subplot(211);
for x0 = X0
    [t, x] = ode45(@prob1c, [0,tend], x0);
    x1 = x(:,1);
    x2 = x(:,2);
    plot(x1, x2);
    hold on
end
title("Stable Equlibrium");
xlabel("x1");
ylabel("x2");
X0 = [[-5;-3] [5;3]];
tend = 3;

subplot(212);
for x0 = X0
    [t, x] = ode45(@prob1c, [0,tend], x0);
    x1 = x(:,1);
    x2 = x(:,2);
    plot(x1, x2);
    hold on
end
title("Unstable Equlibrium");
xlabel("x1");
ylabel("x2");

%% Part 2
%A)
X0 = [[1;1] [-1;-1] [5;5] [-0.1;-0.5]];
xs = [];
for x0=X0
    x = prob2a(x0);
    fprintf("Equlibrium at \nx1=%f\nx2=%f\n",x(1), x(2));
    xs = [xs, x];
end
fprintf('\n');

%B)
for x = xs
    J = zeros(2,2);
    J(1,1) = 2*x(2)*x(1)-exp(x(2)-x(1))/3 + (x(1)*exp(x(2)-x(1)))/3;
    J(1,2) = x(1)^2 - x(1)*exp(x(2)-x(1))/3;
    J(2,1) = 1 + 2*x(2)*x(1);
    J(2,2) = -1 + x(1)^2;
    v = eig(J);
    fprintf("Equlibrium at \nx1=%f\nx2=%f\n ",x(1), x(2));
    fprintf("has eigen values at %f and %f which is ", v(1), v(2));
    if v(1) < 0 && v(2) < 0
        fprintf("stable\n");
    else
        fprintf("unstable\n");
    end
    fprintf('\n');
end

%C)
X0 = [[1;1] [-1;-1]];
tend = 5;
figure();
subplot(211);
for x0 = X0
    [t, x] = ode45(@prob2c, [0,tend], x0);
    x1 = x(:,1);
    x2 = x(:,2);
    plot(x1, x2);
    hold on
end
title("Stable Equlibrium");
xlabel("x1");
ylabel("x2");
X0 = [[-1;1] [0.2;0.1]];
tend = 3;
subplot(212);
for x0 = X0
    [t, x] = ode45(@prob2c, [0,tend], x0);
    x1 = x(:,1);
    x2 = x(:,2);
    plot(x1, x2);
    hold on
end
title("Unstable Equlibrium");
xlabel("x1");
ylabel("x2");

%% Part 3
%A)
X0 = [[-1;1] [1;-1]];
xs = [];
for x0=X0
    x = prob3a(x0);
    fprintf("Equlibrium at \nx1=%f\nx2=%f\n",x(1), x(2));
    xs = [xs, x];
end
fprintf('\n');

%B)
for x = xs
    J = zeros(2,2);
    J(1, 1) = 3*x(1)^2 + 10*x(1)*x(2) + 8*x(2)^2 + 2;
    J(1, 2) = 5*x(1)^2 + 16*x(1)*x(2) + 12*x(2)^2 + 2;
    J(2, 1) = -2^x(1) + 2*x(2);
    J(2, 2) = 2*x(1) + 10*x(2);
    v = eig(J);
    fprintf("Equlibrium at \nx1=%f\nx2=%f\n",x(1), x(2));
    fprintf("has eigen values at %f and %f which is ", v(1), v(2));
    if v(1) < 0 && v(2) < 0
        fprintf("stable\n");
    else
        fprintf("unstable\n");
    end
    fprintf('\n');
end

%% functions
function x = prob1a(x0)
    x = x0;
    F = x^3-2*x^2-6*x+8;
    F_p = 3*x^2-4*x-6;
    dx = F/F_p;
    x = x0 - dx;
    while (abs(dx) > 1e-4)
        F = x^3-2*x^2-6*x+8;
        F_p = 3*x^2-4*x-6;
        dx = F/F_p;
        x = x - dx;
    end
end

function F = prob1c(t,x)
    F = [0;0];
    F(1) = -x(1)+x(2)-3;
    F(2) = x(1)^2+x(1)*x(2)^2-1;
end

function x = prob2a(x0)
    F=[0;0];
    x=x0;
    J = zeros(2,2);
    F(1) = x(1)^2*x(2) - (x(1)*exp(x(2)-x(1)))/3;
    F(2) = x(1) - (1 - x(1)^2)*x(2);
    while (norm(F,inf)>1e-4)
        J(1,1) = 2*x(2)*x(1)-exp(x(2)-x(1))/3 + (x(1)*exp(x(2)-x(1)))/3;
        J(1,2) = x(1)^2 - x(1)*exp(x(2)-x(1))/3;
        J(2,1) = 1 + 2*x(2)*x(1);
        J(2,2) = -1 + x(1)^2;
        dx = J\F;
        x = x - dx;
        F(1) = x(1)^2*x(2) - (x(1)*exp(x(2)-x(1)))/3;
        F(2) = x(1) - (1 - x(1)^2)*x(2);
    end
end

function F = prob2c(t,x)
    F = [0;0];
    F(1) = x(1)^2*x(2) - x(1)*exp(x(2)-x(1))/3;
    F(2) = x(1) - (1 - x(1)^2)*x(2);
end

function x = prob3a(x0)
    F=[0;0];
    x=x0;
    J = zeros(2,2);
    F(1) = x(1)^3 + 5*(x(1)^2)*x(2) + 8*x(1)*x(2)^2 + 4*x(2)^3 + 2*x(1) + 2*x(2) + 1;
    F(2) = -x(1)^2 + 2*x(1)*x(2) + 5*x(2)^2 - 1;
    while (norm(F,inf)>1e-4)
        J(1, 1) = 3*x(1)^2 + 10*x(1)*x(2) + 8*x(2)^2 + 2;
        J(1, 2) = 5*x(1)^2 + 16*x(1)*x(2) + 12*x(2)^2 + 2;
        J(2, 1) = -2^x(1) + 2*x(2);
        J(2, 2) = 2*x(1) + 10*x(2);
        dx = J\F;
        x = x - dx;
        F(1) = x(1)^3 + 5*(x(1)^2)*x(2) + 8*x(1)*x(2)^2 + 4*x(2)^3 + 2*x(1) + 2*x(2) + 1;
        F(2) = -x(1)^2 + 2*x(1)*x(2) + 5*x(2)^2 - 1;
    end
end