%% Part 1
%A)
X0 = [-1.5 0 5];
xs = [ ];

for x0= X0
    x2 = problem1a(x0);
    x1 = x2 - 3;
    fprintf("Equlibria:\n   x1=%f\n   x2=%f\n" , x1 , x2);
    xs = [xs , [x1 ; x2]];
    snapnow()  % take a snapshot image for inclusion in published doc
end

%B)
for x = xs
    J = [-1 1; x(2)^2+2*x(1) 2*x(1).*x(2)];
    v = eig(J); %find the eigenvalue(s) of J
    fprintf("Equlibria:\n   x1=%f\n   x2=%f\n" , x(1) , x(2));      %disp
    fprintf("have eigenvalues at %f and %f so it is" , v(1) , v(2));%disp
    if ((v(1) < 0) && (v(2) < 0))
        fprintf("Stable\n");
    else %if ((v(1) > 0) || (v(2) > 0))
        fprintf("Unstable\n");
    end
    fprintf('\n'); %create new line
end

%C)
figure
X0 = [[0 ; 1] [-3 ; 0] [-1 ; 3] [-4 ; -1] [-2 ; 3] [-2 ; -1]];
num = 5; %End t
subplot(2,1,1);
for x0 = X0
    [t, x] = ode45(@problem1c, [0,num], x0);
    x1 = x(:,1);
    x2 = x(:,2);
    plot(x1, x2);
    hold on
end
num = 3;
X0 = [[-5 ; -3] [5 ; 3]];
title("Stable Equlibria");
xlabel("x1");
ylabel("x2");

subplot(2,1,2);
for x0 = X0
    [t, x] = ode45(@problem1c, [0 , num], x0);
    x1 = x(:,1);
    x2 = x(:,2);
    plot(x1, x2);
    hold on
end
xlabel("x1");
ylabel("x2");
title("Unstable Equlibria");


%% Part 2
%a)
X0 = [1 -1 5];
xs = [ ];
for x0 = X0
    x  = problem2a(x0);
    xs = [xs , x];
    fprintf("Equlibria at \nx1=%f\nx2=%f\n",x(1), x(2));
end
fprintf('\n');

%b)
for x = xs
    J      = zeros(2,2);
    J(1,1) = (x(1).*(exp(x(2)-1)).*(1/3)) + 2*x(2)*x(1)-(1/3)*exp(x(2)-x(1));
    J(1,2) = x(1)^2 - (1/3)*(x(1)*exp(x(2)-x(1)));
    J(2,1) = 2*x(2)*x(1) + 1 ;
    J(2,2) = x(1)^2 - 1;
    v = eig(J); %eigenvalue of J
    fprintf("Equlibria at \nx1=%f \nx2=%f\n",x(1), x(2));
    fprintf("have eigen values at %f and %f which is ", v(1), v(2));
    if v(1) < 0 && v(2) < 0
        fprintf("Stable\n");
    else    %if ((v(1) > 0) || (v(2) > 0))
        fprintf("Unstable\n");
    end
    fprintf('\n');
end

%c)
figure
X0 = [[1 ; 1] [-1 ; -1]];
num = 5;

subplot(2,1,1);
for x0 = X0
    [t, x] = ode45(@problem2c, [0 , num], x0);
    x1 = x(:,1);
    x2 = x(:,2);
    plot(x1, x2);
    hold on
end
xlabel("x1");
ylabel("x2");
title("Stable Equlibria");

subplot(2,1,2);
X0 = [[-1 ; 1] [0.2 ; 0.1]];
num = 3;
for x0 = X0
    [t, x] = ode45(@problem2c, [0,num], x0);
    x1 = x(:,1);
    x2 = x(:,2);
    plot(x1, x2);
    hold on
end
xlabel("x1");
ylabel("x2");
title("Unstable Equlibria");
%% Part 3
%A)
xs = [];
X0 = [[-1 ; 1] [1 ; -1]];
for x0 = X0
    x  = problem3a(x0);
    fprintf("Equlibria at \nx1=%f\nx2=%f\n",x(1), x(2));
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
    
    fprintf("Equlibria at \nx1=%f\nx2=%f\n",x(1), x(2));
    fprintf("have eigenvalues at %f and %f which is ", v(1), v(2));
    
    if v(1) < 0 && v(2) < 0
        fprintf("Stable\n");
    else  %if ((v(1) > 0) || (v(2) > 0))
        fprintf("Unstable\n");
    end
    fprintf('\n');
end

%% functions
function x = problem1a(x0)
    x  = x0;
    f  = 8 + x^3 - 2*x^2 - 6*x;
    fp = 3*x^2 - 4*x-6;
    dx = f/fp;
    x  = x0 - dx;
    
    while (abs(dx) > 0.0001)
        f  = 8 + x^3 - 2*x^2 - 6*x;
        fp = 3*x^2 - 4*x-6;
        dx = f/fp;
        x  = x - dx;
    end
end

function F = problem1c(t,x)
    F    = [0;0];
    F(1) = -x(1)+x(2)-3;
    F(2) = x(1)^2+x(1)*x(2)^2-1;
end

function x = problem2a(x0)
    f    = [0;0];
    x    = [x0;x0];
    J    = zeros(2,2);
    f(1) = x(1)^2*x(2) - (x(1)*exp(x(2)-x(1)))/3;
    f(2) = x(1) - (1 - x(1)^2)*x(2);
    
    while (norm(f,inf)>0.0001)
        J(1,1) = 2*x(2)*x(1) - (1/3)*exp(x(2)-x(1)) + (1/3)*(x(1)*exp(x(2)-x(1)));
        J(1,2) = x(1)^2 - (1/3)*x(1)*exp(x(2)-x(1));
        J(2,1) = 1 + 2*x(2)*x(1);
        J(2,2) = x(1)^2 - 1 ;
        
        dx   = J\f;
        x    = x - dx;
        f(1) = x(1)^2*x(2) - (x(1)*exp(x(2)-x(1)))/3;
        f(2) = x(1) - (1 - x(1)^2)*x(2);
    end
end

function F = problem2c(t,x)
    F    = [0;0];
    F(1) = x(1)^2*x(2) - x(1)*exp(x(2)-x(1))/3;
    F(2) = x(1) - (1 - x(1)^2)*x(2);
end

function x = problem3a(x0)
    F    = [0 ; 0];
    x    = x0;
    J    = zeros(2,2);
    F(1) = x(1)^3 + 5*(x(1)^2)*x(2) + 8*x(1)*x(2)^2 + 4*x(2)^3 + 2*x(1) + 2*x(2) + 1;
    F(2) = -x(1)^2 + 2*x(1)*x(2) + 5*x(2)^2 - 1;
    while (norm(F,inf)>1e-4)
        J(1, 1) = 3*x(1)^2 + 10*x(1)*x(2) + 8*x(2)^2 + 2;
        J(1, 2) = 5*x(1)^2 + 16*x(1)*x(2) + 12*x(2)^2 + 2;
        J(2, 1) = -2^x(1) + 2*x(2);
        J(2, 2) = 2*x(1) + 10*x(2);
        dx   = J\F;
        x    = x - dx;
        F(1) = x(1)^3 + 5*(x(1)^2)*x(2) + 8*x(1)*x(2)^2 + 4*x(2)^3 + 2*x(1) + 2*x(2) + 1;
        F(2) = -x(1)^2 + 2*x(1)*x(2) + 5*x(2)^2 - 1;
    end
end