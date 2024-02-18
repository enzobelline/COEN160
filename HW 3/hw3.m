clear all; 
close all; 
clc;
%% p = -3
x1 = linspace(-2,2,20);
x2 = linspace(-1.5,1.5,20);

for a = 1:length(x1)
    for b = 1:length(x2)
        A(:,a) = [x1(a) ; x2(b)];
    end
end

figure

for a = 1:length(x1)
    x = [x1(a);x2(a)];
    [T,X] = ode45(@odefcn1, [0 10], x);
    hold on
    plot(X(:,1), X(:,2), '-m', x(1), x(2), '-mo');
    grid;
    drawnow;
end

hold on
plot (0, 0, 'X');
axis ([-2.2 2.2 -2 2]);
title ('P = -3');

%% p = -1

for a = 1 : length(x1)
    for b = 1:length(x2)
        A(:,a) = [x1(a) ; x2(b)];
    end
end
figure 
for a = 1:length(x1)
    x = [x1(a);x2(a)];
    [T,X] = ode45(@odefcn2, [0 10], x);
    hold on
    plot(X(:,1), X(:,2), '-r', x(1), x(2), '-ro');
    grid;
    drawnow;
end

plot (0, 0, 'X');
axis ([-2.8 2.8 -3 3]);
title ('P = -1');
%% p = 0.5

for a = 1:length(x1)
    for b = 1:length(x2)
        A(:,a) = [x1(a) ; x2(b)];
    end
end

figure
for a = 1:length(x1)
    x = [x1(a);x2(a)]
    [T,X] = ode45(@odefcn3, [0 10], x);
    hold on
    plot(X(:,1), X(:,2), '-k', x(1), x(2), '-ko');
    grid;
    drawnow;
end


hold on
plot (0, 0, 'X');
axis ([-2.5 2.5 -1.7 1.7]);
title ('P = 0.5');

    
%% p = 3

for a = 1 : length(x1)
    for b = 1:length(x2)
        A(:,a) = [x1(a) ; x2(b)];
    end
end
figure
for a = 1:length(x1)
    x = [x1(a);x2(a)]
    [T,X] = ode45(@odefcn4, [0 10], x);
    hold on
    plot(X(:,1), X(:,2), '-g', x(1), x(2), '-go');
    grid;
    drawnow;
end


plot (0, 0, 'X');
axis ([-2.5 2.5 -1.5 1.5]);
title ('P = 3');


%%
% function definitions
function dxdt = odefcn1(t, x)
    dxdt   = zeros(2,1);
    dxdt(1)= x(2);
    dxdt(2)=-x(1) + (-3).*(1-(x(2)^2))*x(2);
end
function dxdt = odefcn2(t, x)
    dxdt   = zeros(2,1);
    dxdt(1)= x(2);
    dxdt(2)=-x(1) + (-1).*(1-(x(2)^2))*x(2);
end
function dxdt = odefcn3(t, x)
    dxdt   = zeros(2,1);
    dxdt(1)= x(2);
    dxdt(2)=-x(1) + (0.5).*(1-(x(2)^2))*x(2);
end
function dxdt = odefcn4(t,x)
    dxdt   = zeros(2,1);
    dxdt(1)= x(2);
    dxdt(2)=-x(1) + (-3).*(1-(x(2)^2))*x(2);
end

%%
%c
% Hopf Bifurcation because we have two eigenvalues cross the imaginary axis
% because of a variation fo hte system parameters
