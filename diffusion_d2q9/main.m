clear;close all;clc;

% With non-equal weights, things get more complicated in the main solver
% loop, compared to D1Q2 and D2Q4, which both have lattice links of equal
% weight.

% Input parameters.
nodes = [100, 100]; % x nodes, y nodes.
dh = 1; % dh = dx = dy.
dt = 1; % timestep.
timesteps = 400;
alpha = 0.25; % Physical constant.
twall = 1.0; % Left-hand wall temperature.

% Constants.
w = zeros(9,1);
w(1) = 4/9;
w(2:5) = 1/9;
w(6:9) = 1/36;

% Derived inputs.
L = dh*[nodes(1)-1, nodes(2)-1] ; % x , y dimensions of physical domain.
c = dh/dt;
csq = c.^2;
omega = 1 / (3*alpha/(dt*csq)+0.5);

% Initialize.
x = linspace(0,L(1),nodes(1))';
y = linspace(0,L(2),nodes(2))';
rho = zeros(nodes(2),nodes(1));
f = zeros(9,nodes(2),nodes(1));
% BC initialize.
for k = 1:9
    f(k,:,1) = w(k)*twall; 
end

% Main loop.
for iter = 1:timesteps
    % Collision
    rho = sum(f,1);
    feq = w*rho;
    f = omega*feq + (1-omega)*f;
    % Streaming
    f(2,:,2:end) = f(2,:,1:end-1); % East vector.
    f(3,2:end,:) = f(3,1:end-1,:); % North vector.
    f(4,:,1:end-1) = f(4,:,2:end); % West vector.
    f(5,1:end-1,:) = f(5,2:end,:); % South vector.
    f(6,2:end,2:end) = f(6,1:end-1,1:end-1); % Northeast vector.
    f(7,2:end,1:end-1) = f(7,1:end-1,2:end); % Northwest vector.
    f(8,1:end-1,1:end-1) = f(8,2:end,2:end); % Southwest vector.
    f(9,1:end-1,2:end) = f(9,2:end,1:end-1); % Southeast vector.
    % Boundary conditions
    
    % Compute solution
    rho = f1+f2+f3+f4;
end

% Plotting results!
[X, Y] = meshgrid(x,y);
contour(X, Y, rho);
title('Solution');
xlabel('x');
ylabel('y');





