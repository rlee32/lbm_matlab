clear;close all;clc;

% D2Q9 solver
% Simple channel, with a carse-to-fine vertical grid interface in the 
% middle of the channel.
% West: fixed-velocity inlet
% North: wall
% South: wall 
% East: 0th-order-extrapolation outlet.

% We define a two grids; one fine and one coarse.
% The fine and coarse differ by only one level of refinement (4x cells).
% The fine and coarse grids overlap in one layer.
% The fine and coarse grids are both square.

% Physical parameters.
u_p = 0.1;
rho_p = 5;
L_p = 1; % Height of each channel.
nu_p = 1.568e-5; % kinematic viscosity, m^2/s.
% Grid parameters.
nodes_c = 100; % coarse nodes.
dt_c = 1; % coarse timestep.
timesteps = 10;
 
% Derived nondimensional parameters.
Re = u_p*Ly_p/nu_p;
disp(['Reynolds number: ' num2str(Re)]);
% Derived numerical parameters.
nodes_f = 2*nodes_c - 2;
dh_c = 1 / (nodes-1); % coarse spacing.
dh_f = dh_c / 2; % coarse spacing.
dt_f = dt_c / 2;
nu_lb_c = dt_c / dh_c^2 / Re; % coarse viscosity.
nu_lb_f = dt_f / dh_f^2 / Re; % fine viscosity.
tau_c = 3*nu_lb_c + 0.5; % coarse relaxation time.
tau_f = 3*nu_lb_f + 0.5; % fine relaxation time.
omega_c = 1 / ( 3*tau_c + 0.5 );
omega_f = 1 / ( 3*tau_f + 0.5 );
u_lb_c = dh_c / dt_c;
u_lb_f = dh_f / dt_f;
% Lattice link constants.
w = zeros(9,1);
w(1) = 4/9;
w(2:5) = 1/9;
w(6:9) = 1/36;
c = zeros(9,2);
c(1,:) = [0, 0];
c(2,:) = [1, 0];
c(3,:) = [0, 1];
c(4,:) = [-1, 0];
c(5,:) = [0, -1];
c(6,:) = [1, 1];
c(7,:) = [-1, 1];
c(8,:) = [-1, -1];
c(9,:) = [1, -1];

% Initialize.
rho_c = rho_p*ones(nodes_c,nodes_c);
rho_f = rho_p*ones(nodes_f,nodes_f);
u_c = u_lb_c*ones(nodes_c,nodes_c);
u_f = u_lb_f*ones(nodes_f,nodes_f);
v_c = zeros(nodes_c,nodes_c);
v_f = zeros(nodes_f,nodes_f);
f_c = zeros(nodes_c,nodes_c);
f_f = zeros(nodes_f,nodes_f);
% Wall BCs.
u_c(1,:) = 0;
v_c(1,:) = 0;
u_c(end,:) = 0;
v_c(end,:) = 0;
u_f(1,:) = 0;
v_f(1,:) = 0;
u_f(end,:) = 0;
v_f(end,:) = 0;

% Main loop.
reconstruction_time = 0;
collision_time = 0;
streaming_time = 0;
bc_time = 0;
for iter = 1:timesteps
    disp(['Running timestep ' num2str(iter)]);
    % First, we explode coarse cells and map to fine.
    % u
    u_f(1:2:end,1) = u_c(:,end);
    u_f(2:2:end,1) = u_c(2:end,end);
    u_f(1:2:end,2) = u_c(:,end);
    u_f(2:2:end,2) = u_c(2:end,end);
    % v
    v_f(1:2:end,1) = v_c(:,end);
    v_f(2:2:end,1) = v_c(2:end,end);
    v_f(1:2:end,2) = v_c(:,end);
    v_f(2:2:end,2) = v_c(2:end,end);
    % rho
    rho_f(1:2:end,1) = rho_c(:,end);
    rho_f(2:2:end,1) = rho_c(2:end,end);
    rho_f(1:2:end,2) = rho_c(:,end);
    rho_f(2:2:end,2) = rho_c(2:end,end);
    % f
    f_f(1:2:end,1) = f_c(:,end);
    f_f(2:2:end,1) = f_c(2:end,end);
    f_f(1:2:end,2) = f_c(:,end);
    f_f(2:2:end,2) = f_c(2:end,end);
    % Second, we iterate on fine cells.
    for k = 1:2
        % Stream.
        tic;
        f_f = stream(f_f);
        streaming_time = streaming_time + toc;
        % Collide.
        tic;
        f_f = collide(f_f,u_f,v_f,rho_f);
        collision_time = collision_time + toc;
        % BCs.
        tic;
        f_f = outlet_bc(f_f,'east');
        f_f = wall_bc(f_f,'south');
        f_f = wall_bc(f_f,'north');
        bc_time = bc_time + toc;
        % Density and velocity reconstruction.
        tic;
        [u_c, v_c, rho_c] = reconstruct_macro(f_c);
        reconstruction_time = reconstruction_time + toc;
    end
    % Third, we stream on coarse.
    tic;
    f_c = stream(f_c);
    streaming_time = streaming_time + toc;
    % Fourth, coalesce to coarse.
    f_c(1,end) = 0.5 * ( f_f(1,1) + f_f(1,2) );
    f_c(2:end-1,end) = 0.25 * ( f_f(2:2:end-1,1) + f_f(3:2:end-1,1)...
        + f_f(2:2:end-1,2) + f_f(3:2:end-1,2) );
    f_c(end,end) = 0.5 * ( f_f(end,1) + f_f(end,2) );
    % Fifth, collide on coarse.
    tic;
    f_c = collide(f_c,u_c,v_c,rho_c);
    collision_time = collision_time + toc;
    % Coarse BC.
    tic;
    f_c = inlet_bc(f_c, u_lb_c, 'west');
    f_c = wall_bc(f_c,'south');
    f_c = wall_bc(f_c,'north');
    bc_time = bc_time + toc;
    % Density and velocity reconstruction.
    tic;
    [u_c, v_c, rho_c] = reconstruct_macro(f_c);
    reconstruction_time = reconstruction_time + toc;
end

% Timing outputs.
total_time = reconstruction_time + collision_time + streaming_time + bc_time;
disp(['Solution reconstruction time (s): ' num2str(reconstruction_time)]);
disp(['Collision time (s): ' num2str(collision_time)]);
disp(['Streaming time (s): ' num2str(streaming_time)]);
disp(['BC time (s): ' num2str(bc_time)]);
disp(['Solution reconstruction fraction: ' num2str(reconstruction_time/total_time)]);
disp(['Collision fraction: ' num2str(collision_time/total_time)]);
disp(['Streaming fraction: ' num2str(streaming_time/total_time)]);
disp(['BC fraction: ' num2str(bc_time/total_time)]);

% Streamfunction calculation.
strf = zeros(nodes(2),nodes(1));
for i = 2:nodes(1)
    rho_av = 0.5*( rho(1,i-1) + rho(1,i) );
    strf(1,i) = strf(1,i-1) - 0.5*rho_av*( v(1,i-1) + v(1,i) );
    for j = 2:nodes(2)
        rho_m = 0.5 * ( rho(j,i) + rho(j-1,i) );
        strf(j,i) = strf(j-1,i) + 0.5*rho_m*( u(j-1,i) + u(j,i) );
    end
end

% Plotting results!
figure;
L = dh*[nodes(1)-1, nodes(2)-1] ; % x , y dimensions of physical domain.
x = linspace(0,L(1),nodes(1))';
y = linspace(0,L(2),nodes(2))';
[X, Y] = meshgrid(x,y);
contour(X(2:end,2:end), Y(2:end,2:end), strf(2:end,2:end));
title('Solution');
xlabel('x');
ylabel('y');


