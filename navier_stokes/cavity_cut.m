% Lid-driven cavity with a cut corner.
% A Lattice Boltzmann D2Q9 solver.
% This features a non-lattice-aligned wall! 
% Cell centers (nodes) are placed on the boundaries. 
% Author: Robert Lee
% Email: rlee32@gatech.edu

clear;close all;clc;

% Algorithm steps:
% Initialize meso (f)
% Apply meso BCs
% Determine macro variables and apply macro BCs
% Loop:
%   Collide
%   Apply meso BCs
%   Stream
%   Apply meso BCs?
%   Determine macro variables and apply macro BCs

% Physical parameters.
L_p = 0.6;%1.1; % Cavity dimension. 
U_p = 6;%1.1; % Cavity lid velocity.
nu_p = 1.2e-3;%1.586e-5; % Physical kinematic viscosity.
rho0 = 1;
cut_start_y = 0.5; % non-dimensional y-position on the west boundary.
cut_end_x = 0.5; % non-dimensional x-position on the south boundary.
% Discrete/numerical parameters.
nodes = 100;
dt = .002;
timesteps = 10000;

% Derived nondimensional parameters.
Re = L_p * U_p / nu_p;
disp(['Reynolds number: ' num2str(Re)]);
% Derived physical parameters.
t_p = L_p / U_p;
disp(['Physical time scale: ' num2str(t_p) ' s']);
% Derived discrete parameters.
dh = 1/(nodes-1);
nu_lb = dt / dh^2 / Re;
disp(['Lattice viscosity: ' num2str(nu_lb)]);
tau = 3*nu_lb + 0.5;
disp(['Relaxation time: ' num2str(tau)]);
omega = 1 / tau;
disp(['Relaxation parameter: ' num2str(omega)]);
u_lb = dt / dh;
disp(['Lattice speed: ' num2str(u_lb)])

% Determine which lattice vectors are relevant to the cut.
parallel = [-cut_end_x, cut_start_y];
cut_length = norm(parallel);
unit_parallel = parallel / cut_length;
unit_normal = [-parallel(1), parallel(2)] / cut_length;
pgram_height = cut_length * dt;
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
valid = zeros(9,1);
for k = 1:9
    valid(k) = dot(unit_normal,c(k,:)) < 0;
end
c_wall = zeros(sum(valid),2);
counter = 1;
for k = 1:9
    if valid(k)
        c_wall(counter, :) = c(k, :);
        counter = counter + 1;
    end
end
% Pgram defined by pgram_height, cut_length, unit_normal, unit_parallel.
touched = zeros(nodes,nodes,1);
coord_min = dh*(cumsum(ones(nodes,1))-1) - dh/2;
pgram_xmax = cut_end_x + unit_normal(1)*pgram_height;
pgram_ymax = cut_start_y + unit_normal(2)*pgram_height;
max_y_node = ceil((pgram_ymax + dh) / dh + 1);
max_x_node = ceil((pgram_xmax + dh) / dh + 1);
for j = 1:min([max_y_node, nodes])
    for i = 1:min([max_x_node, nodes])
        % Overlap detection between 2 bodies at a time.
        % We detect overlap by projecting the 2 bodies to all of the
        % surface normals (one normal at a time), and check these 1D
        % criteria (all 4 have to meet in order for there to be overlap).
        cell_min = [coord_min(i)-cut_end_x, coord_min(j)];
        cell_max = cell_min+dh;
        % pgram projections
        pgram_projection = [0,pgram_height];
        cell_projection = [dot(unit_normal,cell_min), dot(unit_normal,cell_max)];
        if ( cell_projection(2) <= pgram_projection(1) ) ...
                || ( cell_projection(1) >= pgram_projection )
            continue
        end
        pgram_projection = [0,cut_length];
        cell_projection = [dot(unit_parallel,cell_min), dot(unit_parallel,cell_max)];
        if ( cell_projection(2) <= pgram_projection(1) ) ...
                || ( cell_projection(1) >= pgram_projection )
            continue
        end
        % cell projections
        pgram_projection = [0,pgram_height];
        cell_projection = [dot(unit_normal,cell_min), dot(unit_normal,cell_max)];
        if ( cell_projection(2) <= pgram_projection(1) ) ...
                || ( cell_projection(1) >= pgram_projection )
            continue
        end
        pgram_projection = [0,pgram_height];
        cell_projection = [dot(unit_normal,cell_min), dot(unit_normal,cell_max)];
        if ( cell_projection(2) <= pgram_projection(1) ) ...
                || ( cell_projection(1) >= pgram_projection )
            continue
        end
        % if it makes it here, the cell passes all checks.
        touched(j,i) = 1;
    end
end

% Initialize.
f = ones(nodes,nodes,9);
% Apply meso BCs.
f = moving_wall_bc(f,'north',u_lb);
f = wall_bc(f,'south');
f = wall_bc(f,'east');
f = wall_bc(f,'west');
% Determine macro variables and apply macro BCs
[u,v,rho] = reconstruct_macro_all(f);
u(end,2:end-1) = u_lb;
v(end,2:end-1) = 0;
u(1,:) = 0;
v(1,:) = 0;
u(:,1) = 0;
v(:,1) = 0;
u(:,end) = 0;
v(:,end) = 0;

% Main loop.
disp(['Running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps
    if (mod(iter,timesteps/10)==0)
        disp(['Ran ' num2str(iter) ' iterations']);
    end
    
    % Collision.
    f = collide(f, u, v, rho, omega);
    
    % Apply meso BCs.
    f = moving_wall_bc(f,'north',u_lb);
    f = wall_bc(f,'south');
    f = wall_bc(f,'east');
    f = wall_bc(f,'west');

    % Streaming.
    f = stream(f);
    
    % Apply meso BCs.
    f = moving_wall_bc(f,'north',u_lb);
    f = wall_bc(f,'south');
    f = wall_bc(f,'east');
    f = wall_bc(f,'west');
    
    % Determine macro variables and apply macro BCs
    [u,v,rho] = reconstruct_macro_all(f);
    u(end,2:end-1) = u_lb;
    v(end,2:end-1) = 0;
    u(1,:) = 0;
    v(1,:) = 0;
    u(:,1) = 0;
    v(:,1) = 0;
    u(:,end) = 0;
    v(:,end) = 0;
    
    % VISUALIZATION
    % Modified from Jonas Latt's cavity code on the Palabos website.
    if (mod(iter,10)==0)
        uu = sqrt(u.^2+v.^2) / u_lb;
        imagesc(flipud(uu));
        colorbar
        axis equal off; drawnow
    end
end
disp('Done!');



