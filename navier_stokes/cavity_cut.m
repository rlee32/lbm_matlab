% Lid-driven cavity with a cut corner.
% A Lattice Boltzmann D2Q9 solver.
% This features a non-lattice-aligned wall! 
% Cell centers (nodes) are placed on the boundaries. 
% Author: Robert Lee
% Email: rlee32@gatech.edu

clear;close all;clc;
addpath overlap
addpath freewall
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
cut_start_y = 0.1; % non-dimensional y-position on the west boundary.
cut_end_x = 0.1; % non-dimensional x-position on the south boundary.
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
disp(['Lattice speed: ' num2str(u_lb)]);

% Determine which lattice vectors are relevant to the cut.
parallel = [-cut_end_x, cut_start_y];
cut_length = norm(parallel);
unit_parallel = parallel / cut_length;
unit_normal = [-parallel(1), parallel(2)] / cut_length;
pgram_height = cut_length * dt;
disp(['Ratio of pgram height to dh: ' num2str(pgram_height/dh)]);
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
ci = find(valid==1); % the components relevant to wall.
c_wall = c( valid == 1 , : );
c_lengths = sqrt(sum(c_wall.^2,2));
u_wall = diag(1./c_lengths)*c_wall;
% Let's determine the pgrams.
p0 = [cut_end_x, 0];
v1 = parallel;
v2 = -u_wall * pgram_height; % a v2 for every eligible lattice link.
    
plot([p0(1), p0(1)+v1(1)], [p0(2), p0(2)+v1(2)]);  
hold on;      
plot([p0(1)+v1(1),  p0(1)+v1(1)+v2(1)], [p0(2)+v1(2), p0(2)+v1(2)+v2(2)]);
plot([p0(1), p0(1)+v2(1)], [p0(2), p0(2)+v2(2)]);
plot([p0(1)+v2(1), p0(1)+v2(1)+v1(1)], [p0(2)+v2(2), p0(2)+v2(2)+v1(2)]);

[pgram_links, ~] = size(v2);
areas = zeros(nodes,nodes,pgram_links);
cmin = linspace(0,1,nodes) - dh/2;
j_f = ceil((cut_start_y + pgram_height)/dh + 1);
non_zero_areas = 0;
for j =  1:j_f
    x_line = cut_end_x - (j-1)*(cut_start_y/cut_end_x*dh);
%     disp(num2str(x_line));
    margin = 2;
    i_0 = floor((x_line - pgram_height)/dh - margin);
    i_0 = max([1,i_0]);
    i_f = ceil((x_line + pgram_height)/dh + margin);
    i_f = min([nodes,i_f]);
    for i =  i_0:i_f
        for k = 1:pgram_links
%             disp(['i,j ' num2str(i) ', ' num2str(j)])
            cmin_ = [cmin(i), cmin(j)]';
            areas(j,i,k) = overlap_pgram_cell(p0', v1', v2(k,:)', ...
                cmin_, dh);
            if areas(j,i,k)
                plot([cmin_(1), cmin_(1)+dh],[cmin_(2), cmin_(2)]);
                plot([cmin_(1)+dh, cmin_(1)],[cmin_(2)+dh, cmin_(2)+dh]);
                plot([cmin_(1)+dh, cmin_(1)+dh],[cmin_(2)+dh, cmin_(2)]);
                plot([cmin_(1), cmin_(1)],[cmin_(2), cmin_(2)+dh]);
                areas(j,i,k) = overlap_pgram_cell(p0', v1', v2(k,:)', ...
                    cmin_, dh);
                non_zero_areas = non_zero_areas + 1;
            end
        end
    end
end
% Normalize the areas to the overall area in the pgram.
for k = 1:pgram_links
    areas(:,:,k) = areas(:,:,k) / sum(sum(areas(:,:,k)));
end

% % VISUALIZATION
% % Modified from Jonas Latt's cavity code on the Palabos website.
% uu = areas(:,:,3) / max(max(areas(:,:,3)));
% imagesc(uu);
% colorbar
% axis equal off; drawnow

% Initialize.
f = ones(nodes,nodes,9);
% Apply meso BCs.
f = freewall(f,areas,ci);
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



