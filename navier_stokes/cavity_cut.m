% UNDER CONSTRUCTION

% Lid-driven cavity with a cut corner.
% A Lattice Boltzmann D2Q9 solver.
% This features a non-lattice-aligned wall! 
% Cell centers (nodes) are placed on the boundaries. 
% Author: Robert Lee
% Email: rlee32@gatech.edu

clear;close all;clc;
addpath overlap
addpath freewall
addpath basic
addpath bc

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
surfels = 10; % 'sruface elements', the number of cut surface elements attributed to the cut.

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
[c_wall, ci] = relevant_lattice_speeds(unit_normal);
[considered, ~] = size(c_wall);
% Let's determine the pgrams.
p0 = (0:surfels-1)' * cut_length / surfels * unit_parallel;
p0(:,1) = p0(:,1) + cut_end_x;
v1 = parallel / surfels;
v2 = -c_wall .* repmat(unit_normal,length(c_wall),1) * dt; % a v2 for every eligible lattice link.

figure;
hold on;
plot_lattice_lines(nodes);
for lv = 1:considered
    for k = 1:5:surfels
        plot_surfel(p0(k,:), v1, v2(lv,:));
        [bmin, bmax, imin, imax] = pgram_bounds(p0(k,:), v1, v2(lv,:),dh);
        plot_bounding_box(bmin,bmax);
        weights = pgram_weights(p0(k,:), v1, v2(lv,:),dh);
    end
end

areas = zeros(nodes,nodes,considered);
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
        for k = 1:considered
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
for k = 1:considered
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
% f = freewall(f,areas,ci);
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
        imagesc((uu));
        colorbar
        axis equal off; drawnow
    end
end
disp('Done!');



