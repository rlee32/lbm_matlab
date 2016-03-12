% UNDER CONSTRUCTION
% The dynamic behavior of the distributions and pgrams are hard to capture
%   with MATLAB script... unless I simplify the cuts considerably. Even
%   then, it is still quite a job.
% This may remain defunct (and I may instead implement this in C++ using
%   object-oriented techniques).

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
addpath surfels

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

% To handle the cut cell bc, we:
%   1. Save the wall distributions that would come from inside the wall.
%   2. Advect
%   3. Zero the inner distributions (meso and macro scopically).
%   4. load the saved wall distributions.

% Physical parameters.
L_p = 0.6;%1.1; % Cavity dimension. 
U_p = 6;%1.1; % Cavity lid velocity.
nu_p = 1.2e-3;%1.586e-5; % Physical kinematic viscosity.
rho0 = 1;
cut_start_y = 0.1; % non-dimensional y-position on the west boundary.
cut_end_x = 0.1; % non-dimensional x-position on the south boundary.
% Discrete/numerical parameters.
nodes = 100;
dt = .003;
timesteps = 10000;
surfels = 10; % 'surface elements', the number of cut surface elements attributed to the cut.

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
% Let's determine the pgrams.
p0 = (0:surfels-1)' * cut_length / surfels * unit_parallel;
p0(:,1) = p0(:,1) + cut_end_x;
v1 = parallel / surfels;
% v2 = -c_wall .* repmat(unit_normal,length(c_wall),1) * dt; % a v2 for every eligible lattice link.
v2 = -c_wall * dt; % a v2 for every eligible lattice link.


% Surfel and lattice check.
figure;
hold on;
plot_lattice_lines(nodes);
[considered, ~] = size(c_wall);
for lv = 1:considered
    for k = 1:2:surfels
        plot_surfel(p0(k,:), v1, v2(lv,:));
        [bmin, bmax, imin, imax] = pgram_bounds(p0(k,:), v1, v2(lv,:),dh);
        plot_bounding_box(bmin,bmax);
    end
end
plot([cut_end_x,0],[0,cut_start_y]); 

% get the touched cells, so we can save their prestreaming distributions.
% [tc, lasts] = find_touched_cells([cut_end_x,0;0,cut_start_y],dh, ...
%     cut_start_y,cut_end_x);
disp('Finding cut cells');
tc = find_cut_cells([cut_end_x,0;0,cut_start_y],dh, ...
    cut_start_y,cut_end_x);
disp('Finding fluid areas');
fluid_areas = fluid_areas_table(tc, nodes, dh);
% lets get the surfel objects, which contain pgrams and touched_cells.
disp('Generating surfels');
ss = generate_surfels(p0,v1,dt,dh);

% all-important weights for bc enforcement... 
% weights = surfel_weights(p0,v1,v2,dh,tc); 

% % VISUALIZATION
% % Modified from Jonas Latt's cavity code on the Palabos website.
% uu = areas(:,:,3) / max(max(areas(:,:,3)));
% imagesc(uu);
% colorbar
% axis equal off; drawnow

% Initialize macro, then meso.
rho = rho0*ones(nodes,nodes);
u = zeros(nodes,nodes);
v = zeros(nodes,nodes);
u(end,2:end-1) = u_lb;
% Set f to feq
f = compute_feq(rho,u,v);
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
% Enforce cut corner bc.
% [f,rho,u,v] = zero_out_of_bounds(f,rho,u,v,lasts);

% ad hoc bounced back indices for 45 degree cut.
bounced = [2, 3, 6];

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
    % save cut corner distributions.
%     saved = save_wall_distributions(f,tc,ci);
%     f = preload_inactive_cells(f,lasts);

    f = collect(ss,f,fluid_areas);
    f = stream(f);
    f = area_scale_distributions(f, fluid_areas);
    f = scatter(ss,f,fluid_areas);
    
    % load (the saved) cut corner distributions.
%     f = load_wall_distributions(f,saved,ci);
    % zero inactive cells.
%     [f,~,~,~] = zero_out_of_bounds(f,rho,u,v,lasts);
    
%     % OLD
%     % Now, apply volumetric boundary condition.
%     G = gather(f,weights,ci);
%     f = zero_wall_cells(f,tc,ci);
%     % Apply wall bc here to take care of the cells that touch both wall and
%     %   cut (only 2 cells total should be doing this).
%     f = wall_bc(f,'south');
%     f = wall_bc(f,'west');
%     f = scatter(f,G,weights,ci);

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
    if (mod(iter,1)==0)
        uu = sqrt(u.^2+v.^2) / u_lb;
        imagesc(flipud(uu));
%         imagesc(flipud(f(:,:,2)));
        colorbar
        axis equal off; drawnow
    end
end
disp('Done!');



