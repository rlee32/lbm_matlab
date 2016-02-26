function f = collide_mrt(f, u, v, rho, nu_lb)
% D2Q9 collisions on 2-d matrix.
% Multiple relaxation time formulation.

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

M = [ones(1,9);...
    -4, -ones(1,4), 2*ones(1,4);
    4, -2*ones(1,4), ones(1,4);
    0, 1, 0, -1, 0, 1, -1, -1, 1;
    0, -2, 0, 2, 0, 1, -1, -1, 1;
    zeros(1,4), -1, 1, 1, -1, -1;
    0, 0, -2, 0, 2, 1, 1, -1, -1;
    0, 1, -1, 1, -1, zeros(1,4);
    zeros(1,5), 1, -1, 1 -1];
Minv = 1/36*[4*ones(9,1), M(2,:)', M(3,:)', 6*M(4,:)', 3*M(5,:)',...
    6*M(6,:)', 3*M(7,:)', 9*M(8,:)', 9*M(9,:)'];
nu_term = 2 / ( 1 + 6*nu_lb );
S_vec = [1, 1.2, 1, 1, 1.2, 1, 1.2, nu_term, nu_term]';

% Collide
t1 = u.*u + v.*v;
ru = rho.*u;
rv = rho.*v;
ru2rv2 = ru.^2 + rv.^2;
[rows, cols] = size(rho);
meq = zeros(rows, cols, 9);
meq(:,:,1) = rho;
meq(:,:,2) = -2*rho + 3*ru2rv2;
meq(:,:,3) = rho - 3*ru2rv2;
meq(:,:,4) = ru;
meq(:,:,5) = -ru;
meq(:,:,6) = rv;
meq(:,:,7) = -rv;
meq(:,:,8) = ru.^2-rv.^2;
meq(:,:,9) = ru.*rv;
for j = 1:rows
    for i = 1:cols
        m_vec = M*reshape(f(j,i,:),9,1,1);
        f(j,i,:) = f(j,i,:) ...
            - reshape( ...
                Minv*( ...
                    S_vec.*(m_vec - reshape(meq(j,i,:),9,1,1)) ...
                ),1,1,9 ...
            );
    end
end


