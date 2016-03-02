function [f,rho,u,v] = zero_out_of_bounds(f,rho,u,v,lasts)
% lasts: lasts(j) is the last i position that should be excluded from the 
%   simulation (up to i).

for k = 1:length(lasts)
    f(1:lasts(k), k, :) = 0;
    rho(1:lasts(k), k) = 0;
    u(1:lasts(k), k) = 0;
    v(1:lasts(k), k) = 0;
end