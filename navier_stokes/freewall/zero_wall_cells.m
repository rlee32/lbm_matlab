function f = zero_wall_cells(f, wall_cells, ci)
% Meant to be called before scattering.
% Zeros the relevant distribution function components according to the
%   surfel-incoming components ci.
% Only those cells who touch the wall (and do not get streamed certain
%   components, and therefore have unknown components) are zero'd.

bi = bounceback_components(ci);
for k = 1:size(wall_cells,1)
    f(wall_cells(k,1),wall_cells(k,2),bi) = 0;
end