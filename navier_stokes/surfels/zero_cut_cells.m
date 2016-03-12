function f = zero_cut_cells(cut_cells, bouncedback_indices, f)
% Zeros the distribution functions for the wall-opposite distribution
%   functions, so that wrong distributions are not received from the inner
%   cells of the wall (the non-fluid cells).
% cut_cells: vector of handles to touched_cell objects.
% bouncedback_indices: the indices of the distribution functions that are
%   point into the fluid domain.

% the unknown distributions are those with values >0 when dotted with 
%   surface normal. These should be set to zero after streaming and after
%   gathering, but before scattering!

% steps must take place in this order:
%   gather
%   stream
%   zero
%   scatter

for k = bouncedback_indices
    for h = cut_cells
        f(h.j, h.i, k) = 0;
    end
end