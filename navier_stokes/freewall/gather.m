function gamma = gather(f,areas,ci)
% computes the agglomerated distribution function for the volumetric 
%   boundary condition.
% Assumes one surface described by 'areas' (the relative portion of area of
%   every cell in the pgram).
% Assumes uniform grid; dh denotes the cell dimension.
% f is the distribution function for every cell corresponding to areas.
% ci maps areas to f, by the proper component.
% gamma is the agglomeration, one entry for every z-dimensional entry in
%   areas.

gamma = zeros(length(ci),1);
for k = 1:length(ci)
    gamma(k) = sum(sum(f(:,:,ci(k)).*areas(:,:,k)));
end
