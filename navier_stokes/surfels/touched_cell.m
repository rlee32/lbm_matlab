classdef touched_cell < handle
    properties
        k                       % the lattice velocity index to which this touched_cell pertains.
        i                       % i index of this touched cell.
        j                       % j index of this touched cell.
        pgrams                  % handles of the pgrams that touch this cell.
        overlap_areas           % overlap areas with this cell, corresponding to pgrams.
        non_overlap_area        % area that is not overlapped with the k surfel.
        distributed_particles   % particles distributed to this cell from pgrams.
    end
     methods
        function fa = advected_particles(obj, f)
            fa = f(obj.j,obj.i,obj.k) * obj.non_overlap_area;
        end
        function [fc, pgram_indices] = collected_particles(obj, f, pgrams)
            fc = f(obj.j,obj.i,obj.k) * obj.overlap_areas;
            pgram_indices = pgrams;
        end
    end
end