classdef pgram < handle
    properties
        k                       % the lattice velocity index that pgram corresponds to.
        area                    % area of this pgram.
        touched_cells           % handles of touched_cell objects
        overlap_areas           % the area of surfel-cell overlap corresponding to touched_cells.
        collected_particles     % collected particles in this pgram.
    end
    methods
        function [fd, cell_indices] = distribute_particles(obj)
            for h = obj.touched_cells
                h.
            end
            fd = obj.overlap_areas / obj.area .* obj.collected_particles;
            cell_indices = obj.touched_cells;
        end
    end
end