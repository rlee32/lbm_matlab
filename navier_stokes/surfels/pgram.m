classdef pgram < handle
    properties
        p0                      % one corner of the pgram.
        surface                 % the surface vector associated with this pgram.
        extrusion               % the vector describing the extrusion of 
                                % the surface into the fluid domain.
%         lattice_index           % the lattice velocity index that pgram corresponds to.
        area                    % area of this pgram. Equals magnitude of 
                                % cross(surface, extrusion).
        celli                   % cell i indices touched by this pgram.
        cellj                   % cell j indices touched by this pgram.
%         touched_cells           % handles of touched_cell objects
        overlap_areas           % the area of surfel-cell overlap corresponding to touched_cells.
        weights                 % the normalized overlap_areas
        collected_particles     % collected particles in this pgram.
    end
    methods
        function obj = pgram(p0, surface, extrusion, dh)
            obj.p0 = p0;
            obj.surface = surface;
            obj.extrusion = extrusion;
            
            compute_area(obj, surface, extrusion);
            determine_weights(obj,dh);
            
        end
%         function [fd, cell_indices] = distribute_particles(obj)
%             for k = 1:length(obj.touched_cells)
%                 obj.touched_cells(k).distributed_particles = ...
%                     obj.touched_cells(k).distributed_particles + ...
%                     obj.overlap_areas(k) * obj.collected_particles;
%             end
%             fd = obj.overlap_areas / obj.area .* obj.collected_particles;
%             cell_indices = obj.touched_cells;
%         end
    end
    methods ( Access = private )
        function z = cross2d(v0,v1)
            v3 = cross([v0;0],[v1;0]);
            z = v3(3);
        end
        function compute_area(obj,v0,v1)
            obj.area = abs( cross2d(v0,v1) );
        end
%         function determine_cell_relationships_lattice(obj)
%             obj.touched_cells = [];
%             obj.overlap_areas = [];
%         end
        function [bmin, bmax, imin, imax] = pgram_bounds(p0,v1,v2,dh)
            % returns the bottom left and top right coordinates of the 
            % cell points that make a (rectangular) box that enclose a pgram.
            % p0,v1,v2: describes the pgram.
            % dh: lattice cell dimension.
            % imin: the indices of the cell whose bottom left corner is bmin.
            % imax: the indices of the cell whose top right corner is bmax.

            coordinate_shift = [-dh/2, -dh/2];
            % coordinate_shift = [0,0];

            % First construct all vertices in the pgram.
            p1 = p0 + v1;
            p2 = p0 + v2;
            p3 = p0 + v1 + v2;

            points = [p0;p1;p2;p3];
            mins = min(points,[],1);
            maxs = max(points,[],1);
            bmin = floor((mins-coordinate_shift)/dh)*dh + coordinate_shift;
            bmax = ceil((maxs-coordinate_shift)/dh)*dh + coordinate_shift;

            imin = round((bmin - coordinate_shift) / dh + 1);
            imax = round((bmax - coordinate_shift) / dh);
        end
        function determine_weights(obj,dh)
            % Returns a list of row vectors [i,j,area]
            % the input pgram has contact area 'area' with cell i,j.
            [~, ~, imin, imax] = ...
                pgram_bounds(obj.p0, obj.surface, obj.extrusion, dh);
            x0 = -dh/2;
            y0 = -dh/2;
%             obj.touched_cells = [];
            obj.overlap_areas = [];
            obj.celli = [];
            obj.cellj = [];
            obj.weights = [];
            for i = imin(1):imax(1)
                for j = imin(2):imax(2)
                    minpos = [x0 + dh*(i-1), y0 + dh*(j-1)];
                    overlap_area = overlap_pgram_cell( ...
                        obj.p0, obj.surface, obj.extrusion, minpos, dh);
                    if overlap_area
                        obj.overlap_areas = [obj.overlap_areas, overlap_area];
                        obj.celli = [obj.celli, i];
                        obj.cellj = [obj.cellj, j];
                        obj.weights = [obj.weights, ...
                            overlap_area / obj.area];
%                         obj.touched_cells = [obj.touched_cells, ...
%                             touched_cell( i, j, 0 )];
                    end
                end
            end
        end
    end
end