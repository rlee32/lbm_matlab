function [tc, lasts] = touched_cells(segment,dh)
% gets the indices of the cells touched by segment.
% We assume the segment starts on bottom and ends on left wall.

ray = segment(2,:) - segment(1,:);
nodes = round( 1 / dh );
tc = zeros(nodes*4,2);
lasts = tc(:,1);
[i_start, ~] = cell_point(segment(1,:), dh);
% we assume we have a valid start point.
counter = 1;
finished = 0;
last_j = 1;
for j = 1:nodes
    checking = 1;
    i = i_start;
    while checking
        if i_start < 1
            break
        end
        xmin = -dh/2 + (i-1)*dh;
        ymin = -dh/2 + (j-1)*dh;
        if intersect_segment_cell(segment(1,:), ray, [xmin, ymin], dh)
            if i < 1
                finished = 1;
                break
            end
            tc(counter,:) = [i,j];
            counter = counter + 1;
        elseif i == i_start
            i_start = i_start - 1;
        else
            checking = false;
        end
        i = i-1;
    end
    if tc(counter-1,1)-1 > 0
        lasts(j) = tc(counter-1,1)-1;
        last_j = j;
    end
    if finished
        break
    end
end
tc = tc(1:counter-1,:);
lasts = lasts(1:last_j);
