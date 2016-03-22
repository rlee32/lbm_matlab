function tc = fill_touched_cell_info(tc, ss, nodes)
% from the information in ss, fill in missing info in touched cells (tc).

overlap_areas = zeros(nodes,nodes,9);

for s = ss
    for p = s.pgrams
        for k = 1:length(p.celli)
            i = p.celli(k);
            j = p.cellj(k);
            a = p.overlap_areas(k);
            li = p.lattice_index;
            overlap_areas(j, i, li) = overlap_areas(j, i, li) + a;
        end
    end
end

for tck = 1:length(tc)
    bb = tc(tck).lattice_indices;
    for k = 1:length(tc(tck).lattice_indices)
        bb(k) = bounceback_components(tc(tck).lattice_indices(k));
    end
%     oa = overlap_areas(c.j,c.i,bb);
%     oa = overlap_areas(c.j,c.i,c.lattice_indices);
%     c.overlap_areas = overlap_areas(c.j,c.i,c.lattice_indices);
    for k = bb
        tc(tck).overlap_areas = [tc(tck).overlap_areas, ...
            overlap_areas(tc(tck).j,tc(tck).i,k)]; 
    end
end