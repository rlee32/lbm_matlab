function f = collect(ss, f, fluid_areas)
% Fills each pgram in surfel with particle collections. 

for s = ss
    for p = s.pgrams
        f = collect(p, f, fluid_areas);
    end
end