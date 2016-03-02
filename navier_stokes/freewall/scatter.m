function f = scatter(f, gamma, weights, ci)
% distributes particles to f according to gamma_in.

bi = bounceback_components(ci);
for s = 1:size(weights,1)
  for k = 1:length(bi)
    surfel = weights{s,k};
    touched_cells = size(surfel,1);
    for t = 1:touched_cells
        i = round(surfel(t,1));
        j = round(surfel(t,2));
        w = surfel(t,3);
        f(i,j,bi(k)) = f(i,j,bi(k)) + w*gamma(s,k);
    end
  end
end