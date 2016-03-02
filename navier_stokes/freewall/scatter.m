function f = scatter(f, gamma, weights, ci)
% distributes particles to f according to gamma_in.

bi = bounceback_components(ci);
for s = 1:size(gamma,1)
    for k=1:length(ci)
        f(:,:,bi(k)) = f(:,:,bi(k)) ...
            + weights(:,k).*gamma(:,k);
    end
end
