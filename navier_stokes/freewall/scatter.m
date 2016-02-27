function f = scatter(f, gamma_in, areas, ci)
% distributes particles to f according to gamma_in.

bounceback = zeros(9,1);
bounceback(1) = 1;
bounceback(2) = 4;
bounceback(3) = 5;
bounceback(4) = 2;
bounceback(5) = 3;
bounceback(6) = 8;
bounceback(7) = 9;
bounceback(8) = 6;
bounceback(9) = 7;

for k=1:length(ci)
    f(:,:,bounceback(ci(k))) = f(:,:,bounceback(ci(k))) ...
        + areas(:,:,k)*gamma_in(k);
end
    