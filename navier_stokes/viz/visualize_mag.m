function visualize_mag(u, v, u_lb)
% VISUALIZATION
% Modified from Jonas Latt's cavity code on the Palabos website.

uu = sqrt(u.^2+v.^2) / u_lb;
imagesc(flipud(uu));
%         imagesc(flipud(f(:,:,2)));
colorbar
axis equal off; drawnow
