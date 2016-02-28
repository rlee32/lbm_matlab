function [nut, nutilde] = update_nut(nutilde,nu,dt,d,u,v,ux,vy,uy,vx)
% For the Spalart-Allmaras model in lid-driven cavity.
% input variables are 2d matrices.
% d: wall distances.
% nut is the turbulent viscosity that gets added onto the nominal
%   viscosity to produce the total effective viscosity.

% Model constants.
sigma = 2/3;
cb1 = 0.1355;
cb2 = 0.622;
kappa = 0.41;
cw2 = 0.3;
cw3 = 2;
cv1 = 7.1;

% Derived model constants.
cw1 = cb1 / kappa^2 + ( 1 + cb2 ) / sigma;

% Derived input variables.
S = 2*ux.^2 + 2*vy.^2 + (uy + vx).^2 - 2/3*(ux + vy).^2;
X = nutilde / nu;
fv1 = X.^3 ./ ( X.^3 + cv1^3 );
Stilde = S.^0.5 .* ( 1 ./ X + fv1 );
r = tanh( nutilde ./ ( Stilde*kappa^2*d.^2 ) ) / tanh(1.0);
g = r + cw2 * ( r.^6 - 6 );
fw = g.*( ( 1 + cw3^6 ) ./ ( g.^6 + cw3^6 ) ).^(1/6);

nutildex = ;
nutildey = ;
dnutilde = -( u.*nuildex + v.*nutildey ) + ;

nutilde = nutilde + dt*dnutilde;


