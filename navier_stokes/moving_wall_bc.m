function f = moving_wall_bc(f, side, u_lb)
% D2Q9
% Applies a velocity to enforce a moving wall (as in lid-driven cavity).

if strcmp(side, 'north') % North boundary (moving lid).
    rho_end = f(end,2:end-1,1) + f(end,2:end-1,2) + f(end,2:end-1,4) + ...
        2*( f(end,2:end-1,3) + f(end,2:end-1,7) + f(end,2:end-1,6) );
    f(end,2:end-1,5) = f(end,2:end-1,3);
    f(end,2:end-1,9) = f(end,2:end-1,7) + (u_lb / 6)*rho_end;
    f(end,2:end-1,8) = f(end,2:end-1,6) - (u_lb / 6)*rho_end;
end 
    
    