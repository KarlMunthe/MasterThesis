function xdiff = xdiff(u, hx, nu, yside)

    %calculate diffusion in x direction
    
    nux = 0.5*(nu(1:end-1, 2:end) + nu(1:end-1, 1:end-1));
    ux = (u(1:end-1,  2:end) - u(1:end-1, 1:end-1))./hx;
    
    xdiff = yside(1:end-1).*(nux(:, 2:end).*ux(:, 2:end) - nux(:, 1:end-1).*ux(:, 1:end-1));
    
end