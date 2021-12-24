function ydiff = ydiff(u, hy, nu, xside)
    
    %calculate diffusion in y direction
    
    nuy = 0.5*(nu(2:end, 2:end-1) + nu(1:end-1, 2:end-1));
    nuy = [zeros(1, size(nuy, 2)); nuy];
    uy = (u(2:end, 2:end-1) - u(1:end-1, 2:end-1))./hy;
    uy = [zeros(1, size(uy, 2)); uy];   %Boundary condition

    
    ydiff = xside(2:end-1).*(nuy(2:end, :).*(uy(2:end, :) - nuy(1:end-1, :).*uy(1:end-1, :)));

end