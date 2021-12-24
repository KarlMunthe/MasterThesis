function yflux = yflux(f, dx)
    
    %calculate flux in y direction
    
    yflux = 0.5*(f(1:end-1, 2:end-1) + f(2:end, 2:end-1));
    yflux = [zeros(1, size(f(:, 2:end-1), 2)); yflux];  %boundary condition
    
    yflux = -(yflux(2:end, :) - yflux(1:end-1, :)).*dx(2:end-1);
    
end