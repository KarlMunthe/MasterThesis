function xflux = xflux(f, dy)

    %Calculate flux in x direction
    
    xflux = 0.5*(f(1:end-1, 1:end-1) + f(1:end-1, 2:end));
    
    xflux = -(xflux(:, 2:end) - xflux(:, 1:end-1)).*dy(1:end-1);

end