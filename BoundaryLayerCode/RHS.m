function RHS = RHS(vars, v, hx, hy, xside, yside, gamma, mu)

    rho = vars(:, :, 1);
    m = vars(:, :, 2);
    n = vars(:, :, 3);
    E = vars(:, :, 4);
    nu = mu./rho;
    
    [F1, F2, F3, F4] = F(rho, m, n, E, gamma); 
    [G1, G2, G3, G4] = G(rho, m, n, E, gamma); 
    
    F1flux = xflux(F1, yside);      %Calculate advective terms in x direction for mass equation
    F2flux = xflux(F2, yside);      %Calculate advective terms in x direction for xmomentum equation
    F3flux = xflux(F3, yside);      %Calculate advective terms in x direction for ymomentum equation
    F4flux = xflux(F4, yside);      %Calculate advective terms in x direction for energy equation
    
    G1flux = yflux(G1, xside);      %Calculate advective terms in y direction for mass equation
    G2flux = yflux(G2, xside);      %Calculate advective terms in y direction for xmomentum equation
    G3flux = yflux(G3, xside);      %Calculate advective terms in y direction for ymomentum equation
    G4flux = yflux(G4, xside);      %Calculate advective terms in y direction for energy equation
    
    rhoxdiff = xdiff(rho, hx, nu, yside);   %Calculate diffusive terms in x direction for density equation
    rhoydiff = ydiff(rho, hy, nu, xside);   %Calculate diffusive terms in y direction for density equation
    mxdiff = xdiff(m, hx, nu, yside);       %Calculate diffusive terms in x direction for xmomentum equation
    mydiff = ydiff(m, hy, nu, xside);       %Calculate diffusive terms in y direction for xmomentum equation
    nxdiff = xdiff(n, hx, nu, yside);       %Calculate diffusive terms in x direction for ymomentum equation
    nydiff = ydiff(n, hy, nu, xside);       %Calculate diffusive terms in y direction for ymomentum equation
    Exdiff = xdiff(E, hx, nu, yside);       %Calculate diffusive terms in x direction for energy equation
    Eydiff = ydiff(E, hy, nu, xside);       %Calculate diffusive terms in x direction for energy equation

    RHS = zeros(size(v, 1), size(v,  2), 4);    
    RHS(1:end-1, 2:end-1, 1) = F1flux + G1flux + rhoxdiff + rhoydiff;   %Calculate density
    RHS(1:end-1, 2:end-1, 2) = F2flux + G2flux + mxdiff + mydiff;       %Calculate xmomentum
    RHS(1:end-1, 2:end-1, 3) = F3flux + G3flux + nxdiff + nydiff;       %Calculate ymomentum
    RHS(1:end-1, 2:end-1, 4) = F4flux + G4flux + Exdiff + Eydiff;       %Calculate energy
    
    
end   