function w = RK4(vars, dt, v, hx, hy, xside, yside, gamma, mu, svars)
       
    k1 = RHS(vars, v, hx, hy, xside, yside, gamma, mu);
    
    vars(1, :, 2:3) = svars;    %set boundary conditions on southernn border

    k2 = RHS(vars + 0.5*dt*k1, v, hx, hy, xside, yside, gamma, mu);
    
    vars(1, :, 2:3) = svars;    %set boundary conditions on southernn border
    
    k3 = RHS(vars + 0.5*dt*k2, v, hx, hy, xside, yside, gamma, mu);
    
    vars(1, :, 2:3) = svars;    %set boundary conditions on southernn border
    
    k4 = RHS(vars + dt*k3, v, hx, hy, xside, yside, gamma, mu);
    
    vars(1, :, 2:3) = svars;    %set boundary conditions on southernn border
    
    w = vars + dt./v.*(k1 + 2*k2 + 2*k3 + k4)/6;
    
    w(1, :, 2:3) = svars;       %set boundary conditions on southernn border
    
end