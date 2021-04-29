function [W, CFL, error, h] = CompressibleNSTest(m, epsilon, t_end, x_end, k, c_p, c_V, kappa, rho_sol, mom_sol, E_sol, R1, R2, R3)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);

[~, Q1, ~, ~] = PeriodicD0(m, h);
[~, Q2, ~, ~] = PeriodicD2(m, h);

%Initial conditions and Primitive variables
R = c_p - c_V;
gamma = c_p/c_V;
rho_0 = rho_sol(x(1:end-1), 0)';
mom_0 = mom_sol(x(1:end-1), 0)';
E_0 = E_sol(x(1:end-1), 0)';

RHO = [rho_0, zeros(length(rho_0), length(t)-1)];
MOM = [mom_0, zeros(length(mom_0), length(t)-1)];
E = [E_0, zeros(length(E_0), length(t)-1)];

for i = 1:length(t)-1

    vars = [RHO(:, i), MOM(:, i), E(:, i)];
    
    new_vars = RK4(vars, c_V, R, kappa, k, x, t(i));
    
    RHO(:, i+1) = new_vars(:,1);
    MOM(:, i+1) = new_vars(:,2);
    E(:, i+1) = new_vars(:, 3);

end

RHO(end+1, :) = RHO(1, :);

W = RHO;

error = RHO(:,end) - rho_sol(x, t_end)';
error = sqrt(h*error'*error);

function [w] = RK4(vars, c_V, R, kappa, k, x, t)
    
    k1 = NS(vars, c_V, R, kappa, x, t);
    k2 = NS(vars + 0.5*k*k1, c_V, R, kappa, x, t + 0.5*k);
    k3 = NS(vars + 0.5*k*k2, c_V, R, kappa, x, t + 0.5*k);
    k4 = NS(vars + k*k3, c_V, R, kappa, x, t + k);
    
    w = vars + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = NS(vars, c_V, R, kappa, x, t)
    
    rho = vars(:, 1);
    mom = vars(:, 2);
    E = vars(:, 3);
    p = R/c_V*(E - mom.^2./(2*rho));
    T = p/R;
    
    A1 = -Q1*mom;
    
    B1 = (R1(x(1:end-1), t))';
    
    A2 = -Q1*(mom.^2./rho + p);
    
    B2 = epsilon*4/3*Q2*(mom./rho);

    C2 = (R2(x(1:end-1), t))';
    
    A3 = -Q1*(E.*(mom./rho) + p.*(mom./rho));
    
    B3 = epsilon*4/3*1/2*Q2*((mom./rho).^2) + Q1*(kappa*Q1*T);
    
    C3 = (R3(x(1:end-1), t))';

    flux = [A1 + B1, A2 + B2 + C2, A3 + B3 + C3];

end

end
