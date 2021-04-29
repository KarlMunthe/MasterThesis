function [W, CFL, error, h] = CompressibleNSSTest(m, t_end, x_end, k, c_p, c_V, rho_sol, mom_sol, E_sol, R1, R2, R3)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);

[~, ~, ~, Q1] = PeriodicD0(m, h);

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
    
    new_vars = RK4(vars, c_p, c_V, R, k, x, t(i));
    
    RHO(:, i+1) = new_vars(:,1);
    MOM(:, i+1) = new_vars(:,2);
    E(:, i+1) = new_vars(:, 3);

end

RHO(end+1, :) = RHO(1, :);

W = RHO;

error = RHO(:,end) - rho_sol(x, t_end)';
error = sqrt(h*error'*error);

function [w] = RK4(vars, c_p, c_V, R, k, x, t)
    
    k1 = NSS(vars, c_p, c_V, R, x, t);
    k2 = NSS(vars + 0.5*k*k1, c_p, c_V, R, x, t + 0.5*k);
    k3 = NSS(vars + 0.5*k*k2, c_p, c_V, R, x, t + 0.5*k);
    k4 = NSS(vars + k*k3, c_p, c_V, R, x, t + k);
    
    w = vars + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = NSS(vars, c_p, c_V, R, x, t)
    
    gamma = c_p/c_V;
    rho = vars(:, 1);
    mom = vars(:, 2);
    %p = rho.^gamma;
    E = vars(:, 3);
    
    nu = 1./rho;
    p = R/c_V*(E - mom.^2./(2*rho));
    
    A1 = -Q1*mom;
    
    B1 = Q1*(nu.*(Q1*rho));
    
    C1 = (R1(x(1:end-1), t))';
    
    A2 = -Q1*((mom.^2)./rho + p);
    
    B2 = Q1*(nu.*(Q1*mom));

    C2 = (R2(x(1:end-1), t))';
    
    A3 = -Q1*(E.*(mom./rho) + p.*(mom./rho));
    
    B3 = Q1*(nu.*(Q1*E));
    
    C3 = (R3(x(1:end-1), t))';

    flux = [A1 + B1 + C1, A2 + B2 + C2, A3 + B3 + C3];

end

end
