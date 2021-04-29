function [W, CFL, error, h] = IsentropicNSTest(m, t_end, x_end, k, rho_sol, R1, R2,  mom_sol)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);

[~, ~, ~, D0] = PeriodicD0(m, h);
[~, ~, ~, D2] = PeriodicD2(m, h);

%Initial conditions and Primitive variables
rho_0 = rho_sol(x(1:end-1), 0)';
mom_0 = mom_sol(x(1:end-1), 0)';

RHO = [rho_0, zeros(length(rho_0), length(t)-1)];
MOM = [mom_0, zeros(length(mom_0), length(t)-1)];

for i = 1:length(t)-1

    vars = [RHO(:, i), MOM(:, i)];
    
    new_vars = RK4(vars, k, x, t(i));
    
    RHO(:, i+1) = new_vars(:,1);
    MOM(:, i+1) = new_vars(:,2);
    %p = new_vars(:,3);

end

RHO(end+1, :) = RHO(1, :);

W = RHO;

error = RHO(:,end) - rho_sol(x, t_end)';
error = sqrt(h*error'*error);

function [w] = RK4(vars, k, x, t)
    
    k1 = INS(vars, x, t);
    k2 = INS(vars + 0.5*k*k1, x, t + 0.5*k);
    k3 = INS(vars + 0.5*k*k2, x, t + 0.5*k);
    k4 = INS(vars + k*k3, x, t + k);
    
    w = vars + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = INS(vars, x, t)
    
    rho = vars(:, 1);
    mom = vars(:, 2);
    
    A1 = -D0*mom;
    
    B1 = (R1(x(1:end-1), t))';
    
    A2 = -D0*(mom.^2./rho + rho.^3);
    
    B2 = 4/3*D2*(mom./rho);

    C2 = (R2(x(1:end-1), t))';

    flux = [A1 + B1, A2 + B2 + C2];

end

end
