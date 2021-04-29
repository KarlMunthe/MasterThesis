function [W, CFL, error, h] = PeriodicVariableAdvectionTestWithPade(m, t_end, x_end, k, u, R1,  sol)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/h;

%periodic pade
[LHS1, RHS1, LHS2, RHS2] = PeriodicPade(m,h);
%{
[Q1, ~, ~, ~] = PeriodicD0(m, h);
[Q2, ~, ~, ~] = PeriodicD2(m, h);

%REMEMBER TO USE SECOND ORDER, O(2), DIFFERENCE OPERATORS WHEN USING PADE
%SCHEMES, THE FIRST SIMPLESt OF THE ONES IN THE PERIODICD0/2 FUNCTIONS
PD1 = (eye(m-1) + h^2*Q2/6);
PD2 = (eye(m-1) + h^2*Q2/12);
%}

%Initial conditions and Primitive variables
rho = (sol(x(1:end-1), 0))';

RHO = [rho, zeros(length(rho), length(t)-1)];

for i = 1:length(t)-1

    RHO(:, i+1) = RK4(RHO(:,i), k, t(i));

end

RHO(end+1, :) = RHO(1, :);

W = RHO;

error = RHO(:,end) - sol(x, t_end)';
error = sqrt(h*error'*error);

function [w] = RK4(rho, k, t)
    
    k1 = PVA(rho, t);
    k2 = PVA(rho + 0.5*k*k1, t + 0.5*k);
    k3 = PVA(rho + 0.5*k*k2, t + 0.5*k);
    k4 = PVA(rho + k*k3, t + k);
    
    w = rho + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = PVA(rho, t)
    
    A1 = LHS1\(-RHS1*(u(x(1:end-1)', t).*rho));

    B1 = (R1(x(1:end-1), t))';

    flux = A1 + B1;

end

end
