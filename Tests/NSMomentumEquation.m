function [W, CFL, error, h] = NSMomentumEquation(m, t_end, x_end, k, p, rho, D1,  sol)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);

[LHS1, RHS1, LHS2, RHS2] = PeriodicPade(m,h);

[Q1, ~, ~, ~] = PeriodicD0(m, h);
[Q2, ~, ~, ~] = PeriodicD2(m, h);

%REMEMBER TO USE SECOND ORDER, O(2), DIFFERENCE OPERATORS WHEN USING PADE
%SCHEMES, THE FIRST SIMPLESt OF THE ONES IN THE PERIODICD0/2 FUNCTIONS
PD1 = (eye(m-1) + h^2*Q2/6);
PD2 = (eye(m-1) + h^2*Q2/12);

%Initial conditions and Primitive variables
mom = sol(x(1:end-1), 0)';

MOM = [mom, zeros(length(mom), length(t)-1)];

for i = 1:length(t)-1

    MOM(:, i+1) = RK4(MOM(:,i), rho, p, k, x, t(i));

end

MOM(end+1, :) = MOM(1, :);

W = MOM;

error = MOM(:,end) - sol(x, t_end)';
error = sqrt(h*error'*error);

function [w] = RK4(mom, rho, p, k, x, t)
    
    k1 = PVA(mom, rho, p, x, t);
    k2 = PVA(mom + 0.5*k*k1, rho, p, x, t + 0.5*k);
    k3 = PVA(mom + 0.5*k*k2, rho, p, x, t + 0.5*k);
    k4 = PVA(mom + k*k3, rho, p, x, t + k);
    
    w = mom + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = PVA(mom, rho, p, x, t)
    
    A1 = LHS1\(-RHS1*(mom.^2./(rho(x(1:end-1), t))' + (p(x(1:end-1), t))'));
    
    B1 = LHS2\(4/3*RHS2*(mom./(rho(x(1:end-1), t))'));

    C1 = (D1(x(1:end-1), t))';

    flux = A1 + B1 + C1;

end

end
