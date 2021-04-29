function [W, CFL, error, h] = PeriodicAdvectionTestWithPade(m, k, sol, t_end)

%h = 1/(m-1);       %space step

t = k:k:t_end;

x = linspace(0, 2*pi, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/h;

%initial condition
f = sol(x(1:end-1), 0)';

%periodic pade
%[LHS1, RHS1, LHS2, RHS2] = PeriodicPade(m,h);


[Q1, ~, ~, ~] = PeriodicD0(m, h);
[Q2, ~, ~, ~] = PeriodicD2(m, h);

%REMEMBER TO USE SECOND ORDER, O(2), DIFFERENCE OPERATORS WHEN USING PADE
%SCHEMES, THE FIRST SIMPLESt OF THE ONES IN THE PERIODICD0/2 FUNCTIONS
PD1 = (eye(m-1) + h^2*Q2/6);
PD2 = (eye(m-1) + h^2*Q2/12);

w = [f, zeros(length(f), length(t))];

te = k;

for i = 1:length(t)

    w(:,i+1) = RK4(w(:,i), k);
    
    te = te + k;
    
end

w(end + 1, :) = w(1, :);

error = w(:,end)-sol(x, t_end)';
error = sqrt(h*error'*error);

W = w;

%
%JUST FUNCTIONS BELOW
%

function [w] = RK4(u, k)
    
    k1 = Advection(u);
    k2 = Advection(u + 0.5*k*k1);
    k3 = Advection(u + 0.5*k*k2);
    k4 = Advection(u + k*k3);
    
    w = u + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = Advection(u)
    
    A1 = (Q1*(-u));
    
    %C1 = G(x(1:end-1)', t);

    flux = A1;

end
end