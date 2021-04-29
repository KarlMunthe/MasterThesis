function [W, CFL, error, h] = PeriodicHeatTest(x_end, m, k, sol, t_end, R1)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

h = 2*pi/m;    %individual space points on the grid
x = h*(1:m);

CFL = k/h;

%[LHS1, RHS1, LHS2, RHS2] = PeriodicPade(m,h);

%[~, ~, ~, Q2] = PeriodicD2(m, h);

Q2 = SpectralD2(m, h);

%initial condition
f = sol(x, 0)';


w = [f, zeros(length(f), length(t)-1)];

for i = 1:length(t)-1

    w(:,i+1) = RK4(w(:,i), k, t(i));
    
end

w = [w(end, :) ; w];

ex = [0, x];

error = w(:,end)-sol(ex, t_end)';
error = sqrt(h*error'*error);

W = w;

%
%JUST FUNCTIONS BELOW
%

function [w] = RK4(u, k, t)
    
    k1 = Advection(u, t);
    k2 = Advection(u + 0.5*k*k1, t + 0.5*k);
    k3 = Advection(u + 0.5*k*k2, t + 0.5*k);
    k4 = Advection(u + k*k3, t + k);
    
    w = u + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = Advection(u, t)
    
    A1 = Q2*u;
    
    C1 = R1(x', t);

    flux = A1 + C1;

end
end