function [W, error, h, CFL] = HeatTest(m, k, sol, G, t_end)

%h = 1/(m-1);       %space step

t = k:k:t_end;

x = linspace(0, 2*pi, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);

%initial condition
f = sol(x(1:end-1), 0)';

[~, D2, ~, ~] = PeriodicD2(m, h);


w = [f, zeros(length(f), length(t))];

te = k;

for i = 1:length(t)

    w(:,i+1) = RK4(D2, w(:,i), k, te);
    
    te = te + k;
    
end

w(end + 1, :) = w(1, :);

error = w(:,end)-sol(x, t_end)';
error = sqrt(h*error'*error);

W = w;

%%%
%%%JUST FUNCTIONS BELOW
%%%

function [w] = RK4(D2, u, k, t)
    
    k1 = Advection(D2, u, t);
    k2 = Advection(D2, u + 0.5*k*k1, t + 0.5*k);
    k3 = Advection(D2, u + 0.5*k*k2, t + 0.5*k);
    k4 = Advection(D2, u + k*k3, t + k);
    
    w = u + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = Advection(D2, u, t)
    
    A1 = D2*u;
    
    C1 = cos(x(1:end-1) + t-k)' + sin(x(1:end-1) + t-k)';
    %C1 = G(x(1:end-1), t-k)';

    flux = A1 + C1;

end
end