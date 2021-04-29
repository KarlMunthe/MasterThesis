function [W, error, h, CFL] = MMSTest(m, G, sol, t_end, k)

%h = 1/(m-1);       %space step

t = k:k:t_end;

x = linspace(0, 2*pi, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);

%initial condition
f = sol(x(1:end-1), 0)';

[~, D0, ~, ~] = PeriodicD0(m, h);
[~, D2, ~, ~] = PeriodicD2(m, h);

w = [f, zeros(length(f), length(t))];
    
te = k;

for i = 1:length(t)
    
    w(:,i+1) = RK4(x, G, D0, D2, w(:,i), k, te);

    te = te + k;
    
end

w(end + 1, :) = w(1, :);

W = w;

error = w(:,end) - sol(x, t_end)';
error = sqrt(h*error'*error);

%%%
%%%JUST FUNCTIONS BELOW
%%%

function [w] = RK4(x, G, D0, D2, u, k, te)
    
    k1 = Burgers(x, G, D0, D2, u, te);
    k2 = Burgers(x, G, D0, D2, u + 0.5*k*k1, te + 0.5*k);
    k3 = Burgers(x, G, D0, D2, u + 0.5*k*k2, te + 0.5*k);
    k4 = Burgers(x, G, D0, D2, u + k*k3, te + k);
    
    w = u + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = Burgers(x, G, D0, D2, u, te)
    
    A1 = -D0*(0.5*(u.^2));
    B1 = D2*u;
    C1 = G(x(1:end-1)', te-k);

    flux = A1 + B1 + C1;

end


end
