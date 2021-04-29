function [W, CFL, error, h] = SpectralBurgersDiffusion(m, epsilon, G, sol, t_end, x_end, k)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

h = 2*pi/m;    %individual space points on the grid
x = h*(1:m);

CFL = k/(h^2);

%initial condition
f = sol(x, 0)';

D1 = SpectralD0(m, h);
D2 = SpectralD2(m, h);

w = [f, zeros(length(f), length(t)-1)];

for i = 1:length(t)-1
    
    w(:,i+1) = RK4(x, G, D1, D2, w(:,i), k, t(i));
    
end

w = [w(end, :) ; w];

ex = [0, x];

error = w(:,end)-sol(ex, t_end)';
error = sqrt(h*error'*error);

W = w;


%%%
%%%JUST FUNCTIONS BELOW
%%%

function [w] = RK4(x, G, D0, D2, u, k, t)
    
    k1 = Burgers(x, G, D0, D2, u, t);
    k2 = Burgers(x, G, D0, D2, u + 0.5*k*k1, t + 0.5*k);
    k3 = Burgers(x, G, D0, D2, u + 0.5*k*k2, t + 0.5*k);
    k4 = Burgers(x, G, D0, D2, u + k*k3, t + k);
    
    w = u + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = Burgers(x, G, D0, D2, u, t)
    
    A1 = -D0*(0.5*(u.^2));
    B1 = epsilon*D2*u;
    C1 = G(x', t);

    flux = A1 + B1 + C1;

end


end
