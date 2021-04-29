function [W, CFL, error, h] = PeriodicAdvectionTest(m, k, sol, x_end, t_end)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

h = 2*pi/m; 
x = h*(1:m);    %individual space points on the grid

CFL = k/h;

%initial condition
f = sol(x, 0)';

%[D0, ~, ~, ~] = SpectralD0(m, h);
D0 = SpectralD0(m, h);


w = [f, zeros(length(f), length(t)-1)];


for i = 1:length(t)-1

    w(:,i+1) = RK4(D0, w(:,i), k);
end

w = [w(end, :) ; w];

ex = [0, x];

error = w(:,end)-sol(ex, t_end)';
error = sqrt(h*error'*error);

W = w;

%
%JUST FUNCTIONS BELOW
%

function [w] = RK4(D0, u, k)
    
    k1 = Advection(D0, u);
    k2 = Advection(D0, u + 0.5*k*k1);
    k3 = Advection(D0, u + 0.5*k*k2);
    k4 = Advection(D0, u + k*k3);
    
    w = u + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = Advection(D0, u)
    
    A1 = -D0*u;
    
    %C1 = G(x(1:end-1)', t);

    flux = A1;

end
end