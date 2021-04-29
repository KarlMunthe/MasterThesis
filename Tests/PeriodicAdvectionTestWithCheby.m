function [W, CFL, error, h] = PeriodicAdvectionTestWithCheby(m, k, sol, x_end, t_end)

%h = 1/(m-1);       %space step

t = k:k:t_end;

x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);


%xi = (2*pi*x + sin(2*pi*x))/(2*pi);
xi = (exp(x)-1)/(exp(1)-1);     %transformed grid

%{
cheby_points = 1:m-2;
xi = 0.5*(0 + x_end) + 0.5*(x_end - 0)*cos(pi*(2*cheby_points-1)/(2*(m-2)));
xi = xi(end:-1:1);
xi = [0, xi, x_end];
%}

%J = 1 + cos(2*pi*x);
%J = diag(J(1:end-1));

%Creating matrix norm for non-linear grid
J = 0:h:1;
J = exp(J(1:end-1))/(exp(1)-1);
J = diag(J);


CFL = k/h;

%initial condition
f = sol(xi(1:end-1), 0)';

[D0, ~, ~, ~] = PeriodicD0(m, h);


w = [f, zeros(length(f), length(t))];

te = k;

for i = 1:length(t)

    w(:,i+1) = RK4(D0, w(:,i), k);
    
    te = te + k;
    
end

w(end + 1, :) = w(1, :);

error = w(:,end)-sol(x, t_end)';
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
    
    A1 = (-D0*u);
    
    %C1 = G(x(1:end-1)', t);

    flux = A1;

end
end