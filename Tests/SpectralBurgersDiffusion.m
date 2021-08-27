function [W] = SpectralBurgersDiffusion(m, epsilon, t_end, x_end, k)

%h = 1/(m-1);       %space step

addpath('/Users/karlmunthe/Documents/UiB/UiB Master Oppgave/Code/gitMaster/Difference Operators')

t = 0:k:t_end;

h = 1/m;    %individual space points on the grid
x = h*(1:m);

s = 2*pi/x_end;

CFL = k/(h^2);

%initial condition
f = sin(2*pi*x)';

%D1 = SpectralD0(m, s);
%D2 = SpectralD2(m, s);

D1 = FD0(m, h, 50);
D2 = FD02(m, h, 50);

w = [f, zeros(length(f), length(t)-1)];

for i = 1:length(t)-1
    
    w(:,i+1) = RK4(x, D1, D2, w(:,i), k, t(i));
    
end

w = [w(end, :) ; w];

%ex = [0, x];

%error = w(:,end)-sol(ex, t_end)';
%error = sqrt(h*error'*error);

W = w;


%%%
%%%JUST FUNCTIONS BELOW
%%%

function [w] = RK4(x, D0, D2, u, k, t)
    
    k1 = Burgers(x, D0, D2, u, t);
    k2 = Burgers(x, D0, D2, u + 0.5*k*k1, t + 0.5*k);
    k3 = Burgers(x, D0, D2, u + 0.5*k*k2, t + 0.5*k);
    k4 = Burgers(x, D0, D2, u + k*k3, t + k);
    
    w = u + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = Burgers(x, D0, D2, u, t)
    
    A1 = -D0*(0.5*(u.^2));
    B1 = epsilon*D2*u;
    %C1 = G(x', t);

    flux = A1 + B1;% + C1;

end


end
