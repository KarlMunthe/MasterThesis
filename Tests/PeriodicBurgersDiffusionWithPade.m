function [W, CFL, error, h] = PeriodicBurgersDiffusionWithPade(m, G, sol, t_end, x_end, k)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);


[Q1, ~, ~, ~] = PeriodicD0(m, h);
[Q2, ~, ~, ~] = PeriodicD2(m, h);

%REMEMBER TO USE SECOND ORDER, O(2), DIFFERENCE OPERATORS WHEN USING PADE
%SCHEMES, THE FIRST SIMPLESt OF THE ONES IN THE PERIODICD0/2 FUNCTIONS
PD1 = (eye(m-1) + h^2*Q2/6);% - h^4*Q2^2/30 + h^6*(Q2)^3/140);
PD2 = (eye(m-1) + h^2*Q2/12);% - h^4*(Q2)^2/90+ h^6*(Q2)^3/560);
%}

%initial condition
f = sol(x(1:end-1), 0)';

w = [f, zeros(length(f), length(t)-1)];

for i = 1:length(t)-1
    
    w(:,i+1) = RK4(x, G, w(:,i), k, t(i));
    
end

w(end + 1, :) = w(1, :);

W = w;

error = w(:,end) - sol(x, t_end)';
error = sqrt(h*error'*error);

%%%
%%%JUST FUNCTIONS BELOW
%%%

function [w] = RK4(x, G, u, k, t)
    
    k1 = Burgers(x, G, u, t);
    k2 = Burgers(x, G, u + 0.5*k*k1, t + 0.5*k);
    k3 = Burgers(x, G, u + 0.5*k*k2, t + 0.5*k);
    k4 = Burgers(x, G, u + k*k3, t + k);
    
    w = u + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = Burgers(x, G, u, t)
    
    A1 = PD1\(-Q1*(0.5*(u.^2)));
    B1 = PD2\(Q2*u);
    C1 = G(x(1:end-1)', t);

    flux = A1 + B1 + C1;

end


end
