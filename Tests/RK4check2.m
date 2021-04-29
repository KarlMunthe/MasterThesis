%METHOD OF MANUFACTURED SOLUTION
% manufactured solutions:
%\rho = sin(x+t) + const
%u = cos(x+t)
%nu = cos(x+t) + const
%p = sin(x+t) + const
%A = cos(x+t) + cos(2(x+t)) + 1
%B = 2(cos(x+t) - sin(x+t)) + (cos(x+t) + 3*cos(3(x+t)))/4 + sin(x+t) - 2
%if all constants are set to 1

function [W, CFL, error, h] = RK4check2(m, t_end, x_end, k, c_v, c_p, rho, u, D1, D2, rho_solution)

%constants
%mu = dynamic viscosity
%c_v = heat capacity for a constant volume
%c_p = heat capacity for a constant pressure
%R =  universal gas constant
%K = coefficient in  Fourier's law
%rho =  density initial condition

%h = 1/(m-1);       %space step

mu = 1;

t = k:k:t_end;

x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);


[D0, ~, ~, ~] = PeriodicD0(m, h);
[D02, ~, ~, ~] = PeriodicD2(m, h);

%Initial conditions and Primitive variables
gamma = 2;
rho = (rho(x(1:end-1), 0))';
u = (u(x(1:end-1), 0))';
p = rho.^gamma;

%NSS viscosity coefficient


RHO = [rho, zeros(length(rho), length(t))];

te = k;

for i = 1:length(t)

    vars = [rho, u];

    new_vars = RK4(vars, p, k, te);

    rho = new_vars(:,1);
    u = new_vars(:,2);
    p = rho.^gamma;

    RHO(:, i+1) = rho;
    
    te = te+k;

end

RHO(end+1, :) = RHO(1, :);

W = RHO;

error = (RHO(:,end) - rho_solution(x, t_end)');
error = sqrt(h*error'*error);

function [w] = RK4(vars, p, k, t)
    
    k1 = INS(vars, p, t);
    k2 = INS(vars + 0.5*k*k1, p, t + 0.5*k);
    k3 = INS(vars + 0.5*k*k2, p, t + 0.5*k);
    k4 = INS(vars + k*k3, p, t + k);
    
    w = vars + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = INS(vars, p, t)

    rho = vars(:,1);
    u = vars(:,2);
    
    A1 = -D0*(rho.*u);

    A2 = -D0*(rho.*u.^2 + p);
    B2 = D02*(4*mu*u/3);

    C1 = D1(x(1:end-1)', t-k);
    C2 = D2(x(1:end-1)', t-k);

    flux = [A1 + C1, A2 + B2 + C2];

end

end
