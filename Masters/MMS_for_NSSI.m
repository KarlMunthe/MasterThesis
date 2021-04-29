%METHOD OF MANUFACTURED SOLUTION
% manufactured solutions:
%\rho = sin(x+t) + const
%u = cos(x+t)
%nu = cos(x+t) + const
%p = sin(x+t) + const
%A = cos(x+t) + cos(2(x+t)) + 1
%B = 2(cos(x+t) - sin(x+t)) + (cos(x+t) + 3*cos(3(x+t)))/4 + sin(x+t) - 2
%if all constants are set to 1

function [RHO] = MMS_for_NSSI(m, x, t_end, k, c_v, c_p, D0, rho, mom, C1, C2)


%constants
%mu = dynamic viscosity
%c_v = heat capacity for a constant volume
%c_p = heat capacity for a constant pressure
%R =  universal gas constant
%K = coefficient in  Fourier's law
%rho =  density initial condition

%Initial conditions and Primitive variables
gamma = c_p/c_v;
rho = (rho(x(1:end-1), 0))';
p = rho.^gamma;
mom = (mom(x(1:end-1), 0))';

%Eulerian viscosity coefficient
nu = 1./rho;

%P(:,1) = p;
%P = zeros(
%P(1) = p((end+1)/4);

RHO = [rho, zeros(m-1, round(t_end/k))];

t = 0;
i = 1;

while t <= t_end
    
    A1 = mom;
    B1 = rho;
    
    A2 = (mom.^2)./rho + p;
    B2 = mom;
    
    rho = RK4(rho, A1, B1, C1, k);
    mom = RK4(mom, A2, B2, C2, k);
    
    p = rho.^(gamma);
    %P(:,i+1) = p;
    %P(i+1) = p((end+1)/4);
    RHO(:, i+1) = rho;
    i = i+1;
    t = t+k;
        
end

%P(end+1,:) = P(1,:);

RHO(end+1, :) = RHO(1, :);

function [w] = RK4(v, A, B, C, k)
    
    k1 = -D0*A + D0*(nu.*(D0*B)) + C(x(1:end-1)', t);
    k2 = -D0*(A + 0.5*k*k1) + D0*(nu.*(D0*B + 0.5*k*k1)) + C(x(1:end-1)', t + 0.5*k);
    k3 = -D0*(A + 0.5*k*k2) + D0*(nu.*(D0*B + 0.5*k*k2)) + C(x(1:end-1)', t + 0.5*k);
    k4 = -D0*(A + k*k3) + D0*(nu.*(D0*B + k*k3)) + C(x(1:end-1)', t + k);
    
    w = v + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

end
