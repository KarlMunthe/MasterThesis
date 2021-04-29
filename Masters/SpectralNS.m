function [Rho] = SpectralNS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, Q1, Q2)

%h = 1/(m-1);       %space step

t = 0:k:t_end;

%Initial conditions and Primitive variables
R = c_p - c_V;
gamma = c_p/c_V;

rho = (rho(x))';
p = rho.^(gamma);
v = zeros(m,1);
mom = rho.*v;
T = p./(R*rho);
e = c_V*T;
E = e.*rho+0.5*rho.*(v.^2);

nr_points_save = 100;

%Work = [mom.^2./(2*rho), zeros(m, ceil(t_end/k/nr_points_save)-1)];
Rho = [rho, zeros(m, ceil(t_end/k/nr_points_save)-1)];

j = 1;
for i = 1:length(t)-1

    vars = [rho, mom, E];
    
    new_vars = RK4(vars, k);
    
    rho = new_vars(:, 1);
    mom = new_vars(:, 2);
    E = new_vars(:, 3);
    
    if mod(i, nr_points_save) == 1
        %Work(:, j+1) = mom.^2./(2*rho);
        Rho(:, j+1) = rho;
        %MOM(:, j+1) = mom;
        %ENE(:, j+1) = E;
        j = j+1;
    end
    %}
end

%Work = [Work(end,:) ; Work];
Rho = [Rho ; Rho(1,:)];

%{
RHO = rho;
MOM = mom;
ENE = E;

RHO = [RHO(end, :) ; RHO];
MOM = [MOM(end, :) ; MOM];
ENE = [ENE(end, :) ; ENE];
%}

function [w] = RK4(vars, k)
    
    k1 = NS(vars);
    k2 = NS(vars + 0.5*k*k1);
    k3 = NS(vars + 0.5*k*k2);
    k4 = NS(vars + k*k3);
    
    w = vars + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = NS(vars)
    
    rho = vars(:, 1);
    mom = vars(:, 2);
    E = vars(:, 3);
    p = R/c_V*(E - mom.^2./(2*rho));
    T = p/R;
    
    A1 = -Q1*mom;
    
    A2 = -Q1*(mom.^2./rho + p);
    
    B2 = mu*4/3*Q2*(mom./rho);
    
    A3 = -Q1*(E.*(mom./rho) + p.*(mom./rho));
    
    B3 = mu*4/3*1/2*Q2*((mom./rho).^2) + Q1*(kappa*Q1*T);
    
    flux = [A1, A2 + B2, A3 + B3];

end

end
