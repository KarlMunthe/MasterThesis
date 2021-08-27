clear

m = 10;

k = 1e-6; % time step

t = 0:k:1; %time

x_end = 1;
t_end = 1;

h = x_end/m;    %individual space points on the grid
x = h*(1:m);

y = h*(1:m);
ytrans = (exp(y)-1)/(exp(1)-1);

[xx,yy] = meshgrid(x,y);

mu = 20.64e-06;
kappa = 26.58e-03;
c_p = 915;
c_V = 659;

p_0 = 1e+03;
R = c_p - c_V;
gamma = c_p/c_V;

rho = @(ex, ey) exp(-ex.^2 + ey.^2);

rho_0 = (rho(x, y))';
u_0 = zeros(length(x), 1);
mom_0 = u_0.*rho_0;
p_0 = p_0*rho_0.^(gamma);
T_0 = p_0./(R*rho_0);
e_0 = c_V*T_0;
E_0 = e_0.*rho_0 + 0.5*rho_0.*(u_0.^2);

%renaming for for-loop
rho = rho_0;
mom = mom_0;
E = E_0;

Q1 = FD0(

for i = 1:length(t)-1

    vars = [rho, mom, E];
    
    new_vars = RK4(vars, k);
    
    rho = new_vars(:, 1);
    mom = new_vars(:, 2);
    E = new_vars(:, 3);
        
    if mod(i, nr_points_save) == 0
        Work(:, j+1) = mom.^2./(2*rho);
    %Rho(:, j+1) = rho;
    %Mom(:, j+1) = mom;
    %Ene(:, j+1) = E;
    j = j+1;
    end
    
        
end

function [w] = RK4(vars, k)
    
    k1 = NS(vars);
    k2 = NS(vars + 0.5*k*k1);
    k3 = NS(vars + 0.5*k*k2);
    k4 = NS(vars + k*k3);
    
    w = vars + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = NS(vars)

    mu = 20.64e-06;
    kappa = 26.58e-03;
    c_p = 915;
    c_V = 659;
    R = c_p-c_V;
    
    rho = vars(:, 1);
    mom = vars(:, 2);
    E = vars(:, 3);
    p = R/c_V*(E - mom.^2./(2*rho));
    T = p./(rho*R);
    
    A1 = -Q1*mom;
    
    A2 = -Q1*(mom.^2./rho + p);
    
    B2 = mu*4/3*Q2*(mom./rho);
    
    A3 = -Q1*(E.*(mom./rho) + p.*(mom./rho));
    
    B3 = mu*4/3*1/2*Q2*((mom./rho).^2) + Q1*(kappa*Q1*T);
    
    flux = [A1, A2 + B2, A3 + B3];

end