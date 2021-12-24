function [Work] = SpectralNS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, p_0, Q1, Q2)

%create time vector
t = 0:k:t_end;

%Initial conditions and Primitive variables
R = c_p - c_V;
gamma = c_p/c_V;

%initial conditions
rho_0 = (rho(x))';
u_0 = zeros(length(x), 1);
mom_0 = u_0.*rho_0;
p_0 = p_0*rho_0.^(gamma);
%p = c^2/gamma*rho_0;
T_0 = p_0./(R*rho_0);
e_0 = c_V*T_0;
E_0 = e_0.*rho_0 + 0.5*rho_0.*(u_0.^2);

%used to save every 5th point to save memory
nr_points_save = 5;

%save work variable. Uncomment other lines to save other variables. Also
%creates a list to store the other variables.
Work = [mom_0.^2./(2*rho_0), zeros(m, ceil(t_end/k/nr_points_save)-1)];
%Rho = [rho_0, zeros(m, ceil(t_end/k/nr_points_save)-1)];
%Mom = [mom_0, zeros(m, ceil(t_end/k/nr_points_save)-1)];
%Ene = [E_0, zeros(m, ceil(t_end/k/nr_points_save)-1)];

%renaming initial conditions for the for-loop
rho = rho_0;
mom = mom_0;
E = E_0;

j = 1;
for i = 1:length(t)-1
    
    %The conservative variables.
    vars = [rho, mom, E];
    
    new_vars = RK4(vars, k);
    
    %update conservative variables
    rho = new_vars(:, 1);
    mom = new_vars(:, 2);
    E = new_vars(:, 3);
    
    %saves each n'th point
    if mod(i, nr_points_save) == 0
        Work(:, j+1) = mom.^2./(2*rho);
    j = j+1;
    end
    
end

%Ensure periodicity
Work = [Work(end,:) ; Work];
%Rho = [Rho(end, :) ; Rho];
%Mom = [Mom(end, :) ; Mom];
%Ene = [Ene(end, :) ; Ene];

%4th order Runkge-Kutta method
function w = RK4(vars, k)
    
    k1 = NS(vars);
    k2 = NS(vars + 0.5*k*k1);
    k3 = NS(vars + 0.5*k*k2);
    k4 = NS(vars + k*k3);
    
    w = vars + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function flux = NS(vars)
    
    rho = vars(:, 1);
    mom = vars(:, 2);
    E = vars(:, 3);
    
    p = R/c_V*(E - mom.^2./(2*rho));
    T = p./(rho*R);
    
    A1 = -Q1*mom;               %Advective term in conservation of mass equation
    
    A2 = -Q1*(mom.^2./rho + p); %Advective term in conservation of momentum equation
    
    B2 = mu*4/3*Q2*(mom./rho);  %Advective term in conservation of momentum equation
    
    A3 = -Q1*(E.*(mom./rho) + p.*(mom./rho));   %Advective term in conservation of energy equation
    
    B3 = mu*4/3*1/2*Q2*((mom./rho).^2) + Q1*(kappa*Q1*T); %Diffusive term in conservation of energy equation
    
    flux = [A1, A2 + B2, A3 + B3];  %Vector containing mass, momentum and energy components.

end

end
