function [Work] = SpectralNSS(m, x, t_end, k, c_p, c_V, mu, rho_0, T_0, u, Q1)

t = 0:k:t_end;

%Thermodynamic variables
R = c_p - c_V;
gamma = c_p/c_V;

%rho = (rho(x))';
%p = p_0*rho.^(gamma);
%p = c^2/gamma*rho_0;
%T = p./(R*rho);

%initial conditions
rho_0 = ones(length(x),1).*rho_0;
u_0 = u(x)';
mom_0 = rho_0.*u_0;
e_0 = c_V*T_0;
E_0 = e_0.*rho_0+0.5*rho_0.*(u_0.^2);

%every nth time step is saved
nr_points_save = 10;

%saving memory for variables to save
Work = [mom_0.^2./(2*rho_0), zeros(m, ceil(t_end/k/nr_points_save)-1)];
%Rho = [rho_0, zeros(m, ceil(t_end/k/nr_points_save)-1)];
%Mom = [mom_0, zeros(m, ceil(t_end/k/nr_points_save)-1)];
%Ene = [E_0, zeros(m, ceil(t_end/k/nr_points_save)-1)];

%renaming initial conditions for the for-loop
rho = rho_0;
mom = mom_0;
E = E_0;

%setting j=1 to add points to the variables lists
j = 1;
for i = 1:length(t)-1

    vars = [rho, mom, E];
    
    new_vars = RK4(vars, k);
    
    rho = new_vars(:, 1);
    mom = new_vars(:, 2);
    E = new_vars(:, 3);
        
    if mod(i, nr_points_save) == 1
        Work(:, j+1) = mom.^2./(2*rho);
        %Rho(:, j+1) = rho;
        %Mom(:, j+1) = mom;
        %Ene(:, j+1) = E;
        j = j+1;
    end
    
        
end

%Setting the zeroths values equal to the end values due to periodicity
Work = [Work(end,:) ; Work];
%Rho = [Rho(end,:) ; Rho];
%Mom = [Mom(end,:) ; Mom];
%Ene = [Ene(end,:) ; Ene];

%{
RHO = rho;
MOM = mom;
ENE = E;

RHO = [RHO(end, :) ; RHO];
MOM = [MOM(end, :) ; MOM];
ENE = [ENE(end, :) ; ENE];
%}

function [w] = RK4(vars, k)
    
    k1 = NSS(vars);
    k2 = NSS(vars + 0.5*k*k1);
    k3 = NSS(vars + 0.5*k*k2);
    k4 = NSS(vars + k*k3);
    
    w = vars + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = NSS(vars)
    
    rho = vars(:,1);
    mom = vars(:,2);
    E = vars(:,3);
    
    nu = mu./rho;
    p = R/c_V*(E - mom.^2./(2*rho));
    
    A1 = -Q1*mom;
    
    B1 = Q1*(nu.*(Q1*rho));
        
    A2 = -Q1*((mom.^2)./rho + p);
    
    B2 = Q1*(nu.*(Q1*mom));
    
    A3 = -Q1*(E.*(mom./rho) + p.*(mom./rho));
    
    B3 = Q1*(nu.*(Q1*E));
    
    flux = [A1 + B1, A2 + B2, A3 + B3];

    %flux = [A1, A2, A3];
    
end

end
