function [Work, Heat] = E(m, x, t_end, k, mu, alpha, c_v, c_p, R, D0, rho)
%constants
%mu = dynamic viscosity
%c_v = heat capacity for a constant volume
%c_p = heat capacity for a constant pressure
%R =  universal gas constant
%K = coefficient in  Fourier's law
%rho =  density initial condition

%h = 1/(m-1);
%x = 0:h:1;

%[~, ~, ~, D0] = PeriodicD0(m, h);

%Initial conditions and Primitive variables
gamma = c_p/c_v;
rho = (rho(x(1:end-1)))';
p = rho.^(gamma);
v = zeros(m-1,1);
mom = rho.*v;

T = p./(R*rho);
e = c_v*T;
E = e.*rho+0.5*rho.*(v.^2);

%Eulerian viscosity coefficient
nu = alpha*mu./rho;

%P(:,1) = p;
%P = zeros(length(t),1);
%P = [p(6), zeros(1,ceil((t_end/k)/1000))];
%P = [zeros(ceil(t_end/k)+1, m-1); p'];

%V = [zeros(ceil(t_end/k)+1, m-1); v'];

%E_mech = [E, zeros(m-1, ceil(t_end/k)/t_end)];
Work = [rho.*(v.^2), zeros(m-1, ceil(t_end/k)/t_end)];
%Heat = [E-0.5*rho.*v.^2, zeros(m-1, ceil(t_end/k)/t_end)];

t=0;
i=1;
j=1;

while t < t_end
    
    A1 = mom;
    B1 = rho; 
    RHS1 = -D0*A1 + D0*(nu.*(D0*B1));
    
    A2 = mom.*v + p;
    B2 = mom;
    RHS2 = -D0*A2 + D0*(nu.*(D0*B2));
    
    A3 = E.*v;
    B3 = E;
    RHS3 = -D0*A3 + D0*(nu.*(D0*B3));
    
    %Update primitive variables
    rho = RK4(rho, RHS1, k);
    mom = RK4(mom, RHS2, k);
    E = RK4(E, RHS3, k);
    
    %Primitive variables
    %rho = new_rho;
    %mom = new_mom;
    %E = new_E;
    
    v = mom./rho;
    e = E./rho-0.5*v.^2;
    T = e/c_v;
    p = R*rho.*T;
   
    %P = p';
    %P(:,i+1) = p;
    %P(end-i,:) = p';
    %P(i+1) = p(6);
    %{
    i = i+1;
    if mod(i,t_end) == 1
        P(j+1) = p(6);
        j = j+1;
    end
    %}
    %V(end-i,:) = v';
    %{
    i = i+1;
    if mod(i,t_end) == 1
        E_mech(:, j+1) = E;
        j = j+1;
    end
    %}
    
   
    if mod(i, 10) == 1
        %Heat(:, j+1) = E-0.5*rho.*v.^2;
        Work(:, j+1) = rho.*(v.^2);
        %E_mech(:, j+1) = E;
        j = j+1;
    end
    
    i = i+1;
    
    t = t+k;
    
end

%P(end+1) = P(1);
%P(:,end+1) = P(:,1);
%V(:, end+1) = V(:,1);
%E_mech(end+1, :) = E_mech(1, :);
Work(end+1, :) = Work(1, :);
%Heat(end+1, :) = Heat(1, :);

function [w] = RK4(v, RHS, k)
    
    k1 = RHS;
    k2 = RHS + k*k1/2;
    k3 = RHS + k*k2/2;
    k4 = RHS + k*k3;
    
    w = v + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end
%}
end
