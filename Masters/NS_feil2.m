function [Work] = NS2(m, x, t_end, k, mu, c_v, c_p, K, rho, Q1, Q2)

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
%[~, ~, ~, D2] = PeriodicD2(m, h);

%Initial conditions and Primitive variables
R = c_p - c_V;
gamma = c_p/c_v;

rho = (rho(x(1:end-1)))';
p = rho.^(gamma);
v = zeros(m-1,1);
mom = rho.*v;

T = p./(R*rho);
e = c_v*T;
E = e.*rho+0.5*rho.*(v.^2);

%{
%Creates initial vector 
w1 = rho;
w2 = mom;
w3 = E;
%}

%{
%Determine size of matrices to reduce memory usage
A1 = zeros(m-1,length(t)-1);
A2 = zeros(m-1,length(t)-1);
A3 = zeros(m-1,length(t)-1);

B1 = zeros(m-1,length(t)-1);
B2 = zeros(m-1,length(t)-1);
B3 = zeros(m-1,length(t)-1);
%}

%P = [p, zeros(m-1,length(t)-1)];

%P2 = 1;
%P = [p(6), zeros(1,ceil((t_end/k)/1000))];
%P = [zeros(ceil(t_end/k)+1, m-1); p'];

%E_mech = [E, zeros(m-1, ceil(t_end/k)/t_end)];
Work = [rho.*(v.^2), zeros(m-1, ceil(t_end/k)/t_end)];
%Heat = [E-0.5*rho.*v.^2, zeros(m-1, ceil(t_end/k)/t_end)];

t=0;
i=1;
j=1;
%P2 = [p(3*((end+1)/4)), zeros(1,length(t)-1)];

while t < t_end

    A1 = mom;
    B1 = 0;
    %RHS1 = -D0*A1;
    
    A2 = mom.*v + p;
    B2 = 4/3*mu*v;
    %RHS2 = -D0*A2 + D2*B2;
    
    A3 = E.*v;
    B3 = 4/6*mu*(v.^2) + K*T;
    %RHS3 = -D0*A3 + D2*B3;
    
    %Update primitive variables
    rho = RK4(rho, A1, B1, k);
    mom = RK4(mom, A2, B2, k);
    E = RK4(E, A3, B3, k);
    
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
    %i = i+1;
    %{
    if mod(i,t_end) == 1
        P(j+1) = p(6);
        j = j+1;
    end
    %}
    %{
    i = i+1;
    if mod(i,t_end) == 1
        E_mech(:, j+1) = E;
        j = j+1;
    end
    %}
    
   
    if mod(i, 10) == 1
       %Heat(:, j+1) = E - 0.5*rho.*v.^2;
        Work(:, j+1) = rho.*(v.^2);
       %E_mech(:, j+1) = E;
        j = j+1;
    end
    
    i = i+1;
    
    t = t+k;
    
end

%P(end+1) = P(1);
%P(:,end+1) = P(:,1);

%E_mech(end+1, :) = E_mech(1, :);
Work(end+1, :) = Work(1, :);
%Heat(end+1, :) = Heat(1, :);

function [w] = RK4(v, A, B, k)
    
    k1 = Q1*A + Q2*B;
    k2 = Q1*(A + k*k1/2) + Q2*(B + k*k1/2);
    k3 = Q1*(A + k*k2/2) + Q2*(B + k*k2/2);
    k4 = Q1*(A + k*k3) + Q2*(B + k*k3);
    
    w = v + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

end
