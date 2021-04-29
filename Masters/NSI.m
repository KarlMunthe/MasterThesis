function [P] = NSI(m, x, t, k, mu, c_v, c_p, D0, D2, rho)

%constants
%mu = dynamic viscosity
%c_v = heat capacity for a constant volume
%c_p = heat capacity for a constant pressure
%R =  universal gas constant
%K = coefficient in  Fourier's law
%rho =  density initial condition

%Initial conditions and Primitive variables
gamma = c_p/c_v;
rho = (rho(x(1:end-1)))';
p = rho.^gamma;
v = zeros(m-1,1);
mom = rho.*v;

%P = [p,zeros(m-1,length(t)-1)];
P = [p((end+1)/4), zeros(1,length(t)-1)];

for i = 1:length(t)-1

    A1 = mom;
    B1 = zeros(m-1,1); 
    
    A2 = mom.*v + p;
    B2 = 4/6*mu*v;
    
    rho = RK4(rho, RHS1, k);
    mom = RK4(mom, RHS2, k);
    
    v = mom./rho;
    p = rho.^gamma;
    
    %P(:,i+1) = p;
    P(i+1) = p((end+1)/4);
    
    %{
    if mod(i,10) == 1
        P(end+1) = p((end+1)/4);
    end
    %}
    
end

%P(end+1,:) = P(1,:);

function [w] = RK4(v, A, B, k)
    
    k1 = D0*A + D2*B;
    k2 = D0*(A + k*k1/2) + D2*(B + k*k1/2);
    k3 = D0*(A + k*k2/2) + D2*(B + k*k2/2);
    k4 = D0*(A + k*k3) + D2*(B + k*k3);
    
    w = v + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

end
