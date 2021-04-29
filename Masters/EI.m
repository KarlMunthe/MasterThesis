function [P] = EI(m, x, t, k, mu, alpha, c_v, c_p, D0, rho)


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
p = rho.^(gamma);
v = zeros(m-1,1);
mom = rho.*v;

%Eulerian viscosity coefficient
nu = alpha*mu./rho;

%Creates initial vector
w1 = rho;
w2 = mom;

%P(:,1) = p;
%P = zeros(
P(1) = p((end+1)/4);

for i = 1:length(t)-1
    
    A1(:,i) = -mom;
    B1(:,i) = rho; 
    
    A2(:,i) = -mom.*v - p;
    B2(:,i) = mom;
   
    rho = RK4(w1(:,i), A1, B1, k);
    mom = RK4(w2(:,i), A2, B2, k);
    
    v = mom./rho;
    p = rho.^gamma;
    
    %P(:,i+1) = p;
    P(i+1) = p((end+1)/4);
        
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
