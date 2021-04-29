clear

m=10; %  points
h=1/(m-1); % grid size.

k = 1e-4;        %time step
    
x = 0:h:1;    %individual space points on the grid
t = 0:k:80;    %individual time points on the grid

%constants
mu = 1.85e-4;
c_v = 717;
c_p = 1004;
gamma = c_p/c_v;
R = 287;
K = 0.026;

%Initial conditions and Primitive variables
rho = (2e-04*(sin(2*pi*x(1:end-1)))+1.2)';
p = rho.^(gamma);
v = zeros(m-1,1);
mom = rho.*v;
T = p./(R*rho);
e = c_v*T;
E = e.*rho+0.5*rho.*(v.^2);

%Creates initial vector
w1 = rho;
w2 = mom;
w3 = E;

%Creates first order derivative with periodic boundary conditions
D1 = diag(zeros(m-1, 1), 0) - diag(0.5*ones(m-2, 1), -1) + diag(+0.5*ones(m-2, 1), 1);
D1(end, 1) = 0.5; D1(1, end) = -0.5;
D1 = D1/h;

%Creates second order derivative with periodic boundary conditions
Dp = diag(-ones(m-1, 1), 0) + diag(ones(m-2, 1), 1);
Dm = diag(ones(m-1, 1), 0) + diag(-ones(m-2, 1), -1);
D2 = Dp*Dm;
D2(end,end-1) = 1; D2(end,end) = -2; D2(end,1) = 1; D2(1,end) = 1;
D2 = 2*D2/(h*h);

P(:,1) = p;
RHO(:,1) = rho;
MOM(:,1) = mom;

for i = 1:length(t)-1

    A1(:,i) = mom;
    B1(:,i) = zeros(m-1,1); 
    RHS1 = -D1*A1(:,i) + D2*B1(:,i);
    
    A2(:,i) = mom.*v + p;
    B2(:,i) = 4/3*mu*v;
    RHS2 = -D1*A2(:,i) + D2*B2(:,i);
    
    A3(:,i) = E.*v;
    B3(:,i) = 4/6*mu*(v.^2) + K*T;
    RHS3 = -D1*A3(:,i) + D2*B3(:,i);
    
    w1(:,i+1) = RK4(w1(:,i), RHS1, k);
    w2(:,i+1) = RK4(w2(:,i), RHS2, k);
    w3(:,i+1) = RK4(w3(:,i), RHS3, k);
    
    %Primitive variables
    rho = w1(:,end);
    mom = w2(:,end);
    E = w3(:,end);
    v = mom./rho;
    e = E./rho-0.5*v.^2;
    T = e/c_v;
    p = R*rho.*T;
    P(:,i+1) = p;
    RHO(:, i+1) = rho;
    MOM(:, i+1) = mom;
        
end

tiledlayout(2,1) 
nexttile
plot(t,P(2,:))
title('left side'), xlabel('t'), ylabel('p')

nexttile
plot(t,P(round(end),:))
title('right side'), xlabel('t'), ylabel('p')

function [w] = RK4(v, RHS, k)
    
    k1 = RHS;
    k2 = RHS + k*k1/2;
    k3 = RHS + k*k2/2;
    k4 = RHS + k*k3;
    
    w = v + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end