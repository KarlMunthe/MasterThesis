clear

%MU ER  STØRRELSESORDEN STØRRE ENN DEN SKAL VÆRE

m=20; %  points
h=1/(m-1); % grid size.

k = 1e-4;        %time step
    
x = 0:h:1;    %individual space points on the grid
t = 0:k:5;    %individual time points on the grid

%constants
mu = 1.85e-4;
c_v = 717;
c_p = 1004;
gamma = c_p/c_v;

%Initial conditions and Primitive variables
rho = (2e-04*(sin(2*pi*x(1:end-1)))+1.2)';
%rho = 1.2*ones(m-1,1);
p = rho.^(gamma);
v = zeros(m-1,1);
mom = rho.*v;

%Eulerian viscosity coefficient alpha can equal 1 or 4/3
alpha = 1;
nu = alpha*mu./rho;

w1 = rho;
w2 = mom;

%Creates central first order derivative with periodic boundary conditions
D0 = diag(zeros(m-1, 1), 0) - diag(0.5*ones(m-2, 1), -1) + diag(+0.5*ones(m-2, 1), 1);
D0(end, 1) = 0.5; D0(1, end) = -0.5;
D0 = D0/h;

%Creates second order derivative with periodic boundary conditions
Dp = diag(-ones(m-1, 1), 0) + diag(ones(m-2, 1), 1);
Dm = diag(ones(m-1, 1), 0) + diag(-ones(m-2, 1), -1);
D2 = Dp*Dm;
D2(end,end-1) = 1; D2(end,end) = -2; D2(end,1) = 1; D2(1,end) = 1;
D2 = 2*D2/(h*h); 

P(:,1) = p;
RHO(:,1) = rho;
V(:,1) = v;
MOM(:,1) = mom;


for i = 1:length(t)-1
    
    A1(:,i) = mom;
    B1(:,i) = rho; 
    RHS1 = -D0*A1(:,i) + D0*(nu.*(D0*B1(:,i)));
    
    A2(:,i) = mom.*v + p;
    B2(:,i) = mom;
    RHS2 = -D0*A2(:,i) + D0*(nu.*(D0*B2(:,i)));
    
    w1(:,i+1) = RK4(w1(:,i), RHS1, k);
    w2(:,i+1) = RK4(w2(:,i), RHS2, k);
    
    
    %Primitive variables
    rho = w1(:,end);
    mom = w2(:,end);
    v = mom./rho;
    p = rho.^gamma;
    P(:,i+1) = p;
    RHO(:, i+1) = rho;
    MOM(:, i+1) = mom;
    V(:, i+1) = v;
    
  
end

%EulerianIsentropicPressure = P;

%{
f1 = figure('Name','Eulerian Model');
tiledlayout(2,1) 
nexttile
plot(t,EulerianIsentropicPressure(2,:))
title('left side'), xlabel('t'), ylabel('p')
%}
%nexttile

%%%%%%%%%
%LOOK HERE, SUPPOSED TO BE LIKE THIS CUZ FIRST AND LAST VALUE ARE "DEFINED"
%AS THE SAME??


%mesh(x(1:end),t,P')
%title('Pressure'), xlabel('x'), ylabel('t'), zlabel('P')

P(end+1,:) = P(1,:);
plot(t,P(end/4,:))
title('right side'), xlabel('t'), ylabel('p')

function [w] = RK4(v, RHS, k)
    
    k1 = RHS;
    k2 = RHS + k*k1/2;
    k3 = RHS + k*k2/2;
    k4 = RHS + k*k3;
    
    w = v + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end


