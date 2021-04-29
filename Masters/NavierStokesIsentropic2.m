clear

%MU ER  STØRRELSESORDEN STØRRE ENN DEN SKAL VÆRE

m=20; %  points
h=1/(m-1); % grid size.

k = 1e-4;        %time step
    
x = 0:h:1;    %individual space points on the grid
t = 0:k:160;    %individual time points on the grid

%constants
mu = 1.85e-04;
c_p = 1004;
c_v = 717;
gamma = c_p/c_v;

%Initial conditions and Primitive variables
rho = (2e-04*(sin(2*pi*x(1:end-1)))+1.2)';
%rho = 1.2*ones(m-1,1);
p = rho.^(gamma);
v = zeros(m-1,1);
mom = rho.*v;

w1 = rho;
w2 = mom;

%Creates first order derivative with periodic boundary conditions
%{
D0 = diag(zeros(m-1, 1), 0) - diag(0.5*ones(m-2, 1), -1) + diag(+0.5*ones(m-2, 1), 1);
D0(end, 1) = 0.5; D0(1, end) = -0.5;
D0 = D0/h;
%}


%8th order central difference method first derivative
D0 = -diag(4/5*ones(m-2, 1), -1) + diag(4/5*ones(m-2, 1), 1)...
    + diag(1/5*ones(m-3, 1), -2) - diag(1/5*ones(m-3, 1), 2) - diag(4/105*ones(m-4, 1), -3)...
    + diag(4/105*ones(m-4, 1), 3) + diag(1/280*ones(m-5, 1), -4) - diag(1/280*ones(m-5, 1), 4);
D0(1:4,end-3:end) = [1/280, -4/105, 1/5, -4/5; 0, 1/280, -4/105, 1/5; 0, 0,1/280, -4/105; 0, 0, 0, 1/280];
D0(end-3:end,1:4)=fliplr(flipud(-D0(1:4,end-3:end)));
D0 = D0/h;
%}

%Creates second order derivative with periodic boundary conditions
%{
Dp = diag(-ones(m-1, 1), 0) + diag(ones(m-2, 1), 1);
Dm = diag(ones(m-1, 1), 0) + diag(-ones(m-2, 1), -1);
D2 = Dp*Dm;
D2(end,end-1) = 1; D2(end,end) = -2; D2(end,1) = 1; D2(1,end) = 1;
D2 = 2*D2/(h*h); 
%}

%8th order second derivative
D2 = -205/72*diag(ones(m-1,1),0) + 8/5*diag(ones(m-2,1),1) - 1/5*diag(ones(m-3,1),2) + 8/315*diag(ones(m-4,1),3) - 1/560*diag(ones(m-5,1), 4)...
    + 8/5*diag(ones(m-2,1),-1) - 1/5*diag(ones(m-3,1),-2) + 8/315*diag(ones(m-4,1),-3) - 1/560*diag(ones(m-5,1), -4);
D2(1:4, end-3:end) = [-1/560, 8/315, -1/5, 8/5; 0, -1/560, 8/315, -1/5; 0, 0, -1/560, 8/315;0, 0, 0, -1/560];
D2(end-3:end,1:4)=fliplr(flipud(D2(1:4,end-3:end)));
D2 = D2/(h*h);
%}

%P = p((end+1)/4);
P(:,1) = p;
%RHO(:,1) = rho;
%V(:,1) = v;
%MOM(:,1) = mom;

for i = 1:length(t)-1
    
    A1(:,i) = mom;
    B1(:,i) = zeros(m-1,1); 
    RHS1 = -D0*A1(:,i) + D2*B1(:,i);
    
    A2(:,i) = mom.*v + p;
    B2(:,i) = 4/3*mu*v;
    RHS2 = -D0*A2(:,i) + D2*B2(:,i);
    
    w1(:,i+1) = RK4(w1(:,i), RHS1, k);
    w2(:,i+1) = RK4(w2(:,i), RHS2, k);
    
    %Primitive variables
    rho = w1(:,end);
    mom = w2(:,end);
    
    v = mom./rho;
    p = rho.^(gamma);
    %P(i+1) = p((end+1)/4);
    P(:,i+1) = p;
    %RHO(:, i+1) = rho;
    %MOM(:, i+1) = mom;
    %V(:, i+1) = v;
    
  
end

%NavierStokesIsentropicPressure = P;
%{
f1 = figure('Name','Navier Stokes Model');
tiledlayout(2,1) 
nexttile
plot(t,NavierStokesIsentropicPressure(2,:))
title('left side'), xlabel('t'), ylabel('p')

nexttile
%}

P(end+1,:) = P(1,:);

plot(t,P((end)/4,:))
%title('right side'), xlabel('t'), ylabel('p')

%mesh(x(1:end-1),t,P')
%title('Pressure'), xlabel('x'), ylabel('t'), zlabel('P')

function [w] = RK4(v, RHS, k)
    
    k1 = RHS;
    k2 = RHS + k*k1/2;
    k3 = RHS + k*k2/2;
    k4 = RHS + k*k3;
    
    w = v + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end


