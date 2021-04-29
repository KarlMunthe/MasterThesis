clear

m=10; %  points
h=1/(m-1); % grid size.

k = 1e-4;        %time step
    
x = 0:h:1;    %individual space points on the grid
t = 0:k:k;    %individual time points on the grid

%constants
mu = 1.85e-4;
c_v = 717;
c_p = 1004;
gamma = c_p/c_v;

%Initial conditions and Primitive variables
%rho = (2e-04*(sin(2*pi*x(1:end-1)))+1.2)';
rho = 1.2*ones(m-1,1);
p = rho.^(gamma);
v = zeros(m-1,1);
mom = rho.*v;

w = [rho; mom];

I = [1,0;0,1];

%Creates first order derivative with periodic boundary conditions
D1 = diag(zeros(m-1, 1), 0) - diag(0.5*ones(m-2, 1), -1) + diag(+0.5*ones(m-2, 1), 1);
D1(end, 1) = 0.5; D1(1, end) = -0.5;
D1 = D1/h;
D1 = kron(I, D1);

%Creates second order derivative with periodic boundary conditions
Dp = diag(-ones(m-1, 1), 0) + diag(ones(m-2, 1), 1);
Dm = diag(ones(m-1, 1), 0) + diag(-ones(m-2, 1), -1);
D2 = Dp*Dm;
D2(end,end-1) = 1; D2(end,end) = -2; D2(end,1) = 1; D2(1,end) = 1;
D2 = 2*D2/(h*h); 
D2 = kron(I, D2);

P(:,1) = p;
RHO(:,1) = rho;
V(:,1) = v;
MOM(:,1) = mom;
%z = animatedline;
%axis([0 1, 0 0.0035])

for i = 1:length(t)-1
    
    A(:,i) = [mom; rho.*(v.^2) + p];
    B(:,i) = [zeros(m-1,1); 4/3*mu*v];
    
    k1 = -D1*(A(:,i)) + D2*(B(:,i));
    k2 = -D1*(A(:,i) + k*k1/2) + D2*(B(:,i) + k*k1/2);
    k3 = -D1*(A(:,i) + k*k2/2) + D2*(B(:,i) + k*k2/2);
    k4 = -D1*(A(:,i) + k*k3) + D2*(B(:,i) + k*k3);
    
    w(:,i+1) = w(:,i) + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
    %Primitive variables
    rho = w(1:end/2, end);
    mom = w(end/2+1:end, end);
    p = rho.^gamma;
    P(:,i+1) = p;
    RHO(:, i+1) = rho;
    MOM(:, i+1) = mom;
    
    %{
    clearpoints(z)
    addpoints(z, x, rho)
    drawnow limitrate
    %}
end
%drawnow

%{
mesh(x,t,u1')
hold on 
mesh(x,t,u2')
title('Isentropic Sound Attenuation'), xlabel('x'), ylabel('t'), zlabel('u')
legend('density [\rho]', 'momentum [\rho \cdot v]')
%}


tiledlayout(2,1) 
nexttile
mesh(x(1:end-1),t,RHO')
title('Density'), xlabel('x'), ylabel('t'), zlabel('\rho')

nexttile
mesh(x(1:end-1),t,MOM')
title('Momentum'), xlabel('x'), ylabel('t'), zlabel('\rho v')


%{
tiledlayout(2,1) 
nexttile
hold on
plot(x(1:end-1),RHO(:,1))
plot(x(1:end-1),RHO(:,end))
legend('1st', 'last')
title('Density'), xlabel('x')

nexttile
plot(x(1:end-1),V(:,1))
hold on
plot(x(1:end-1),V(:,end))
legend('1st', 'last')
title('Velocity'), xlabel('x')
%}

%{
plot(x(1:end-1),P(:,1))
hold on
plot(x(1:end-1),P(:,end))
legend('first', 'last')
%}