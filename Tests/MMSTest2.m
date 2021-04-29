clear

m = 32;         %number of grid points

h = 1/(m-1);        %space step
k = 1e-6;        %time step

t_end = 128*k;

x = linspace(0,1,m);    %individual space points on the grid

alpha = 1;

%CFL = k/h

syms ex te;
sol = cos(2*pi*ex + 2*pi*te);
A = diff(sol, te);
B = sol*diff(sol, ex);
C = diff(sol, ex, 2);
D = A + B - alpha*C;

G = matlabFunction(D);
sol = matlabFunction(sol);

%Source term
%{
G = @(x, t) alpha*4*pi^2*cos(2*pi*x+2*pi*t) - pi*sin(4*pi*(x+t)) ...
    - 2*pi*sin(2*pi*x+2*pi*t);
%}

%initial condition
f = sol(x(1:end-1), 0)';

[~, D0, ~, ~] = PeriodicD0(m, h);
[~, ~, ~, D2] = PeriodicD2(m, h);

w = [f, zeros(length(f), round(t_end/k))];

%Solve Pade Schemes with RK4
t = k;
i = 1;
while t < t_end
       
    k1 = -D0*(0.5*(w(:,i).^2)) + alpha*D2*w(:,i) + G(x(1:end-1)', t);
    k2 = -D0*(0.5*(w(:,i) + k*k1/2).^2) + alpha*D2*(w(:,i) + k*k1/2)...
        + G(x(1:end-1)' + k*k1/2, t + k/2);
    k3 = -D0*(0.5*(w(:,i) + k*k2/2).^2) + alpha*D2*(w(:,i) + k*k2/2)...
        + G(x(1:end-1)' + k*k2/2, t + k/2);
    k4 = -D0*(0.5*(w(:,i) + k*k3).^2) + alpha*D2*(w(:,i) + k*k3)...
        + G(x(1:end-1)' + k*k3, t + k);
    
    w(:,i+1) = w(:,i) + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
    t = t + k;
    i = i+1;
end
%}

w(end + 1, :) = w(1, :);
w = w(:,1:end); 

time = linspace(0, t_end, length(w));

%mesh(x,time,w'), xlabel('x'), ylabel('t'), zlabel('u')

w = w(:,1:end-1);
plot(x, w(:,end))
hold on
plot(x, cos(2*pi*x + 2*pi*t_end))
legend('simulation', 'analytical')
error = norm(w(:,end)-sol(x, t_end)')

%}


