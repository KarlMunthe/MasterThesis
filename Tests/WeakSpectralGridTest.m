k = 1e-3; % time step

x_end = 1;
t_end = 1;
syms ex te

u = cos(2*pi*ex - 2*pi*te);

sol = matlabFunction(u);

m = 7;
    
t = 0:k:t_end;

h = x_end/m; 
x = h*(1:m);    %individual space points on the grid
CFL = k/h;

s = 2*pi/(2*x_end);


%initial condition
f = sol(x, 0)';

%D0 = FD0(m, h, 40);
%D0 = SpectralD0(m, s);

TD0 = testoperator(m, h);
D0 = FD0(m, h, 2);


w = [f, zeros(length(f), length(t)-1)];


for i = 1:length(t)-1

    w(:,i+1) = RK4(D0, w(:,i), k);
end

w = [w(end, :) ; w];

ex = [0, x];

error = w(:,end)-sol(ex, t_end)';
error = sqrt(h*error'*error)

   
%{
for i = 1:N-1
    p_list(i) = log(error_list(i+1)/error_list(i))/log(h_list(i+1)/h_list(i));
end
%}
%{
for i = 1:N-1
    p_list(i) = log2(error_list(i)/error_list(i+1));
end
%}
%p_list


x = linspace(0, x_end, m+1);
plot(x, w(:,end))
hold on
plot(x, sol(x,t_end))
legend('numerical', 'analytical')


f1 = figure;
time1 = 0:k:t_end;
mesh(x, time1, w'), xlabel('x'), ylabel('t'), zlabel('u')

%
%JUST FUNCTIONS BELOW
%

function [w] = RK4(D0, u, k)
    
    k1 = Advection(D0, u);
    k2 = Advection(D0, u + 0.5*k*k1);
    k3 = Advection(D0, u + 0.5*k*k2);
    k4 = Advection(D0, u + k*k3);
    
    w = u + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function [flux] = Advection(D0, u)
    
    A1 = -D0*u;
    
    %C1 = G(x(1:end-1)', t);

    flux = A1;

end
