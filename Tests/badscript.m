clear

m = 128;

%h = 1/(m-1); % grid size.
k = 1e-4; % time step

x_end = 2*pi;
t_end = 0.5*pi;

syms ex te;
%sol = sin(ex - te);
sol = cos(ex + te);
A = diff(sol, te);
B = diff(sol^2/2, ex);
C = diff(sol, ex, 2);
D = A + B - C;

G = matlabFunction(D);
sol = matlabFunction(sol);

%[W1, error1, h1, CFL1] = HeatTest(m, k, sol, G, t_end);
%[W2, error2, h2, CFL2] = HeatTest(2*m-1, k, sol, G, t_end);

%[W1, error1, h1] = PeriodicAdvectionTest(m, k, sol, t_end);
%[W2, error2, h2] = PeriodicAdvectionTest(2*m - 1, k, sol, t_end);

[W1, error1, h1, CFL1] = PeriodicBurgersDiffusion(m, G, sol, t_end, k);
[W2, error2, h2, CFL2] = PeriodicBurgersDiffusion(2*m-1, G, sol, t_end, k);

%w1 = W1(:,end);
%w2 = W2(:,end);

%W2_vals = reshape(w2(1:end-1,:), 2, []);
%W2_vals = [W2_vals(1,:), w2(end)];


p1 = log(error2/error1)/log(h2/h1)

%p2 = log2(diff1/diff2)

%p1_2 = log2(err1/err2)
f1 = figure;
time1 = linspace(0, t_end, size(W1,2));
x = linspace(0, 2*pi, m);
mesh(x,time1,W1'), xlabel('x'), ylabel('t'), zlabel('u')

f2 = figure;
plot(x, W1(:,end))
hold on
plot(x, sol(x,t_end))
legend('simulation', 'analytical')
%}

