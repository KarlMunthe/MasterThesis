clear

addpath('/Users/karlmunthe/Documents/UiB/UiB Master Oppgave/Code/gitMaster/Difference Operators')

%h = 1/(m-1); % grid size.
k = 1e-6; % time step

x_end = 1;
t_end = 0.01;
 

c_V = 1;     %heat capacity at constant volume
c_p = 2;     %heat capacity at constant pressure


%%%
%%%ADVECTION EQUATION WITH
%%%
%{

syms ex te

u = cos(2*pi*ex - 2*pi*te);

sol = matlabFunction(u);

N = ;
CFL_list = zeros(1,N);
error_list = zeros(1,N);
h_list = zeros(1,N);
p_list = zeros(1,N-1);
m = 10;

for i = 0:N-1
    
    [W, CFL, error, h] = PeriodicAdvectionTest(m, k, sol, x_end, t_end);
    
    Wi = W;
    CFL_list(i+1) = CFL;
    error_list(i+1) = error;
    h_list(i+1) = h;
    
    m = 2*m-1;
    
end

%{
for i = 1:N-1
    p_list(i) = log(error_list(i+1)/error_list(i))/log(h_list(i+1)/h_list(i));
end
%}
for i = 1:N-1
    p_list(i) = log2(error_list(i)/error_list(i+1));
end
p_list
error_list


x = linspace(0, x_end, (m+1)/2 + 1);
plot(x, W(:,end))
hold on
plot(x, sol(x,t_end))
legend('numerical', 'analytical')

%{
f1 = figure;
time1 = 0:k:t_end;
x = linspace(0, x_end, (m+1)/2 + 1);
mesh(x, time1, W'), xlabel('x'), ylabel('t'), zlabel('u')
%}
%}

%%%
%%%HEAT EQUATION
%%%
%{
%MANUFACTURED SOLUTION
syms ex te
rho = sin(2*pi*(ex + te));
u = cos(2*pi*(ex + te));

A1 = diff(rho, te);
B1 = diff(rho, ex, 2);
R1 = A1 - B1;

sol = matlabFunction(rho);

rho = matlabFunction(rho);
u = matlabFunction(u);
R1 = matlabFunction(R1);

N = 3;
CFL_list = zeros(1,N);
error_list = zeros(1,N);
h_list = zeros(1,N);
p_list = zeros(1,N-1);
m = 12;

for i = 0:N-1
    
    [W, CFL, error, h] = PeriodicHeatTest(x_end, m, k, sol, t_end, R1);
    
    Wi = W;
    CFL_list(i+1) = CFL;
    error_list(i+1) = error;
    h_list(i+1) = h;
    
    m = m + 10;
    
end

for i = 1:N-1
    p_list(i) = log(error_list(i+1)/error_list(i))/log(h_list(i+1)/h_list(i));
end
p_list
error_list

f1 = figure;
time1 = 0:k:t_end;
x = linspace(0, x_end, m-10+1);
mesh(x, time1, W'), xlabel('x'), ylabel('t'), zlabel('u')

%}

%%%
%%%BURGERS EQUATION WITH DIFFUSION
%%%
%{
%MANUFACTURED SOLUTION
syms ex te
rho = sin(2*pi*(ex + te)) + 2;
u = cos(2*pi*(ex + te)) + 2;
epsilon = 1e-6;

A1 = diff(u, te);
B1 = diff(u^2/2, ex);
C1 = epsilon*diff(u, ex, 2);
source_term = A1 + B1 - C1;

sol = matlabFunction(u);
source_term = matlabFunction(source_term);

%Normal FDM schemes
N = 3;
CFL_list = zeros(1,N);
error_list = zeros(1,N);
h_list = zeros(1,N);
p_list = zeros(1,N-1);
m = 11;

for i = 0:N-1
    
    [W, CFL, error, h] = SpectralBurgersDiffusion(m, epsilon, source_term, sol, t_end, x_end, k);
    
    Wi = W;
    CFL_list(i+1) = CFL;
    error_list(i+1) = error;
    h_list(i+1) = h;
    
    m = m + 10;
    
end

for i = 1:N-1
    p_list(i) = log(error_list(i+1)/error_list(i))/log(h_list(i+1)/h_list(i));
end
p_list
error_list

f1 = figure;
time1 = 0:k:t_end;
x = linspace(0, x_end, m-10+1);
mesh(x, time1, W'), xlabel('x'), ylabel('t'), zlabel('u')

%}

%%%
%%%VARIABLE COEFFICIENT ADVECTION EQUATION (NS CONTINUITY EQUATION)
%%%
%{
%MANUFACTURED SOLUTION
syms ex te
rho = sin(2*pi*(ex + te));
u = cos(2*pi*(ex + te));

A1 = diff(rho, te);
B1 = diff(rho*u, ex);
R1 = A1 + B1;

sol = matlabFunction(rho);

rho = matlabFunction(rho);
u = matlabFunction(u);
R1 = matlabFunction(R1);

N = 3;
CFL_list = zeros(1,N);
error_list = zeros(1,N);
h_list = zeros(1,N);
p_list = zeros(1,N-1);
m = 10;

for i = 0:N-1
    
    [W, CFL, error, h] = PeriodicAdvectionTest(m, k, sol, x_end, t_end);
    
    Wi = W;
    CFL_list(i+1) = CFL;
    error_list(i+1) = error;
    h_list(i+1) = h;
    
    m = 2*m-2^i;
    
end

for i = 1:N-1
    p_list(i) = log(error_list(i+1)/error_list(i))/log(h_list(i+1)/h_list(i));
end
p_list
error_list
%}

%%%
%%%NS MOMENTUM EQUATION EQUATION
%%%
%{
%MANUFACTURED SOLUTION
syms ex te
rho = sin(2*pi*ex + 2*pi*te) + 2;
mom = (sin(2*pi*ex + 2*pi*te) + 2)^2;
p = rho^3;

A1 = diff(mom, te);
B1 = diff(mom^2/rho + p, ex);
C1 = diff(4/3*mom/rho, ex, 2);
D1 = A1 + B1 - C1;

sol = matlabFunction(mom);

rho = matlabFunction(rho);
p = matlabFunction(p);
D1 = matlabFunction(D1);

%{
xxx = 0:1/(m-1):x_end;
plot(xxx, D1(xxx, 0))
%}

N = 6;
CFL_list = zeros(1,N);
error_list = zeros(1,N);
h_list = zeros(1,N);
p_list = zeros(1,N-1);
m = 16;

for i = 0:N-1
    
    [W, CFL, error, h] = NSMomentumEquation(m, t_end, x_end, k, p, rho, D1, sol);
    
    CFL_list(i+1) = CFL;
    error_list(i+1) = error;
    h_list(i+1) = h;
    
    m = 2*m-2^i;
    
end

for i = 1:N-1
    p_list(i) = log(error_list(i+1)/error_list(i))/log(h_list(i+1)/h_list(i));
end
p_list

%{
f1 = figure;
time1 = 0:k:t_end;
x = linspace(0, x_end, m);
mesh(x, time1, W1'), xlabel('x'), ylabel('t'), zlabel('u')

f2 = figure;
plot(x, W1(:,end))
hold on
plot(x, sol(x,t_end))
legend('simulation', 'analytical')
%}
%}

%%%
%%%COUPLED ISENTROPIC NS EQUATIONS
%%%
%{
%MANUFACTURED SOLUTION
syms ex te
rho = sin(ex + te) + 2;
mom = (sin(ex + te) + 2)^2;
p = rho^3;

A1 = diff(rho, te);
B1 = diff(mom, ex);
C1 = 0;
R1 = A1 + B1 - C1;

A2 = diff(mom, te);
B2 = diff(mom^2/rho + p, ex);
C2 = diff(4/3*mom/rho, ex, 2);
R2 = A2 + B2 - C2;

rho_sol = matlabFunction(rho);
mom_sol = matlabFunction(mom);
p = matlabFunction(p);

R1 = matlabFunction(R1);
R2 = matlabFunction(R2);

%{
xxx = 0:1/(m-1):x_end;
plot(xxx, D1(xxx, 0))
%}

[W1, CFL1, error1, h1] = IsentropicNSTest(m, t_end, x_end, k, rho_sol, R1, R2, mom_sol);
[W2, CFL2, error2, h2] = IsentropicNSTest(2*m-1, t_end, x_end, k, rho_sol, R1, R2, mom_sol);

p = log(error2/error1)/log(h2/h1)
%{
f1 = figure;
time1 = 0:k:t_end;
x = linspace(0, x_end, m);
mesh(x, time1, W1'), xlabel('x'), ylabel('t'), zlabel('u')

f2 = figure;
plot(x, W1(:,end))
hold on
plot(x, rho_sol(x,t_end))
legend('simulation', 'analytical')
%}
%}

%%%
%%%COUPLED ISENTROPIC NSS EQUATIONS
%%%
%{
%MANUFACTURED SOLUTION
gamma = c_p/c_V;
syms ex te
rho = sin(ex + te) + 2;
mom = (sin(ex + te) + 2)^2;

nu = 1/rho;
p = rho^gamma;

A1 = diff(rho, te);
B1 = diff(mom, ex);
C1 = diff(nu*diff(rho, ex), ex);
R1 = A1 + B1 - C1;

A2 = diff(mom, te);
B2 = diff(mom^2/rho + p, ex);
C2 = diff(nu*diff(mom, ex), ex);
R2 = A2 + B2 - C2;

rho_sol = matlabFunction(rho);
mom_sol = matlabFunction(mom);
p = matlabFunction(p);

R1 = matlabFunction(R1);
R2 = matlabFunction(R2);

%{
xxx = 0:1/(m-1):x_end;
plot(xxx, D1(xxx, 0))
%}

W = [zeros(1, 6)];
CFL = [zeros(1, 6)];
error = [zeros(1, 6)];
h = [zeros(1, 6)];

%{
for i = 1:6;
    [W(i), CFL(i), error(i), h(i)] = IsentropicNSSTest((2^i)*m, t_end, x_end, k, nu, rho_sol, R1, R2, mom_sol);
end
%}

[W1, CFL1, error1, h1] = IsentropicNSSTest(m, t_end, x_end, k, c_p, c_V, nu, rho_sol, R1, R2, mom_sol);
[W2, CFL2, error2, h2] = IsentropicNSSTest(2*m-1, t_end, x_end, k, c_p, c_V, nu, rho_sol, R1, R2, mom_sol);

p = log(error2/error1)/log(h2/h1)
%{
f1 = figure;
time1 = 0:k:t_end;
x = linspace(0, x_end, m);
mesh(x, time1, W1'), xlabel('x'), ylabel('t'), zlabel('u')

f2 = figure;
plot(x, W1(:,end))
hold on
plot(x, rho_sol(x,t_end))
legend('simulation', 'analytical')
%}
%}

%%%
%%%FULL COUPLED COMPRESSIBLE NS EQUATIONS
%%%

%MANUFACTURED SOLUTION
R = c_p - c_V;
gamma = c_p/c_V;
kappa = 1;
syms ex te
rho = cos(2*pi*(ex + te)) + sin(2*pi*(ex + te)) + 2;
mom = cos(2*pi*(ex + te)) - sin(2*pi*(ex + te));
E = -cos(2*pi*(ex + te)) + sin(2*pi*(ex + te));
p = R/c_V*(E - mom.^2./(2*rho));
T = p/(rho*R);

epsilon = 1;

A1 = diff(rho, te);
B1 = diff(mom, ex);
C1 = 0;
R1 = A1 + B1 - C1;

A2 = diff(mom, te);
B2 = diff(mom^2/rho + p, ex);
C2 = epsilon*diff(4/3*mom/rho, ex, 2);
R2 = A2 + B2 - C2;

A3 = diff(E, te);
B3 = diff(E*mom/rho + p*mom/rho, ex);
C3 = epsilon*4/3*diff(mom/rho*diff(mom/rho, ex), ex) + diff(kappa*diff(T, ex), ex);
R3 = A3 + B3 - C3;

rho_sol = matlabFunction(rho);
mom_sol = matlabFunction(mom);
E_sol = matlabFunction(E);

R1 = matlabFunction(R1);
R2 = matlabFunction(R2);
R3 = matlabFunction(R3);

%{
xxx = 0:1/(m-1):x_end;
plot(xxx, D1(xxx, 0))
%}

N = 5;
CFL_list = zeros(1,N);
RHO_error_list = zeros(1,N);
MOM_error_list = zeros(1,N);
E_error_list = zeros(1,N);
h_list = zeros(1,N);

RHO_p_list = zeros(1,N-1);
MOM_p_list = zeros(1,N-1);
E_p_list = zeros(1,N-1);
m = 11;

for i = 0:N-1
    
    [CFL, RHO_error,  MOM_error , E_error, h] = SpectralCompressibleNSTest(m, epsilon, t_end, x_end, k, c_p, c_V, kappa, rho_sol, mom_sol, E_sol, R1, R2, R3);
    
    CFL_list(i+1) = CFL;
    RHO_error_list(i+1) = RHO_error;
    MOM_error_list(i+1) = MOM_error;
    E_error_list(i+1) = E_error;
    h_list(i+1) = h;
    
    m = 2*m - 1;
    
end


for i = 1:N-1
    RHO_p_list(i) = log(RHO_error_list(i+1)/RHO_error_list(i))/log(h_list(i+1)/h_list(i));
    MOM_p_list(i) = log(MOM_error_list(i+1)/MOM_error_list(i))/log(h_list(i+1)/h_list(i));
    E_p_list(i) = log(E_error_list(i+1)/E_error_list(i))/log(h_list(i+1)/h_list(i));

end

%RHO_p_list
%MOM_p_list
%E_p_list

NS_RHO_error_list = RHO_error_list
NS_MOM_error_list = MOM_error_list
NS_E_error_list = E_error_list
%CFL_list
%%
%{
f1 = figure;
time1 = 0:k:t_end;
x = linspace(0, x_end, (m+1)/2+1);
mesh(x, time1, W'), xlabel('x'), ylabel('t'), zlabel('u')
%}

%{
f2 = figure;
plot(x, W1(:,end))
hold on
plot(x, rho_sol(x,t_end))
legend('simulation', 'analytical')
%}
%}

%%%
%%%FULL COUPLED COMPRESSIBLE NSS EQUATIONS
%%%

%MANUFACTURED SOLUTION
R = c_p - c_V;
gamma = c_p/c_V;
syms ex te
rho = cos(2*pi*(ex + te)) + sin(2*pi*(ex + te)) + 2;
mom = cos(2*pi*(ex + te)) - sin(2*pi*(ex + te));
E = -cos(2*pi*(ex + te)) + sin(2*pi*(ex + te));
p = R/c_V*(E-mom^2/(2*rho));

nu = 1/rho;

A1 = diff(rho, te);
B1 = diff(mom, ex);
C1 = diff(nu*diff(rho, ex), ex);
R1 = A1 + B1 - C1;

A2 = diff(mom, te);
B2 = diff(mom^2/rho + p, ex);
C2 = diff(nu*diff(mom, ex), ex);
R2 = A2 + B2 - C2;

A3 = diff(E, te);
B3 = diff(E*(mom/rho) + p*(mom/rho), ex);
C3 = diff(nu*diff(E, ex), ex);
R3 = A3 + B3 - C3;

rho_sol = matlabFunction(rho);
mom_sol = matlabFunction(mom);
E_sol = matlabFunction(E);

R1 = matlabFunction(R1);
R2 = matlabFunction(R2);
R3 = matlabFunction(R3);

N = 5;
CFL_list = zeros(1,N);
RHO_error_list = zeros(1,N);
MOM_error_list = zeros(1,N);
E_error_list = zeros(1,N);
h_list = zeros(1,N);
p_list = zeros(1,N-1);
m = 11;

for i = 0:N-1
    
    [CFL, RHO_error,  MOM_error , E_error, h] = SpectralCompressibleNSSTest(m, t_end, x_end, k, c_p, c_V, rho_sol, mom_sol, E_sol, R1, R2, R3);
    
    CFL_list(i+1) = CFL;
    RHO_error_list(i+1) = RHO_error;
    MOM_error_list(i+1) = MOM_error;
    E_error_list(i+1) = E_error;
    h_list(i+1) = h;
    
    m = 2*m-1;
    
end

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

NSS_RHO_error_list = RHO_error_list
NSS_MOM_error_list = MOM_error_list
NSS_E_error_list = E_error_list
%CFL_list

%{
f1 = figure;
time1 = 0:k:t_end;
x = linspace(0, x_end, (m+1)/2 + 1);
mesh(x, time1, W'), xlabel('x'), ylabel('t'), zlabel('u')
%}
%}

%%%
%%%FULL COUPLED COMPRESSBLE NSS EQUATIONS USING COMPACT PADE SCHEMES
%%%
%{
%MANUFACTURED SOLUTION
R = c_p - c_V;
gamma = c_p/c_V;
syms ex te
rho = sin(2*pi*ex + 2*pi*te) + 2;
mom = (sin(2*pi*ex + 2*pi*te) + 2)^2;
E = cos(2*pi*ex + 2*pi*te) + 2; 
p = R/c_V*(E-mom^2/(2*rho));

nu = 1/rho;

A1 = diff(rho, te);
B1 = diff(mom, ex);
C1 = diff(nu*diff(rho, ex), ex);
R1 = A1 + B1 - C1;

A2 = diff(mom, te);
B2 = diff(mom^2/rho + p, ex);
C2 = diff(nu*diff(mom, ex), ex);
R2 = A2 + B2 - C2;

A3 = diff(E, te);
B3 = diff(E*(mom/rho) + p*(mom/rho), ex);
C3 = diff(nu*diff(E, ex), ex);
R3 = A3 + B3 - C3;

rho_sol = matlabFunction(rho);
mom_sol = matlabFunction(mom);
E_sol = matlabFunction(E);

R1 = matlabFunction(R1);
R2 = matlabFunction(R2);
R3 = matlabFunction(R3);

N = 6;
CFL_list = zeros(1,N);
error_list = zeros(1,N);
h_list = zeros(1,N);
p_list = zeros(1,N-1);
m = 4;

for i = 0:N-1
    
    [W, CFL, error, h] = CompressibleNSSTestWithPade(m, t_end, x_end, k, c_p, c_V, rho_sol, mom_sol, E_sol, R1, R2, R3);
    
    CFL_list(i+1) = CFL;
    error_list(i+1) = error;
    h_list(i+1) = h;
    
    m = 2*m-2^i;
    
end

for i = 1:N-1
    p_list(i) = log(error_list(i+1)/error_list(i))/log(h_list(i+1)/h_list(i));
end
p_list
error_list
%}