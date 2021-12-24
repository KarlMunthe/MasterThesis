clear

%Dimensionless temperature and pressure
T0 = 1;
p0 = 1;

%fluid properties for argon
%{
mu = 22.61e-06;
kappa = 0.0178;
c_V = 313;
c_p = 520;
R = c_p - c_V;
gamma = c_p/c_V;
%}

%fluid properties for a dimensionless fluid
mu = 1;
c_V = 1;
c_p = 2;
R = c_p - c_V;
gamma = c_p/c_V;

%set up time
dt = 1e-5;
tend = 10;
time = dt:dt:tend;

%set up grid in x direction
Nx = 33;
xstart = 9;
xend = 11;
x = linspace(xstart, xend, Nx);

%set up grid in y direction
Ny = 128;
ystart = 0;
yend = 200;
y = linspace(ystart, yend, Ny);

%transform grid in ydirection
q = 30;
ytrans = ((ystart*exp(yend/q) - yend*exp(ystart/q) + (yend-ystart)*exp(y/q))/(exp(yend/q)-exp(ystart/q)))';
hy = ytrans(2:end)-ytrans(1:end-1);

%calculate volume size in y-direction
yside = 0.5*(ytrans(3:end) - ytrans(1:end-2));
nyside = 0.5*(ytrans(2) - ytrans(1));
syside = 0.5*(ytrans(end) - ytrans(end - 1));
yside = [nyside; yside; syside];

%calculate volume size in x-direction
hx = x(2:end)-x(1:end-1);
xside = 0.5*(x(3:end) - x(1:end-2));
wxside = 0.5*(x(2) - x(1));
exside = 0.5*(x(end) - x(end-1));
xside = [wxside, xside, exside];

%calculate size of each volume
vol = (yside*xside);
%display size of smallest volume for CFL information
min(min(vol))

%inlet velocity and background density
U = ones(length(y), 1);
rho0 = p0./(R*T0)*ones(size(vol));

%calculate initical conditions for velocity in x and y direction, u annd v.
nu = mu./rho0;
[u, v] = InitialCondition(x, ytrans, U, nu);

m0 = rho0.*u;
n0 = rho0.*v;
e0 = c_V*T0;
E0 = e0.*rho0 + 0.5*(m0.^2 + n0.^2)./rho0;

%rename variables for the for loop
rho = rho0;
m = m0;
n = n0;
E = E0;
p = (gamma-1)*(E - 0.5*(m0.^2 + n0.^2)./rho);

vars = zeros(size(vol, 1), size(vol, 2), 4);
vars(:, :, 1) = rho;
vars(:, :, 2) = m;
vars(:, :, 3) = n;
vars(:, :, 4) = E;

svars = vars(1, :, 2:3);
    
t=0;
l2error = Inf;
while l2error > 1e-7
    t = t+1;
    
    vars = RK4(vars, dt, vol, hx, hy, xside, yside, gamma, mu, svars);
   
    %save density values every 499 time steps 
    if mod(t, 500) == 499
        testvars = vars(:, :, 1);
    end

    if mod(t, 500) == 0
        
        %calculate error between current time step and the last. 
        l2error = sqrt(sum(sum(vol.*(vars(:, :, 1) - testvars).^2)))
        
        %plot every 500 time steps to see if is stable or not. lets you end
        %simulation prematurely if it isn't stable.
        mesh(x, ytrans, vars(:, :, 1))
        title('density')
        drawnow
    end
    
end

%plot initial density
mesh(x, ytrans, rho0)
xlabel('x')
ylabel('y')

%plot current density
f2 = figure;
mesh(x, ytrans, vars(:, :, 1))
xlabel('x')
ylabel('y')
    