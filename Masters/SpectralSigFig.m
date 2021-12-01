%Oxygen values
mu = 20.64e-06;
kappa = 26.58e-03;
c_p = 915;
c_V = 659;

%Gas constant and heat capacity
R = c_p - c_V;
gamma = c_p/c_V;

%Standard temperature and pressure (STP)
T_0 = 293.15;
p_0 = 1e+05;

rho_0 = p_0/(R*T_0);  %background density at STP
rho = @(x) rho_0 + 1e-02*(sin(2*pi*x));   %reynolds decomposed pressure

nu_0 = mu/rho_0;
c = sqrt(gamma*R*T_0);    %speed of sound
lambda = 1;        %wave length
omega = 2*pi*c;
eta = 0;

m = 11;

k = 1e-06; % time step

x_end = 0.001;
t_end = 0.13; %0.63; %0.33 % 1.3; %c ~ 1000

%Eulerian viscosity coefficient alpha can equal 1 or 4/3
alpha = 1;

h = x_end/m;    %individual space points on the grid
x = h*(1:m);

s = 2*pi;


Q1 = SpectralD0(m, s);
Q2 = SpectralD2(m, s);

N = 1;

RHO_p_list = zeros(1, 1);
MOM_p_list = zeros(1, 1);
ENE_p_list = zeros(1, 1);

RHO_list = zeros(1, N);
MOM_list = zeros(1, N);
ENE_list = zeros(1, N);

h_list = zeros(1, N);
CFL_list = zeros(1, N);

M = 89;
%X = linspace(0, x_end, M);    %individual space points on the grid
%H= X(2) - X(1);


H = 2*pi/M;    %individual space points on the grid
X = H*(1:M);

CFL = k/(H^2);

Q1 = SpectralD0(M, H);
Q2 = SpectralD2(M, H);

%Q1 = PeriodicD0(M, H);
%Q2 = PeriodicD2(M, H);

[RHO_sol, MOM_sol, ENE_sol] = SpectralNS(M, X, t_end, k, c_p, c_V, kappa, mu, rho, p_0, Q1, Q2);

for i = 1:N
    
    h = 2*pi/m;    %individual space points on the grid
    x = h*(1:m);
    
    %x = linspace(0, x_end, m);
    %h = x(2)-x(1)
    %CFL = k/(h^2);

    Q1 = SpectralD0(m, h);
    Q2 = SpectralD2(m, h);
    
    CFL = k/(h^2);
    CFL_list(i) = CFL;

    [RHO, MOM, ENE] = SpectralNS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, Q1, Q2);
    rho_error = RHO - RHO_sol(1:2^(N-i+2):end);
    rho_error = sqrt(h*(rho_error'*rho_error));
    RHO_list(i) = rho_error;
    
    mom_error = MOM - MOM_sol(1:2^(N-i+2):end);
    mom_error = sqrt(h*(mom_error'*mom_error));
    MOM_list(i) = mom_error;
    
    ene_error = ENE - ENE_sol(1:2^(N-i+2):end);
    ene_error = sqrt(h*(ene_error'*ene_error));
    ENE_list(i) = ene_error;
    
    m = 2*m;
    
end