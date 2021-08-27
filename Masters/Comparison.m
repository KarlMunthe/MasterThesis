clear

addpath('/Users/karlmunthe/Documents/UiB/UiB Master Oppgave/Code/gitMaster/Difference Operators')

%Gas properties
%Hydrogen values
%{
mu = 8.90e-06;
kappa = 0.1819;
c_p = 14.32e+03;
c_V = 10.16e+03;
%}

%Oxygen values

mu = 20.64e-06;
kappa = 26.58e-03;
c_p = 915;
c_V = 659;
%}

%Helium values
%{
mu = 20.64e-06;   %coefficient of viscosity
kappa = 0.1357; %thermal conductivity 
c_V = 3157;     %heat capacity at constant volume
c_p = 5240;     %heat capacity at constant pressure
%}


%Gas constant and heat capacity
R = c_p - c_V;
gamma = c_p/c_V;

%Standard temperature and pressure (STP)
T_0 = 293.15;
p_0 = 1e+05;

rho_0 = p_0/(R*T_0);  %background density at STP
rho = @(x) rho_0 + 1e-02*(sin(2*pi*x));   %reynolds decomposed pressure
%u = @(x) 1e-01*(sin(2*pi*x));

nu_0 = mu/rho_0;
c = sqrt(gamma*R*T_0);    %speed of sound
lambda = 1;        %wave length
omega = 2*pi*c;
eta = 0;

%Coefficient of absorption as a function of time

PropCoeffNS = omega^2/(2*rho_0*c^2)*((4/3*mu + eta) + (kappa/c_V)*(1/c_V-1/c_p)); %page 301 in Landou & Lipschitz
PropCoeffNSS = omega^2/(2*rho_0*c^2)*(mu + c_V*mu*(gamma-1)*(T_0/c^2 + 1/c_p));
%}

%Coefficient of absorption as a function of space
%{
PropCoeffNS = omega^2/(2*rho_0*c^3)*((4/3*mu + eta) + (kappa/c_V)*(1/c_V-1/c_p)); %page 301 in Landou & Lipschitz
PropCoeffNSS = omega^2/(2*rho_0*c^3)*(mu + c_V*mu*(gamma-1)*(T_0/c^2 + 1/c_p));
%}

article = PropCoeffNS*lambda;
article2 = PropCoeffNSS*lambda;

m = 11;

k = 1e-06; % time step

x_end = 1;
t_end = 0.13; %0.63; %0.33 % 1.3; %c ~ 1000

%x = 0:h:x_end;    %individual space points on the grid

%Eulerian viscosity coefficient alpha can equal 1 or 4/3
alpha = 1;
alpha2 = 4/3;

%Creates Pade scheme
%{
Dp = diag(-ones(m-1, 1), 0) + diag(ones(m-2, 1), 1);
Dp(end, 1) = 1;
Dp = Dp/h;
Dm = diag(ones(m-1, 1), 0) + diag(-ones(m-2, 1), -1);
Dm(1, end) = -1;
Dm = Dm/h;
I = diag(ones(m-1,1));

%REMEMBER TO USE SECOND ORDER, O(2), DIFFERENCE OPERATORS WHEN USING PADE
%SCHEMES, THE FIRST SIMPLESt OF THE ONES IN THE PERIODICD0/2 FUNCTIONS
PD1 = (I + h^2*Dp*Dm/6 - h^4*(Dp*Dm)^2/30 + h^6*(Dp*Dm)^3/140);
PD2 = (I + h^2*Dp*Dm/12 - h^4*(Dp*Dm)^2/90 + h^6*(Dp*Dm)^3/560);
%}

%Pade schemes
%{
%[NS_Work] = Pade_NS(m, x, t_end, k, mu, c_v, c_p, R, K, D0, D2, PD1, PD2, rho);

%[NSS_Work] = Pade_NSS(m, x, t_end, k, mu, alpha, c_v, c_p, R, D0, PD1, PD2, rho);
%}

%Normal difference schemes
%{
x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);

[~, ~, ~, Q1] = PeriodicD0(m, h);
[~, ~, ~, Q2] = PeriodicD2(m, h);

[NS_Work] = NS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, Q1, Q2);
%[NS_Work] = NS_feil(m, x, t_end, k, mu, c_V, c_p, R, kappa, Q1, Q2, rho);


[NSS_Work] = NSS(m, x, t_end, k, c_p, c_V, mu, rho, Q1);
%[NSS_Work] = E2(m, x, t_end, k, mu, alpha, c_V, c_p, R, Q1, rho);
%}

%Spectral difference schemes

h = x_end/m;    %individual space points on the grid
x = h*(1:m);

s = 2*pi;

%[Q1, ~,~,~] = PeriodicD0(m,h)

Q1 = SpectralD0(m, s);
Q2 = SpectralD2(m, s);

%[NSS_Work] = SpectralNSS(m, x, t_end, k, c_p, c_V, mu, rho_0, T_0, u, Q1);
[NSS_Work] = SpectralNSStrash(m, x, t_end, k, c_p, c_V, mu, rho, p_0, Q1);

[NS_Work] = SpectralNS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, p_0, Q1, Q2);

x = [0, x];

%Conservative check
%{
x = [0, x];
%{
NS_Rho  = NS_Rho(1:end-1,:);
NS_Mom  = NS_Mom(1:end-1,:);
NS_Ene  = NS_Ene(1:end-1,:);
%}

NS_Rho2 = NS_Rho(:,1:round(end/2));
NS_Mom2 = NS_Mom(:,1:round(end/2));
NS_Ene2 = NS_Ene(:,1:round(end/2));

NS_int_rho1 = trapz(x, NS_Rho(:,end));
NS_int_mom1 = trapz(x, NS_Mom(:,end));
NS_int_ene1 = trapz(x, NS_Ene(:,end));

NS_int_rho2 = trapz(x, NS_Rho2(:,end));
NS_int_mom2 = trapz(x, NS_Mom2(:,end));
NS_int_ene2 = trapz(x, NS_Ene2(:,end));

NS_diff_int_rho1 = abs(NS_int_rho1 - NS_int_rho2);
NS_diff_int_mom1 = abs(NS_int_mom1 - NS_int_mom2);
NS_diff_int_ene1 = abs(NS_int_ene1 - NS_int_ene2);

%{
NSS_Rho  = NSS_Rho(1:end-1,:);
NSS_Mom  = NSS_Mom(1:end-1,:);
NSS_Ene  = NSS_Ene(1:end-1,:);
%}

NSS_Rho2 = NSS_Rho(:,1:round(end/2));
NSS_Mom2 = NSS_Mom(:,1:round(end/2));
NSS_Ene2 = NSS_Ene(:,1:round(end/2));

NSS_int_rho1 = trapz(x, NSS_Rho(:,end));
NSS_int_mom1 = trapz(x, NSS_Mom(:,end));
NSS_int_ene1 = trapz(x, NSS_Ene(:,end));

NSS_int_rho2 = trapz(x, NSS_Rho2(:,end));
NSS_int_mom2 = trapz(x, NSS_Mom2(:,end));
NSS_int_ene2 = trapz(x, NSS_Ene2(:,end));

NSS_diff_int_rho1 = abs(NSS_int_rho1 - NSS_int_rho2);
NSS_diff_int_mom1 = abs(NSS_int_mom1 - NSS_int_mom2);
NSS_diff_int_ene1 = abs(NSS_int_ene1 - NSS_int_ene2);
%}

%Error check
%{
NSS_Work = 0.5*NSS_Mom.^2./NSS_Rho;
NS_Work = 0.5*NS_Mom.^2./NS_Rho;

N = 4;
h1 = h;
NSS_Rho_list = [NSS_Rho(1:end-1,end), zeros(m,N-1)];
NSS_Mom_list = [NSS_Mom(1:end-1,end), zeros(m,N-1)];
NSS_Ene_list = [NSS_Ene(1:end-1,end), zeros(m,N-1)];
NSS_Work_list = [NSS_Work(1:end-1,end), zeros(m, N-1)];

NS_Rho_list = [NS_Rho(1:end-1,end), zeros(m,N-1)];
NS_Mom_list = [NS_Mom(1:end-1,end), zeros(m,N-1)];
NS_Ene_list = [NS_Ene(1:end-1,end), zeros(m,N-1)];
NS_Work_list = [NS_Work(1:end-1,end), zeros(m, N-1)];

for i = 1:N-1
    
m = 2*m;

h = x_end/m;    %individual space points on the grid
x = h*(1:m);

Q1 = SpectralD0(m, s);
Q2 = SpectralD2(m, s);
    
[NSS_Rho, NSS_Mom, NSS_Ene] = SpectralNSStrash(m, x, t_end, k, c_p, c_V, mu, rho, p_0, Q1);
NSS_Rho_list(:, i+1) = NSS_Rho(1:2^i:end-1, end);%; NSS_Rho(end,end)];
NSS_Mom_list(:, i+1) = NSS_Mom(1:2^i:end-1, end);%; NSS_Mom(end,end)];
NSS_Ene_list(:, i+1) = NSS_Ene(1:2^i:end-1, end);%; NSS_Ene(end,end)];
NSS_Work_list(:, i+1) = 0.5*NSS_Mom(1:2^i:end-1, end).^2./NSS_Rho(1:2^i:end-1, end);

[NS_Rho, NS_Mom, NS_Ene] = SpectralNS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, p_0, Q1, Q2);
NS_Rho_list(:, i+1) = NS_Rho(1:2^i:end-1, end);%; NS_Rho(end,end)];
NS_Mom_list(:, i+1) = NS_Mom(1:2^i:end-1, end);%; NS_Mom(end,end)];
NS_Ene_list(:, i+1) = NS_Ene(1:2^i:end-1, end);%; NS_Ene(end,end)];
NS_Work_list(:, i+1) = 0.5*NS_Mom(1:2^i:end-1, end).^2./NS_Rho(1:2^i:end-1, end);

end

NSS_Rho_error = zeros(1, N-1);
NSS_Mom_error = zeros(1, N-1);
NSS_Ene_error = zeros(1, N-1);
NSS_Work_error = zeros(1, N-1);

NS_Rho_error = zeros(1, N-1);
NS_Mom_error = zeros(1, N-1);
NS_Ene_error = zeros(1, N-1);
NS_Work_error = zeros(1, N-1);


for i = 1:N-1
    NSS_Rho_error(i) = sum(sqrt(h1*(NSS_Rho_list(:,end) - NSS_Rho_list(:,i)).^2));
    NSS_Mom_error(i) = sum(sqrt(h1*(NSS_Mom_list(:,end) - NSS_Mom_list(:,i)).^2));
    NSS_Ene_error(i) = sum(sqrt(h1*(NSS_Ene_list(:,end) - NSS_Ene_list(:,i)).^2));
    NSS_Work_error(i) = sum(sqrt(h1*(NSS_Work_list(:,end) - NSS_Work_list(:,i)).^2));
    
    NS_Rho_error(i) = sum(sqrt(h1*(NS_Rho_list(:,end) - NS_Rho_list(:,i)).^2));
    NS_Mom_error(i) = sum(sqrt(h1*(NS_Mom_list(:,end) - NS_Mom_list(:,i)).^2));
    NS_Ene_error(i) = sum(sqrt(h1*(NS_Ene_list(:,end) - NS_Ene_list(:,i)).^2));
    NS_Work_error(i) = sum(sqrt(h1*(NS_Work_list(:,end) - NS_Work_list(:,i)).^2));

end

NSS_Rho_error = NSS_Rho_error/mean(abs(NSS_Rho(:, end)));
NSS_Mom_error = NSS_Mom_error/mean(abs(NSS_Mom(:, end)));
NSS_Ene_error = NSS_Ene_error/mean(abs(NSS_Ene(:, end)));
NSS_Work_error = NSS_Work_error/mean(abs(NSS_Work(:, end)));

NS_Rho_error = NS_Rho_error/mean(abs(NS_Rho(:, end)));
NS_Mom_error = NS_Mom_error/mean(abs(NS_Mom(:, end)));
NS_Ene_error = NS_Ene_error/mean(abs(NS_Ene(:, end)));
NS_Work_error = NS_Work_error/mean(abs(NS_Work(:, end)));

NSS_Work_error
NS_Work_error

%[NS_Work] = SpectralNS(m, x, t_end, k, c_p, c_V, kappa, mu, rho_0, T_0, u, Q1, Q2);

%{
[Q1,~,~,~] = PeriodicD0(m,h);
[Q2,~,~,~] = PeriodicD2(m,h);

[NSS_Work2] = SpectralNSS(m, x, t_end, k, c_p, c_V, mu, rho, Q1);

[~,Q1,~,~] = PeriodicD0(m,h);
[~,Q2,~,~] = PeriodicD2(m,h);

[NSS_Work3] = SpectralNSS(m, x, t_end, k, c_p, c_V, mu, rho, Q1);

[~,~,Q1,~] = PeriodicD0(m,h);
[~,~,Q2,~] = PeriodicD2(m,h);

[NSS_Work4] = SpectralNSS(m, x, t_end, k, c_p, c_V, mu, rho, Q1);

[~,~,~,Q1] = PeriodicD0(m,h);
[~,~,~,~] = PeriodicD2(m,h);

[NSS_Work5] = SpectralNSS(m, x, t_end, k, c_p, c_V, mu, rho, Q1);
%[NS_Work] = NS_feil(m, x, t_end, k, mu, c_V, c_p, R, kappa, Q1, Q2, rho)
%}

%{
time = linspace(0, t_end, size(NS_Rho,2));

f1 = figure;
plot(time, NS_Rho(1,:))
title('density at boundary')
f2 = figure;
plot(time, NS_Mom(1,:))
title('momentum at boundary')
f3 = figure;
plot(time, NS_Ene(1,:))
title('energy at boundary')
f4 = figure;
mesh(x, time(round(9*end/10):end), NS_Rho(:,round(9*end/10):end)')
title('last 10% of density')
f5 = figure;
mesh(x, time(round(9*end/10):end), NS_Mom(:,round(9*end/10):end)')
title('last 10% of momentum')
f4 = figure;
mesh(x, time(round(9*end/10):end), NS_Ene(:,round(9*end/10):end)')
title('last 10% of energy')

%plot(x, Mom(:,end))
%}
%}

%%


%NSS_Work = 0.5*NSS_Mom.^2./NSS_Rho;
%NS_Work = 0.5*NS_Mom.^2./NS_Rho;

%[NSS_Work] = SpectralNSS(m, x, t_end, k, c_p, c_V, mu, rho, Q1);

%[NSS_Work] = SpectralNSS(m, x, t_end, k, c_p, c_V, mu, rho, Q1);
%[NSS_Work] = E2(m, x, t_end, k, mu, alpha, c_V, c_p, R, Q1, rho);
%}

%Normal differnce schemes for isentropic (constant entropy)
%{
NSIP = NSI(m, x, t, k, mu, c_v, c_p, D0, D2, rho);

EIP = EI(m, x, t, k, mu, alpha1, c_v, c_p, D0, rho);
%}

%time-step convergence error thing
%{
NS1 = NS(m, t_end, k, mu, c_v, c_p, R, K, rho);
NS2 = NS(2*m-1, t_end, k, mu, c_v, c_p, R, K, rho);
NS2vals = reshape(NS2(1:end-1), 2, []);
NS2vals = [NS2vals(1,:), NS2(end)];


E2 = E(2*m-1, t_end, k, mu, alpha1, c_v, c_p, R, rho);
E2vals = reshape(E2(1:end-1), 2, []);
E2vals = [E2vals(1,:), E2(end)];
E1 = E(m, t_end, k, mu, alpha1, c_v, c_p, R, rho);

E_error = norm(E2vals-E1);
NS_error = norm(NS2vals-NS1);
%}

%l^2 convergence with spectral schemes
%{
N = 1;

RHO_p_list = zeros(1, N-1);
MOM_p_list = zeros(1, N-1);
ENE_p_list = zeros(1, N-1);

RHO_list = zeros(1, N);
MOM_list = zeros(1, N);
ENE_list = zeros(1, N);

h_list = zeros(1, N);
CFL_list = zeros(1, N);

M = 40;
%X = linspace(0, x_end, M);    %individual space points on the grid
%H= X(2) - X(1);


H = 2*pi/M;    %individual space points on the grid
X = H*(1:M);

CFL = k/(H^2);

Q1 = SpectralD0(M, H);
Q2 = SpectralD2(M, H);

%Q1 = PeriodicD0(M, H);
%Q2 = PeriodicD2(M, H);

[RHO_sol, MOM_sol, ENE_sol] = SpectralNS(M, X, t_end, k, c_p, c_V, kappa, mu, rho, Q1, Q2);

for i = 1:N
    
    h = 2*pi/m;    %individual space points on the grid
    x = h*(1:m);
    
    %x = linspace(0, x_end, m);
    %h = x(2)-x(1)
    %CFL = k/(h^2);

    Q1 = SpectralD0(m, h);
    Q2 = SpectralD2(m, h);
    
    %Q1 = PeriodicD0(m, h);
    %Q2 = PeriodicD2(m, h);

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


%{
for i = 1:N-1
    RHO_p_list(i) = log(RHO_list(i+1)/RHO_list(i))/log(h_list(i+1)/h_list(i));
end
RHO_p_list
for i = 1:N-1
    MOM_p_list(i) = log(MOM_list(i+1)/MOM_list(i))/log(h_list(i+1)/h_list(i));
end
MOM_p_list
for i = 1:N-1
    ENE_p_list(i) = log(ENE_list(i+1)/ENE_list(i))/log(h_list(i+1)/h_list(i));
end
ENE_p_list
%}
RHO_list
MOM_list
ENE_list
%}

%l^2 convergence with FD schemes
%{
N = 1;

RHO_p_list = zeros(1, N-1);
MOM_p_list = zeros(1, N-1);
ENE_p_list = zeros(1, N-1);

RHO_list = zeros(1, N);
MOM_list = zeros(1, N);
ENE_list = zeros(1, N);

h_list = zeros(1, N);
CFL_list = zeros(1, N);

M = 41;
X = linspace(0, x_end, M);    %individual space points on the grid
H= X(2) - X(1);

CFL = k/(H^2);

[~, ~, ~, Q1] = PeriodicD0(M, H);
[~, ~, ~, Q2] = PeriodicD2(M, H);

[RHO_sol, MOM_sol, ENE_sol] = NS(M, X, t_end, k, c_p, c_V, kappa, mu, rho, Q1, Q2);

for i = 1:N
    
    x = linspace(0, x_end, m);
    
    h = x(2)-x(1);
    h_list(i) = h;
    
    CFL = k/(h^2);

    [~, ~, ~, Q1] = PeriodicD0(m, h);
    [~, ~, ~, Q2] = PeriodicD2(m, h);

    CFL = k/(h^2);
    CFL_list(i) = CFL;

    [RHO, MOM, ENE] = NS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, Q1, Q2);
    
    rho_error = RHO - RHO_sol(1:2^(N-i+2):end);
    rho_error = sqrt(h*(rho_error'*rho_error));
    RHO_list(i) = rho_error;
    
    mom_error = MOM - MOM_sol(1:2^(N-i+2):end);
    mom_error = sqrt(h*(mom_error'*mom_error));
    MOM_list(i) = mom_error;
    
    ene_error = ENE - ENE_sol(1:2^(N-i+2):end);
    ene_error = sqrt(h*(ene_error'*ene_error));
    ENE_list(i) = ene_error;
    
    m = 2*m-1;
    
end


%{
for i = 1:N-1
    RHO_p_list(i) = log(RHO_list(i+1)/RHO_list(i))/log(h_list(i+1)/h_list(i));
end
RHO_p_list
for i = 1:N-1
    MOM_p_list(i) = log(MOM_list(i+1)/MOM_list(i))/log(h_list(i+1)/h_list(i));
end
MOM_p_list
for i = 1:N-1
    ENE_p_list(i) = log(ENE_list(i+1)/ENE_list(i))/log(h_list(i+1)/h_list(i));
end
ENE_p_list
%}
RHO_list
MOM_list
ENE_list
%}

%Round off error check
%{
x = linspace(0, x_end, m);    %individual space points on the grid
h = x(2) - x(1);

CFL = k/(h^2);

[~, ~, ~, Q1] = PeriodicD0(m, h);
[~, ~, ~, Q2] = PeriodicD2(m, h);

iter = 2;

NS_Work = zeros(length(x),iter);
NSS_Work = zeros(length(x),iter);

for i = 1:iter
    Navier_Stokes = NS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, Q1, Q2);
    NS_Work(:, i) = Navier_Stokes(:, end-1);
    Navier_Stokes_Svard = NSS(m, x, t_end, k, c_p, c_V, mu, rho, Q1);
    NSS_Work(:, i) = Navier_Stokes_Svard(:, end-1);
    k = 2*k;
end
%}

%mesh plot
%{
mesh(x,t, NSIP')
xlabel('t [sec]')
ylabel('Pressure')
title('Anisentropic Navier-Stokes')
%}

%Individual plot
%{
t = linspace(0,t_end,ceil(ceil(t_end/k)/100)+1);
plot(t, NS)
xlabel('t [sec]')
ylabel('Pressure')
title('Isentropic Navier-Stokes, 8th order operators')
%}

%Comparrison plots

%Comparison of decay of pressure
%{
%Getting max pressure amplitude data from Navier-Stokes-Svärd 
N_S_S = nonzeros(N_S_S);
N_S_S = N_S_S(1:round(length(N_S_S),-1));
N_S_S_plot = reshape(N_S_S, t_end, []);
N_S_S_plot = N_S_S_plot(1,:);
N_S_S_flux_plot = N_S_S_plot - p_0;         %subtracting background pressure to only show change in pressure

%Getting max pressure amplitude data from Navier-Stokes
N_S = nonzeros(N_S);
N_S = N_S(1:round(length(N_S),-1));
N_S_plot = reshape(N_S, t_end, []);
N_S_plot = N_S_plot(1,:);
N_S_flux_plot = N_S_plot - p_0;       %subtracting background pressure to only show change in pressure


tN_S_S = linspace(0,t_end,length(N_S_S_flux_plot));
tN_S = linspace(0, t_end, length(N_S_flux_plot));

%finding exponential coeficients
N_S_S_amp = max(reshape(N_S_S, [], t_end));
N_S_S_amp_flux = N_S_S_amp - p_0;
N_S_amp = max(reshape(N_S, [], t_end));
N_S_amp_flux = N_S_amp - p_0;
t2 = linspace(0,t_end,length(N_S_S_amp_flux)); %E_amp_flux and NS_amp_flux are same length

%least squares approximation to find y = exp(ax+b) for Navier-Stokes-Svärd
logN_S_S_amp_flux = log(N_S_S_amp_flux);
linearN_S_S_amp_flux = polyfit(t2, logN_S_S_amp_flux, 1);

%function giving least squares approximation of the data for
%Navier-Stokes-Svärd
LeastSquared_N_S_S = exp(linearN_S_S_amp_flux(2))*exp(linearN_S_S_amp_flux(1)*t2);

%least squares approximation to find y = exp(ax + b) for Navier-Stokes
logN_S_amp_flux = log(N_S_amp_flux);
linearN_S_amp_flux = polyfit(t2, logN_S_amp_flux, 1);

%function giving least squares approximation of the data for Navier-Stokes
LeastSquared_N_S = exp(linearN_S_amp_flux(2))*exp(linearN_S_amp_flux(1)*t2);

%initial condition
decayAmp = (LeastSquared_N_S_S(1) + LeastSquared_N_S(1))/2;

%Theoretical expectation for decay of Navier-Stokes sound waves
decayN_S = decayAmp*exp(-PropCoeffNS*t2);

%Theoretical expectation for decay pf Navier-Stokes-Svärd sound waves
decayN_S_S = decayAmp*exp(-PropCoeffNSS*t2);

%plots
tiledlayout(2,1)
nexttile
hold on
plot(tN_S, N_S_flux_plot);% 'color', 'blue')
plot(t2, LeastSquared_N_S, 'color', 'cyan')
plot(tN_S_S, N_S_S_flux_plot)%, 'color', 'yellow')
plot(t2, LeastSquared_N_S_S, 'color', 'red')
legend('Navier-Stokes', 'Navier-Stokes best fit', 'Eulerian', 'Eulerian best fit')
xlabel('t [sec]')
ylabel('Pressure')
title('Navier-Stokes vs. Eulerian, step size = 1e-6, ')
grid on

nexttile
hold on
plot(t2, LeastSquared_N_S, 'color', 'cyan')
plot(t2, decayN_S);
plot(t2, LeastSquared_N_S_S, 'color', 'red')
plot(t2, decayN_S_S);
title('Navier-Stokes and Eulerian model with amplitude approximation')
legend('Navier-Stokes', 'Navier-Stokes approximation', 'Navier-Stokes-Svärd', 'Navier-Stokes-Svärd approximation')
xlabel('t [sec]')
ylabel('Pressure')
grid on

%}


%Comparison of decay of energy
%Getting max work amplitude data from Navier-Stokes-Svärd 
int_NSS_Work = trapz(x, NSS_Work);
int_NSS_Work = nonzeros(int_NSS_Work);
int_NSS_Work = int_NSS_Work(1:round(length(int_NSS_Work), -1));

%t_end*1e+02/1.3 for p_0 = 10^5, t_end*1e+02/3.3 for p_0 = 10^3, t_end*1e+02/6.3
int_NSS_Work_plot = reshape(int_NSS_Work, [], 100);

[~, max_pos_NSS] = max(int_NSS_Work_plot);
col_size = size(int_NSS_Work_plot, 1);
for i = 0:length(max_pos_NSS)-1
    max_pos_NSS(i+1) = col_size*i*k + max_pos_NSS(i+1)*k;
end
max_pos_NSS = max_pos_NSS*2;

int_NSS_Work_plot = max(int_NSS_Work_plot);
logint_NSS_Work_plot = log(int_NSS_Work_plot);
t2NSS = linspace(0, t_end, length(logint_NSS_Work_plot));
Linear_NSS_Work = polyfit(t2NSS, logint_NSS_Work_plot, 1);
LeastSquared_NSS_Work = exp(Linear_NSS_Work(2))*exp(Linear_NSS_Work(1)*t2NSS);
%}
%Getting max heat amplitude data from Navier-Stokes-Svärd
%{
NSS_Heat = trapz(x, N_S_SHeat);
NSS_Heat = nonzeros(NSS_Heat);
NSS_Heat = NSS_Heat(1:round(length(NSS_Heat),-1));
NSS_Heat_plot = reshape(NSS_Heat, [], t_end);
NSS_Heat_plot = min(NSS_Heat_plot);
log_NSS_Heat_plot = log(NSS_Heat_plot);
t2NSS = linspace(0, t_end, length(log_NSS_Heat_plot));
Linear_NSS_Heat = polyfit(t2NSS, log_NSS_Heat_plot, 1);
LeastSquared_NSS_Heat = exp(Linear_NSS_Heat(2))*exp(Linear_NSS_Heat(1)*t2NSS);
%}
%Getting max work amplitude data from Navier-Stokes
int_NS_Work = trapz(x, NS_Work);
int_NS_Work = nonzeros(int_NS_Work);
int_NS_Work = int_NS_Work(1:round(length(int_NS_Work), -1));

%t_end*1e+02/1.3 for p_0 = 10^5, t_end*1e+02/3.3 for p_0 = 10^3,
%t_end*1e+02/6.3 for t_end = 10^4
int_NS_Work_plot = reshape(int_NS_Work, [], 100);

[~, max_pos_NS] = max(int_NS_Work_plot);
col_size = size(int_NS_Work_plot, 1);
for i = 0:length(max_pos_NS)-1
    max_pos_NS(i+1) = col_size*i*k + max_pos_NS(i+1)*k;
end
max_pos_NS = max_pos_NS*2;

int_NS_Work_plot = max(int_NS_Work_plot);
logint_NS_Work_plot = log(int_NS_Work_plot);
t2NS = linspace(0, t_end, length(logint_NS_Work_plot));
Linear_NS_Work = polyfit(t2NS, logint_NS_Work_plot, 1);
LeastSquared_NS_Work = exp(Linear_NS_Work(2))*exp(Linear_NS_Work(1)*t2NS);

%initial condition
%EnergyDecayAmp = (LeastSquared_NSS_Work(1) + LeastSquared_NS_Work(1))/2;

%Theoretical expectation for decay of Navier-Stokes sound waves
EnergyDecay_NS = LeastSquared_NS_Work(1)*exp(-2*PropCoeffNS*t2NS);

%Theoretical expectation for decay pf Navier-Stokes-Svärd sound waves
EnergyDecay_NSS = LeastSquared_NSS_Work(1)*exp(-2*PropCoeffNSS*t2NSS);

t2 = linspace(0, t_end, length(int_NSS_Work));

%plots
f1 = figure;
hold on
plot(t2, int_NS_Work);% 'color', 'blue')
plot(t2NS, LeastSquared_NS_Work, 'color', 'cyan')
plot(t2, int_NSS_Work)%, 'color', 'yellow')
plot(t2NSS, LeastSquared_NSS_Work, 'color', 'red')
t3 = t2(1:t_end*1e+03/1.3:end);
plot(max_pos_NS, int_NS_Work_plot, '*')
plot(max_pos_NSS, int_NSS_Work_plot, '*')
legend('NS', 'NS least squared', 'NSS', 'NSS least squared', 'NS data point', 'NSS data point')
xlabel('t [sec]')
ylabel('Work')
title('Work for NS and NSS, step size = 1e-6, t_end = 1.3 sec, ~ 1k waves')
grid on


f2 = figure;
hold on
ylim([0, inf])
plot(t2NS, LeastSquared_NS_Work, 'color', 'cyan')
plot(t2NS, EnergyDecay_NS);
plot(t2NSS, LeastSquared_NSS_Work, 'color', 'red')
plot(t2NSS, EnergyDecay_NSS);
title('NS vs. NSS with Theoretical decay')
legend('NS', 'NS theoretical', 'NSS', 'NSS theoretical')
xlabel('t [sec]')
ylabel('Work')
grid on
%}

%Theoretical_NS_coefficient = PropCoeffNS*2
Numerical_NS_coefficient = -Linear_NS_Work(1)

%Theoretical_NSS_coefficient = PropCoeffNSS*2
Numerical_NSS_coefficient = -Linear_NSS_Work(1)
%}


%change in frequency
%{
t_freq = linspace(0, t_end, length(E));
freq_dataE = [];
freq_indexE = [];
for i = 2:length(E)
    if E(i-1) > p_0 && E(i) < p_0 || E(i-1) < p_0 && E(i) > p_0
        freq_dataE = [freq_dataE, (E(i)+E(i-1))/2];
        freq_indexE = [freq_indexE, t_freq(i)];
    end
end

freq_dataNS = [];
freq_indexNS = [];
for i = 2:length(NS)
    if NS(i-1) > p_0 && NS(i) < p_0 || NS(i-1) < p_0 && NS(i) > p_0
        freq_dataNS = [freq_dataNS, (NS(i)+NS(i-1))/2];
        freq_indexNS = [freq_indexNS, i];
    end
end
        

%{

%exponential approximation values
%gamma = c_p/c_v;
%background_pressure = rho_bar^gamma;
%initial_amplitude = (rho_bar + 2e-04)^gamma - background_pressure;

E2_amplitude = E2(1) - background_pressure;
slope_E2 = log((E2(end)-background_pressure)/(E2_amplitude*background_pressure))/(tE2(end)-tE2(1));

NS2_amplitude = NS(1) - background_pressure;
slope_NS2 = log((NS2(end)-background_pressure)/(NS2_amplitude*background_pressure))/(tNS2(end)-tNS2(1));
%slope_E2 = log(E2(end)/E2(1))/(tE2(end)-tE2(1));
%slope_NS2 = log(NS2(end)/NS2(1))/(tNS2(end)-tNS2(1));
f = initial_amplitude*exp(slope_E2*tE2) + background_pressure;
g = initial_amplitude*exp(slope_NS2*tNS2) + background_pressure;
nexttile
hold on
plot(tNS2, NS2)
plot(tNS2,g)
plot(tE2, E2)
plot(tE2,f)
%semilogy(tNS2, NS2)
%semilogy(tE2, E2)
%}
%}

%}

%Model plot with decay
%{
decay = EP1(1)*exp(-PropCoeff*t);

t = t(1:100:end);
EP1 = EP1(1:100:end);
plot(t,EP1)
hold on
plot(t,decay)
%}

%Tiled plots
%{
tiledlayout(2,1)
nexttile
plot(t,NSP)
xlabel('t [sec]')
ylabel('Pressure')
title('Navier-Stokes at x = 1/4')

nexttile
plot(t,NSP2)
xlabel('t [sec]')
ylabel('Pressure')
title('Navier-Stokes at x = 3/4')
%}