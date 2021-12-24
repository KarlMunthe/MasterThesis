clear
%This script simulates and compares the acoustic attenuation for the
%Navier-Stokes and Magnus Svärd modified Navier-Stokes equations. 
addpath('/Users/karlmunthe/Documents/UiB/UiBMasterOppgave/Code/MasterThesis/DifferenceMatrices')

%Standard temperature and pressure (STP)
T_0 = 273.15;
p_0 = 1e+05;    

%Oxygen values

mu = 20.64e-06;
kappa = 26.58e-03;
c_p = 915;
c_V = 659;

%p_0 changes the speed of sound, and to obtain 10^4 waves the simulation
%has to run a certain amount of time. Empirically we have found these
%values. If p_0 = 10^3 we only extimate 10^3 waves because otherwise the
%simulation would take too long.
%t_end_list = [0.63, 6.3, 63]; % for p_0 = 10^4
t_end_list = [0.32, 3.2, 32]; %for p_0 = 10^3
%t_end_list = [0.13, 1.3, 13]; %for p_0 = 10^5
num_wave_list = [10^2, 10^3, 10^4];
num_wave_list = num_wave_list/10; %for p_0 = 10^3
%}

%Argon values
%{
mu = 22.61e-06;
kappa = 0.0178;
c_V = 313;
c_p = 520;

t_end_list = [0.104, 1.04, 10.4]; %for p_0 = 10^5
%t_end_list = [0.69, 6.9, 69]; %for p_0 = 10^4
%t_end_list = [0.47, 4.7, 47]; %for p_0 = 10^3;
num_wave_list = [10^2, 10^3, 10^4]; 
%num_wave_list = num_wave_list/10; %for p_0 = 10^3
%}

R = c_p - c_V;          %heat capacity difference
gamma = c_p/c_V;        %heat capacity ratio

rho_0 = p_0/(R*T_0);  %background density
rho_fluctuation = 1e-6;

nu_0 = mu/rho_0;        %background kinematic viscosity
c = sqrt(gamma*R*T_0);    %speed of sound
lambda = 1;        %wave length
omega = 2*pi*c;     %angular frequency
eta = 0;

PropCoeffNS = omega^2/(rho_0*c^2)*((4/3*mu + eta) + (kappa)*(1/c_V-1/c_p));

PropCoeffNSS = omega^2/(rho_0*c^2)*mu*(R*T_0/c^2 + 2 - 1/gamma);


m = 11;     %number of grid points

k = 1e-06;  %time step

x_end = 1;  %end of domain


%Spectral difference schemes

h = x_end/m;    %individual space points on the grid

s = 2*pi;       %scaling variable for spectral methods

Q1 = SpectralD0(m, s);  %difference matrices
Q2 = SpectralD2(m, s);  %difference matrices

%determine if you want to simulate the wave fluctuate 10^2, 10^3, or 10^4
%times.
t_end = t_end_list(1);
num_wave = num_wave_list(1);

%decomposes density into background density and density fluctuation
rho = @(x) rho_0 + rho_fluctuation*(sin(2*pi*x));   

%grid pointd
x = h*(1:m);

%computes the work for the NS equations
[NS_Work] = SpectralNS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, p_0, Q1, Q2);

%computes the work for the NSS equations
[NSS_Work] = SpectralNSS(m, x, t_end, k, c_p, c_V, mu, rho, p_0, Q1);

%add a grid point at the boundary for periodicity
x = [0, x];

%calculates the absorption coefficient of the NS and NSS equations by
%by means of least squared regression.
[LSNSW, LSNSSW, time] = WorkAnalysis(x, t_end, NS_Work, NSS_Work, num_wave);

%Analytical solution to linearized NS and NSS equations
TheoNS = LSNSW(3, 1).*exp(-PropCoeffNS*time(3, :));
TheoNSS = LSNSSW(3, 1).*exp(-PropCoeffNSS*time(3, :));


%plots the last third of the data. Plots both numerical and analytical
%solution fo linearized equations for both the NS and NSS equations
fig(1) = semilogy(time(3, :), TheoNS, 'color', 'black');
hold on
fig(2) = semilogy(time(3, :), TheoNSS, 'color', 'magenta');
fig(3) = semilogy(time(3, :), LSNSW(3, :), 'color', 'red');
fig(4) = semilogy(time(3, :), LSNSSW(3, :), 'color', 'blue');
legend([fig(1), fig(2), fig(3), fig(4)], 'Theoretical NS', 'Theoretical modNS', 'Numerical NS', 'Numerical modNS')
grid on
xlabel('time')
ylabel('work')