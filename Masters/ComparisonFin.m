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
p_0 = 1e+03;
%p_0_list = [1e+05, 1e+04, 1e+03];

rho_0 = p_0/(R*T_0);  %background density at STP
rho_list = [10^-8, 10^-7, 10^-6];
%rho = @(x) rho_0 + 1e-08*(sin(2*pi*x));   %reynolds decomposed pressure
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
%t_end = 0.13; %0.33 % 1.3; %c ~ 1000
%t_end_list = [0.63, 6.3, 63]; % for p_0 = 10^4
t_end_list = [0.32, 3.2, 32]; %for p_0 = 10^3
%t_end_list = [0.13, 1.3, 13]; %for p_0 = 10^5
num_wave_list = [10^2, 10^3, 10^4];
num_wave_list = num_wave_list/10;

%x = 0:h:x_end;    %individual space points on the grid

%Eulerian viscosity coefficient alpha can equal 1 or 4/3
alpha = 1;
alpha2 = 4/3;

%Spectral difference schemes

h = x_end/m;    %individual space points on the grid


s = 2*pi;

Q1 = SpectralD0(m, s);
Q2 = SpectralD2(m, s);

NS_Work_list = zeros(3,3);
NSS_Work_list = zeros(3,3);


for i = 1:3
    
    i
        
    t_end = t_end_list(i);
    num_wave = num_wave_list(i);

    for j = 1:3
        
        j

        rho = @(x) rho_0 + rho_list(j)*(sin(2*pi*x));   %reynolds decomposed pressure

        x = h*(1:m);

        [NSS_Work] = SpectralNSStrash(m, x, t_end, k, c_p, c_V, mu, rho, p_0, Q1);

        [NS_Work] = SpectralNS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, p_0, Q1, Q2);

        x = [0, x];

        [NumNSCoeff, NumNSSCoeff] = WorkAnalysis(x, k, t_end, PropCoeffNS, PropCoeffNSS, NS_Work, NSS_Work, num_wave);
        NS_Work_list(i, j) = NumNSCoeff;
        NSS_Work_list(i, j) = NumNSSCoeff;
        
    end
end

%REMEMBER TO TRANSPOSE THE NS(S)_WORK_LIST WHEN YOU COPY THEM INTO EXCEL


