clear

addpath('/Users/karlmunthe/Documents/UiB/UiBMasterOppgave/Code/gitMaster/Difference Operators')

%Standard temperature and pressure (STP)
T_0 = 273.15;
p_0 = 1e+05;
%p_0_list = [1e+05, 1e+04, 1e+03];

%Oxygen values

mu = 20.64e-06;
kappa = 26.58e-03;
c_p = 915;
c_V = 659;

%t_end = 0.13; %0.33 % 1.3; %c ~ 1000
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

%Gas constant and heat capacity
R = c_p - c_V;
gamma = c_p/c_V;

rho_0 = p_0/(R*T_0);  %background density at STP
rho_list = [10^-8, 10^-7, 10^-3];

nu_0 = mu/rho_0;
c = sqrt(gamma*R*T_0);    %speed of sound
lambda = 1;        %wave length
omega = 2*pi*c;
eta = 0;

%Coefficient of absorption as a function of time
%%%DENNE BRUKTE JEG OPPRINNELIG, MEN HVOR KOMMER C_V FRA?? DET ER IKKE FRA
%%%LANDAU OG LIPSCHITZ HVERTFALL. KAPPA/C_V SKAL BARE VÆRE KAPPA OG C_V*MU SKAL BARE VÆRE MU. 
%{
PropCoeffNS = omega^2/(2*rho_0*c^2)*((4/3*mu + eta) + (kappa/c_V)*(1/c_V-1/c_p)); %page 301 in Landou & Lipschitz
PropCoeffNSS = omega^2/(2*rho_0*c^2)*(mu + c_V*mu*(gamma-1)*(T_0/c^2 + 1/c_p));
%}

PropCoeffNS = omega^2/(2*rho_0*c^2)*((4/3*mu + eta) + (kappa)*(1/c_V-1/c_p));
PropCoeffNSS = c_V*omega^2/(2*rho_0*c^2)*(mu/c_V + mu*(gamma-1)*(T_0/c^2 + 1/c_p));

PropCoeffNSS = omega^2/(2*rho_0*c^2)*mu*(R*T_0/c^2 + 2 - 1/gamma);

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

%x = 0:h:x_end;    %individual space points on the grid

%Eulerian viscosity coefficient alpha can equal 1 or 4/3
alpha = 1;
alpha2 = 4/3;

%Spectral difference schemes

h = x_end/m;    %individual space points on the grid


s = 2*pi;

Q1 = SpectralD0(m, s);
Q2 = SpectralD2(m, s);

NS_Work_list = zeros(1,2);
NSS_Work_list = zeros(1,2);

for i = 1:1
        
    t_end = t_end_list(i);
    num_wave = num_wave_list(i);

    for j = 3:3
        
        rho = @(x) rho_0 + rho_list(j)*(sin(2*pi*x));   %reynolds decomposed pressure

        x = h*(1:m);

        [NSS_Work] = SpectralNSS2(m, x, t_end, k, c_p, c_V, mu, rho, p_0, Q1);

        [NS_Work] = SpectralNS(m, x, t_end, k, c_p, c_V, kappa, mu, rho, p_0, Q1, Q2);

        x = [0, x];

        %[LNSW, LNSSW] = WorkAnalysis(x, k, t_end, PropCoeffNS, PropCoeffNSS, NS_Work, NSS_Work, num_wave);
        [LNSW, LSNSW, LNSSW, LSNSSW, time] = WorkAnalysis2(x, k, t_end, PropCoeffNS, PropCoeffNSS, NS_Work, NSS_Work, num_wave);
        %LNSW = LNSW';
        %LNSSW = LNSSW';
    end
end

%NS_Work_list = [[-2*PropCoeffNS, mean(NS_Work_list(:, 2))]; NS_Work_list];
%NSS_Work_list = [[-2*PropCoeffNSS, mean(NSS_Work_list(:, 2))]; NSS_Work_list];
%timestart = t_end/10;
%timestep = (t_end - t_end/10)/3;
%time = linspace(timestart, timestart + timestep, length(LSNSW));
%ime = 

TheoNS = LSNSW(3, 1).*exp(-2*PropCoeffNS*time(3, :));
TheoNSS = LSNSSW(3, 1).*exp(-2*PropCoeffNSS*time(3, :));

NSPlotList = LSNSW(3, :);
NSSPlotList = LSNSSW(3, :);

%NSPlotList = exp(LNSW(:, 2)).*exp(LNSW(:, 1)*time);
%NSSPlotList = exp(LNSSW(:, 2)).*exp(LNSSW(:, 1)*time);

fig(1) = semilogy(time(3, :), TheoNS, 'color', 'black');
hold on
fig(2) = semilogy(time(3, :), TheoNSS, 'color', 'magenta');
fig(3) = semilogy(time(3, :), NSPlotList, 'color', 'red');
fig(4) = semilogy(time(3, :), NSSPlotList, 'color', 'blue');
legend([fig(1), fig(2), fig(3), fig(4)], 'Theoretical NS', 'Theoretical modNS', 'Numerical NS', 'Numerical modNS')
grid on
xlabel('time')
ylabel('work')
ylim([min(TheoNS), max(NSSPlotList)])
%REMEMBER TO TRANSPOSE THE NS(S)_WORK_LIST WHEN YOU COPY THEM INTO EXCEL


