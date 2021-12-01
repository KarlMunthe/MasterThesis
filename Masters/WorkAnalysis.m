function [LNSW, LNSSW] = WorkAnalysis(x, k, t_end, PropCoeffNS, PropCoeffNSS, NS_Work, NSS_Work, num_wave)

nr_points_save = 5;

int_NSS_Work = trapz(x, NSS_Work);
int_NSS_Work = nonzeros(int_NSS_Work);
int_NSS_Work = int_NSS_Work(1:round(length(int_NSS_Work), -1));

%t_end*1e+02/1.3 for p_0 = 10^5, t_end*1e+02/3.3 for p_0 = 10^3, t_end*1e+02/6.3
int_NSS_Work_plot = reshape(int_NSS_Work, [], num_wave);

[~, max_pos_NSS] = max(int_NSS_Work_plot);
col_size = size(int_NSS_Work_plot, 1);
for i = 0:length(max_pos_NSS)-1
    max_pos_NSS(i+1) = col_size*i*k + max_pos_NSS(i+1)*k;
end
max_pos_NSS = max_pos_NSS*nr_points_save;

int_NSS_Work_plot = max(int_NSS_Work_plot);
logint_NSS_Work_plot = log(int_NSS_Work_plot);
t2NSS = linspace(0, t_end, length(logint_NSS_Work_plot));
Linear_NSS_Work = polyfit(t2NSS, logint_NSS_Work_plot, 1);
LeastSquared_NSS_Work = exp(Linear_NSS_Work(2))*exp(Linear_NSS_Work(1)*t2NSS);

%Getting max work amplitude data from Navier-Stokes
int_NS_Work = trapz(x, NS_Work);
int_NS_Work = nonzeros(int_NS_Work);
int_NS_Work = int_NS_Work(1:round(length(int_NS_Work), -1));

%t_end*1e+02/1.3 for p_0 = 10^5, t_end*1e+02/3.3 for p_0 = 10^3, t_end*1e+02/6.3
int_NS_Work_plot = reshape(int_NS_Work, [], num_wave);

[~, max_pos_NS] = max(int_NS_Work_plot);
col_size = size(int_NS_Work_plot, 1);
for i = 0:length(max_pos_NS)-1
    max_pos_NS(i+1) = col_size*i*k + max_pos_NS(i+1)*k;
end
max_pos_NS = max_pos_NS*nr_points_save;

int_NS_Work_plot = max(int_NS_Work_plot);
logint_NS_Work_plot = log(int_NS_Work_plot);
t2NS = linspace(0, t_end, length(logint_NS_Work_plot));
Linear_NS_Work = polyfit(t2NS, logint_NS_Work_plot, 1);
LeastSquared_NS_Work = exp(Linear_NS_Work(2))*exp(Linear_NS_Work(1)*t2NS);

%initial condition
%EnergyDecayAmp = (LeastSquared_NSS_Work(1) + LeastSquared_NS_Work(1))/2;

%Theoretical expectation for decay of Navier-Stokes sound waves
EnergyDecay_NS = LeastSquared_NS_Work(1)*exp(-PropCoeffNS*t2NS);

%Theoretical expectation for decay pf Navier-Stokes-Svärd sound waves
EnergyDecay_NSS = LeastSquared_NSS_Work(1)*exp(-PropCoeffNSS*t2NSS);

t2 = linspace(0, t_end, length(int_NSS_Work));

%plots

f1 = figure;
hold on
%{
loglog(t2, int_NS_Work);% 'color', 'blue')
loglog(t2NS, LeastSquared_NS_Work, 'color', 'cyan')
loglog(t2, int_NSS_Work)%, 'color', 'yellow')
loglog(t2NSS, LeastSquared_NSS_Work, 'color', 'red')
loglog(max_pos_NS, int_NS_Work_plot, '*')
loglog(max_pos_NSS, int_NSS_Work_plot, '*')
%}

plot(t2, int_NS_Work);% 'color', 'blue')
plot(t2NS, LeastSquared_NS_Work, 'color', 'cyan')
plot(t2, int_NSS_Work)%, 'color', 'yellow')
plot(t2NSS, LeastSquared_NSS_Work, 'color', 'red')
plot(max_pos_NS, int_NS_Work_plot, '*')
plot(max_pos_NSS, int_NSS_Work_plot, '*')
%}
legend('NS', 'NS least squared', 'modNS', 'modNS least squared', 'NS data point', 'modNS data point')
xlabel('Time')
ylabel('Work')
%title('Work for NS and NSS, step size = 1e-6, t_{end} = 3.2 sec, 100 waves')
grid on


%Theoretical_NS_coefficient = PropCoeffNS*2
NumNSCoeff = -Linear_NS_Work(1);

%Theoretical_NSS_coefficient = PropCoeffNSS*2
NumNSSCoeff = -Linear_NSS_Work(1);

lNumNSCoeff = num2str(NumNSCoeff)
lPropCoeffNS = num2str(2*PropCoeffNS)
lNumNSSCoeff = num2str(NumNSSCoeff)
lPropCoeffNSS = num2str(2*PropCoeffNSS)

f2 = figure;
semilogy(t2NS, LeastSquared_NS_Work, 'color', 'cyan')
hold on
semilogy(t2NS, EnergyDecay_NS);
semilogy(t2NSS, LeastSquared_NSS_Work, 'color', 'red')
semilogy(t2NSS, EnergyDecay_NSS);

title('NS vs. NSS with Theoretical decay')
legend('NS', 'NS theoretical', 'NSS', 'NSS theoretical')
xlabel('t [sec]')
ylabel('Work')
grid on
%}
LNSW = Linear_NS_Work;
LNSSW = Linear_NSS_Work;
%}
