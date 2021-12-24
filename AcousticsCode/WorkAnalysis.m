function [LSNSW, LSNSSW, t] = WorkAnalysis(x, t_end, NS_Work, NSS_Work, num_wave)
%calculates the absorption coefficient of the NS and NSS equations by
%by means of least squared regression.

%Getting max work amplitude data from Navier-Stokes
int_NS_Work = trapz(x, NS_Work);    %integrate work in space
int_NS_Work = nonzeros(int_NS_Work);    %delete zero entries
int_NS_Work = int_NS_Work(1:round(length(int_NS_Work), -1));    %round to nearest 10
int_NS_Work_plot = reshape(int_NS_Work, [], num_wave);  %reshape such that one periode is approximately one column
int_NS_Work_plot = int_NS_Work_plot(:, size(int_NS_Work_plot, 2)/10+1:end);     %delete first 10% to eliminate initial errors
max_int_NS_Work_plot = max(int_NS_Work_plot);               %calculate the maximum value for each entry

NS_mech(1, :) = max_int_NS_Work_plot(1:size(max_int_NS_Work_plot, 2)/3);    %first third of work
NS_mech(2, :) = max_int_NS_Work_plot(size(max_int_NS_Work_plot, 2)/3+1:2*size(max_int_NS_Work_plot, 2)/3);   %second third of work
NS_mech(3, :) = max_int_NS_Work_plot(2*size(max_int_NS_Work_plot, 2)/3+1:end);       %last third of work
logNS_mech = log(NS_mech);  %ln of work

%Getting max work amplitude data from Navier-Stokes-Svärd
int_NSS_Work = trapz(x, NSS_Work);  %integrate work in space
int_NSS_Work = nonzeros(int_NSS_Work);  %delete zero entries
int_NSS_Work = int_NSS_Work(1:round(length(int_NSS_Work), -1)); %round to nearest 10
int_NSS_Work_plot = reshape(int_NSS_Work, [], num_wave);    %reshape such that one periode is approximately one column
int_NSS_Work_plot = int_NSS_Work_plot(:, size(int_NSS_Work_plot, 2)/10+1:end); %delete first 10% to eliminate initial errors
max_int_NSS_Work_plot = max(int_NSS_Work_plot); %calculate the maximum value for each entry

NSS_mech(1, :) = max_int_NSS_Work_plot(1:size(max_int_NSS_Work_plot, 2)/3); %first third of work
NSS_mech(2, :) = max_int_NSS_Work_plot(size(max_int_NSS_Work_plot, 2)/3+1:2*size(max_int_NSS_Work_plot, 2)/3);  %second third of work
NSS_mech(3, :) = max_int_NSS_Work_plot(2*size(max_int_NSS_Work_plot, 2)/3+1:end);   %last third of work
logNSS_mech = log(NSS_mech);    %ln of work

%create time vectors for first, second, and last third of the data
%excluding the first 10%
t2start = t_end/10;     %remove first 10% of time
t2step = (t_end - t_end/10)/3;  %time step
t = zeros(3, length(logNSS_mech));
t(1, :) = linspace(t2start, t2start + t2step, length(logNSS_mech));     %first third chunk of time
t(2, :) = linspace(t2start + t2step, t2start + 2*t2step, length(logNSS_mech));  %second third chunk of time
t(3, :) = linspace(t2start + 2*t2step, t_end, length(logNSS_mech));     %last third chunk of time

%Least squared regression of Navier-Stokes work points for each time vector
for i = 1:3
    Linear_NS_Work(i, :) = polyfit(t(i, :), logNS_mech(i, :), 1);
end
LeastSquared_NS_Work = exp(Linear_NS_Work(:, 2)).*exp(Linear_NS_Work(:, 1).*t);

%Least squared regression of Navier-Stokes-Svärd work points for eaech time
%vector
for i = 1:3
    Linear_NSS_Work(i, :) = polyfit(t(i, :), logNSS_mech(i, :), 1);
end

LeastSquared_NSS_Work = exp(Linear_NSS_Work(:, 2)).*exp(Linear_NSS_Work(:, 1).*t);

LSNSW = LeastSquared_NS_Work;
LSNSSW = LeastSquared_NSS_Work;
