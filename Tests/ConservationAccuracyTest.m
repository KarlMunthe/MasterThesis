clear

addpath('/Users/karlmunthe/Documents/UiB/UiB Master Oppgave/Code/gitMaster/Difference Operators')

m = 7;
s = 1;
x_end = 2*pi;

h = x_end/m;    %individual space points on the grid
x = h*(1:m);

v_0 = sin(x)';
v = v_0;

Q1 = SpectralD0(m, s);
corr = Q1*v-cos(x)';

%%

for num = 1:10
    v = Q1*v;
end
v
%{
plot(x, v_0)
hold on
plot(x, v)
legend('first', 'not first')
%}
