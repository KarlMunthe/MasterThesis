function [D] = FDp(m, h, o)

%m = number of points
%h = space step
%o = order of accuracy

%creates D_+
Dp = diag(-ones(m, 1), 0) + diag(ones(m-1, 1), 1);

%comment out if you don't want periodic operators
Dp(end, 1) = 1;

%finds coefficints of taylor series
x = sym('x');
f = log(1+x);
vector = double(coeffs(taylor(f, x, 'Order', o+1)));

%first value
q = (vector(1))*eye(m);

for n = 1:o-1
    q = q + vector(n+1)*Dp^n;
end

D = Dp/h*q;