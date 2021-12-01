function [D] = FDm2(m, h, o)

%m = number of points
%h = space step
%o = order of accuracy

%creates D_-
Dm = diag(ones(m, 1), 0) + diag(-ones(m-1, 1), -1);

%comment out if you don't want periodic operators
Dm(1, end) = -1;

%finds coefficints of taylor series
x = sym('x');
f = (log(1-x))^2;
vector = double(coeffs(taylor(f, x, 'Order', o + 2)));

%first value
q = (vector(1))*eye(m);

for n = 1:o-1
    q = q + vector(n+1)*Dm^n;
end

D = (Dm/h)^2*q;