function [D] = FD02(m, h, o)

%m = number of points
%h = space step
%o = order of accuracy

%creates D_+ and D_-
Dp = diag(-ones(m, 1), 0) + diag(ones(m-1, 1), 1);
Dm = diag(ones(m, 1), 0) + diag(-ones(m-1, 1), -1);
D0 = - diag(0.5*ones(m-1, 1), -1) + diag(0.5*ones(m-1, 1), 1);

%comment out if you don't want periodic operators
Dp(end, 1) = 1;
Dm(1, end) = -1;
D0(end, 1) = 0.5; D0(1, end) = -0.5;

%finds coefficints of taylor series
x = sym('x');
f = (asin(x))^2;
vector = double(coeffs(taylor(f, x, 'Order', o + 2)));

%first value
q = (vector(1))*eye(m);

for n = 1:o/2-1
    q = q + (-1)^n*vector(n+1)/(2^(2*n))*(Dp*Dm)^n;
end

D = Dp*Dm/(h)^2*q;