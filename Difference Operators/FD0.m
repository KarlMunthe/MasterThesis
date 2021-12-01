function [D] = FD0(m, h, o)

%m = number of points
%h = space step
%o = order of accuracy

%creates D_+ and D_-
Dp = diag(-ones(m, 1), 0) + diag(ones(m-1, 1), 1);
Dm = diag(ones(m, 1), 0) + diag(-ones(m-1, 1), -1);
D0 = - diag(0.5*ones(m-1, 1), -1) + diag(+0.5*ones(m-1, 1), 1);

%comment out if you don't want periodic operators
Dp(end, 1) = 1;
Dm(1, end) = -1;
D0(end, 1) = 0.5; D0(1, end) = -0.5;


%finds coefficints of taylor series
x = sym('x');
f = asin(x)/sqrt(1-x^2);
vector = double(coeffs(taylor(f, x, 'Order', o + 1)));

%first value
q = eye(m);

for n = 1:o/2-1
    q = q + (-1)^n*vector(n+1)/(2^(2*n))*(Dp*Dm)^n;
end

D = D0/h*q;

%{
C = [B(1:round(end/2),1); zeros(1,1); B(round(end/2)+1:end,1)];

R = [B(1, 1:round(end/2)), zeros(1,1), B(1, round(end/2)+1:end)];

D = toeplitz(C,R);

%}

%making the diagonal entries = 0
for i = 1:m
    D(i,i) = 0;
end