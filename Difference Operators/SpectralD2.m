function [D2] = SpectralD2(m, s)

if rem(m, 2) == 0

    h = 2*pi/m; 

    column = [-pi^2/(3*h^2)-1/6, -0.5*(-1).^(1:m-1)./sin(h*(1:m-1)/2).^2];

    D2 = s^2*toeplitz(column);%,column([1, m:-1:2]));
    
end

if rem(m, 2) == 1
    
    D2 = SpectralD0(m, s);
    D2 = D2^2;
    
end

%{
h = 2*pi/m; 
    
column = [1/12-pi^2/(3*h^2), -0.25*(-1).^(1:m-1).*cos(h*(1:m-1)/2)./sin(h*(1:m-1)/2) - (-1).^(1:m-1).*cos(h*(1:m-1)/2)/(sin(h*(1:m-1)/2).^2)];

D2 = s^2*toeplitz(column);
%}