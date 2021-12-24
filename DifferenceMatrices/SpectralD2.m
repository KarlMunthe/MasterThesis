function [D2] = SpectralD2(m, s)
%This function computes the spectral difference matrix approximating the
%second derivative for an equispaced grid. 
%m is the number of points. s is a scaling factor. 
%Parts of this code is inspired form Nick Trefethens book "Spectral Methods
%in Matlab".

%computes matrix for even number of grid points.
if rem(m, 2) == 0

    h = 2*pi/m;     %stepsize

    column = [-pi^2/(3*h^2)-1/6, -0.5*(-1).^(1:m-1)./sin(h*(1:m-1)/2).^2];  %computes column

    D2 = s^2*toeplitz(column);  %puts together the matrix.
    
end

%computes matrix for odd number of grid points.
if rem(m, 2) == 1
    
    %for odd number of grid points one can simply square the first
    %derivative matrix.
    D2 = SpectralD0(m, s);  
    D2 = D2^2;
    
end

%{
h = 2*pi/m; 
    
column = [1/12-pi^2/(3*h^2), -0.25*(-1).^(1:m-1).*cos(h*(1:m-1)/2)./sin(h*(1:m-1)/2) - (-1).^(1:m-1).*cos(h*(1:m-1)/2)/(sin(h*(1:m-1)/2).^2)];

D2 = s^2*toeplitz(column);
%}