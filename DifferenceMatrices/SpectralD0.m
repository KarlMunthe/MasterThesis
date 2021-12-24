function [D0] = SpectralD0(m, s)
%This function computes the spectral difference matrix approximating the
%first derivative for an equispaced grid. 
%m is the number of points. s is a scaling factor. 
%Parts of this code is inspired form Nick Trefethens book "Spectral Methods
%in Matlab".

%computes matrix for even number of grid points.
if rem(m,2) == 0
    
    h = 2*pi/m;     %step size

    column = s*[0, 0.5*(-1).^(1:m-1).*cot((1:m-1)*h/2)];    %computes first column

    D0 = toeplitz(column,column([1, m:-1:2]));  %creates difference matrix
    
end

%computes matrix for odd number of grid points.
if rem(m,2) == 1
    
    h = 2*pi/m;     %step size

    column = s*[0, 0.5*(-1).^(1:m-1)./sin((1:m-1)*h/2)];    %computes first column

    D0 = toeplitz(column,column([1, m:-1:2]));  %creates difference matrix
  
end
