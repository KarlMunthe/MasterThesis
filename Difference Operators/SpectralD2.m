function [D2] = SpectralD2(m, s)

h = 2*pi/m; 

column = s^2*[-pi^2/(3*h^2)-1/6, -0.5*(-1).^(1:m-1)./sin(h*(1:m-1)/2).^2]';

D2 = toeplitz(column,column([1, m:-1:2]));
  
end