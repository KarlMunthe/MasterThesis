function [D0] = SpectralD0(m, s)

h = 2*pi/m;

column = s*[0, 0.5*(-1).^(1:m-1).*cot((1:m-1)*h/2)]';

D0 = toeplitz(column,column([1, m:-1:2]));
  
end