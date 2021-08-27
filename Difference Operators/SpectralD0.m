function [D0] = SpectralD0(m, s)

if rem(m,2) == 0
    
    h = 2*pi/m;

    column = s*[0, 0.5*(-1).^(1:m-1).*cot((1:m-1)*h/2)];

    D0 = toeplitz(column,column([1, m:-1:2]));
    
end


if rem(m,2) == 1
    
    h = 2*pi/m;

    column = s*[0, 0.5*(-1).^(1:m-1)./sin((1:m-1)*h/2)];

    D0 = toeplitz(column,column([1, m:-1:2]));
  
end
