clear

H = zeros(1000, 1e+04);

for i = 1:1e+05-1
    
    H(:,i+1) = H(:,i) -1 + 2*rand(1000,1);
    
end

a = 0;
for i = 1:1000
    if H(i,end) <= (1e+04)^(3/4)
        a = a + 1;
    end
end

%{
x = 0:1e+04;
y = sqrt(x);
    
plot(H')
hold on
plot(x,y)
%}