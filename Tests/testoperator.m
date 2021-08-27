function D = testoperator(m, h)

B = FD0(3, 1, 100);

C = [B(1:round(end/2),1); zeros(m-3,1); B(round(end/2)+1:end,1)];

R = [B(1, 1:round(end/2)), zeros(1,m-3), B(1, round(end/2)+1:end)];

D = toeplitz(C,R);

D = D/h;

end
