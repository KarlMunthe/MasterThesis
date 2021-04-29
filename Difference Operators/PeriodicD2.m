function [D22, D24, D26, D28] = PeriodicD2(m,h)


%Difference Operators
%Creates second order derivative with periodic boundary conditions
%D02 = 2nd order
%D04 = 4th order
%D06 = 6th order
%D08 = 8th order

m = m+1;

%Creates second order derivative with periodic boundary conditions
Dp = diag(-ones(m-1, 1), 0) + diag(ones(m-2, 1), 1);
Dp(end, 1) = 1;
Dp = Dp/h;
Dm = diag(ones(m-1, 1), 0) + diag(-ones(m-2, 1), -1);
Dm(1, end) = -1;
Dm = Dm/h;

%determine accuracy
I = diag(ones(m-1,1));
q4 = (I - h^2*Dp*Dm/12);
q6 = (I - h^2*Dp*Dm/12 + h^4*(Dp*Dm)^2/90);
q8 = (I - h^2*Dp*Dm/12 + h^4*(Dp*Dm)^2/90 - h^6*(Dp*Dm)^3/560);

D22 = Dp*Dm;
D24 = Dp*Dm*q4;
D26 = Dp*Dm*q6;
D28 = Dp*Dm*q8;

%{
D22 = D22/(h*h);
%her hadde jeg skrevet D22 = 2*D22/(h*h);
%%%%%
%%%%%
%%%%%
D24 = -5/2*diag(ones(m-1,1),0) + 4/3*diag(ones(m-2,1),-1)...
    +4/3*diag(ones(m-2,1),1) - 1/12*diag(ones(m-3,1), -2) ...
    - 1/12*diag(ones(m-3,1), 2);
D24(1:2,end-1:end) = [-1/12, 4/3; 0, -1/12];
D24(end-1:end,1:2)=fliplr(flipud(D24(1:2,end-1:end)));
D24 = D24/(h*h);
%%%%%
%%%%%
%%%%%
D26 = -49/18*diag(ones(m-1,1),0) + 3/2*diag(ones(m-2,1),1) - 3/20*diag(ones(m-3,1),2) + 1/90*diag(ones(m-4,1),3)...
    + 3/2*diag(ones(m-2,1),-1) - 3/20*diag(ones(m-3,1),-2) + 1/90*diag(ones(m-4,1),-3);
D26(1:3,end-2:end) = [1/90, -3/20, 3/2; 0, 1/90, -3/20; 0, 0, 1/90];
D26(end-2:end,1:3)=fliplr(flipud(D26(1:3,end-2:end)));
D26 = D26/(h*h);
%%%%%
%%%%%
%%%%%
D28 = -205/72*diag(ones(m-1,1),0) + 8/5*diag(ones(m-2,1),1) - 1/5*diag(ones(m-3,1),2) + 8/315*diag(ones(m-4,1),3) - 1/560*diag(ones(m-5,1), 4)...
    + 8/5*diag(ones(m-2,1),-1) - 1/5*diag(ones(m-3,1),-2) + 8/315*diag(ones(m-4,1),-3) - 1/560*diag(ones(m-5,1), -4);
D28(1:4, end-3:end) = [-1/560, 8/315, -1/5, 8/5; 0, -1/560, 8/315, -1/5; 0, 0, -1/560, 8/315;0, 0, 0, -1/560];
D28(end-3:end,1:4)=fliplr(flipud(D28(1:4,end-3:end)));
D28 = D28/(h*h);

%}


