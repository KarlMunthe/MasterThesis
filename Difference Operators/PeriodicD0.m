function [D02, D04, D06, D08] = PeriodicD0(m,h)

%Difference Operators
%Creates first order derivative with periodic boundary conditions
%D02 = 2nd order
%D04 = 4th order
%D06 = 6th order
%D08 = 8th order

m = m+1;

Dp = diag(-ones(m-1, 1), 0) + diag(ones(m-2, 1), 1);
Dp(end, 1) = 1;
Dp = Dp/h;
Dm = diag(ones(m-1, 1), 0) + diag(-ones(m-2, 1), -1);
Dm(1, end) = -1;
Dm = Dm/h;

%determines accuracy
q4 = (eye(m-1) - h^2*Dp*Dm/6);
q6 = (eye(m-1) - h^2*Dp*Dm/6 + h^4*(Dp*Dm)^2/30);
q8 = (eye(m-1) - h^2*Dp*Dm/6 + h^4*(Dp*Dm)^2/30 - h^6*(Dp*Dm)^3/140);

D02 = - diag(0.5*ones(m-2, 1), -1) + diag(+0.5*ones(m-2, 1), 1);
D02(end, 1) = 0.5; D02(1, end) = -0.5;
D02 = D02/h;

D04 = D02*q4;
D04 = D04.*~eye(size(D04));

D06 = D02*q6;
D06 = D06.*~eye(size(D06));

D08 = D02*q8;
D08 = D08.*~eye(size(D08));

%{

D04= -1/12*diag(ones(m-3,1),2)+8/12*diag(ones(m-2,1),1)- ...
    8/12*diag(ones(m-2,1),-1)+1/12*diag(ones(m-3,1),-2);
D04(1:2,end-1:end) = [1/12, -2/3;0,1/12];
D04(end-1:end,1:2)=fliplr(flipud(-D04(1:2,end-1:end)));
D04 = D04/h;

D06 = -diag(3/4*ones(m-2, 1), -1) + diag(3/4*ones(m-2, 1), 1)...
    + diag(3/20*ones(m-3, 1), -2) - diag(3/20*ones(m-3, 1), 2) - diag(1/60*ones(m-4, 1), -3)...
    + diag(1/60*ones(m-4, 1), 3);
D06(1:3,end-2:end) = [-1/60, 3/20, -3/4; 0, 3/20, -3/4; 0, 0, -3/4];
D06(end-2:end,1:3)=fliplr(flipud(-D06(1:3,end-2:end)));
D06 = D06/h;

D08 = -diag(4/5*ones(m-2, 1), -1) + diag(4/5*ones(m-2, 1), 1)...
    + diag(1/5*ones(m-3, 1), -2) - diag(1/5*ones(m-3, 1), 2) - diag(4/105*ones(m-4, 1), -3)...
    + diag(4/105*ones(m-4, 1), 3) + diag(1/280*ones(m-5, 1), -4) - diag(1/280*ones(m-5, 1), 4);
D08(1:4,end-3:end) = [1/280, -4/105, 1/5, -4/5; 0, 1/280, -4/105, 1/5; 0, 0,1/280, -4/105; 0, 0, 0, 1/280];
D08(end-3:end,1:4)=fliplr(flipud(-D08(1:4,end-3:end)));
D08 = D08/h;
%}



