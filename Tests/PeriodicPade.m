function [LHS1, RHS1, LHS2, RHS2] = PeriodicPade(m,h)

%10TH ORDER
%FIRST DERIVATIVE

alpha = 1/2;
beta = 1/20;
a = 17/12;
b = 101/150;
c = 1/100;
%}
%SPECTRAL LIKE
%{
alpha = 0.5771439;
beta = 0.0896406;
a = 1.3025166;
b = 0.9935500;
c = 0.03750245;
%}

%FIRST DERIVATIVE
LHS1 = eye(m-1) + beta*diag(ones(m-3,1), -2) + alpha*diag(ones(m-2,1), -1)...
    + beta*diag(ones(m-3,1), 2) + alpha*diag(ones(m-2,1), 1);

ABottomCorner = [beta, 0;alpha, beta];
ATopCorner = rot90(ABottomCorner, 2);

LHS1(1:2, end-1:end) = ATopCorner;
LHS1(end-1:end, 1:2) = ABottomCorner;


RHS1 = - c/6*diag(ones(m-4,1), -3) - b/4*diag(ones(m-3,1), -2) - a/2*diag(ones(m-2,1), -1)...
    + c/6*diag(ones(m-4,1), 3) + b/4*diag(ones(m-3,1), 2) + a/2*diag(ones(m-2,1), 1);

BBottomCorner = [c/6, 0, 0; b/4, c/6, 0; a/2, b/4, c/6];
BTopCorner = rot90(-BBottomCorner, 2);

RHS1(1:3, end-2:end) = BTopCorner;
RHS1(end-2:end, 1:3) = BBottomCorner;

RHS1 = RHS1/h;

%SECOND DERIVATIVE
alpha = 334/899;
beta = 43/1798;
a = 1065/1798;
b = 1038/899;
c = 79/1798;

LHS2 = eye(m-1) + beta*diag(ones(m-3,1), -2) + alpha*diag(ones(m-2,1), -1)...
    + beta*diag(ones(m-3,1), 2) + alpha*diag(ones(m-2,1), 1);

CBottomCorner = [beta, 0;alpha, beta];
CTopCorner = rot90(CBottomCorner, 2);

LHS2(1:2, end-1:end) = CTopCorner;
LHS2(end-1:end, 1:2) = CBottomCorner;


RHS2 = -2*(a+b/4+c/9)*diag(ones(m-1, 1)) + c/9*diag(ones(m-4,1), -3)...
    + b/4*diag(ones(m-3,1), -2) + a*diag(ones(m-2,1), -1)...
    + c/9*diag(ones(m-4,1), 3) + b/4*diag(ones(m-3,1), 2)...
    + a*diag(ones(m-2,1), 1);

DBottomCorner = [c/9, 0, 0; b/4, c/9, 0; a, b/4, c/9];
DTopCorner = rot90(DBottomCorner, 2);

RHS2(1:3, end-2:end) = DTopCorner;
RHS2(end-2:end, 1:3) = DBottomCorner;

RHS2 = RHS2/(h^2);
%}
