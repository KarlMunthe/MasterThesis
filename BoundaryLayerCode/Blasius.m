function [BlasProf, intBlasProf, vBlasProf] = Blasius(eta, fasit, tol)

f1 = 0;
f2 = 0;

output = [[f1;f2;0], zeros(3, length(eta)-1)];
%set initial max and min values for the initial conditions to Blasius ODE.
%The actual initial value must be betwen vunder and vover.
vunder = 0;
vover = 1;
vnew = (vover + vunder)/2;  %calculates first guess
hitnew = BlasiusProfile(f1, f2, vnew, eta, output); %calculates if this is over or under the true solution
q = 0;

while abs(hitnew - fasit) > tol
    q = q + 1;
    if hitnew == fasit      %if new initial condition is the right one stop loop.
    end
    if hitnew > fasit       %if new intial condition is to high create new upper bound
        vover = vnew;
        vnew = (vnew+vunder)/2;
        hitnew = BlasiusProfile(f1, f2, vnew, eta, output);
    end
    if hitnew < fasit       %if new intial condition is to high create new lower bound
        vunder = vnew;
        vnew = (vnew+vover)/2; 
        hitnew = BlasiusProfile(f1, f2, vnew, eta, output);
    end
    
end

output(3,1) = vnew;     %calculates Blasius solution with initial condition that satisfies the given tolerance
[~, output] = BlasiusProfile(f1, f2, vnew, eta, output);
BlasProf = output(2, :);    %Boundary layer profile
intBlasProf = output(1, :);     %Stream function
vBlasProf = eta.*output(2, :) - output(1, :);   %calculates y velocity component (it is not scaled coorrectly here)

%Calculates f, df/d\eta, d^2f/d\eta^2
function [hitnew, output] = BlasiusProfile(f1, f2, f3, eta, output)

    for i = 1:length(eta)-1
        
        k = eta(i+1)-eta(i);
        LHS = [f1, f2, f3];
        RHS = [f2, f3, -0.5*f1*f3];

        w = RK4Blaus(LHS, RHS, k);

        f1 = w(1);
        f2 = w(2);
        f3 = w(3);

        hitnew = f2;
        
        output(:,i+1) = [f1; f2; f3];

    end
end

%fourth order Runge Kutta method
function w = RK4Blaus(LHS, RHS, k)
    
    k1 = RHS;
    k2 = RHS + 0.5*k*k1;
    k3 = RHS + 0.5*k*k2;
    k4 = RHS + k*k3;
    
    w = LHS + k*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end
end
