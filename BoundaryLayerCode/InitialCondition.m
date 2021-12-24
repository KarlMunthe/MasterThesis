function [u, v] = InitialCondition(x, y, U, nu)

%make matrix where u and v are to be saved
u = zeros(length(y), length(x));
v = zeros(length(y), length(x));

delta = sqrt(nu.*x./U);
gridVal = y./delta;

eta = linspace(0, max(max(gridVal)), length(y)*10);

%calculate various derivatives of the Blasius solution
[BlasProf, intBlasProf, vBlasProf] = Blasius(eta, 1, 1e-15);

for col = 1:size(gridVal, 1)
    for row = 1:size(gridVal, 2)
        
        [~, closest_index] = min(abs(gridVal(col, row) - eta));
        u(col, row) = BlasProf(closest_index);
        
        [~, closest_index] = min(abs(gridVal(col, row) - eta));
        v(col, row) = vBlasProf(closest_index);
        
    end
end

u = U.*u;
v = 0.5*nu./delta.*v;