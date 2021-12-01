clear

%Standard Temperature and Pressure
T0 = 273.15;
p0 = 10^5;

%fluid properties for argon
%
mu = 22.61e-1;
%kappa = 0.0178;
c_V = 313;
c_p = 520;
R = c_p - c_V;
gamma = c_p/c_V;
%}

%fluid properties for a nice fluid
%{
mu = 1e-1;
%kappa = 0.0178;
c_V = 300;
c_p = 500;
R = c_p - c_V;
gamma = c_p/c_V;
%}

%set up time
dt = 1e-7;
tend = 1e+6*dt;
time = dt:dt:tend;

%set up grid
Nx = 64;
xstart = 10;
xend = 12;
x = linspace(xstart, xend, Nx);
%}

Ny = 128;
ystart = 0;
yend = 256;
y = linspace(ystart, yend, Ny);

%ytrans = y';
%ytrans = ((ystart*exp(yend) - yend*exp(ystart) + (yend-ystart)*exp(y))/(exp(yend)-exp(ystart)))';
q = 25;
ytrans = ((ystart*exp(yend/q) - yend*exp(ystart/q) + (yend-ystart)*exp(y/q))/(exp(yend/q)-exp(ystart/q)))';
hy = ytrans(2:end)-ytrans(1:end-1);
yside = 0.5*(ytrans(3:end) - ytrans(1:end-2));
nyside = 0.5*(ytrans(2) - ytrans(1));
syside = 0.5*(ytrans(end) - ytrans(end - 1));
yside = [nyside; yside; syside];

hx = x(2:end)-x(1:end-1);
xside = 0.5*(x(3:end) - x(1:end-2));
wxside = 0.5*(x(2) - x(1));
exside = 0.5*(x(end) - x(end-1));
xside = [wxside, xside, exside];

%dy = 1;
v = (yside*xside);
min(min(v))

[BlasProfField, intBlasProfField] = InitialCondition2(x, ytrans);
U = 0.001*ones(length(y),1);
BLthick = sqrt(mu*x./U);
%initial conditions
%exp(-100*((x-0.5).^2 + (ytrans-0.5).^2)) + 1;
%ones(length(y), length(x));
m0 = BlasProfField;
n0 = 0.5*mu./BLthick.*(BlasProfField - ytrans.*intBlasProfField);
rho0 = p0./(R*T0)*ones(size(v));
nu0 = mu./rho0;
%p0 = rho0.^gamma;
%T0 = p0./(rho0*R);
e0 = c_V*T0;
E0 = e0.*rho0 + 0.5*(m0.^2 + n0.^2)./rho0;
%}
rho = rho0;
m = m0;
n = n0;
E = E0;
p = (gamma-1)*(E - 0.5*(m0.^2 + n0.^2)./rho);

vars = zeros(size(v, 1), size(v, 2), 4);
vars(:, :, 1) = rho;
vars(:, :, 2) = m;
vars(:, :, 3) = n;
vars(:, :, 4) = E;

wBC = [rho(:, 1), m(:, 1), n(:, 1), E(:, 1)];   %Standard BC
eBC = [rho(:, end), m(:, end), n(:, end), E(:, end)];   %Standard BC
nBC = [rho(1, :); m(1, :); n(1, :); E(1, :)];   %Standard BC
sBC = zeros(4, size(v, 1));     %Flux BC

wvars = vars(:, 1, :);      %Injection BC
evars = vars(:, end, :);    %Injection BC
nvars = vars(1, :, 2:3);
svars = vars(end, :, 2:3);
%[wBC, eBC, nBC, sBC] = BC(rho0, m0, n0, E0, x, y);
    
for t = 1:length(time)
    
    
    vars = RK4(vars, dt, v, hx, hy, xside, yside, gamma, mu, wvars, evars, nvars, svars, wBC, eBC, nBC, sBC, rho0);
   
    
    if mod(t, 500) == 499
        testvars = vars(:, :, 1);
    end
    if mod(t, 500) == 0
        
        %TT = (vars(:, :, 4) - (vars(:, :, 2).^2 + vars(:, :, 3).^2)./(2*vars(:, :, 1)))./(vars(:, :, 1).*c_V);
        l2error = sqrt(sum(sum(v.*(vars(:, :, 1) - testvars).^2)))
        plot(ytrans, m0(:, end/2))
        hold on
        plot(ytrans, vars(:, end/2, 2))
        hold off
        %mesh(x, ytrans, vars(:, :, 1))
        %{
        tiledlayout(2, 1)
        nexttile
        plot(ytrans, vars(:, end/2, 2))
        hold on
        plot(ytrans, m0(:, end/2))
        xlim([0, 10])
        ylim([0, 1.5])
        hold off
        nexttile
        plot(ytrans, vars(:, end/2, 2))
        hold on
        plot(ytrans, m0(:, end/2))
        xlim([90, 100])
        ylim([0, 1.5])
        %}
        drawnow
    end
    %}
    
end

%error = vars(:, :, 1) - rho0(x, tend, ytrans);
%l2error = sqrt(sum(sum(v.*error.*error)))
mesh(x, ytrans, rho0)
xlabel('x')
ylabel('y')
f2 = figure;
mesh(x, ytrans, vars(:, :, 1))
xlabel('x')
ylabel('y')

function w = RK4(vars, dt, v, hx, hy, xside, yside, gamma, mu, wvars, evars, nvars, svars, wBC, eBC, nBC, sBC, rho0)
    
    vars(:, 1, :) = wvars;
    vars(:, end, :) = evars;
    vars(1, :, 2:3) = nvars;
    %vars(end, :, 2:3) = svars;
    
    k1 = RHS(vars, v, hx, hy, xside, yside, gamma, mu, wBC, eBC, nBC, sBC);
    
    vars(:, 1, :) = wvars;% + 0.5*dt*k1(:, 1, :);
    vars(:, end, :) = evars;% + 0.5*dt*k1(:, end, :);
    vars(1, :, 2:3) = nvars;% + 0.5*dt*k1(end, :, 2:3);
    %vars(end, :, 2:3) = svars;

    k2 = RHS(vars + 0.5*dt*k1, v, hx, hy, xside, yside, gamma, mu, wBC, eBC, nBC, sBC);
    
    vars(:, 1, :) = wvars;% + 0.5*dt*k1(:, 1, :);
    vars(:, end, :) = evars;% + 0.5*dt*k1(:, end, :);
    vars(1, :, 2:3) = nvars;% + 0.5*dt*k1(end, :, 2:3);
    %vars(end, :, 2:3) = svars;
    
    k3 = RHS(vars + 0.5*dt*k2, v, hx, hy, xside, yside, gamma, mu, wBC, eBC, nBC, sBC);
    
    vars(:, 1, :) = wvars;% + dt*k1(:, 1, :);
    vars(:, end, :) = evars;% + dt*k1(:, end, :);
    vars(1, :, 2:3) = nvars;% + dt*k1(end, :, 2:3);
    %vars(end, :, 2:3) = svars;
    
    k4 = RHS(vars + dt*k3, v, hx, hy, xside, yside, gamma, mu, wBC, eBC, nBC, sBC);
    
    w = vars + dt./v.*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

function RHS = RHS(vars, v, hx, hy, xside, yside, gamma, mu, wBC, eBC, nBC, sBC)

    rho = vars(:, :, 1);
    m = vars(:, :, 2);
    n = vars(:, :, 3);
    E = vars(:, :, 4);
    nu = mu./rho;
    
    [F1, F2, F3, F4] = F(rho, m, n, E, gamma);
    [G1, G2, G3, G4] = G(rho, m, n, E, gamma);
    
    F1flux = xflux(F1, yside, wBC(:, 1), eBC(:, 1));
    F2flux = xflux(F2, yside, wBC(:, 2), eBC(:, 2));
    F3flux = xflux(F3, yside, wBC(:, 3), eBC(:, 3));
    F4flux = xflux(F4, yside, wBC(:, 4), eBC(:, 4));
    
    G1flux = yflux(G1, xside, nBC(1, :), sBC(1, :));
    G2flux = yflux(G2, xside, nBC(2, :), sBC(2, :));
    G3flux = yflux(G3, xside, nBC(3, :), sBC(3, :));
    G4flux = yflux(G4, xside, nBC(4, :), sBC(4, :));
    
    rhoxdiff = xdiff(rho, hx, nu, yside, wBC(:, 1), eBC(:, 1));
    rhoydiff = ydiff(rho, hy, nu, xside, nBC(1, :), sBC(1, :));
    mxdiff = xdiff(m, hx, nu, yside, wBC(:, 2), eBC(:, 2));
    mydiff = ydiff(m, hy, nu, xside, nBC(2, :), sBC(2, :));
    nxdiff = xdiff(n, hx, nu, yside, wBC(:, 3), eBC(:, 3));
    nydiff = ydiff(n, hy, nu, xside, nBC(3, :), sBC(3, :));
    Exdiff = xdiff(E, hx, nu, yside, wBC(:, 4), eBC(:, 4));
    Eydiff = ydiff(E, hy, nu, xside, nBC(4, :), sBC(4, :));

    RHS = zeros(size(v, 1), size(v,  2), 4);
    %{
    RHS(:, :, 1) = F1flux + G1flux; 
    RHS(:, :, 2) = F2flux + G2flux;
    RHS(:, :, 3) = F3flux + G3flux;
    RHS(:, :, 4) = F4flux + G4flux;
    %}
    %{
    RHS(:, :, 1) = rhoxdiff + rhoydiff;
    RHS(:, :, 2) = mxdiff + mydiff;
    RHS(:, :, 3) = nxdiff + nydiff;
    RHS(:, :, 4) = Exdiff + Eydiff;
    %}
    
    RHS(:, :, 1) = F1flux + G1flux + rhoxdiff + rhoydiff; 
    RHS(:, :, 2) = F2flux + G2flux + mxdiff + mydiff;
    RHS(:, :, 3) = F3flux + G3flux + nxdiff + nydiff;
    RHS(:, :, 4) = F4flux + G4flux + Exdiff + Eydiff;
    %}
    
end   

function [F1, F2, F3, F4] = F(rho, m, n, E, gamma)

    p = (gamma - 1)*(E - 0.5*(m.^2 + n.^2)./rho);

    F1 = m; 
    F2 = m.^2./rho + p;
    F3 = m.*n./rho;
    F4 = (E + p).*m./rho;
       
end

function [G1, G2, G3, G4] = G(rho, m, n, E, gamma)

    p = (gamma - 1)*(E - 0.5*(m.^2 + n.^2)./rho);
    
    G1 = n;
    G2 = m.*n./rho;
    G3 = n.^2./rho + p;
    G4 = (E + p).*n./rho;
       
end

function xflux = xflux(f, dy, wBC, eBC)

    xflux = 0.5*(f(:, 1:end-1) + f(:, 2:end));
    %xflux = [wBC, xflux, eBC];
    %xflux = [zeros(size(xflux, 1), 1), xflux, zeros(size(xflux, 1), 1)];
    %xflux = [f(:, 1), xflux, f(:, end)];
    xflux = [xflux(:, 1), xflux, xflux(:, end)];

    
    xflux = -(xflux(:, 2:end) - xflux(:, 1:end-1)).*dy;

end
%}


function yflux = yflux(f, dx, nBC, sBC)
    
    yflux = 0.5*(f(1:end-1, :) + f(2:end, :));
    %yflux = [zeros(1, size(yflux, 2));  yflux; zeros(1, size(yflux, 2))];
    %yflux = [nBC; yflux; sBC];
    %yflux = [f(1, :); yflux; f(end, :)];
    yflux = [yflux(1, :); yflux; yflux(end, :)];
    
    yflux = -(yflux(2:end, :) - yflux(1:end-1, :)).*dx;
    
end

function xdiff = xdiff(u, hx, nu, yside, wBC, eBC)
    
    nux = 0.5*(nu(:, 2:end) + nu(:, 1:end-1));
    nux = [nu(:, 1), nux, nu(:, end)];
    %nux = [zeros(size(nux, 1), 1), nux, zeros(size(nux, 1), 1)];
    ux = (u(:,  2:end) - u(:, 1:end-1))./hx;
    %ux = [zeros(size(ux, 1), 1), ux, zeros(size(ux, 1), 1)];
    %ux = [wBC, ux, eBC];
    %ux = [ux(:, end), ux, ux(:, 1)];
    ux = [ux(:, 1), ux, ux(:, end)];
    
    xdiff = yside.*(nux(:, 2:end).*ux(:, 2:end) - nux(:, 1:end-1).*ux(:, 1:end-1));
    
end

function ydiff = ydiff(u, hy, nu, xside, nBC, sBC)
    
    nuy = 0.5*(nu(2:end, :) + nu(1:end-1, :));
    %nuy = [nu(1, :); nuy; nu(end, :)];
    nuy = [zeros(1, size(nuy, 2)); nuy; zeros(1, size(nuy, 2))];
    uy = (u(2:end, :) - u(1:end-1, :))./hy;
    uy = [zeros(1, size(uy, 2));  uy; zeros(1, size(uy, 2))];
    %uy = [uy(end, :); uy; uy(1, :)];
    %uy = [nBC; uy; sBC];
    %uy = [uy(1, :); uy; uy(end, :)];
    
    ydiff = xside.*(nuy(2:end, :).*(uy(2:end, :) - nuy(1:end-1, :).*uy(1:end-1, :)));

end

function northBC = northBC(oldnBC, rho_nBC, E_nBC)

    %rho_nBC = vars(1, :, 1);
    %E_nBC = vars(1, :, 4);
    oldnBC(1, :) = rho_nBC;
    oldnBC(4, :) = E_nBC;
    northBC = oldnBC;

end

function [wBC, eBC, nBC, sBC] = BC(rho0, m0, n0, E0, x, y)

    wBC = zeros(length(x), 4);
    wBC(:, 1) = rho0(:, 1);
    wBC(:, 2) = m0(:, 1);
    wBC(:, 3) = n0(:, 1);
    wBC(:, 4) = E0(:, 1);
    
    eBC = zeros(length(x), 4);
    eBC(:, 1) = rho0(:, end);
    eBC(:, 2) = m0(:, end);
    eBC(:, 3) = n0(:, end);
    eBC(:, 4) = E0(:, end);
    
    nBC = zeros(4, length(y));
    nBC(1, :) = rho0(1, :);
    nBC(2, :) = m0(1, :);
    nBC(3, :) = n0(1, :);
    nBC(4, :) = E0(1, :);

    sBC = zeros(4, length(y));
    sBC(1, :) = zeros(1, length(y));
    sBC(2, :) = zeros(1, length(y));
    sBC(3, :) = zeros(1, length(y));
    sBC(4, :) = zeros(1, length(y));

end

    