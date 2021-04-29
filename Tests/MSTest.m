clear; clf;

N   = 200;
CFL = 1;

T   = 1;

x=linspace(0,1,N)';

x=x(1:end-1);

h=x(2)-x(1)

% remember to change error and final time (T=1). 
%udisc = 0*x;
%udisc(floor(N/2):end) = 1 ;

u=sin(2*pi*x(1:end));

%u=udisc


k1=0*u;
k2=0*u;
k3=0*u;
k4=0*u;

%% time step
k = CFL*h;

nsteps = ceil(T/k);

k= T/nsteps;

%% Diff op

D=diag(ones(N-2,1),1)-diag(ones(N-2,1),-1);

D(1,end)=-1;
D(end,1)= 1;

D=D/(2*h);


%% Check derivative:

%norm(D*u-2*pi*cos(2*pi*x(1:end)))*h

%pause



for ii=1:nsteps

  % update solution
 
  k1=D*u;
  k2=D*(u+k/2*k1);
  k3=D*(u+k/2*k2);
  k4=D*(u+k  *k3);
  
  u = u + k/6*(k1 + 2*k2 +2*k3 + k4);

  
  t=ii*k;

end

error = u-sin(2*pi*(x(1:end)+t));
%error = u-udisc;

error = h*norm(error)

plot(x,u);
hold on
plot(x, sin(2*pi*(x(1:end)+t)))