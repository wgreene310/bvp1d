function stringVibration
T=10; % tension in the string
rho=.2; % mass/length of the string
L=20; % length of the string
analOmega = omegaAnal(L, T, rho);
n=5;
solinit.x=linspace(0,L,n);
solinit.y=[ones(1,n); zeros(1,n)];
solinit.parameters = 0;
odeFunc = @(x, u, lambda) stringODE(x, u, lambda,T,rho);
sol=bvp1d(odeFunc, @stringBC, solinit);
figure; plot(sol.x, sol.y(1,:), 'x-');
fprintf('Analytical frequency = %7.5f, bvp1d frequency = %7.5f\n',
   analOmega, sqrt(sol.parameters));
end

function dudx=stringODE(x, u, lambda,T,rho)
c2 = T/rho;
dudx = [u(2) -lambda/c2*u(1)]';
end

function g=stringBC(ya, yb, lambda)
g=[ya(1) yb(1) ya(2)-.1]';
end

function omega = omegaAnal(L, T, rho)
% lowest vibration frequency of the string
c = sqrt(T/rho);
omega = c*pi/L;
end