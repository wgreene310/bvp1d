% bvp1d Numerical solution of systems of boundary value differential equations with a single
% independent variable. Typically the single independent variable is a spatial
% dimension and boundary conditions are prescribed at both ends of the interval.
%
% Usage:
% solution = bvp1d(bvpFunc, bcFunc, initSolution)
% where bvpFunc and bcFunc are function handles to user-defined
% functions defining the differential equations and boundary conditions
% for the boundary value problem. initSolution is a structure that contains
% fields defining the initial spatial mesh and a guess for the solution.
% solution is a structure that contains several fields defining the solution to
% the problem.
%
% bvp1d expects problems to be defined in the following form: y' = f(x, y) where the
% independent spatial variable x ranges from x=a to x=b. The dependent variable, 
% y, is the solution to be obtained (a vector of dependent variables when more than one
% equation is included in the system). The "prime" symbol (') indicates the derivative
% of y with respect to x. The boundary conditions are defined as
% bc(y(a), y(b))=0. The number of boundary conditions (i.e. the number of rows
% in the bc variable) must equal the number of differential equations in the system.
%
% initSolution is a structure than must contain the following fields: the initial
% mesh, x, and a guess of the solution on that mesh, y. The solution structure
% returned by bvp1d contains the following fields: x, y, and yp. bvp1d refines
% the mesh so that the solution meets a prescribed accuracy so the points in
% the returned x field typically are different from the initial x variable
% defined in the initSolution structure. The field y is the solution defined
% on this mesh and the yp field is the derivative of this solution with
% respect to x (i.e. x').
%
% The usage of bvp1d is best explained with a concrete example.
% Example:
% Solve u'' + |u| = 0
% where '' indicates the second derivative of y with respect to x and || indicates
% absolute value. The range of x is from zero to four. The boundary conditions at
% x=0 are u=0 and at x=4 are u=-2. To convert the differential equation into the
% form required by bvp1d, a new dependent variable and differential equation are
% defined: v = u'. The original equation can be rewritten as two first order
% equations:
%     u' = v
%     v' = -|u|
% The variable, y, is a vector containing u and v. The function, bvpFunc, is
% defined as:
%    function dydx = example1Func(x, y)
%      dydx = [y(2); -abs(y(1))];
%    end
% The boundary condition function, bcFunc, is defined as
%    function bc=example1BCFunc(ya, yb)
%      bc = [ya(1); yb(2)+2];
%    end
% We will define five evenly-spaced points in the initial mesh.
%    initSolution.x = linspace(0,4,5);
% This differential equation has two solutions. The particular solution returned
% by bvp1d depends on the initial guess. The initial guess
%    initSolution.y = [ones(1,5); zeros(1,5)];
% returns one solution and the guess
%    initSolution.y = [-ones(1,5); zeros(1,5)];
% returns the second solution. bvp1d is called as follows:
%    sol = bvp1d(@example1Func, @example1BCFunc, initSolution);
% An x-y plot of the solution can be produced with the following statement.
%    plot(sol.x, sol.y(1,:));
%
% In addition to solving differential equations like those described above, bvp1d
% has the capability to solve differential equations that include one or more
% unknown "parameters" to be determined. Generally, both the differential equations
% and the boundary conditions depend on these unknown parameters and the number
% of "boundary conditions" returned by the bcFunc must equal the number of
% differential equations plus the number of parameters. This capability can
% best be explained by an example.
%
% Example: Vibration of a String With Tension
% In this example we compute the vibration frequency and the mode shape for
% a string stretched between two points that has been pre-tensioned. The equation
% describing this behavior is:
%    u'' + lambda*rho/T u = 0
% where u is the lateral deflection of the string, rho is the mass per unit
% length, T is the tension force in the string, and lambda is the square of the
% vibration frequency. Both the deflection, u, and the parameter lambda are
% unknowns to be determined. The boundary conditions are u=0 at both ends of
% the string. Because there is a single parameter in this system, one additional
% constraint must be imposed. We, somewhat arbitrarily, will require the slope
% of the mode shape at the left end of the string to equal 0.1. The input for
% this example is:
%    function stringVibration
%    T=10; % tension in the string
%    rho=.2; % mass/length of the string
%    L=20; % length of the string
%    analOmega = omegaAnal(L, T, rho);
%    n=5;
%    solinit.x=linspace(0,L,n);
%    solinit.y=[ones(1,n); zeros(1,n)];
%    solinit.parameters = 0;
%    odeFunc = @(x, u, lambda) stringODE(x, u, lambda,T,rho);
%    sol=bvp1d(odeFunc, @stringBC, solinit);
%    figure; plot(sol.x, sol.y(1,:), 'x-');
%    fprintf('Analytical frequency = %7.5f, bvp1d frequency = %7.5f\n',
%         analOmega, sqrt(sol.parameters));
%    end
% 
%    function dudx=stringODE(x, u, lambda,T,rho)
%    c2 = T/rho;
%    dudx = [u(2) -lambda/c2*u(1)]';
%    end
%    
%    function g=stringBC(ya, yb, lambda)
%    g=[ya(1) yb(1) ya(2)-.1]';
%    end
%    
%    function omega = omegaAnal(L, T, rho)
%    % lowest vibration frequency of the string
%    c = sqrt(T/rho);
%    omega = c*pi/L;
%    end
%
% As can be seen, the parameters (in this case lambda) are passed to the
% BVP and boundary functions as a third argument. In addition, the
% initSolution structure must contain a field, parameters, with an initial
% guess for the value of lambda. Five points are used in the initial mesh
% but the mesh is refined by bvp1d to satisfy the accuracy requirements.
%
% An optional fourth argument may be passed to bvp1d. This argument is a
% structure that contains fields that affect the behavior of bvp1d. The
% currently supported options are:
%    AbsTol - absolute accuracy tolerance (default=1e-6)
%    RelTol - relative  accuracy tolerance (default=1e-3)
%             (how AbsTol and RelTol are used is described below)
%    NMax   - bvp1d will not increase the number mesh points to be larger than this value
%    Stats  - if set to 'on', statistics about the solution process are printed
%
% Solution Accuracy
% On each subinterval of the mesh, a residual error, r(i), is computed. This
% residual error is scaled by the right hand side, f(i), so that
%    r_scaled(i) = r(i)/max(f(i),AbsTol/RelTol)
% The Euclidean norm of r_scaled is required to be <= RelTol.
