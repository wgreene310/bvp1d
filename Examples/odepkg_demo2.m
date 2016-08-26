function odepkg_demo2
 a            = 0;
 b            = 4;
 Nint         = 2;
 Nvar         = 2;
 s            = 3;
 t            = linspace(a,b,Nint+1);
 u_1          = -ones(1, Nint+1);
 u_2          = 0*u_1;
 u_0          = [u_1 ; u_2];
 f            = @(t,u) [ u(2); -abs(u(1)) ];
 g            = @(ya,yb) [ya(1); yb(1)+2];
 solinit.x = t; solinit.y=u_0;
 sol = bvp1d(f,g,solinit);
 plot (sol.x,sol.y,'x-')
end