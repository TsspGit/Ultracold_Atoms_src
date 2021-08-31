function [ivp]=deivpdp1f1(ar,ai,br,bi,xf,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using         %
% the Dormand-Prince method to solve the initial value problem (3.2), as  %
% in Appendix G.6.                                                        %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         xf=Point where computation is desired                           %
%         n=Number of mesh points used                                    %
% Output: ivp=Computed value of 1F1(a;b;z)                                %
%-------------------------------------------------------------------------%

% Compute a,b in terms of ar,ai,br,bi
a=ar+ai*1i;
b=br+bi*1i;

% Define coefficients of differential operators
x = linspace(0,xf,n)';

A1 = @(x,y) x ;
A2 = @(x,y) b-x ;
A3 = @(x,y) -a ;

f1 = @(x,y,z) z ;
f2 = @(x,y,z) -1/A1(x,y)*(A2(x,y)*z+A3(x,y)*y) ;

alpha=1;x0=0 ;
y0=alpha ;

dydx1 = a/b ;

% Define mesh points
h=(xf-x0)/n ;

for i = 1:n+1
    x(i)=x0+(i-1)*h;
end

% Define ICs using another method or a Taylor series expansion
y1(1)=y0;
z1(1)=dydx1;
y1(2)=taylora1f1(ar,ai,br,bi,h,0,1e-15);
z1(2)=a/b*taylora1f1(ar+1,ai,br+1,bi,h,0,1e-15);

% RK4 method
for i = 2:n
    k1y=f1(x(i),y1(i),z1(i));
    k1z=f2(x(i),y1(i),z1(i));
    k2y=f1(x(i)+h/5,y1(i)+k1y*h/5,z1(i)+k1z*h/5);
    k2z=f2(x(i)+h/5,y1(i)+k1y*h/5,z1(i)+k1z*h/5);
    k3y=f1(x(i)+3/10*h,y1(i)+3/40*k2y*h,z1(i)+3/40*k2z*h);
    k3z=f2(x(i)+3/10*h,y1(i)+3/40*k2y*h,z1(i)+3/40*k2z*h);
    k4y=f1(x(i)+4/5*h,y1(i)+44/45*k1y*h-56/15*k2y*h+32/9*k3y*h,...
        z1(i)+44/45*k1z*h-56/15*k2z*h+32/9*k3z*h);
    k4z=f2(x(i)+4/5*h,y1(i)+44/45*k1y*h-56/15*k2y*h+32/9*k3y*h,...
        z1(i)+44/45*k1z*h-56/15*k2z*h+32/9*k3z*h);
    k5y=f1(x(i)+8/9*h,y1(i)+19372/6561*k1y*h-25360/2187*k2y*h...
        +64448/6561*k3y*h-212/729*k4y*h,z1(i)+19372/6561*k1z*h...
        -25360/2187*k2z*h+64448/6561*k3z*h-212/729*k4z*h);
    k5z=f2(x(i)+8/9*h,y1(i)+19372/6561*k1y*h-25360/2187*k2y*h...
        +64448/6561*k3y*h-212/729*k4y*h,z1(i)+19372/6561*k1z*h...
        -25360/2187*k2z*h+64448/6561*k3z*h-212/729*k4z*h);
    k6y=f1(x(i)+h,y1(i)+9017/3168*k1y*h-355/33*k2y*h+46732/5247*k3y*h...
        +49/176*k4y*h-5103/18656*k5y*h,z1(i)+9017/3168*k1z*h...
        -355/33*k2z*h+46732/5247*k3z*h+49/176*k4z*h-5103/18656*k5z*h);
    k6z=f2(x(i)+h,y1(i)+9017/3168*k1y*h-355/33*k2y*h+46732/5247*k3y*h...
        +49/176*k4y*h-5103/18656*k5y*h,z1(i)+9017/3168*k1z*h...
        -355/33*k2z*h+46732/5247*k3z*h+49/176*k4z*h-5103/18656*k5z*h);
    y1(i+1)=y1(i)+h*(35/384*k1y+500/1113*k3y+125/192*k4y...
        -2187/6784*k5y+11/84*k6y);
    z1(i+1)=z1(i)+h*(35/384*k1z+500/1113*k3z+125/192*k4z...
        -2187/6784*k5z+11/84*k6z);
end

% Return solution
ivp=y1(end);