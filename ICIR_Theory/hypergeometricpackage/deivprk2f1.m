function [ivp]=deivprk2f1(ar,ai,br,bi,cr,ci,xf,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the Gauss hypergeometric function 2F1(a,b;c;z) using the RK4   %
% method to solve the initial value problem (3.2), as in Section 4.5.     %                                                                    %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         xf=Point where computation is desired                           %
%         n=Number of mesh points used                                    %
% Output: ivp=Computed value of 2F1(a,b;c;z)                              %
%-------------------------------------------------------------------------%

% Compute a,b,c using ar,ai,br,bi,cr,ci
a=ar+ai*1i;
b=br+bi*1i;
c=cr+ci*1i;

% Define coefficients of differential operators
x = linspace(0,xf,n)';

A1 = @(x,y) x*(1-x) ;
A2 = @(x,y) c-(a+b+1)*x ;
A3 = @(x,y) -a*b ;

f1 = @(x,y,z) z ;
f2 = @(x,y,z) -1/A1(x,y)*(A2(x,y)*z+A3(x,y)*y) ;

alpha=1;
x0=0 ;
y0=alpha ;

dydx1 = a*b/c ;

% Define mesh points
h=(xf-x0)/n ;

for i = 1:n+1
    x(i)=x0+(i-1)*h;
end

% Define ICs using other method or Taylor series expansion
y1(1)=y0;
z1(1)=dydx1;
y1(2)=taylora2f1(ar,ai,br,bi,cr,ci,h,0,1e-15);
z1(2)=a*b/c*taylora2f1(ar+1,ai,br+1,bi,cr+1,ci,h,0,1e-15);

% RK4 method
for i = 2:n
    k1y=f1(x(i),y1(i),z1(i));
    k1z=f2(x(i),y1(i),z1(i));
    k2y=f1(x(i)+0.5*h,y1(i)+0.5*k1y*h,z1(i)+0.5*k1z*h);
    k2z=f2(x(i)+0.5*h,y1(i)+0.5*k1y*h,z1(i)+0.5*k1z*h);
    k3y=f1(x(i)+0.5*h,y1(i)+0.5*k2y*h,z1(i)+0.5*k2z*h);
    k3z=f2(x(i)+0.5*h,y1(i)+0.5*k2y*h,z1(i)+0.5*k2z*h);
    k4y=f1(x(i)+h,y1(i)+k3y*h,z1(i)+k3z*h);
    k4z=f2(x(i)+h,y1(i)+k3y*h,z1(i)+k3z*h);
    y1(i+1)=y1(i)+h/6*(k1y+2*k2y+2*k3y+k4y);
    z1(i+1)=z1(i)+h/6*(k1z+2*k2z+2*k3z+k4z);
end

% Return solution
ivp=y1(end);