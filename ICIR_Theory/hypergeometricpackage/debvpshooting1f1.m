function [sh]=debvpshooting1f1(ar,ai,br,bi,np)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) on the range  %
% [-1,1] using the shooting method applied to the differential equation   %
% (3.2).                                                                  %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         np=Number of mesh points used                                   %
% Output: sh=Profile of 1F1(a;b;z) on [-1,1]                              %
%-------------------------------------------------------------------------%

% Compute a,b in terms of ar,ai,br,bi
a=ar+ai*1i;
b=br+bi*1i;

% Define coefficients of differential operators
A1 = @(x,y) x ;
A2 = @(x,y) b-x ;
A3 = @(x,y) -a ;

f1 = @(x,y,z) z ;
f2 = @(x,y,z) -1/A1(x,y)*(A2(x,y)*z+A3(x,y)*y) ;

% Set BCs on [-1,0]
alpha=hypergeom(a,b,-1);
beta=hypergeom(a,b,0);

x0=-1 ;
y0=alpha ;

xf=0 ;
yf=beta ;

dydx1 = 0 ;
dydx2 = 1 ;

% Set up grid
n=ceil(np/2) ;
h=(xf-x0)/n ;

for i = 1:n+1
    x(i)=x0+(i-1)*h;
end

y1(1)=y0;
z1(1)=dydx1;

% RK4 method
for i = 1:n
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

y2(1)=0;
z2(1)=dydx2;

for i = 1:n
    k1y=f1(x(i),y2(i),z2(i));
    k1z=f2(x(i),y2(i),z2(i));
    k2y=f1(x(i)+0.5*h,y2(i)+0.5*k1y*h,z2(i)+0.5*k1z*h);
    k2z=f2(x(i)+0.5*h,y2(i)+0.5*k1y*h,z2(i)+0.5*k1z*h);
    k3y=f1(x(i)+0.5*h,y2(i)+0.5*k2y*h,z2(i)+0.5*k2z*h);
    k3z=f2(x(i)+0.5*h,y2(i)+0.5*k2y*h,z2(i)+0.5*k2z*h);
    k4y=f1(x(i)+h,y2(i)+k3y*h,z2(i)+k3z*h);
    k4z=f2(x(i)+h,y2(i)+k3y*h,z2(i)+k3z*h);
    y2(i+1)=y2(i)+h/6*(k1y+2*k2y+2*k3y+k4y);
    z2(i+1)=z2(i)+h/6*(k1z+2*k2z+2*k3z+k4z);
end

% Apply shooting method formaula on [-1,0]
y3=y1+((beta-y1(end))/y2(end))*y2;

%%

% BCs on [0,1]
alpha=hypergeom(a,b,1);
beta=hypergeom(a,b,0);

x0=0 ;
y0=alpha ;

xf=1 ;
yf=beta ;

% Define differential operators with transformation as described in [73]
A1 = @(x,y) (1-x) ;
A2 = @(x,y) -(b-(1-x)) ;
A3 = @(x,y) -a ;

f1 = @(x,y,z) z ;
f2 = @(x,y,z) -1/A1(x,y)*(A2(x,y)*z+A3(x,y)*y) ;

dydx1 = 0 ;
dydx2 = 1 ;

% Set up grid
n=ceil(np/2) ;
h=(xf-x0)/n ;

for i = 1:n+1
    x3(i)=x0+(i-1)*h;
    X(i)=xf-x3(i);
end

y1(1)=y0;
z1(1)=dydx1;

% RK4 method
for i = 1:n
    k1y=f1(x3(i),y1(i),z1(i));
    k1z=f2(x3(i),y1(i),z1(i));
    k2y=f1(x3(i)+0.5*h,y1(i)+0.5*k1y*h,z1(i)+0.5*k1z*h);
    k2z=f2(x3(i)+0.5*h,y1(i)+0.5*k1y*h,z1(i)+0.5*k1z*h);
    k3y=f1(x3(i)+0.5*h,y1(i)+0.5*k2y*h,z1(i)+0.5*k2z*h);
    k3z=f2(x3(i)+0.5*h,y1(i)+0.5*k2y*h,z1(i)+0.5*k2z*h);
    k4y=f1(x3(i)+h,y1(i)+k3y*h,z1(i)+k3z*h);
    k4z=f2(x3(i)+h,y1(i)+k3y*h,z1(i)+k3z*h);
    y1(i+1)=y1(i)+h/6*(k1y+2*k2y+2*k3y+k4y);
    z1(i+1)=z1(i)+h/6*(k1z+2*k2z+2*k3z+k4z);
end  

y2(1)=0;
z2(1)=dydx2;

for i = 1:n
    k1y=f1(x3(i),y2(i),z2(i));
    k1z=f2(x3(i),y2(i),z2(i));
    k2y=f1(x3(i)+0.5*h,y2(i)+0.5*k1y*h,z2(i)+0.5*k1z*h);
    k2z=f2(x3(i)+0.5*h,y2(i)+0.5*k1y*h,z2(i)+0.5*k1z*h);
    k3y=f1(x3(i)+0.5*h,y2(i)+0.5*k2y*h,z2(i)+0.5*k2z*h);
    k3z=f2(x3(i)+0.5*h,y2(i)+0.5*k2y*h,z2(i)+0.5*k2z*h);
    k4y=f1(x3(i)+h,y2(i)+k3y*h,z2(i)+k3z*h);
    k4z=f2(x3(i)+h,y2(i)+k3y*h,z2(i)+k3z*h);
    y2(i+1)=y2(i)+h/6*(k1y+2*k2y+2*k3y+k4y);
    z2(i+1)=z2(i)+h/6*(k1z+2*k2z+2*k3z+k4z);
end

% Apply shooting method rule on [0,1]
y3new=y1+((beta-y1(end))/y2(end))*y2;

x1=X(1:end-1);
y31=y3new(1:end-1);

%%

% Piece together complete solution
for l=1:length(x1)
    x1new(l)=x1(length(x1)-l+1);
    y31new(l)=y31(length(x1)-l+1);
end

%%
x2=[x,x1new];
y32=[y3,y31new];
plot(x2,y32,'r')

sh=y32;