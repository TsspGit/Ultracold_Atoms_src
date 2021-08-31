function [sh]=debvpshooting2f1(ar,ai,br,bi,cr,ci,xend,np)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the Gauss hypergeometric function 2F1(a,b;c;z) on the range    %
% [-1,xend], 0<xend<1, using the shooting method applied to the           %
% differential equation (4.2).                                                                  %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         xend=Right-hand end point of interval of computation            %
%         np=Number of mesh points used                                   %
% Output: sh=Profile of 2F1(a,b;c;z) on [-1,xend]                         %
%-------------------------------------------------------------------------%

% Compute a,b,c in terms of ar,ai,br,bi,cr,ci
a=ar+ai*1i;
b=br+bi*1i;
c=cr+ci*1i;

% Define coefficients of differential operators
A1 = @(x,y) x*(1-x) ;
A2 = @(x,y) c-(a+b+1)*x ;
A3 = @(x,y) -a*b ;

f1 = @(x,y,z) z ;
f2 = @(x,y,z) -1/A1(x,y)*(A2(x,y)*z+A3(x,y)*y) ;

% BCs for [-1,0]
alpha=hypergeom([a,b],c,-1);
beta=hypergeom([a,b],c,0);

x0=-1 ;
y0=alpha ;

xf=0 ;
yf=beta ;

dydx1 = 0 ;
dydx2 = 1 ;

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

% Apply shooting method formula on [-1,0]
y3=y1+((beta-y1(end))/y2(end))*y2;

%%

% BCs on [0,xend]
alpha=hypergeom([a,b],c,xend);
beta=hypergeom([a,b],c,0);

x0=0 ;
y0=alpha ;

xf=xend ;
yf=beta ;

% Define coefficients of differential operators
A1 = @(x,y) (xend-x)*(1-xend+x) ;
A2 = @(x,y) c-(a+b+1)*(xend-x) ;
A3 = @(x,y) -a*b ;

f1 = @(x,y,z) z ;
f2 = @(x,y,z) -1/A1(x,y)*(A2(x,y)*z+A3(x,y)*y) ;

dydx1 = 0 ;
dydx2 = 1 ;

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

% Apply shooting method formula on [0,xend]
y3new=y1+((beta-y1(end))/y2(end))*y2;

x1=X(1:end-1);
y31=y3new(1:end-1);

%%

% Piece together solution on [-1,xend]
for l=1:length(x1)
    x1new(l)=x1(length(x1)-l+1);
    y31new(l)=y31(length(x1)-l+1);
end

%%
x2=[x,x1new];
y32=[y3,y31new];
plot(x2,y32,'r')
%hold on
%plot(x2,hypergeom(a,b,x2))
%hold off

sh=y32;
%e1 = norm(hypergeom(a,b,x2)-y32,1);
%e2 = norm(hypergeom(a,b,x2)-y32,2);
%einf = norm(hypergeom(a,b,x2)-y32,inf);