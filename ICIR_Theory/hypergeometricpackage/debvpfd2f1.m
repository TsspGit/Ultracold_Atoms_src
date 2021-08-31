function [fd]=debvpfd2f1(ar,ai,br,bi,cr,ci,alpha,beta,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the Gauss hypergeometric function 2F1(a,b;c;z) on the range    %
% [alpha,beta] using a finite difference method of solving the boundary   %
% value problem from the differential equation (4.2).                     %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         alpha=Left point of solution interval                           %
%         beta=Right point of solution interval                           %
%         n =Number of mesh points used                                   %
% Output: fd=Profile of 2F1(a,b;c;z) on [alpha,beta]                      %
%-------------------------------------------------------------------------%

% Compute a,b,c in terms of ar,ai,br,bi,cr,ci
a=ar+ai*1i;
b=br+bi*1i;
c=cr+ci*1i;

% Set up grid
x = linspace(alpha,beta,n)';
h = x(2)-x(1);

% Define differentiation matrices
v1 = ones(n,1);
D2 = (diag(-2*v1) + diag(v1(2:end),1)+diag(v1(2:end),-1))/h^2;
D1 = (-diag(v1(1:end-1),-1) + diag(v1(2:end),1))/(2*h);

% Set up matrices corresponding to coefficients of derivatives in (4.2)
B=zeros(n,n);
for j=1:n
    B(:,j)=x;
end
C=zeros(n,n);
for j=1:n
    C(:,j)=x.*(1-x);
end

% Define matrix corresponding to differential equation
Mfd = C.*D2+(c-(a+b+1)*B).*D1-a*b*eye(n);
Mfd(1,:) = [1 zeros(1,n-1)];
Mfd(n,:) = [zeros(1,n-1) 1];

% Define RHS matrix corresponding to boundary conditions
d = zeros(n,1);
d(1)=hypergeom([a,b],c,alpha);
d(end)=hypergeom([a,b],c,beta);

% Solve system to obtain solution
fd = Mfd\d;
plot(x,fd)