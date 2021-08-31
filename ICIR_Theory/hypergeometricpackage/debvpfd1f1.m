function [fd]=debvpfd1f1(ar,ai,br,bi,alpha,beta,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) on the range  %
% [alpha,beta] using a finite difference method of solving the boundary   %
% value problem from the differential equation (3.2).                     %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         alpha=Left point of solution interval                           %
%         beta=Right point of solution interval                           %
%         n =Number of mesh points used                                   %
% Output: fd=Profile of 1F1(a;b;z) on [alpha,beta]                        %
%-------------------------------------------------------------------------%

% Compute a,b in terms of ar,ai,br,bi
a=ar+ai*1i;
b=br+bi*1i;

% Set up grid
x = linspace(alpha,beta,n)';
h = x(2)-x(1);

% Define differentiation matrices
v1 = ones(n,1);
D2 = (diag(-2*v1) + diag(v1(2:end),1)+diag(v1(2:end),-1))/h^2;
D1 = (-diag(v1(1:end-1),-1) + diag(v1(2:end),1))/(2*h);

% Set up matrices corresponding to coefficients of derivatives in (3.2)
B=zeros(n,n);
for j=1:n
    B(:,j)=x;
end

% Define matrix corresponding to differential equation
Mfd = B.*D2+b*D1-B.*D1-a*eye(n);
Mfd(1,:) = [1 zeros(1,n-1)];
Mfd(n,:) = [zeros(1,n-1) 1];

% Define RHS matrix corresponding to boundary conditions
d = zeros(n,1);
d(1)=hypergeom(a,b,alpha);
d(end)=hypergeom(a,b,beta);

% Solve system to obtain solution
fd = Mfd\d;
plot(x,fd)