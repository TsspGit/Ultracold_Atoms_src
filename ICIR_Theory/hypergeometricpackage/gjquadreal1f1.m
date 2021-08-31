function [h]=gjquadreal1f1(a,b,z,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) for real a,b,z%
% by applying Gauss-Jacobi quadrature on (3.25), as in Section 3.6.       %

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         z=Variable value z                                              %
%         n=Number of mesh points, denoted as Nmesh in Section 3.6        %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Apply code qrule.m [72] to integral (3.25) to find nodes and weights
[x,w]=qrule(n,7,b-a-1,a-1);

% Compute relevant Gamma functions
e1=gamma(a);e2=gamma(b);e3=gamma(b-a);

% Apply x,w,e1,e2,e3 to (3.25)
h=e2/e1/e3/(2^(b-1))*exp(0.5*z)*sum(w.*(exp(0.5*(zr+zi*1i)*x)));