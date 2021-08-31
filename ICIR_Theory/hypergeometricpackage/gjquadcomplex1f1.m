function [h]=gjquadcomplex1f1(ar,ai,br,bi,zr,zi,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) for complex   %
% a,b,z by applying Gauss-Jacobi quadrature on (3.25), as in Section 3.6. %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         n=Number of mesh points, denoted as Nmesh in Section 3.6        %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Apply code qrule.m [72] to integral (3.25) to find nodes and weights
[x,w]=qrule(n,7,b-a-1,a-1);

% Compute relevant Gamma functions using cgama.m [71]
[gr1,gi1]=cgama(ar,ai,1);
[gr2,gi2]=cgama(br,bi,1);
[gr3,gi3]=cgama(br-ar,bi-ai,1);
e1=gr1+gi1*1i;e2=gr2+gi2*1i;e3=gr3+gi3*1i;

% Apply x,w,e1,e2,e3 to (3.25)
h=e2/e1/e3/(2^(b-1))*exp(0.5*z)*sum(w.*(exp(0.5*(zr+zi*1i)*x)));