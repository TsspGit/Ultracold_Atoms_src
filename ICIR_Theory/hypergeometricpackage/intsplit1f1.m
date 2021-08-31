function [h]=intsplit1f1(a,b,zr,zi,c,nout)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) for parameters%
% a,b and variable z by applying the method of splitting the integral as  %
% detailed in Appendix G.5.                                               %

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         c=Length of outer integrals we wish to take                     %
%         nout=Number of mesh points we use on each outer integral        %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute z in terms of zr,zi
z=zr+zi*1i;

% Compute first integral in (G.10) using Gauss-Jacobi quadrature with the
% program qrule.m [72]
[xend1,wend1]=qrule(nout,7,0,a-1);
end1=(c/2)^(b-1)*exp(0.5*z*c)*sum(wend1.*(exp(0.5*z*c*xend1)...
    .*((2-c)/c-xend1).^(b-a-1)));

% Compute third integral in (G.10) using Gauss-Jacobi quadrature with the
% program qrule.m [72]
[xend2,wend2]=qrule(nout,7,b-a-1,0);
end2=(c/2)^(b-1)*exp(z*(1-0.5*c))*sum(wend2.*(exp(0.5*z*c*xend2)...
    .*((2-c)/c+xend2).^(a-1)));

% Compute the middle integral of (G.10) using 'quad'
f=@(y) exp(z*y).*y.^(a-1).*(1-y).^(b-a-1);
mid=quad(f,c,1-c,1e-15);

% Compute solution using (G.10)
h=gamma(b)/gamma(a)/gamma(b-a)*(end1+mid+end2);