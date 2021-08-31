function [h]=intsplit2f1(a,b,c,zr,zi,d,nout)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the Gauss hypergeometric function 2F1(a,b;c;z) for parameters  %
% a,b,c and variable z by applying the method of splitting the integral as%
% detailed in Appendix H.1.                                               %

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         c=Parameter value c                                             %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         d=Length of outer integrals we wish to take                     %
%         nout=Number of mesh points we use on each outer integral        %
% Output: h=Computed value of 2F1(a,b;c;z)                                %
%-------------------------------------------------------------------------%

% Compute z in terms of zr,zi
z=zr+zi*1i;

% Compute first integral in (H.1) using Gauss-Jacobi quadrature with the
% program qrule.m [72]
[xend1,wend1]=qrule(nout,7,0,b-1);
end1=(d/2)^(c-1)*sum(wend1.*((1-d/2*z*(xend1+1)).^(-a)...
    .*((2-d)/d-xend1).^(c-b-1)));

% Compute third integral in (H.1) using Gauss-Jacobi quadrature with the
% program qrule.m [72]
[xend2,wend2]=qrule(nout,7,c-b-1,0);
end2=(d/2)^(c-1)*sum(wend2.*((1-z*(d/2*(xend2-1)+1)).^(-a)...
    .*((2-d)/d+xend2).^(b-1)));

% Compute the middle integral of (H.1) using 'quad'
f=@(y) (1-z*y).^(-a).*y.^(b-1).*(1-y).^(c-b-1);
mid=quad(f,d,1-d,1e-15);

% Compute solution using (H.1)
h=gamma(c)/gamma(b)/gamma(c-b)*(end1+mid+end2);