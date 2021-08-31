function [h]=bacvas1f1(ar,ai,br,bi,zr,zi,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) for parameters%
% a,b and variable z by applying the method of computing oscillatory      %
% integrals as detailed in Appendix G.5.                                  %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         n=Number of mesh points to be used                              %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Find Gauss quadrature nodes and weights using qrule.m [72]
omega=zi/2;
[x,w]=qrule(n,1);

P=zeros(n,2);
a1=zeros(n,1);

% Set up integral in (G.19}
r=zeros(n,1);
for i=1:n
    r(i)=exp(zr*x(i)/2)*(1+x(i))^(a-1)*(1-x(i))^(b-a-1);
end

% Compute Legendre polynomials
for l=1:n
    for i=1:n
        for k=1:n
            lv=legendre(k-1,x(i));
            P(k,i)=lv(1);
        end
    end
    a1(l)=sum(w'.*P(l,:)'.*r);
end

b1=zeros(n,1);
c1=b1;
d1=c1;

% Compute b1, i^(k-1) term in (G.18), c1, 2k-1 term, and d1, Bessel
% function term
for l=1:n
    b1(l)=1i^(l-1);
    c1(l)=2*l-1;
    d1(l)=besselj(l-1/2,omega);
end

% Compute sum in (G.18)
e1=sqrt(pi/2/omega)*sum(b1.*c1.*d1.*a1);

% Multiply by appropriate Gamma functions in (3.25) and return solution
h=gamma(b)/gamma(a)/gamma(b-a)*exp(z/2)/2^(b-1)*e1;