function [h]=romberg1f1(a,b,z,n,m)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using Romberg %
% integration as described in Appendix G.5.                               %

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         z=Variable value z                                              %
%         n=Value n as denoted in (G.15)                                  %
%         m=Value m as denoted in (G.15)                                  %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Define integrand
f=@(y) exp(z*y)*y^(a-1)*(1-y)^(b-a-1);

% Define recursively S(n,m) as in (G.14)-(G.16)
S=zeros(n+1,m+1);
S(1,1)=(f(0)+f(1))/2;
for i=2:n+1
    v=zeros(2^(i-1),1);
    for k=1:2^(i-2)
        v(k)=f((2*k-1)/2^(i-1));
    end
    S(i,1)=S(i-1,1)/2+sum(v)/2^(i-1);
end
for i=2:n+1
    for j=2:m
        S(i,j)=(4^m*S(i,j-1)-S(i-1,j-1))/(4^m-1);
    end
end

% Return solution to (3.25)
h=S(n,m)*gamma(b)/gamma(a)/gamma(b-a);