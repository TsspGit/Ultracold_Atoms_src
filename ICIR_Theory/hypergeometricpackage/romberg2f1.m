function [r]=romberg2f1(a,b,c,z,n,m)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 2F1(a,b;c;z) using       %
% integration as described in Appendix H.1.                               %

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         c=Parameter value c                                             %
%         z=Variable value z                                              %
%         n=Value n as denoted in (G.15)                                  %
%         m=Value m as denoted in (G.15)                                  %
% Output: h=Computed value of 2F1(a,b;c;z)                                  %
%-------------------------------------------------------------------------%

% Define integrand
f=@(y) (1-z*y)^(-a)*y^(b-1)*(1-y)^(c-b-1);

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

% Return solution to (4.8)
r=S(n,m)*gamma(c)/gamma(b)/gamma(c-b);