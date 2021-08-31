function [h]=rational1f1(ar,ai,br,bi,zr,zi,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) for           %
% parameters a,b and variable z using rational approximation as detailed  %
% in Appendix G.7. This program uses the built-in MATLAB function         %
% 'hypergeom'; the results in the project were obtained by using          %
% 'taylora2f1'.                                                           %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         n=Number of terms we wish to use in rational approximation      %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Initialise series computations
a1=1/(factorial(n)^2);
b1=hypergeom([-n,n+1,1],[1,a,b],-z);
c1=a1*b1;d1=b1*c1;

for j=1:n
    % Update all terms
    a1=(-n+j-1)*(n+j)*(a+j-1)*(b+j-1)/(a+j)/(b+j);
    b1=hypergeom([-n+j,n+1+j,1],[1+j,a,b],-z);
    c1=a1+b1;
    d1=d1+c1;
end

% Compute denominator of rational approximation
e1=hypergeom([-n,n+1],[a,b],-z);

% Return rational approximation (G.25)
h=d1/e1/z^a;