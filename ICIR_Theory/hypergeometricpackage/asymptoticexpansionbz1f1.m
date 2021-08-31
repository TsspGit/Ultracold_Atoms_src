function [h]=asymptoticexpansionbz1f1(ar,ai,br,bi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using         %
% asymptotic expansion for large |b| and |z| of Section G.3, for          %
% Re(b-z-a-1)>0.                                                          %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance as detailed in Section 3.2              %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Initialise c1, coefficient of a2n(z) and a2n-1(z), d1, coefficient of
% a2n(z), e1, product of c1 and d1, and f1, sum of all terms computed thus
% far
c1=gamma(a-1)/(b-a-z-1)^a;
d1=a-1;
e1=c1*d1;
f1=e1;

for n=1:500
    % Update c1, d2 (current term) and d3 (sum of all terms for a2n(z))
    c1=c1*(2*n+a-3)*(2*n+a-2)/(b-a-z-1)^2;
    d2=(b-a-1)^(2*n)/factorial(2*n);
    d3=d2;
    for k=1:2*n
        d2=d2*(a-b+k)/(b-a-1)/k*(2*n-k+1);
        d3=d3+d2;
    end
    % Update d4 (current term) and d5 (sum of all terms for a2n-1(z))
    d4=(b-a-1)^(2*n-1)/factorial(2*n-1);
    d5=d4;
    for k=1:2*n-1
        d4=d4*(a-b+k)/(b-a-1)/k*(2*n-k);
        d5=d5+d4;
    end
    % Update d1, combination of a2n(z) and a2n-1(z), e1, multiplication of
    % that with coefficient, and f1, sum of all terms computed so far
    d1=(2*n+a-1)*d3+(b-a-z-1)*d5;
    e1=c1*d1;
    f1=f1+e1;
    % Terminate summation if stopping criterion is satisfied
    if abs(e1)/abs(f1)<tol
        break
    end
    % If 500 terms have been computed before stopping criterion is
    % satisfied, state this
    if (n==500)
        ['Reached ' num2str(n) ' limit']
        return
    end
end

% Return solution
h=gamma(b)/gamma(a)/gamma(b-a)*f1;