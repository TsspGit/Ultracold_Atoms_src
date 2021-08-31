function [h]=incgammaexpansion1f1(a,b,z,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;-z) using the    %
% incomplete gamma expansion of Section G.2.                              %

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         z=Variable value z                                              %
%         tol=Specified tolerance as detailed in Section 3.2              %
% Output: h=Computed value of 1F1(a;b;-z)                                 %
%-------------------------------------------------------------------------%

% Initialise a1, coefficient of incomplete gamma function in series in 
% (G.3), b1, incomplete gamma function term in (G.3), and c1, their product
a1=zeros(1,1);b1=zeros(1,1);
a1(1)=1;b1(1)=gammainc(z,a)*gamma(a);c1=b1(1);

for j=1:500
    % Update a1,b1,c1
    a1(j+1)=(a-b+j)/j/z*a1(j);
    b1(j+1)=a1(j+1)*gammainc(z,j+a)*gamma(j+a);
    c1=c1+b1(j+1);
    % Terminate summation if stopping criterion is satisfied
    if abs(b1(j))/abs(c1)<tol && abs(b1(j+1))/abs(c1)<tol
        break
    end
    % State that 500 terms were computed before stopping criterion is
    % satisfied if this is the case
    if (j==500)
        ['Reached ' num2str(j) ' limit']
        return
    end
end

% Return solution as in (G.3)
h=real(gamma(b)/gamma(b-a)/gamma(a)*z^(-a)*c1);