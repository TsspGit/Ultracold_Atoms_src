function [h]=buchholzbreal1f1(a,b,z,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using method 2%
% of Section 3.4 for real a,b,z.                                          %

%-------------------------------------------------------------------------%
% Input:  a  =Parameter value a                                           %
%         b  =Parameter value b                                           %
%         z  =Variable value z                                            %
%         tol=Specified tolerance, defined as the relative size of one    %
%             term in series (3.18) in terms of sum of all previous terms %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Initialise vector A, which corresponds to entries Dj on (3.20)
A=zeros(3,1);
A(1)=1;A(2)=0;A(3)=b/2;

% Initialise a1, individual term of (3.20), and b1, sum of all terms
% computed thus far
a1=besselj(b-1,sqrt(z*(2*b-4*a)))/(sqrt(z*(2*b-4*a))^(b-1))...
    +b/2*z^2*besselj(b+1,sqrt(z*(2*b-4*a)))/(sqrt(z*(2*b-4*a))^(b+1));
b1=a1;

for j=3:500
    % Compute A(j+1) in terms of A(j), and a1 and b1 in terms of previous
    % terms
    A(j+1)=((j-2+b)*A(j-1)+(2*a-b)*A(j-2))/j;
    a1=A(j+1)*z^(j)*besselj(b-1+j,sqrt(z*(2*b-4*a)))/(sqrt(z*(2*b-4*a))^(b-1+j));
    b1=b1+a1;
    % If stopping criterion is satisfied, terminate computation
    if abs(a1)/abs(b1)<tol
        break
    end
    % If 500 terms are computed without stopping criterion being satusfied,
    % say so
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Return RHS of (3.20) using sum of series
h=gamma(b)*exp(z/2)*2^(b-1)*b1(end);