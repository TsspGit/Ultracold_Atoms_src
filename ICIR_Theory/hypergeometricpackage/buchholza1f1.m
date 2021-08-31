function [h]=buchholza1f1(ar,ai,br,bi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using method 1%
% of Section 3.4.                                                         %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance, defined as the relative size of one    %
%             term in series (3.18) in terms of sum of all previous terms %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Compute Gamma(b) using routine cgama.m from MathWorks website [71]
[gr1,gi1]=cgama(br,bi,1);

% Initialise computation of individual term, a1, and sum of all terms
% computed so far, b1
a1=buchholzpoly(br,bi,zr,zi,0)*besselj(b-1,sqrt(z*(2*b-4*a)))...
    /(sqrt(z*(2*b-4*a))^(b-1));
b1=a1;

for j=1:170 % Sum up to 170, for at this point Bernoulli numbers for
            % Buchholz polynomial computation are infinite 
    [gr11,gi11]=cgama(real(b+j),imag(b+j),1);
    % Compute latest term
    a1=buchholzpoly(br,bi,zr,zi,j)*(sqrt(z*(2*b-4*a))/2)^(b-1+j)/(gr11+gi11*1i)...
        *hgf0f1(br+j,bi,-real(z*(2*b-4*a))/4,-imag(z*(2*b-4*a))/4,1e-15)...
        /(z*(2*b-4*a))^(0.5*(b-1+j));
    % Sum this with sum of all previous terms
    b1=b1+a1;
    % Stopping criterion
    if abs(a1)/abs(b1)<tol
        break
    end
    % If 170 terms have been computed without stopping criterion being
    % satisfied, say so
    if (j==170)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Return correct value of 1F1 using the series and its coefficient in (3.18)
h=(gr1+gi1*1i)*exp(z/2)*2^(b-1)*b1;