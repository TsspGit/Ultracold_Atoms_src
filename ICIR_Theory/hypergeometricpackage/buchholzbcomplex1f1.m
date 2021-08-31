function [h]=buchholzbcomplex1f1(ar,ai,br,bi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using method 2%
% of Section 3.4 for complex a,b,z.                                       %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance, defined as the relative size of one    %
%             term in series (3.20) in terms of sum of all previous terms %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Initialise vector of coefficientd Dj
A=zeros(3,1);
A(1)=1;A(2)=0;A(3)=b/2;

% Compute relevant gamma functions using cgama.m from MathWorks website
% [71]
[gr100,gi100]=cgama(real(b),imag(b),1);[gr101,gi101]=cgama(real(b)+2,imag(b),1);

% Initialise a1, first 2 terms, and b1, sum of all terms computed so far
a1=(sqrt(z*(2*b-4*a))/2)^(b-1)/(gr100+gi100*1i)...
    *hgf0f1(real(b-1)+1,imag(b-1),-real(z*(2*b-4*a))/4,-imag(z*(2*b-4*a))/4,1e-15)...
    /(sqrt(z*(2*b-4*a))^(b-1))+b/2*z^2*(sqrt(z*(2*b-4*a))/2)^(b+1)/(gr101+gi101*1i)...
    *hgf0f1(real(b+1)+1,imag(b+1),-real(z*(2*b-4*a))/4,-imag(z*(2*b-4*a))/4,1e-15)...
    /(sqrt(z*(2*b-4*a))^(b+1));
b1=a1;

for j=3:500
    % Update coefficient Dj
    A(j+1)=((j-2+b)*A(j-1)+(2*a-b)*A(j-2))/j;
    % Compute relevant gamma function
    [gr11,gi11]=cgama(real(b+j),imag(b+j),1);
    % Update latest term to be computed
    a1=A(j+1)*z^(j)*(sqrt(z*(2*b-4*a))/2)^(b-1+j)/(gr11+gi11*1i)...
        *hgf0f1(br+j,bi,-real(z*(2*b-4*a))/4,-imag(z*(2*b-4*a))/4,1e-15)...
        /(z*(2*b-4*a))^(0.5*(b-1+j));
    % Update sum of all terms
    b1=b1+a1;
    % If stopping criterion is satisfied, terminate computations
    if abs(a1)/abs(b1)<tol
        break
    end
    % If 500 terms have been computed before stopping criterion has been
    % satisfied, say so
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Compute 1F1 in terms of series
h=(gr100+gi100*1i)*exp(z/2)*2^(b-1)*b1(end);