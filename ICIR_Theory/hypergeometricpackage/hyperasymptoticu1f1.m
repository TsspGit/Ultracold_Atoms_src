function [h]=hyperasymptoticu1f1(ar,ai,br,bi,zr,zi,tol,prob)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Model code for computing the hypergeometric function 1F1(a;b;z) using   %
% the hyperasymptotic expansion (G.5). This method cannot work effectively%
% in MATLAB due to its inability to compute the incomplete gamma function %
% for all parameter and variable values within the complex plane.         %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance as detailed in Section 3.5              %
%         prob=1 if asymptotic expansion for U(a;b;z) is wanted, and 2 if %
%              hyperasymptotic expansion is wanted                        %
% Output: h=Computed value of U(a;b;z)                                    %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Initialise a1, vector of current values, and b1, sum of all values
% computed so far
a1=zeros(1,1);
a1(1)=1;b1=1;

% Sum up to ncount terms
for j=1:500
    % Update count of terms
    ncount=j+1;
    % Update current term and sum of terms
    a1(j+1)=(a+j-1)*(a-b+j)/j/(-z)*a1(j);
    b1=b1+a1(j+1);
    % Terminate computation if stopping criterion is satisfied
    if abs(a1(j))/abs(b1)<tol && abs(a1(j+1))/abs(b1)<tol
        break
    end
    % If 500 terms have been computed before stopping criterion has been
    % satisfied, state this
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Asymptotic expansion for U(a;b;z)
if prob==1
    h=z^(-a)*b1;
end

% Initialise c1, coefficient of Gnu(z) in (G.6), d1, equal to Gnu(z), e1,
% equal to individual term, and f1, equal to the sum of all terms computed
% so far
c1=zeros(1,1);d1=zeros(1,1);e1=zeros(1,1);
c1(1)=1;
d1(1)=gamma(ncount+2*a-b)*gammainc(z,1-ncount-2*a+b)*gamma(1-ncount-2*a+b);
e1(1)=c1*d1;
f1=e1;

% Sum up to ncount/2 terms
for j=1:fix(ncount/2)-1
    % Update c1,d1,e1,f1
    c1(j+1)=(-a+j)*(b-a+j-1)/j/(-z)*c1(j);
    d1(j+1)=exp(z)/2/pi*gamma(ncount+2*a-b-j)*gammainc(z,1-ncount-2*a+b+j);
    e1=c1(j+1)*d1(j+1);
    f1=f1+e1;
end

% Hyperasymptotic expansion for U(a;b;z)
if prob==2
    h=z^(-a)*b1+(-1)^(ncount)*2*pi*z^(a-b)/gamma(a)/gamma(a-b+1)*f1;
end