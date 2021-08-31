function [h]=asymptoticb1f1(ar,ai,br,bi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using         %
% asymptotic series method (b) of Section 3.5.                            %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance as detailed in Section 3.5              %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Initialise r(j), which represents current term, and A(j), which
% represents sum of all terms computed thus far, for first series in (3.23)
% and (3.24)
r=zeros(2,1);
A=zeros(2,1);
r(1)=(b-a)*(1-a);
r(2)=(b-a+1)*(2-a)/2;
A(1)=1+r(1)/z;
A(2)=A(1)+r(1)*r(2)/z^2;

for j=3:500
    % Update r(j) and A(j) in terms of r(j-1), A(j-1) and A(j-2)
    r(j)=(b-a+j-1)*(j-a)/j;
    A(j)=A(j-1)+(A(j-1)-A(j-2))*r(j)/z;
    % Terminate summation if stopping criterion is satisfied
    if abs(A(j)-A(j-1))/abs(A(j-1))<tol && abs(A(j-1)-A(j-2))/abs(A(j-2))<tol
        break
    end
    % If 500 terms computed before stopping criterion is satisfied, state
    % this
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Initialise s(j), which represents current term, and B(j), which
% represents sum of all terms computed thus far, for second series in
% (3.23) and (3.24)
s=zeros(2,1);
B=zeros(2,1);
s(1)=a*(a-b+1);
s(2)=(a+1)*(a-b+2)/2;
B(1)=1+s(1)/(-z);
B(2)=B(1)+s(1)*s(2)/z^2;

for j=3:500
    % Update s(j) and B(j) in terms of s(j), B(j-1) and B(j-2)
    s(j)=(a+j-1)*(a-b+j)/j;
    B(j)=B(j-1)+(B(j-1)-B(j-2))*s(j)/(-z);
    % Stopping criterion
    if abs(B(j)-B(j-1))/abs(B(j-1))<tol && abs(B(j-1)-B(j-2))/abs(B(j-2))<tol
        break
    end
    % If 500 terms computed, state this
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Compute relevant Gamma functions using cgama.m function [71]
[gr1,gi1]=cgama(ar,ai,1);[gr2,gi2]=cgama(br,bi,1);[gr3,gi3]=cgama(br-ar,bi-ai,1);

% Compute (3.23) or (3.24)
if zr<0 && zi<0
    d=A(end)*exp(z)*z^(a-b)/(gr1+gi1*1i);
    e=B(end)*exp(-pi*1i*a)*z^(-a)/(gr3+gi3*1i);
    c=(gr2+gi2*1i)*(d+e);
else
    d=A(end)*exp(z)*z^(a-b)/(gr1+gi1*1i);
    e=B(end)*exp(pi*1i*a)*z^(-a)/(gr3+gi3*1i);
    c=(gr2+gi2*1i)*(d+e);
end

% Take real part if Im(a)=Im(b)=Im(z)=0, otherwise simply take what has
% been computed
if imag(a)==0 && imag(b)==0 && imag(z)==0
    h=real(c);
else
    h=c;
end