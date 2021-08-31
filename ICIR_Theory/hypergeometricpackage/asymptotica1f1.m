function [h]=asymptotica1f1(ar,ai,br,bi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using         %
% asymptotic series method (a) of Section 3.5.                            %

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

% Initialise a1, which represents current term, and b1, which represents
% sum of all terms computed thus far, for first series in (3.23) and (3.24)
a1=zeros(1,1);
a1(1)=1;b1=1;

for j=1:500
    % Update a1(j) and b1 in terms of previously computed a1(j-1) and b1
    a1(j+1)=(b-a+j-1)*(-a+j)/j/z*a1(j);
    b1=b1+a1(j+1);
    % If stopping criterion is satisfied
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

% Initialise c1, which represents current term, and d1, which represents
% sum of all terms computed thus far, for second series in (3.23) and
% (3.24)
c1=zeros(1,1);
c1(1)=1;d1=1;

for k=1:500
    % Update c1 and d1
    c1(k+1)=(a+k-1)*(a-b+k)/k/(-z)*c1(k);
    d1=d1+c1(k+1);
    % Stopping criterion
    if abs(c1(k))/abs(d1)<tol && abs(c1(k+1))/abs(d1)<tol
        break
    end
    % Specify if 500 terms computed
    if (k==500)
        [' ' num2str(k) ' terms computed']
        return
    end
end

% Take last terms computed
h1=b1;h2=d1;

% Compute relevant Gamma functions using cgama.m [71]
[gr1,gi1]=cgama(br,bi,1);[gr2,gi2]=cgama(ar,ai,1);[gr3,gi3]=cgama(br-ar,bi-ai,1);

% Compute (3.23) or (3.24)
if (zr>0)
    h3=(gr1+gi1*1i)*exp(z)*z^(a-b)/(gr2+gi2*1i)*h1...
        +exp(pi*1i*a)*z^(-a)/(gr3+gi3*1i)*h2;
elseif (zr<0)
    h3=(gr1+gi1*1i)*(exp(z)*z^(a-b)/(gr2+gi2*1i)*h1...
        +exp(-pi*1i*a)*z^(-a)/(gr3+gi3*1i)*h2);
end

% Take real part if Im(a)=Im(b)=Im(z)=0, otherwise simply take what has
% been computed
if imag(a)==0 && imag(b)==0 && imag(z)==0
    h=real(h3);
else
    h=h3;
end