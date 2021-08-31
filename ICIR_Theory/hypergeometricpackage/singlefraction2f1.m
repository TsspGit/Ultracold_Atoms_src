function [h]=singlefraction2f1(ar,ai,br,bi,cr,ci,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 2F1(a,b;c;z) using the   %
% single fraction method detailed in Section 4.3.                         %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance as detailed in Section 4.3              %
% Output: h=Computed value of 2F1(a,b;c;z)                                %
%-------------------------------------------------------------------------%

% Compute a,b,c,z in terms of ar,ai,br,bi,cr,ci,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
c=cr+ci*1i;
z=zr+zi*1i;

% Initialise vectors of terms denoted as alpha(j), beta(j), gamma(j) and
% xi(j) in Section 4.3
a1=[0,c];b1=[1,a*b*z];c1=[1,c];d1=[1,(c+a*b*z)/c];

for j=3:500
    % Update alpha(j), beta(j), gamma(j) in terms of previous entries
    a1(j)=(a1(j-1)+b1(j-1))*(j-1)*(c+j-2);
    b1(j)=b1(j-1)*(a+j-2)*(b+j-2)*z;
    c1(j)=c1(j-1)*(j-1)*(c+j-2);
    % Stop if any terms become infinitely large
    if (a1(j)==Inf)||(b1(j)==Inf)||(c1(j)==Inf)
        break
    end
    % Apply (4.7) in text
    d1(j)=(a1(j)+b1(j))/c1(j);
    % Terminate computation if stopping criterion is satisfied
    if abs(d1(j)-d1(j-1))/abs(d1(j-1))<tol && abs(d1(j-1)-d1(j-2))/abs(d1(j-2))<tol
        break
    end
    % If 500 terms have been computed before stopping criterion has been
    % satisfied, state this
    if (j==1000)
        ['Reached ' num2str(j) ' limit']
        return
    end
end

% Return latest single fraction as solution
h=d1(end);