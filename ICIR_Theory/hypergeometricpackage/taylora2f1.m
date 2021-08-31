function [h]=taylora2f1(ar,ai,br,bi,cr,ci,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 2F1(a,b;c;z) using Taylor%
% series method (a) of Section 4.2.                                       %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance as detailed in Section 4.2              %
% Output: h=Computed value of 2F1(a,b;c;z)                                %
%-------------------------------------------------------------------------%

% Compute a,b,c,z in terms of ar,ai,br,bi,cr,ci,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
c=cr+ci*1i;
z=zr+zi*1i;

% Initialise vector of individual terms and sum of terms computed thus far
a1=zeros(1,1);
a1(1)=1;b1=1;

for j=1:500
    % Update value of a1, current term, and b1, sum of all terms in terms
    % of previous values
    a1(j+1)=(a+j-1)*(b+j-1)/(c+j-1)*z/j*a1(j);
    b1=b1+a1(j+1);
    % Terminate summation if stopping criterion is satisfied
    if abs(a1(j))/abs(b1)<tol && abs(a1(j+1))/abs(b1)<tol
        break
    end
    % If 500 terms have been computed, say so
    if (j==1000)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Return sum as solution
h=b1;