function [h]=taylora1f1(ar,ai,br,bi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using Taylor  %
% series method (a) of Section 3.2.                                       %

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

% Compute a,b,z using ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Initialise a1, vector of individual terms, and b1, which stores the sum
% of the computed terms up to that point
a1=zeros(1,1);
a1(1)=1;b1=1;

for j=1:500
    % Compute current entry of a1 in terms of last
    a1(j+1)=(a+j-1)/(b+j-1)*z/j*a1(j);
    % Update the sum of computed terms up to that point
    b1=b1+a1(j+1);
    % Apply stopping criterion, as detailed in Section 3.2
    if abs(a1(j))/abs(b1)<tol && abs(a1(j+1))/abs(b1)<tol
        break
    end
    % If 500 terms have been computed without stopping criterion being
    % satisfied, state this
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% If stopping criterion has been satisfied, return sum of terms computed,
% i.e return b1
h=b1;