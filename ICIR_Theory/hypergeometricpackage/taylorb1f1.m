function [h]=taylorb1f1(ar,ai,br,bi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using Taylor  %
% series method (b) of Section 3.2.                                       %

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

% Compute a,b,c using ar,ai,br,bi,z,r,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Initialise vector r, which stores rj's as details in Section 3.2
r=zeros(2,1);
r(1)=a/b;
r(2)=(a+1)/2/(b+1);

% Initialise vector A, which stores sum of first j+1 terms, denoted as Sj
% in Section 3.2
A=zeros(2,1);
A(1)=1+z*r(1);
A(2)=A(1)+z^2*a/b*r(2);

for j=3:500
    % Compute current rj
    r(j)=(a+j-1)/j/(b+j-1);
    % Compute Aj in terms of Aj-1 and Aj-2, as detailed in Section 3.2
    A(j)=A(j-1)+(A(j-1)-A(j-2))*r(j)*z;
    % If stopping criterion is satisfied, terminate summation
    if abs(A(j)-A(j-1))/abs(A(j-1))<tol && abs(A(j-1)-A(j-2))/abs(A(j-2))<tol
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
% i.e return Aend
h=A(end);