function [h]=taylorb2f1(ar,ai,br,bi,cr,ci,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 2F1(a,b;c;z) using Taylor%
% series method (b) of Section 4.2.                                       %

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

% Initialise r(j) as detailed in Section 4.2
r=zeros(2,1);
r(1)=a*b/c;
r(2)=(a+1)*(b+1)/2/(c+1);

% Initialise A(j) as detailed in Section 4.2
A=zeros(2,1);
A(1)=1+z*r(1);
A(2)=A(1)+z^2*a*b/c*r(2);

for j=3:500
    % Update r(j) and A(j) in terms of previous values
    r(j)=(a+j-1)*(b+j-1)/j/(c+j-1);
    A(j)=A(j-1)+(A(j-1)-A(j-2))*r(j)*z;
    % If stopping criterion is satisfied, terminate computation
    if abs(A(j)-A(j-1))/abs(A(j-1))<tol && abs(A(j-1)-A(j-2))/abs(A(j-2))<tol
        break
    end
    % If 500 terms have been computed before stopping criterion has been
    % satisfied, state this
    if (j==1000)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Return sum as solution
h=A(end);