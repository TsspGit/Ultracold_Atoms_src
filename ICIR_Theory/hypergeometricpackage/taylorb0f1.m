function [h]=taylorb0f1(ar,ai,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric limit function 0F1(;a;z) using    %
% method (b) as described in Appendix I.                                  %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance, defined as in Appendix I               %
% Output: h=Computed value of 0F1(;a;z)                                   %
%-------------------------------------------------------------------------%

% Compute a,z in terms of ar,ai,zr,zi
a=ar+ai*1i;
z=zr+zi*1i;

% Initialise r(j) as described in Appendix I
r=zeros(2,1);
r(1)=1/a;
r(2)=1/2/(a+1);

% Initialise vector of individual terms
A=zeros(2,1);
A(1)=1+z*r(1);
A(2)=A(1)+z^2/a*r(2);

for j=3:500
    % Update r(j) and A(j) in terms of previous entries
    r(j)=1/j/(a+j-1);
    A(j)=A(j-1)+(A(j-1)-A(j-2))*r(j)*z;
    % Terminate summation if stopping criterion is satisfied
    if abs(A(j)-A(j-1))/abs(A(j-1))<tol
        break
    end
    % If 500 terms computed before stopping criterion is satisfied, state
    % this
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Return sum as solution
h=A(end);