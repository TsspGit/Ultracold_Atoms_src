function [h]=taylora0f1(ar,ai,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric limit function 0F1(;a;z) using    %
% method (a) as described in Appendix I.                                  %

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

% Initialise a1, current term, and b1, sum of all terms computed so far
a1=1;b1=1;

for j=1:500
    % Compute next term from previous one
    a1=1/(a+j-1)*z/j*a1;
    % Update sum
    b1=b1+a1;
    % If stopping criterion is fulfilled, terminate summation
    if abs(a1)/abs(b1)<tol
        break
    end
    % If 500 terms are computed before stopping criterion is satisfied, say
    % so
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Return sum as solution
h=b1;