function [h]=buchholzc1f1(ar,ai,br,bi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using method 3%
% of Section 3.4.                                                         %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance, defined as the relative size of one    %
%             term in series (3.22) in terms of sum of all previous terms %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Compute t as shown in (3.22)
t=z*(a-b/2);
tr=real(t);
ti=imag(t);

% Initialise a1, value of 1/(2^j*(b)j), b1, value of
% pj(b,z)*0F1( ;b+j;xi)/a1, and c1, sum of all b1's computed so far
a1=1;
b1=buchholzpoly(tr,ti,zr,zi,0)*hgf0f1(br,bi,tr,ti,tol);
c1=b1;

for j=1:170
    % Update a1,b1,c1 in terms of last
    a1=a1*2*(b+j-1);
    b1=buchholzpoly(br,bi,zr,zi,j)*hgf0f1(br+j,bi,tr,ti,tol2)/a1;
    c1=c1+b1;
    % If stopping criterion is satisfied, terminate summation
    if abs(b1)/abs(c1)<tol
        break
    end
    % If 170 terms have been computed without stopping criterion being
    % satisfied, say so
    if (j==170)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Return h, which is exp(z/2) multiplied by the series in (3.22)
h=exp(z/2)*c1;