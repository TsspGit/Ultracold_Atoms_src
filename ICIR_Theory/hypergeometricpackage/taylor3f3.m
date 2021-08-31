function [h]=taylor3f3(ar,ai,br,bi,cr,ci,dr,di,er,ei,fr,fi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the hypergeometric function 3F3(a,b,c;d,e,f;z) using a Taylor  %
% method, equivalent to method (a) from Sections 3.2 and 4.2.             %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         dr=Re(d)                                                        %
%         di=Im(d)                                                        %
%         er=Re(e)                                                        %
%         ei=Im(e)                                                        %
%         fr=Re(f)                                                        %
%         fi=Im(f)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance, discussed in Sections 3.2 and 4.2      %
% Output: h=Computed value of 3F3(a,b,c;d,e,f;z)                          %
%-------------------------------------------------------------------------%

% Compute a,b,c,d,e,f,z from ar,ai,br,bi,cr,ci,dr,di,er,ei,fr,fi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
c=cr+ci*1i;
d=dr+di*1i;
e=er+ei*1i;
f=fr+fi*1i;
z=zr+zi*1i;

% Initialise individual term computed, and sum of all terms computed so far
a1=1;b1=1;

for j=1:500
    % Update current term in terms of previous one
    a1=(a+j-1)*(b+j-1)*(c+j-1)/(d+j-1)/(e+j-1)/(f+j-1)*z/j*a1;
    % Update sum of terms computed so far
    b1=b1+a1;
    % Terminate summation if stopping criterion is satisfied
    if abs(a1)/abs(b1)<tol
        break
    end
    % If 500 terms have been computed before stopping criterion is
    % satisfied, state this
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Return sum as solution
h=b1;