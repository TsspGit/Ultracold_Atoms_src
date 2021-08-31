function [h]=buhringa2f1(ar,ai,br,bi,cr,ci,zr,zi,z0,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the Gauss hypergeometric function 2F1(a,b;c;z) using the       %
% analytic continuation formula (4.21) of Section 4.7. This is the code   %
% used for real a,b,c; code for complex a,b,c would be obtained by        %
% replacing gamma.m with cgama.m [71].                                    %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         z0=Value at which expansion is taken (usually 1/2)              %
%         tol=Specified tolerance as detailed in Section 4.2              %
% Output: h=Computed value of 2F1(a,b;c;z)                                %
%-------------------------------------------------------------------------%

% Compute a,b,c,z in terms of ar,ai,br,bi,cr,ci,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
c=cr+ci*1i;
z=zr+zi*1i;

% Initialise vector for coefficients dj(a,z0) in (4.21)
d=zeros(2,1);
d(1)=(1+a-b)^(-1)*a*((a+1)*(1-2*z0)+(a+b+1)*z0-c);
d(2)=(2*(2+a-b))^(-1)*(a+1)*(((a+2)*(1-2*z0)+(a+b+1)*z0-c)...
    *(1+a-b)^(-1)*a*((a+1)*(1-2*z0)+(a+b+1)*z0-c)+z0*(1-z0)*a);

% Initialise value of individual term and sum of all terms in first series
% of (4.21) computed so far
a1=1+d(1)*(z-z0)^(-1)+d(2)*(z-z0)^(-2);
b1=a1;
for n=3:500
    % Update d(n), a1 and b1
    d(n)=(n*(n+a-b))^(-1)*(n+a-1)*...
        (((n+a)*(1-2*z0)+(a+b+1)*z0-c)*d(n-1)+z0*(1-z0)*(n+a-2)*d(n-2));
    a1=d(n)*(z-z0)^(-n);
    b1=b1+a1;
    % Terminate summation if stopping criterion is satisfied
    if abs(a1)/abs(b1)<tol
        break
    end
    % If 500 terms have been computed before stopping criterion is
    % satisfied, state this
    if (n==500)
        [' ' num2str(n) ' terms computed']
        return
    end
end

% Initialise vector for dj(b,z0) as in (4.21)
e=zeros(2,1);
e(1)=(1-a+b)^(-1)*b*((b+1)*(1-2*z0)+(a+b+1)*z0-c);
e(2)=(2*(2-a+b))^(-1)*(b+1)*(((b+2)*(1-2*z0)+(a+b+1)*z0-c)...
    *(1-a+b)^(-1)*b*((b+1)*(1-2*z0)+(a+b+1)*z0-c)+z0*(1-z0)*b);

% Initialise values of current terms and sum of all terms in second series
% of (4.21) computed so far
a2=1+e(1)*(z-z0)^(-1)+e(2)*(z-z0)^(-2);
b2=a2;
for n=3:500
    % Update e(n), a2 and b2
    e(n)=(n*(n-a+b))^(-1)*(n+b-1)*...
        (((n+b)*(1-2*z0)+(a+b+1)*z0-c)*e(n-1)+z0*(1-z0)*(n+b-2)*e(n-2));
    a2=e(n)*(z-z0)^(-n);
    b2=b2+a2;
    % Terminate summation if stopping criterion is satisfied
    if abs(a2)/abs(b2)<tol
        break
    end
    % If 500 terms have been computed before stopping criterion is
    % satisfied, state this
    if (n==500)
        [' ' num2str(n) ' terms computed']
        return
    end
end

% Apply formula (4.21) in terms of two series evaluated
h=gamma(c)*(gamma(b-a)/gamma(b)/gamma(c-a)*(z0-z)^(-a)*b1...
    +gamma(a-b)/gamma(a)/gamma(c-b)*(z0-z)^(-b)*b2);