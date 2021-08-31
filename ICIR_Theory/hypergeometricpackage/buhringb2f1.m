function [h]=buhringb2f1(ar,ai,cr,ci,zr,zi,z0,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the Gauss hypergeometric function 2F1(a,a;c;z) using the       %
% analytic continuation formula (4.23) of Section 4.7. This is the code   %
% used for real a,c; code for complex a,c would be obtained by replacing  %
% gamma.m with cgama.m [71], and s14ae.m with s14af from the NAG Toolbox. %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         z0=Value at which expansion is taken (usually 1/2)              %
%         tol=Specified tolerance as detailed in Section 4.2              %
% Output: h=Computed value of 2F1(a,b;c;z)                                %
%-------------------------------------------------------------------------%

% Compute a,c,z in terms of ar,ai,cr,ci,zr,zi
a=ar+ai*1i;
c=cr+ci*1i;
z=zr+zi*1i;

% Initialise ej and fj
e=zeros(2,1);
e(1)=(a+1)*(1-2*z0)+(2*a+1)*z0-c;
e(2)=1/2*(((a+2)*(1-2*z0)+(2*a+1)*z0-c)*((a+1)*(1-2*z0)+(2*a+1)*z0-c)...
    +z0*(1-z0));
f=zeros(2,1);
f(1)=1-2*z0;
f(2)=1/2*(((a+2)*(1-2*z0)+(2*a+1)*z0-c)*(1-2*z0)+(1-2*z0)...
    *((a+1)*(1-2*z0)+(2*a+1)*z0-c)+2*z0*(1-z0));

% Initialise values of current term and value of the sum of all terms
a1=(2*s14ae(1,int32(0))-s14ae(a,int32(0))...
    -s14ae(c-a,int32(0))+log(z0-z))...
    +a*(e(1)*(2*s14ae(2,int32(0))-s14ae(a+1,int32(0))...
    -s14ae(c-a,int32(0))+log(z0-z))-f(1))*(z0-z)^(-1)...
    +a*(a+1)/2*(e(2)*(2*s14ae(3,int32(0))...
    -s14ae(a+2,int32(0))...
    -s14ae(c-a,int32(0))+log(z0-z))-f(2))*(z0-z)^(-2);
b1=a1;
pochgam=a*(a+1)/2;
for n=3:500
    % Update all terms
    pochgam=pochgam*(a+n-1)/n;
    e(n)=1/n*(((a+n)*(1-2*z0)+(2*a+1)*z0-c)*e(n-1)+z0*(1-z0)*(n-1)*e(n-2));
    f(n)=1/n*(((a+n)*(1-2*z0)+(2*a+1)*z0-c)*f(n-1)+z0*(1-z0)*(n-1)*f(n-2)...
        +(1-2*z0)*e(n-1)+2*z0*(1-z0)*e(n-2));
    a1=pochgam*(e(n)*(2*s14ae(1+n,int32(0))-s14ae(a+n,int32(0))...
        -s14ae(c-a,int32(0))+log(z0-z))-f(n))*(z-z0)^(-n);
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

% Compute expansion
[kr1,ki1]=cgama(ar,ai,1);k1=kr1+ki1*1i;
[kr2,ki2]=cgama(cr,ci,1);k2=kr2+ki2*1i;
[kr3,ki3]=cgama(cr-ar,ci-ai,1);k3=kr3+ki3*1i;
h=gamma(c)/gamma(a)/gamma(c-a)*(z0-z)^(-a)*b1;