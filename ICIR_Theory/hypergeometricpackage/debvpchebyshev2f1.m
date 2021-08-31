function [ch]=debvpchebyshev2f1(ar,ai,br,bi,cr,ci,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the Gauss hypergeometric function 2F1(a,b;c;z) on the range    %
% [-1,1] using a Chebyshev differentiation matrix method of solving the   %
% boundary value problem from the differential equation (4.2).            %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         n =Number of mesh points used                                   %
% Output: ch=Profile of 2F1(a,b;c;z) on [-1,1]                            %
%-------------------------------------------------------------------------%

% Compute a,b,c using ar,ai,br,bi,cr,ci
a=ar+ai*1i;
b=br+bi*1i;
c=cr+ci*1i;

% Set up mesh points
x = -cos((0:n-1)*pi/(n-1))';
A = zeros(n,n);
B = A;

% Set up matrices of Chebyshev points and their derivatives
for k = 0:n-1
    A(:,k+1) = cos(acos(x)*k);
    B(:,k+1) = k*sin(acos(x)*k)./sqrt(1-x.^2);
    B(1,k+1) = (-1)^(k+1)*k^2;
    B(end,k+1) = k^2;
end

% Solve matrix system to find differentiation matrix
D = B/A;

E=zeros(n,n);
F=zeros(n,n);

% Set up coefficient matrices
for j=1:n
    E(:,j)=x;
end
for j=1:n
    F(:,j)=x.*(1-x);
end

% Set up Mch for differential equation
Mch=F.*D^2+(c-(a+b+1)*E).*D-a*b*eye(n);
Mch(1,:) = [1 zeros(1,n-1)];
Mch(n,:) = [zeros(1,n-1) 1];

% Set up RHS matrix corresponding to boundary conditions
d = zeros(n,1);
d(1)=hypergeom([a,b],c,-1);
d(end)=hypergeom([a,b],c,1);

% Solve system to obtain profile
ch=Mch\d;
plot(x,ch)