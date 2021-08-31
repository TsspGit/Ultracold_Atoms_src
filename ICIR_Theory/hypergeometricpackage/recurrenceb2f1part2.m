function [r]=recurrenceb2f1part2(a,b,c,z,n,k)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the minimal solution of recurrence relation 2 of Section 4.8,  %
% forwards, using Miller's algorithm, outside the curve shown in Figure 6.%

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         z=Variable value z                                              %
%         k=Number of digits forward we wish to compute the recurrence    %
%           relation                                                      %
%         n=Number >k sufficiently large to generate the minimal solution %
% Output: r=Computed value of 2F1(a+k,b+k;1+a+b-c+3k;1-z)                 %
%-------------------------------------------------------------------------%

% Initialise f and v in Miller's algorithm in Section 4.8
f=zeros(k+1,1);
v=zeros(n+1,1);
a1=zeros(n-1,1);
b1=zeros(n-1,1);

% Define coefficients in recurrence relation
for i1=1:n-1
    a1(i1)=-(a-c+2*i1)*(a-c+2*i1-1)*(b-c+2*i1-1)*(b-c+2*i1)*z...
        *(z*(a+b-c+3*i1+1)*(a+b-c+3*i1+2)+(a+i1)*(b+i1)*(1-z))...
        /((a+i1)*(b+i1)*(c-i1)*(c-i1-1)...
        *(1-z)^3*((1-z)*(1-a-b-2*i1+(a+i1)*(b+i1))+z*(a+b-c+3*i1-1)*(a+b-c+3*i1-2)));
    b1(i1)=(c-i1)*(((1-z)*(b-c+2*i1)*(b+i1-1)*(a+i1-1+z*(b-c+2*i1-1)))...
        *(z*(a+b-c+3*i1+1)*(a+b-c+3*i1+2)+(a+i1)*(b+i1)*(1-z))...
        +((b+i1)*(b+1-c+2*i1)*(1-z)*(a+i1+z*(b-c+2*i1+2)))...
        *((1-z)*(1-a-b-2*i1+(a+i1)*(b+i1))+z*(a+b-c+3*i1-1)*(a+b-c+3*i1-2))...
        +(c-2*b-3*i1-(a-b)*z)*(z*(a+b-c+3*i1+1)*(a+b-c+3*i1+2)+(a+i1)*(b+i1)*(1-z))...
        *((1-z)*(1-a-b-2*i1+(a+i1)*(b+i1))+z*(a+b-c+3*i1-1)*(a+b-c+3*i1-2)))...
        /((a+i1)*(b+i1)*(c-i1)*(c-i1-1)*(1-z)^3*...
        ((1-z)*(1-a-b-2*i1+(a+i1)*(b+i1))+z*(a+b-c+3*i1-1)*(a+b-c+3*i1-2)));
end

% Input minimal solution with k=0
f1=gamma(1+b-c)*gamma(1+a-c)/gamma(1-c)/gamma(1+a+b-c)*hypergeom([a,b],1+a+b-c,1-z);%gausshypergeometricseries(a,0,b,0,c,0,z,0,1e-15);
v(end)=0;v(end-1)=1;

% Compute recurrence backwards
for i2=2:n
    v(n+1-i2)=-(v(n+3-i2)+b1(n+1-i2)*v(n+2-i2))/a1(n+1-i2);
end

% Apply last line of Miller's algorithm and return f(end)
for i3=1:k+1
    f(i3)=f1/v(1)*v(i3)*gamma(-c+i3)*gamma(1+a+b-c+3*(i3-1))/gamma(1+b-c+2*(i3-1))/gamma(1+a-c+2*(i3-1));
end

r=f(end);