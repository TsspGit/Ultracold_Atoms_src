function [r]=recurrencec2f1part3(a,b,c,z,n,k)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the minimal solution of recurrence relation 3 of Section 4.8,  %
% forwards, using Miller's algorithm, outside curves shown in Figure 6.   %

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         z=Variable value z                                              %
%         k=Number of digits forward we wish to compute the recurrence    %
%           relation                                                      %
%         n=Number >k sufficiently large to generate the minimal solution %
% Output: r=Computed value of 2F1(1-b,-b+c-k;1+a-b+k;1/z)                 %
%-------------------------------------------------------------------------%

% Initialise f and v in Miller's algorithm in Section 4.8
f=zeros(k+1,1);
v=zeros(n+1,1);
a1=zeros(n-1,1);
b1=zeros(n-1,1);

% Define coefficients in recurrence relation
for i1=1:n-1
    a1(i1)=-z*(a-c+2*i1)*(a-c+2*i1-1)*(b-c+i1)*(a+i1+z*(b+1-c+i1))...
        /(a+i1)/(c-i1)/(c-i1-1)/(a+i1-1+z*(b-c+i1))/(1-z)^2;
    b1(i1)=-(c-i1)*((a+i1)*(a+i1-1)*(c-i1-1)+(a+i1)*(a+i1-1)*(a+3*b-4*c+5*i1+2)*z...
        +(b-c+i1)*(b-c+i1+1)*(4*a-c+5*i1-1)*z^2-(a-b+i1)*(b-c+i1)*(b-c+i1+1)*z^3)...
        /(a+i1)/(c-i1)/(c-i1-1)/(a+i1-1+z*(b-c+i1))/(1-z)^2;
end

% Input minimal solution with k=0
f1=gamma(1+a-c)/gamma(1-c)/gamma(1+a-b)*hypergeom([1-b,-b+c],1+a-b,1/z);
v(end)=0;v(end-1)=1;

% Compute recurrence backwards
for i2=2:n
    v(n+1-i2)=-(v(n+3-i2)+b1(n+1-i2)*v(n+2-i2))/a1(n+1-i2);
end

% Apply last line of Miller's algorithm and return f(end)
for i3=1:k+1
    f(i3)=f1/v(1)*v(i3)/(-z/(1-z)^2)^(i3-1)*gamma(-c+i3)*gamma(a-b+i3)...
        /gamma(1+a-c+2*(i3-1));
end

r=f(end);