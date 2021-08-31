function [r]=recurrenced2f1part1(a,b,c,z,n,k)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the minimal solution of recurrence relation 4 of Section 4.8,  %
% forwards, using Miller's algorithm, in the region Re(z)<1/2.            %

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         z=Variable value z                                              %
%         k=Number of digits forward we wish to compute the recurrence    %
%           relation                                                      %
%         n=Number >k sufficiently large to generate the minimal solution %
% Output: r=Computed value of 2F1(a,b;c+k;z)                              %
%-------------------------------------------------------------------------%

% Initialise f and v in Miller's algorithm in Section 4.8
f=zeros(k+1,1);
v=zeros(n+1,1);
a1=zeros(n-1,1);
b1=zeros(n-1,1);

% Define coefficients in recurrence relation
for i1=1:n-1
    a1(i1)=(c+i1)*(c+i1-1)*(z-1)/(c+i1-a)/(c+i1-b)/z;
    b1(i1)=(c+i1)*(c+i1-1-(2*(c+i1)-a-b-1)*z)/(c+i1-a)/(c+i1-b)/z;
end

% Input minimal solution with k=0
f1=hypergeom([a,b],c,z);
v(end)=0;v(end-1)=1;

% Compute recurrence backwards
for i2=2:n
    v(n+1-i2)=-(v(n+3-i2)+b1(n+1-i2)*v(n+2-i2))/a1(n+1-i2);
end

% Apply last line of Miller's algorithm and return f(end)
for i3=1:k+1
    f(i3)=f1/v(1)*v(i3);
end

r=f(end);